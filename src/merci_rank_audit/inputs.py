"""Validate explicit score, matched thresholds, optional nuclear inputs, and identity."""
import hashlib
import io
import json
import math
import re
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_distribution, binom

from .sampler import CUTOFFS, GateSampler, exact_integers

INTEGER_CONTRACT = 'exact-integer-v2'


def fraction_decimal(value):
    with localcontext() as ctx:
        ctx.prec = 40
        return Decimal(value.numerator) / Decimal(value.denominator)


def ordinary_interval(hits, n):
    low = 0.0 if hits == 0 else float(beta_distribution.ppf(.025, hits, n-hits+1))
    high = 1.0 if hits == n else float(beta_distribution.ppf(.975, hits+1, n-hits))
    return low, high


def percentile_with_MC_interval(values, q):
    x = np.sort(values[np.isfinite(values)]); n = len(x)
    if n == 0:
        return (None, None, None)
    lo = int(binom.ppf(.025, n, q)); hi = int(binom.ppf(.975, n, q)) + 1
    return (float(np.quantile(x, q)), float(x[lo-1]) if lo >= 1 else None,
            float(x[hi-1]) if hi <= n else None)


def required_integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return value


def input_path(base, name):
    if not isinstance(name, str) or not name:
        raise ValueError('Every input must be a nonempty explicit path')
    path = (base / name).resolve()
    if not path.is_file():
        raise ValueError(f'Input file does not exist: {path}')
    return path


def legacy_safe_integers(values):
    """Unversioned samples are reusable only for digit-only small integer tokens."""
    return all(re.fullmatch(r'[0-9]+', value) and Decimal(value) <= 2**31-1 for value in values)


def validate_nuclear(path, sampler, score):
    with np.load(path, allow_pickle=False) as source:
        needed = {'ids', 'observable', 'X', 'original', 'alternative', 'names'}
        if not needed.issubset(source.files):
            raise ValueError(f'Nuclear input missing fields: {sorted(needed-set(source.files))}')
        data = {key: source[key] for key in needed}
    ids = data['ids']; observed = data['observable']; names = data['names']
    if ids.ndim != 1 or ids.dtype.kind not in 'US' or len(ids) != sampler.n:
        raise ValueError('Nuclear ids must be a one-dimensional string array of length N')
    if np.any(ids == '') or len(np.unique(ids)) != len(ids):
        raise ValueError('Nuclear IDs must be nonempty and unique')
    index = pd.Index(sampler.ids).get_indexer(ids)
    if (index < 0).any() or len(np.unique(index)) != sampler.n:
        raise ValueError('Nuclear IDs must match all score IDs exactly')
    if observed.dtype != np.bool_ or observed.shape != (sampler.n,):
        raise ValueError('Nuclear observable must be a Boolean array of length N')
    if 'observable' in score:
        flags = score.set_index('cell').loc[ids, 'observable'].astype(str).str.lower()
        if not flags.isin(['true','false','1','0']).all():
            raise ValueError('Score observable must contain Boolean values')
        if not np.array_equal(flags.isin(['true','1']).to_numpy(), observed):
            raise ValueError('Score and nuclear observable flags differ')
    if names.ndim != 1 or names.dtype.kind not in 'US' or not len(names) or np.any(names == '') or len(np.unique(names)) != len(names):
        raise ValueError('Nuclear outcome names must be unique nonempty strings')
    nobs = int(observed.sum())
    if nobs < 1:
        raise ValueError('Nuclear observable domain is empty')
    for key in ('X','original','alternative'):
        value = data[key]
        if value.ndim != 2 or value.shape[0] != nobs or value.shape[1] < 1 or value.dtype.kind not in 'fiub' or not np.isfinite(value).all():
            raise ValueError(f'Nuclear {key} must be a finite numeric matrix with observable rows')
    if data['original'].shape[1] != len(names) or data['alternative'].shape[1] != len(names):
        raise ValueError('Nuclear outcome columns must match names')
    data['index'] = index
    return data


def validate_manifest(path):
    manifest = json.loads(path.read_text())
    budget = required_integer(manifest.get('budget'), 'budget', 1)
    cases = manifest.get('cases')
    if not isinstance(cases, list) or not cases:
        raise ValueError('cases must be a nonempty array')
    validated = []; names = set()
    for case in cases:
        if not isinstance(case, dict):
            raise ValueError('Each case must be an object')
        name = case.get('name')
        if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9_-]+', name) or name in names:
            raise ValueError('Case names must be unique and contain only letters, digits, _ or -')
        names.add(name)
        score_path = input_path(path.parent, case.get('score'))
        threshold_path = input_path(path.parent, case.get('thresholds'))
        score_bytes = score_path.read_bytes()
        threshold_bytes = threshold_path.read_bytes()
        score = pd.read_csv(io.BytesIO(score_bytes), sep='\t', dtype={'cell':str,'Neff':str,'RNA_rank':str,'DNA_rank':str}, keep_default_na=False,
                            compression=pd.io.common.infer_compression(str(score_path), 'infer'))
        if not {'cell','Neff','RNA_rank'}.issubset(score.columns):
            raise ValueError(f'{name}: score requires cell, Neff and RNA_rank')
        if score.cell.eq('').any():
            raise ValueError(f'{name}: cell IDs must not be empty')
        integer_fields = ['Neff','RNA_rank'] + (['DNA_rank'] if 'DNA_rank' in score else [])
        legacy_safe = all(legacy_safe_integers(score[col]) for col in integer_fields)
        for col in integer_fields:
            score[col] = exact_integers(score[col], f'{name}: {col}')
        thresholds = pd.read_csv(io.BytesIO(threshold_bytes), dtype={'N':str,'cutoff':str,'maximum':str}, keep_default_na=False,
                                 compression=pd.io.common.infer_compression(str(threshold_path), 'infer'))
        if not {'N','cutoff','maximum'}.issubset(thresholds.columns):
            raise ValueError(f'{name}: thresholds require N, cutoff and maximum')
        legacy_safe = legacy_safe and all(legacy_safe_integers(thresholds[col]) for col in ('N','cutoff','maximum'))
        for col in ('N','cutoff','maximum'):
            thresholds[col] = exact_integers(thresholds[col], f'{name}: threshold {col}')
        if len(thresholds) != 10 or thresholds.cutoff.duplicated().any() or set(thresholds.cutoff) != set(CUTOFFS):
            raise ValueError(f'{name}: exactly the ten distinct cutoffs 100..10 are required')
        if not thresholds.N.eq(len(score)).all():
            raise ValueError(f'{name}: every threshold N must equal the score row count')
        maxima = thresholds.set_index('cutoff').loc[CUTOFFS, 'maximum'].to_numpy()
        sampler = GateSampler(score.cell.to_numpy(str), score.Neff.to_numpy(), score.RNA_rank.to_numpy(), maxima)
        original_counts = None
        if 'DNA_rank' in score:
            rank = score.DNA_rank.to_numpy()[sampler.order]
            if not sampler.validate_order(rank):
                raise ValueError(f'{name}: DNA_rank is not a legal rank for Neff')
            original_counts = sampler.counts(rank)
        observable = None
        if 'observable' in score:
            flags = score.observable.astype(str).str.lower()
            if not flags.isin(['true','false','1','0']).all():
                raise ValueError(f'{name}: score observable must contain Boolean values')
            observable = flags.isin(['true','1']).to_numpy()[sampler.order]
        nuclear_bytes = input_path(path.parent, case['nuclear']).read_bytes() if 'nuclear' in case else None
        nuclear = validate_nuclear(io.BytesIO(nuclear_bytes), sampler, score) if nuclear_bytes is not None else None
        identity = dict(scheme='sha256-file-bytes-v1', integer_contract=INTEGER_CONTRACT, score=hashlib.sha256(score_bytes).hexdigest(),
                        thresholds=hashlib.sha256(threshold_bytes).hexdigest(),
                        nuclear=hashlib.sha256(nuclear_bytes).hexdigest() if nuclear_bytes is not None else None)
        if nuclear is not None:
            observable = np.empty(sampler.n, dtype=bool)
            observable[nuclear['index']] = nuclear['observable']
        validated.append(dict(name=name, score=score, sampler=sampler, nuclear=nuclear, observable=observable, input_identity=identity, legacy_integer_safe=legacy_safe,
            original_counts=original_counts, seed_uniform=required_integer(case.get('seed_uniform'), 'seed_uniform'),
            seed_conditional=required_integer(case.get('seed_conditional'), 'seed_conditional')))
    return budget, validated


def require_identity(case, metadata, require_nuclear=False):
    """Reject sample reuse unless identity came from the bytes parsed at sampling time."""
    saved = metadata.get('input_identity')
    current = case['input_identity']
    keys = ('score', 'thresholds', 'nuclear') if require_nuclear else ('score', 'thresholds')
    if not isinstance(saved, dict) or saved.get('scheme') != current['scheme']:
        raise ValueError(f"{case['name']}: saved sample has no supported input identity")
    if any(not isinstance(saved.get(key), str) or saved[key] != current[key] for key in keys):
        raise ValueError(f"{case['name']}: saved sample input identity does not match {', '.join(keys)}")
    version = saved.get('integer_contract')
    if version != INTEGER_CONTRACT and not (version is None and case['legacy_integer_safe']):
        raise ValueError(f"{case['name']}: saved sample integer parsing is not known to match")


