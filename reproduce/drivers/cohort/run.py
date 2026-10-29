"""Run original exact-family annotation and fixed-N/nested cohort comparisons.

The manifest binds private cell-level inputs explicitly. Outputs containing
member IDs must remain private. Run from any directory.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd

from core import confusion, digest_ids, model_case, ratio, sha, write_members
from design import LEVELS, condition_rows, compare_levels, make_condition


def write_rows(path, rows):
    if rows:
        with path.open('x', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys(), lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)


def task_arrays(path):
    if path.suffix == '.npz':
        with np.load(path, allow_pickle=False) as src:
            return {k: src[k] for k in src.files}
    df = pd.read_csv(path, sep='\t', compression='infer', dtype={'cell': str})
    required = {'cell', 'Neff', 'RNA_rank', 'DNA_rank', 'truth', 'observed_source_features'}
    if not required.issubset(df.columns):
        raise ValueError(f'missing score fields: {sorted(required-set(df.columns))}')
    truth = df.truth.astype(str).str.lower()
    if not truth.isin(['true', 'false', '1', '0']).all():
        raise ValueError('truth must be true/false or 1/0')
    return {'ids': df.cell.to_numpy(str), 'neff': df.Neff.to_numpy(),
            'rna_rank': df.RNA_rank.to_numpy(), 'dna_rank': df.DNA_rank.to_numpy(),
            'truth': truth.isin(['true', '1']).to_numpy(),
            'observed_source_features': pd.to_numeric(df.observed_source_features, errors='coerce').to_numpy()}


def run_task(spec, base, members):
    path = (base / spec['source']).resolve()
    if spec.get('score_sha256') and sha(path) != spec['score_sha256']:
        raise ValueError(f'original score digest differs for {path}')
    arrays = task_arrays(path)
    if spec.get('ordered_ID_sha256') and digest_ids(arrays['ids']) != spec['ordered_ID_sha256']:
        raise ValueError(f'ordered ID digest differs for {path}')
    model = model_case(arrays, np.asarray(spec['maxima'], dtype=int))
    if spec.get('N') is not None and model['N'] != spec['N']:
        raise ValueError('case N differs from frozen input')
    if spec.get('original_gate') is not None and model['gate'] != spec['original_gate']:
        raise ValueError('supplied-order gate differs from frozen input')
    truth = model['truth']
    identity = {key: spec[key] for key in ('experiment', 'stage', 'split', 'case')}
    common = dict(**identity, N=model['N'], CC=int(truth.sum()), MC=int((~truth).sum()),
                  original_gate=model['gate'], gate_can_pass=model['can'], gate_must_pass=model['must'])
    policies, partitions, errors = [], [], []
    for name, selected in (('original_raw', model['raw_call']), ('original_default', model['default_call']),
                           ('raw_core', model['raw_state'] == 0), ('default_core', model['default_state'] == 0),
                           ('conditional_core', None if model['conditional_state'] is None else model['conditional_state'] == 0)):
        policies.append(dict(**common, policy=name, status='NA_G_IMPOSSIBLE' if selected is None else 'DEFINED',
                             **confusion(selected, truth)))
    for family in ('raw', 'default', 'conditional'):
        state = model[family + '_state']
        for i, label in enumerate(('core', 'variable', 'outside')):
            for truth_label, truth_mask in (('CC', truth), ('MC', ~truth)):
                if state is None:
                    counts = dict(count=None, original_raw_selected=None, original_default_selected=None, Neff_zero=None,
                                  observed=None, unobserved=None, observation_unknown=None)
                else:
                    mask = (state == i) & truth_mask
                    obs = model['observed']; finite = np.isfinite(obs)
                    counts = dict(count=int(mask.sum()), original_raw_selected=int((mask & model['raw_call']).sum()),
                                  original_default_selected=int((mask & model['default_call']).sum()),
                                  Neff_zero=int((mask & (model['neff'] == 0)).sum()),
                                  observed=int((mask & finite & (obs > 0)).sum()),
                                  unobserved=int((mask & finite & (obs == 0)).sum()),
                                  observation_unknown=int((mask & ~finite).sum()))
                partitions.append(dict(**common, family=family, partition=label, truth=truth_label,
                                       status='NA_G_IMPOSSIBLE' if state is None else 'DEFINED', **counts))
        for baseline in ('raw', 'default'):
            call = model[baseline + '_call']; fp = call & ~truth; fn = ~call & truth
            core_fp = None if state is None else int((fp & (state == 0)).sum())
            outside_fn = None if state is None else int((fn & (state == 2)).sum())
            errors.append(dict(**common, family=family, baseline='original_' + baseline,
                               status='NA_G_IMPOSSIBLE' if state is None else 'DEFINED',
                               original_FP=int(fp.sum()), original_FN=int(fn.sum()), FP_in_core=core_fp,
                               FN_outside_union=outside_fn, FP_core_fraction=ratio(core_fp, int(fp.sum())) if core_fp is not None else None,
                               FN_outside_fraction=ratio(outside_fn, int(fn.sum())) if outside_fn is not None else None))
    conditional, raw = model['conditional_state'], model['raw_state']
    diff = dict(**common, status='NA_G_IMPOSSIBLE' if conditional is None else 'DEFINED',
                core_same=None, union_same=None, core_added=None, core_added_CC=None, core_added_MC=None,
                union_removed=None, union_removed_CC=None, union_removed_MC=None)
    if conditional is not None:
        added = (raw != 0) & (conditional == 0); removed = (raw != 2) & (conditional == 2)
        diff.update(core_same=not bool(added.any()), union_same=not bool(removed.any()), core_added=int(added.sum()),
                    core_added_CC=int((added & truth).sum()), core_added_MC=int((added & ~truth).sum()),
                    union_removed=int(removed.sum()), union_removed_CC=int((removed & truth).sum()),
                    union_removed_MC=int((removed & ~truth).sum()))
    member = members / f"{spec['experiment']}_{spec['split']:03d}_{spec['case']}.npz"
    digest = write_members(model, member)
    record = dict(**common, member_output_sha256=digest, original_counts=json.dumps(model['original_counts'].tolist()),
                  lower_counts=json.dumps(model['lower_counts'].tolist()), upper_counts=json.dumps(model['upper_counts'].tolist()),
                  maxima=json.dumps(model['maxima'].tolist()))
    return record, policies, partitions, errors, diff


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); data = json.loads(a.manifest.read_text()); base = a.manifest.resolve().parent
    a.output.mkdir(parents=True, exist_ok=False); members = a.output / 'private_members'; members.mkdir(mode=0o700)
    results = [run_task(s, base, members) for s in data['cases']]
    for i, name in enumerate(('task_summary', 'policy_metrics', 'partition_truth', 'error_fractions', 'family_differences')):
        write_rows(a.output / (name + '.csv'), [row for result in results for row in (result[i] if isinstance(result[i], list) else [result[i]])])
    null = pd.read_csv((base / data['null_maxima']).resolve())
    null_map = {int(n): group.set_index('cutoff').loc[list(range(100, 0, -10)), 'maximum'].to_numpy(int)
                for n, group in null.groupby('N')}
    tables = {key: [] for key in ('conditions', 'transitions', 'channels', 'joint', 'calls')}
    for split in data['splits']:
        prepared = (base / split['prepared']).resolve()
        if split.get('prepared_sha256') and sha(prepared) != split['prepared_sha256']:
            raise ValueError('prepared split digest differs; regenerate and rebind it with prepare_split.py')
        with np.load(prepared, allow_pickle=False) as src:
            arrays = {k: src[k] for k in src.files}
        for label, field in (('CC_order','cc_order'), ('anchor','anchor'), ('background','background')):
            expected = split.get(label + '_ID_sha256')
            if expected and digest_ids(arrays['ids'][arrays[field]]) != expected:
                raise ValueError('frozen selection ID digest differs: ' + label)
        m = split['M']; cc = len(arrays['cc_order']); background = len(arrays['background'])
        if (m != len(arrays['anchor']) + background or len(arrays['anchor']) != split['anchor_N']
                or background != split['background_slots'] or cc != split['CC_available']):
            raise ValueError('composition counts differ from frozen manifest')
        indices = np.concatenate((arrays['anchor'], arrays['background'], arrays['cc_order']))
        if len(np.unique(indices)) != len(arrays['ids']) or len(indices) != len(arrays['ids']):
            raise ValueError('composition pools must partition the parent')
        for mode in ('nested', 'fixedN'):
            models = {}
            for level in LEVELS:
                n = m + cc * level // 100 if mode == 'nested' else m
                model = make_condition(arrays, split, mode, level, null_map[n], members)
                models[level] = model
                tables['conditions'].extend(condition_rows(split, mode, level, model))
            for left, right in ((0,25), (25,50), (50,75), (75,100), (0,100)):
                for key, rows in zip(('transitions','channels','joint','calls'),
                                     compare_levels(split, mode, left, right, models[left], models[right])):
                    tables[key].extend(rows)
    for name, rows in tables.items():
        write_rows(a.output / (name + '.csv'), rows)
    receipt = {'tasks': len(results), 'splits': len(data['splits']), 'conditions': len(data['splits']) * 10,
               'data_scope': 'user-provided; member outputs contain IDs', 'status': 'COMPLETE'}
    (a.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))

if __name__ == '__main__':
    main()
