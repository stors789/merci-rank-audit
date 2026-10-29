"""Regenerate a frozen composition split from its strict mixed and negative scores."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from core import digest_ids, sha
from run import task_arrays


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--experiment', required=True)
    p.add_argument('--split', type=int, required=True)
    p.add_argument('--mixed', type=Path, required=True)
    p.add_argument('--negative', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--output-manifest', type=Path, help='Write a private manifest rebound to the regenerated file')
    a = p.parse_args(); config = json.loads(a.manifest.read_text())
    spec = next(s for s in config['splits'] if s['experiment'] == a.experiment and s['split'] == a.split)
    def read(path, target):
        item = next(c for c in config['cases'] if c['experiment'] == a.experiment and c['split'] == a.split
                    and c['case'] == 'strict_regular_' + target)
        if sha(path) != item['score_sha256']:
            raise ValueError('strict score differs from frozen input: ' + str(path))
        raw = task_arrays(path)
        frame = pd.read_csv(path, sep='\t', compression='infer', dtype={'cell': str}, float_precision='round_trip')
        raw['rna_score'] = pd.to_numeric(frame.RNA_score, errors='raise').to_numpy(float)
        raw['rna_tied'] = pd.Series(raw['rna_score']).duplicated(keep=False).to_numpy()
        if digest_ids(raw['ids']) != item['ordered_ID_sha256']:
            raise ValueError('frozen ordered IDs differ')
        return raw
    mixed, negative = read(a.mixed, 'mixed'), read(a.negative, 'negative')
    idx = pd.Index(mixed['ids']).get_indexer(negative['ids'])
    if (negative['truth'].any() or (idx < 0).any() or len(negative['ids']) != int((~mixed['truth']).sum())
            or set(negative['ids']) != set(mixed['ids'][~mixed['truth']])
            or not np.array_equal(mixed['neff'][idx], negative['neff'])
            or not np.array_equal(mixed['rna_score'][idx], negative['rna_score'])):
        raise ValueError('strict MC input IDs or saved scores differ')
    mc = np.flatnonzero(~mixed['truth']); cc = np.flatnonzero(mixed['truth'])
    mc = mc[np.argsort(mixed['ids'][mc], kind='stable')]; cc = cc[np.argsort(mixed['ids'][cc], kind='stable')]
    m = len(mc); anchor_n = m // 2
    cc_order = cc[np.random.default_rng(spec['CC_seed']).permutation(len(cc))]
    anchor_perm = np.random.default_rng(spec['anchor_seed']).permutation(m)
    anchor = mc[anchor_perm[:anchor_n]]; remaining = np.sort(mc[anchor_perm[anchor_n:]])
    remaining = remaining[np.argsort(mixed['ids'][remaining], kind='stable')]
    background = remaining[np.random.default_rng(spec['background_seed']).permutation(len(remaining))]
    untied = np.zeros(len(mixed['ids']), bool); untied[idx] = ~mixed['rna_tied'][idx] & ~negative['rna_tied']
    for label, values in (('CC_order', cc_order), ('anchor', anchor), ('background', background)):
        if digest_ids(mixed['ids'][values]) != spec[label + '_ID_sha256']:
            raise ValueError('frozen selection ID digest differs: ' + label)
    if int(untied.sum()) != spec['untied_MC_N']:
        raise ValueError('untied endpoint count differs')
    if a.output.exists() or (a.output_manifest and a.output_manifest.exists()):
        raise FileExistsError('output already exists')
    a.output.parent.mkdir(parents=True, exist_ok=True)
    keys = ('ids','neff','rna_rank','dna_rank','rna_score','truth','observed_source_features')
    np.savez_compressed(a.output, **{k: mixed[k] for k in keys}, cc_order=cc_order,
                        anchor=anchor, background=background, untied_endpoints=untied)
    actual = sha(a.output)
    if a.output_manifest:
        # These paths are private local runtime paths; do not redistribute this generated manifest.
        base = a.manifest.resolve().parent
        for item in config['cases']:
            item['source'] = str((base / item['source']).resolve())
        for item in config['splits']:
            item['prepared'] = str((base / item['prepared']).resolve())
        config['null_maxima'] = str((base / config['null_maxima']).resolve())
        spec['historical_prepared_sha256'] = spec.get('historical_prepared_sha256', spec.get('prepared_sha256'))
        spec['prepared_sha256'] = actual; spec['prepared'] = str(a.output.resolve())
        a.output_manifest.parent.mkdir(parents=True, exist_ok=True)
        a.output_manifest.write_text(json.dumps(config, indent=2) + '\n')
    print(json.dumps({'status':'COMPLETE','N':len(mixed['ids']),'M':m,'sha256':actual}))

if __name__ == '__main__':
    main()
