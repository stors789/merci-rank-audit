"""Bind an already prepared real score/threshold pair into a benchmark NPZ."""
import argparse
import json
from pathlib import Path
import numpy as np
from merci_rank_audit.inputs import validate_manifest
from case_inputs import digest_arrays, sha_file


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--case', required=True)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--output-npz', type=Path, required=True)
    p.add_argument('--output-config', type=Path, required=True)
    a = p.parse_args(); config = json.loads(a.config.read_text())
    spec = next(x for x in config['cases'] if x['case'] == a.case)
    _, cases = validate_manifest(a.manifest)
    if len(cases) != 1 or spec['scope'] != 'real':
        raise ValueError('provide one real case in the input manifest')
    case = cases[0]; sampler = case['sampler']; identity = case['input_identity']
    if identity['score'] != spec['score_sha256'] or identity['thresholds'] != spec['thresholds_sha256']:
        raise ValueError('score/threshold file bytes differ from the historical frozen inputs')
    ids, neff, rna, maxima = sampler.ids, sampler.neff, sampler.rna, sampler.maxima
    if digest_arrays(ids, neff, rna, maxima) != spec['canonical_array_sha256']:
        raise ValueError('canonical arrays differ from the historical input')
    if a.output_npz.exists() or a.output_config.exists():
        raise FileExistsError('output already exists')
    a.output_npz.parent.mkdir(parents=True, exist_ok=True); a.output_config.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.output_npz, ids=ids, neff=neff, rna=rna, maxima=maxima)
    # Absolute paths here are local runtime inputs; do not redistribute this generated config.
    spec['prepared_input'] = str(a.output_npz.resolve()); spec['prepared_sha256'] = sha_file(a.output_npz)
    inputs = json.loads(a.manifest.read_text())['cases'][0]
    spec['score'] = str((a.manifest.parent / inputs['score']).resolve())
    spec['thresholds'] = str((a.manifest.parent / inputs['thresholds']).resolve())
    # Keep other case paths relative to the source config, not the new config.
    for item in config['cases']:
        if item is not spec:
            for field in ('prepared_input','score','thresholds'):
                if field in item: item[field] = str((a.config.parent / item[field]).resolve())
    a.output_config.write_text(json.dumps(config, indent=2) + '\n')
    print(json.dumps({'status':'COMPLETE','case':a.case,'prepared_sha256':spec['prepared_sha256']}))

if __name__ == '__main__':
    main()
