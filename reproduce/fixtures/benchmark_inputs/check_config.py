"""Check all historical split/case selectors without running biological inputs."""
import importlib.util
import json
from pathlib import Path

repro = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('benchmark_inputs', repro / 'drivers/benchmark_inputs/run.py')
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
config = json.loads((repro / 'configs/cohort_formal_portable.json').read_text())
assert len(config['cases']) == 720 and len(config['splits']) == 120
seeds = set()
for experiment, end in (('bench1', 100), ('bench2', 20)):
    for split in range(1, end + 1):
        cases, splits = module.select_config(config, experiment, split)
        parent = splits[0]
        m, cc = (3786, 2904) if experiment == 'bench1' else (1272, 1151)
        assert parent['M'] == m and parent['CC_available'] == cc
        for field in ('CC_seed', 'anchor_seed', 'background_seed'):
            assert parent[field] not in seeds
            seeds.add(parent[field])
        for case in cases:
            expected_cc = cc if case['case'].endswith('_mixed') else 0
            assert (case['N'], case['CC'], case['MC']) == (m + expected_cc, expected_cc, m)
            assert len(case['ordered_ID_sha256']) == len(case['score_sha256']) == 64
        for suffix in ('mixed', 'negative'):
            maxima = {tuple(x['maxima']) for x in cases if x['case'].endswith('_' + suffix)}
            assert len(maxima) == 1
for experiment, split in (('bench1', 0), ('bench1', 101), ('bench2', 0), ('bench2', 21), ('other', 1)):
    try:
        module.select_config(config, experiment, split)
    except ValueError:
        pass
    else:
        raise AssertionError((experiment, split))
broken = dict(config); broken['cases'] = [dict(x) for x in config['cases']]
broken['cases'][0]['stage'] = 'wrong'
try:
    module.select_config(broken, broken['cases'][0]['experiment'], broken['cases'][0]['split'])
except ValueError:
    pass
else:
    raise AssertionError('wrong stage accepted')
print('PASS: 720 case selectors, 120 split identities, 360 composition seeds, stages and invalid inputs')
