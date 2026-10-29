"""Execute four small input-to-output checks and redraw every aggregated figure."""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pypdf import PdfReader
from scipy.stats import hypergeom
from merci_rank_audit.inputs import validate_manifest


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    a.output.mkdir(parents=True); out=a.output.resolve(); root=Path(__file__).resolve().parent
    env=os.environ.copy()
    for var in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):env[var]='1'
    env['MPLCONFIGDIR']=str(out/'matplotlib-config');env['PYTHONDONTWRITEBYTECODE']='1'
    def run(path,*args):
        command=[sys.executable,str(root/path),*[str(x) for x in args]]
        subprocess.run(command,check=True,env=env,cwd=out)
    fig=out/'figures'
    for script,data in (('main_figures.py','figures'),('supplement_figures.py','figures'),('bio06_figure.py','bio06')):
        run('figures/'+script,'--data',root/'data'/data,'--output',fig)
    name='formal_n100_tie50_single50_rare_target'
    run('drivers/benchmark/run.py','--config',root/'configs/benchmark_smoke.json','--case',name,'--output',out/'benchmark')
    run('drivers/cohort/run.py','--manifest',root/'configs/cohort_smoke.json','--output',out/'cohort')
    run('fixtures/b1/check_prepare.py','--output',out/'b1-preparation-fixture')
    run('fixtures/benchmark_inputs/check_config.py')
    run('fixtures/gbml/check_thresholds.py')
    formal=json.loads((root/'configs/cohort_formal_portable.json').read_text())
    for split in range(21,101):
        cases=[x for x in formal['cases'] if x['experiment']=='bench1' and x['split']==split
               and x['case'] in ('strict_regular_mixed','strict_regular_negative')]
        selections=[x for x in formal['splits'] if x['experiment']=='bench1' and x['split']==split]
        assert len(cases)==2 and len(selections)==1 and sorted(x['N'] for x in cases)==[3786,6690]
        assert (selections[0]['M'],selections[0]['CC_available'],selections[0]['anchor_N'],
                selections[0]['background_slots'])==(3786,2904,1893,1893)
    run('drivers/nuclear/run.py','--manifest',root/'configs/nuclear_smoke_manifest.json','--config',root/'configs/nuclear_smoke.json',
        '--stream','reference','--output',out/'nuclear')
    run('drivers/sensitivity/run.py','--config',root/'configs/sensitivity_smoke.json','--output',out/'sensitivity')
    log=out/'benchmark'/'results'/name/'rep_01.jsonl'; rows=[json.loads(x) for x in log.read_text().splitlines()]
    def checkpoint(method):
        return next(x for x in rows if x.get('method')==method and x.get('target')==100 and x.get('proposals')==100)
    uniform=checkpoint('uniform');weighted=checkpoint('conditional_weighted');accepted=checkpoint('conditional_accepted')
    assert uniform['accepted']==0 and uniform['member_RMSE'] is None
    for v in (weighted,accepted):
        assert v['accepted']==100 and abs(v['member_RMSE']-0.0300521323095101)<1e-12
    cohort=out/'cohort';receipt=json.loads((cohort/'receipt.json').read_text())
    assert (receipt['tasks'],receipt['splits'],receipt['conditions'])==(1,1,10)
    policies=pd.read_csv(cohort/'policy_metrics.csv').set_index('policy')
    assert tuple(int(policies.loc['original_raw',k]) for k in ('selected','TP','FP','TN','FN'))==(11,8,3,9,0)
    diff=pd.read_csv(cohort/'family_differences.csv').iloc[0];assert bool(diff.core_same) and bool(diff.union_same)
    trans=pd.read_csv(cohort/'transitions.csv')
    assert (trans.groupby(['design','from_level','to_level','group','family'])['count'].sum()
            == trans.groupby(['design','from_level','to_level','group','family'])['group_N'].first()).all()
    _,cases=validate_manifest(root/'configs/nuclear_smoke_manifest.json');case=cases[0];s=case['sampler'];d=case['nuclear']
    with np.load(out/'nuclear'/'samples.npz',allow_pickle=False) as x:
        assert x['betas'].shape==(100,4)
        rank=x['first_ranks'][0];z=((rank>s.n-s.ks[5])&(s.rna>s.n-s.ks[5]))[d['index']][d['observable']].astype(float)
        y=np.column_stack([d['original'],d['alternative']]);ols=np.linalg.lstsq(np.column_stack([d['X'],z]),y,rcond=None)[0][-1]
        assert int(z.sum())==8 and np.allclose(x['betas'][0],ols,atol=1e-12,rtol=0)
        assert np.isfinite(x['betas']).sum(axis=0).tolist()==[96]*4
    nuclear=pd.read_csv(out/'nuclear'/'quantiles.csv');assert len(nuclear)==12 and nuclear.estimate.isna().all() and (nuclear.status=='LOW_VALID_N').all()
    r=pd.read_csv(out/'sensitivity'/'rarity.csv').set_index('null_id')
    assert abs(r.loc['original','P_G']-hypergeom.sf(8,20,11,11))<1e-14
    assert abs(r.loc['null_01','P_G']-hypergeom.sf(9,20,11,11))<1e-14
    assert len(pd.read_csv(out/'sensitivity'/'scenarios.csv'))==4
    transitions=pd.read_csv(out/'sensitivity'/'transitions.csv');assert len(transitions)==12 and (transitions.changed==0).all()
    expected=['Figure1_workflow','Figure2_sampling','Figure3_annotation','Figure4_composition','Figure5_nuclear',
              'FigureS1_gate','FigureS2_cohorts','FigureS3_precision','FigureS4_union','FigureS5_bio06']
    for name in expected:
        assert (fig/(name+'.png')).stat().st_size>1000
        assert len(PdfReader(fig/(name+'.pdf')).pages)==1
    result={'status':'PASS','benchmark':'known original-seed 100-proposal prefix',
            'cohort':'one synthetic task, ten design conditions and exhaustive transition totals',
            'b1_preparation':'two synthetic score tasks, frozen selection seeds and float64/NA prepared representation',
            'b1_frozen_split_configs_checked':80,
            'benchmark_input_configs_checked':720,
            'gbml_frozen_threshold_integrity':'original table accepted; altered maximum rejected',
            'nuclear':'100 proposals, first beta matches independent OLS; low-valid quantiles withheld',
            'sensitivity':'four synthetic scenarios; two single-gate probabilities equal hypergeometric truth',
            'redrawn_figures':expected,'scope':'engineering checks; no full raw-data rerun'}
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
