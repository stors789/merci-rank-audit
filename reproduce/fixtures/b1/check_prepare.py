"""Check B1 cohort preparation and output identity using synthetic scores."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=False)
    fixture=Path(__file__).resolve().parent;repro=fixture.parents[1];driver=repro/'drivers/cohort'
    manifest=output/'manifest.json';prepared=output/'prepared.npz'
    subprocess.run([sys.executable,str(driver/'prepare_split.py'),'--manifest',str(fixture/'manifest.json'),
                    '--experiment','fixture','--split','1','--mixed',str(fixture/'mixed.tsv'),
                    '--negative',str(fixture/'negative.tsv'),'--output',str(prepared),'--output-manifest',str(manifest)],check=True)
    with np.load(prepared,allow_pickle=False) as data:
        assert data['observed_source_features'].dtype==np.dtype('float64')
        assert np.isnan(data['observed_source_features'][0]) and np.isfinite(data['observed_source_features'][1:]).all()
        assert len(data['ids'])==20 and len(data['anchor'])==6 and len(data['background'])==6 and len(data['cc_order'])==8
    subprocess.run([sys.executable,str(driver/'run.py'),'--manifest',str(manifest),'--output',str(output/'cohort')],check=True)
    receipt=json.loads((output/'cohort/receipt.json').read_text())
    assert (receipt['tasks'],receipt['splits'],receipt['conditions'])==(2,1,10)
    table=pd.read_csv(output/'cohort/policy_metrics.csv')
    raw=table[(table['case']=='strict_regular_mixed')&(table['policy']=='original_raw')].iloc[0]
    assert tuple(int(raw[x]) for x in ('selected','TP','FP','TN','FN'))==(11,8,3,9,0)
    partitions=pd.read_csv(output/'cohort/partition_truth.csv')
    mixed=partitions[(partitions['case']=='strict_regular_mixed')&(partitions['family']=='raw')]
    assert int(mixed['observation_unknown'].sum())==1
    transitions=pd.read_csv(output/'cohort/transitions.csv')
    keys=['design','from_level','to_level','group','family']
    assert (transitions.groupby(keys)['count'].sum()==transitions.groupby(keys)['group_N'].first()).all()
    print('PASS: synthetic score binding, float64/NA preservation, two tasks and ten conditions')

if __name__=='__main__':main()
