"""Bind a newly regenerated P2 score, fixed maxima and nuclear input."""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

TABLE_SHA='7d75d473b3861ec69ed5b3368305730fdcb8ee190d517329014b2381dbdeb0d6'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('cells','threshold-table','nuclear','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if sha(a.threshold_table)!=TABLE_SHA:raise ValueError('Published fixed-threshold table bytes differ')
    data=pd.read_csv(a.cells,sep='\t',dtype={'barcode':str,'DMTvar_count':str,'MTvar_rank':str,'W_Donor_rank':str},keep_default_na=False)
    cuts=pd.read_csv(a.threshold_table);cuts=cuts[(cuts.case=='P2')&(cuts.null_id=='original')].copy()
    cuts=cuts.set_index('cutoff').loc[list(range(100,0,-10))].reset_index()
    if len(data)!=3265 or data.barcode.duplicated().any() or len(cuts)!=10 or not cuts.N.eq(len(data)).all():
        raise ValueError('P2 scored axis or matched-N maxima differ')
    a.output.mkdir(parents=True,mode=0o700,exist_ok=False)
    score=a.output/'score.tsv';thresholds=a.output/'thresholds.csv'
    data[['barcode','DMTvar_count','W_Donor_rank','MTvar_rank']].rename(columns={
        'barcode':'cell','DMTvar_count':'Neff','W_Donor_rank':'RNA_rank','MTvar_rank':'DNA_rank'}).to_csv(score,sep='\t',index=False)
    cuts[['N','cutoff','maximum']].to_csv(thresholds,index=False)
    manifest={'budget':100,'cases':[{'name':'P2','score':str(score.resolve()),'thresholds':str(thresholds.resolve()),
                'nuclear':str(a.nuclear.resolve()),'seed_uniform':0,'seed_conditional':0}]}
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    result={'status':'COMPLETE','N':len(data),'score_sha256':sha(score),'threshold_sha256':sha(thresholds),
            'nuclear_sha256':sha(a.nuclear),'no_old_sample_reuse':True}
    (a.output/'binding_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
