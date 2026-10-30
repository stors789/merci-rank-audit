"""Check the P2 scorer rejects fixed thresholds for a different synthetic N."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('upstream','r-library','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--rscript',default='Rscript');p.add_argument('--full-size',action='store_true')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    repro=Path(__file__).resolve().parents[2];folder=a.output/'inputs';folder.mkdir()
    donors=[f'd{i:02d}' for i in range(20)];targets=[f'r{i:04d}' for i in range(3265 if a.full_size else 40)];cells=donors+targets
    features=['MT_1_A-C','MT_2_C-G','MT_3_G-T']
    m=np.zeros((3,len(cells)));m[0,:20]=.5;m[0,20:25]=.2;m[1,:10]=.3;m[1,20:30]=.1;m[2,10:20]=.2;m[2,30:40]=.4
    pd.DataFrame(m,index=features,columns=cells).to_csv(folder/'VAF.tsv',sep='\t')
    coverage=pd.DataFrame(np.full((3,len(cells)),10),index=[1,2,3],columns=cells)
    coverage.to_csv(folder/'Coverage.tsv',sep='\t')
    rows=[]
    for i,j in np.argwhere(m>0):rows.append({'Cell':cells[j],'ID':features[i],'pos':i+1,'s_reads':int(round(10*m[i,j])),
                                             't_reads':10,'AF':m[i,j],'Cell_reads':20})
    pd.DataFrame(rows).to_csv(folder/'MT_variants.tsv',sep='\t',index=False)
    genes=sorted(['MT-ATP6','MT-ATP8','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6'])
    values=np.array([[1+(i*7+j*3)%19 for j in range(len(cells))] for i in range(len(genes))],dtype=float)
    pd.DataFrame(values,index=genes,columns=cells).to_csv(folder/'MT_expression.tsv',sep='\t')
    for name,items in (('donors.txt',donors),('targets.txt',targets),('features.txt',features)):
        (folder/name).write_text(''.join(x+'\n' for x in items))
    (folder/'ready.json').write_text(json.dumps({'status':'ready','parameters':{}})+'\n')
    out=a.output/'score';script=repro/'drivers/p2/score_merci.R';table=repro/'data/bio06/all_thresholds.csv'
    result=subprocess.run([a.rscript,'--vanilla',str(script),str(folder.resolve()),str(out.resolve()),
                           str(a.upstream.resolve()),str(a.r_library.resolve()),str(table.resolve())],capture_output=True,text=True)
    (a.output/'score.log').write_text(result.stdout+result.stderr)
    receipt=json.loads((out/'summary.json').read_text())
    if a.full_size:
        if result.returncode or receipt['result']['status']=='run_failed':
            raise AssertionError('Full-size artificial scorer did not complete; inspect score.log')
        score=pd.read_csv(out/'cells.tsv',sep='\t');curve=pd.read_csv(out/'gate_curve.tsv',sep='\t')
        if len(score)!=3265 or score.barcode.tolist()!=targets or len(curve)!=10 or not curve.N.eq(3265).all():
            raise AssertionError('Full-size artificial score or threshold axes differ')
        selected=score.MTvar_rank.gt(3265-1633)&score.W_Donor_rank.gt(3265-1633)
        if not np.array_equal(selected,score.raw_prediction.eq('Receiver')):
            raise AssertionError('Original 50 percent selection and rank intersection differ')
        print('PASS: the P2 scorer writes the complete artificial matched-N score and ten fixed gates')
    else:
        if result.returncode==0 or receipt['result']['status']!='run_failed' or 'supplied$N' not in receipt['result']['error'] or (out/'cells.tsv').exists():
            raise AssertionError('Different-N fixed thresholds were not rejected after scoring')
        print('PASS: P2 scorer reaches and rejects the original fixed thresholds for a different synthetic N')


if __name__=='__main__':main()
