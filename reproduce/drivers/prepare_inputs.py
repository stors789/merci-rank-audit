"""Convert explicit saved cell scores or the original own P2 prepared NPZ to audit inputs."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='kind',required=True)
    score=sub.add_parser('score',help='Convert a saved per-cell table and ten selected cutoff rows')
    score.add_argument('--cells',type=Path,required=True);score.add_argument('--cutoffs',type=Path,required=True)
    score.add_argument('--dna-column',default='DNA_rank');score.add_argument('--maximum-column',default='random_max')
    score.add_argument('--include-rna-score',action='store_true');score.add_argument('--output',type=Path,required=True)
    nuclear=sub.add_parser('nuclear',help='Convert P2 prepared NPZ using its original field names')
    nuclear.add_argument('--prepared',type=Path,required=True);nuclear.add_argument('--score',type=Path,required=True)
    nuclear.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.kind=='score':
        if a.output.exists():raise FileExistsError(a.output)
        cells=pd.read_csv(a.cells,sep='\t',compression='infer',dtype={'cell':str,'Neff':str,'RNA_rank':str,a.dna_column:str},keep_default_na=False)
        cuts=pd.read_csv(a.cutoffs,compression='infer',dtype={'N':str,'cutoff':str,a.maximum_column:str},keep_default_na=False)
        fields=['cell','Neff','RNA_rank',a.dna_column]+(['RNA_score'] if a.include_rna_score else [])
        if not set(fields).issubset(cells.columns):raise ValueError('missing cell fields')
        if not {'N','cutoff',a.maximum_column}.issubset(cuts.columns):raise ValueError('missing cutoff fields')
        if len(cuts)!=10 or set(cuts.cutoff.astype(int))!=set(range(100,0,-10)) or not cuts.N.astype(int).eq(len(cells)).all():
            raise ValueError('provide exactly ten matched-N rows for one case')
        a.output.mkdir(parents=True);cells[fields].rename(columns={a.dna_column:'DNA_rank'}).to_csv(a.output/'score.tsv',sep='\t',index=False)
        cuts[['N','cutoff',a.maximum_column]].rename(columns={a.maximum_column:'maximum'}).to_csv(a.output/'thresholds.csv',index=False)
        print(f'Generated score.tsv and thresholds.csv for N={len(cells)}; verify historical byte digests separately.')
    else:
        if a.output.exists():raise FileExistsError(a.output)
        with np.load(a.prepared,allow_pickle=False) as src:
            mapping={'ids':'barcodes','observable':'observable','X':'covariates','original':'original','alternative':'alternative','names':'names'}
            if not set(mapping.values()).issubset(src.files):raise ValueError('P2 prepared NPZ missing required arrays')
            data={key:src[value] for key,value in mapping.items()}
        scores=pd.read_csv(a.score,sep='\t',dtype={'cell':str},keep_default_na=False)
        if set(scores.cell)!=set(data['ids'].astype(str)) or len(scores)!=len(data['ids']):raise ValueError('P2 score and nuclear IDs differ')
        nobs=int(data['observable'].sum())
        if any(data[key].shape[0]!=nobs for key in ('X','original','alternative')):raise ValueError('observable nuclear row axis differs')
        if data['original'].shape[1]!=len(data['names']) or data['alternative'].shape[1]!=len(data['names']):raise ValueError('nuclear feature axis differs')
        a.output.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(a.output,**data)
        print(f'Generated nuclear NPZ for N={len(scores)}; observable rows={nobs}; features={len(data["names"])}.')

if __name__=='__main__':main()
