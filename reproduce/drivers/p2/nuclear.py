"""Recreate the fixed P2 52×2 nuclear matrices and original-label models."""
import argparse
import json
from pathlib import Path
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

from regression import fit,bh

CONTACT=['RAC1','CDC42','RHOB','MYH9','VIM','FLNA']
HUMAN='MKI67 PCNA MCM2 MCM5 TYMS TOP2A CCNB1 CDK1 UBE2C STMN1'.split()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('counts','metadata','cells','gmt','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();started=time.monotonic()
    meta=pd.read_csv(a.metadata,sep='\t',dtype={'barcode':str}).set_index('barcode')
    target=meta[meta.working_type.eq('hepatic')].copy()
    score=pd.read_csv(a.cells,sep='\t',dtype={'barcode':str},float_precision='round_trip').set_index('barcode')
    if len(target)!=3428 or len(score)!=3265 or not score.index.is_unique or not score.index.is_monotonic_increasing:
        raise ValueError('P2 target or scored axis differs')
    ix=pd.Index(target.index).get_indexer(score.index);observed=score.source_marker_observable.to_numpy(bool)
    if (ix<0).any() or observed.sum()!=2741:raise ValueError('P2 source-observable axis differs')
    sets={line.split('\t')[0]:line.rstrip().split('\t')[2:] for line in a.gmt.read_text().splitlines()}
    if len(sets)!=50:raise ValueError('Expected 50 fixed Hallmark sets')
    data=ad.read_h5ad(a.counts,backed='r');symbols=data.var.gene_symbol.to_numpy(dtype=str)
    nuclear=~np.char.startswith(symbols,'MT-');ids=target.cell_id.tolist()
    x=data[ids,:].X.tocsr();nuc=target.nuclear_UMI.to_numpy()
    if not np.array_equal(x[:,nuclear].sum(1).A1,nuc) or not np.array_equal(x.sum(1).A1,target.UMI.to_numpy()):
        raise ValueError('New GEO count denominators differ')
    hallmarks=[]
    for members in sets.values():
        mask=nuclear & np.isin(symbols,members)
        hallmarks.append(np.log1p(10000*x[:,mask].sum(1).A1/nuc))
    if not all((symbols==g).any() for g in HUMAN):raise ValueError('Human10 gene missing')
    h=np.column_stack([np.log1p(10000*x[:,symbols==g].sum(1).A1/nuc) for g in HUMAN])
    proliferation=np.column_stack([(rankdata(h[:,j],method='average')-1)/(len(h)-1) for j in range(len(HUMAN))]).mean(1)
    names=list(sets)+['P2_EXISTING_CONTACT_TRANSPORT','HUMAN10_PROLIFERATION_RANK']
    original=np.column_stack(hallmarks+[target.contact_transport.to_numpy(),proliferation])
    del x
    sets['P2_EXISTING_CONTACT_TRANSPORT']=CONTACT;sets['HUMAN10_PROLIFERATION_RANK']=HUMAN
    wanted={g for members in sets.values() for g in members if not g.startswith('MT-')}
    columns=np.flatnonzero(np.isin(symbols,list(wanted)) & nuclear);unique=np.unique(symbols[columns])
    mapping=sparse.csr_matrix((np.ones(len(columns)),(np.arange(len(columns)),np.searchsorted(unique,symbols[columns]))),
                              shape=(len(columns),len(unique)))
    gene_counts=np.empty((len(target),len(unique)),float)
    for start in range(0,len(target),128):
        chunk=data[ids[start:start+128],:].X.tocsr()[:,columns]
        gene_counts[start:start+len(chunk.indptr)-1]=(chunk@mapping).toarray()
    data.file.close();log_gene=np.log1p(10000*gene_counts/nuc[:,None])
    percentile=(rankdata(log_gene,method='average',axis=0)-1)/(len(target)-1)
    alternative=[]
    for name in names:
        member_index=np.flatnonzero(np.isin(unique,sets[name]))
        if not len(member_index):raise ValueError('A fixed nuclear feature has no present members: '+name)
        if name=='HUMAN10_PROLIFERATION_RANK':
            if len(member_index)!=10:raise ValueError('Human10 gene count differs')
            alt=log_gene[:,member_index].mean(axis=1)
        else:alt=percentile[:,member_index].mean(axis=1)
        alternative.append(alt)
    alternative=np.column_stack(alternative)
    identity=pd.get_dummies(meta.loc[score.index,'identity_cluster'].astype(str),drop_first=True,dtype=float)
    full=np.column_stack([np.log1p(meta.loc[score.index,'full_library_accepted_MT_reads']),
                          np.log1p(meta.loc[score.index,'nuclear_UMI']),identity.to_numpy(),
                          score.observed_source_feature_sites])
    variable=full[observed].std(0)>0;covariates=np.column_stack([np.ones(int(observed.sum())),full[observed][:,variable]])
    raw=score.raw_prediction.eq('Receiver').to_numpy()
    # The original column selection yielded Fortran-order matrices. Preserve its
    # reduction order for standard deviations and standardized coefficients.
    original_scored=np.asfortranarray(original[ix][observed])
    alternative_scored=np.asfortranarray(alternative[ix][observed])
    a.output.mkdir(parents=True,mode=0o700,exist_ok=False)
    np.savez_compressed(a.output/'nuclear.npz',ids=score.index.to_numpy(dtype=str),observable=observed,X=covariates,
                        original=original_scored,alternative=alternative_scored,names=np.array(names,dtype=str))
    for representation,y in (('original',original_scored),('alternative',alternative_scored)):
        table=bh(fit(y,covariates,raw[observed],names),52)
        table['representation']=representation;table['analysis']='raw_sensitivity_analysis';table['selection_state']='original_fixed_MT_RNA_label'
        table['full_scored_domain']=len(score);table['full_raw_candidates']=int(raw.sum());table['full_missing_candidates']=int((raw&~observed).sum())
        table.to_csv(a.output/('P2_primary_'+representation+'.csv'),index=False)
    result={'status':'COMPLETE','all_hepatic_target_domain':len(target),'scored_domain':len(score),'observable_domain':int(observed.sum()),
            'nuclear_features':len(names),'representations':2,'raw_selected':int(raw.sum()),'raw_observable':int(raw[observed].sum()),
            'raw_missing':int((raw&~observed).sum()),'covariate_columns':covariates.shape[1],
            'percentile_rank_domain':len(target),'duplicate_old_OXPHOS_excluded':True,'seconds':round(time.monotonic()-started,3)}
    (a.output/'nuclear_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
