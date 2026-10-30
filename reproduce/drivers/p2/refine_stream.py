"""Reconstruct the original C5 refinement using the new joint counts and roles."""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc

from prepare_counts import CONTACT


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--counts',type=Path,required=True);p.add_argument('--global-roles',type=Path,required=True)
    p.add_argument('--gmt',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,mode=0o700,exist_ok=False)
    obs=pd.read_csv(args.global_roles,sep='\t',dtype={'cell_id':str,'identity_cluster':str,'GSM':str,'barcode':str},float_precision='round_trip').set_index('cell_id')
    data=ad.read_h5ad(args.counts,backed='r');genes=data.var.gene_symbol.to_numpy().astype(str)
    if obs.index.tolist()!=data.obs_names.tolist():raise ValueError('New joint count and cluster axes differ')
    original=obs.identity_cluster.astype(str);members=original.eq('5').to_numpy()
    if members.sum()!=1216:raise ValueError('Original C5 refinement domain differs')
    X=data.X[np.flatnonzero(members),:].tocsr();data.file.close()
    n=obs.loc[members,'nuclear_UMI'].to_numpy();logn=(sparse.diags(1e4/n) @ X.astype(np.float32)).tocsr()
    logn.data=np.log1p(logn.data)
    hallmark={line.split('\t')[0]:line.rstrip().split('\t')[2:] for line in args.gmt.read_text().splitlines()}
    mt=np.char.startswith(genes,'MT-');eligible=~mt & ~np.isin(genes,hallmark['HALLMARK_OXIDATIVE_PHOSPHORYLATION']+CONTACT+['ACTA2'])
    var=pd.DataFrame({'gene_symbol':genes},index=data.var_names)
    b=ad.AnnData(X=logn[:,eligible].copy(),obs=obs.loc[members].copy(),var=var.loc[eligible].copy())
    sc.pp.highly_variable_genes(b,n_top_genes=2000,flavor='seurat');b=b[:,b.var.highly_variable].copy()
    sc.pp.scale(b,max_value=10);sc.tl.pca(b,n_comps=30,svd_solver='arpack',random_state=0)
    sc.pp.neighbors(b,n_neighbors=15,n_pcs=30,random_state=0)
    sc.tl.leiden(b,resolution=.5,random_state=0,flavor='igraph',n_iterations=2)
    sub=b.obs['leiden'].astype(str).to_numpy();clusters=sorted(set(sub),key=int)
    labels={c:str(max(map(int,original.unique()))+1+i) for i,c in enumerate(clusters)}
    decisions=pd.read_csv(Path(__file__).resolve().parent/'source_identity_decisions.csv',dtype={'identity_cluster':str}).set_index('identity_cluster')
    if set(labels.values())!=set(decisions.index):raise ValueError('Refined cluster IDs differ from frozen decisions')
    obs['global_identity_cluster']=original;obs['identity_cluster']=original.astype(object);obs['working_type']=obs.working_type.astype(str).astype(object)
    pos=np.where(members)[0]
    for c in clusters:
        cells=pos[sub==c];label=labels[c]
        obs.iloc[cells,obs.columns.get_loc('identity_cluster')]=label
        obs.iloc[cells,obs.columns.get_loc('working_type')]=decisions.loc[label,'working_type']
    obs['compatible_identity_set']=obs.identity_cluster.map(decisions['compatible_identity_set'])
    obs.to_csv(args.output/'refined_roles.tsv.gz',sep='\t')
    result={'status':'COMPLETE','global_C5_cells':int(members.sum()),'subclusters':len(clusters),
            'cluster_counts':{labels[c]:int((sub==c).sum()) for c in clusters},'fixed_source_decisions_applied':True}
    (args.output/'refinement_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
