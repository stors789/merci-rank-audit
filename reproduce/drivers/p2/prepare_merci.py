"""Build P2 MERCI inputs from newly regenerated source and GEO count matrices."""
import argparse
import csv
import json
from pathlib import Path
import time

import anndata as ad
import numpy as np
import pandas as pd

GENES=sorted(['MT-ATP6','MT-ATP8','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ND1',
              'MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('qualification','counts','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();started=time.monotonic();q=a.qualification
    receipt=json.loads((q/'qualification_receipt.json').read_text())
    if receipt['status']!='COMPLETE' or not receipt['donor_enriched']:
        raise ValueError('The new strict source qualification is not eligible for scoring')
    meta=pd.read_csv(q/'full_metadata.tsv.gz',sep='\t',dtype={'barcode':str}).set_index('barcode')
    donors=sorted((q/'donors.txt').read_text().splitlines());targets=sorted((q/'targets.txt').read_text().splitlines());cells=donors+targets
    if set(donors)&set(targets) or not meta.loc[donors,'working_type'].eq('myeloid').all() or not meta.loc[targets,'working_type'].eq('hepatic').all():
        raise ValueError('Source or target role mismatch')
    if not meta.loc[cells,'full_library_accepted_MT_reads'].gt(1000).all():raise ValueError('Strict read gate mismatch')
    stats=pd.read_csv(q/'feature_statistics.csv');features=stats.loc[stats.catogary.eq('donor_enrich'),'ID'].tolist()
    vaf=pd.read_csv(q/'VAF.tsv',sep='\t',index_col=0)
    if not set(features)<=set(vaf.index) or set(cells)!=set(vaf.columns):raise ValueError('New VAF axis differs')
    a.output.mkdir(parents=True,mode=0o700,exist_ok=False)
    vaf.loc[:,cells].to_csv(a.output/'VAF.tsv',sep='\t',na_rep='NA')
    for name,values in (('donors.txt',donors),('targets.txt',targets),('features.txt',features)):
        (a.output/name).write_text(''.join(x+'\n' for x in values))
    with np.load(q/'qualification_inputs.npz',allow_pickle=False) as data:
        barcodes=list(data['barcodes']);ids=list(data['IDs']);positions=[int(x.split('_')[1]) for x in features]
        ci=np.array([barcodes.index(c) for c in cells]);ri=[ids.index(x) for x in features];coverage={}
        with (a.output/'MT_variants.tsv').open('x',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=['Cell','ID','pos','s_reads','t_reads','AF','Cell_reads'],delimiter='\t');writer.writeheader()
            for feature,index,pos in zip(features,ri,positions):
                depth=data['depth'][index,ci];alt=data['ALT'][index,ci]
                if pos in coverage:
                    if not np.array_equal(coverage[pos],depth):raise ValueError('Feature position coverage differs')
                else:coverage[pos]=depth
                for j in np.flatnonzero(alt>0):
                    if depth[j]<=0:raise ValueError('ALT with zero depth')
                    writer.writerow({'Cell':cells[j],'ID':feature,'pos':pos,'s_reads':int(alt[j]),'t_reads':int(depth[j]),
                                     'AF':float(alt[j]/depth[j]),'Cell_reads':int(meta.loc[cells[j],'full_library_accepted_MT_reads'])})
    pd.DataFrame.from_dict(coverage,orient='index',columns=cells).sort_index().to_csv(a.output/'Coverage.tsv',sep='\t',index_label='native_position1')
    counts=ad.read_h5ad(a.counts,backed='r');symbols=counts.var.gene_symbol.astype(str).to_numpy();wanted=np.flatnonzero(np.isin(symbols,GENES))
    if not set(GENES)<=set(symbols[wanted]):raise ValueError('Original GEO MT13 genes missing')
    x=counts[meta.loc[cells,'cell_id'].tolist(),wanted].X.tocsr();counts.file.close()
    values=np.vstack([x[:,symbols[wanted]==g].sum(1).A1 for g in GENES]);total=meta.loc[cells,'UMI'].to_numpy()
    expression=np.log1p(10000*values/total[None,:])
    if not np.isfinite(expression).all() or (expression<0).any():raise ValueError('Invalid MT expression')
    pd.DataFrame(expression,index=GENES,columns=cells).to_csv(a.output/'MT_expression.tsv',sep='\t',index_label='gene')
    meta.loc[cells].to_csv(a.output/'cell_roles.tsv',sep='\t')
    result={'status':'ready','donor_cells':len(donors),'candidate_target_cells':len(targets),'source_features':features,
            'coverage_positions_native1based':len(coverage),'RNA_genes':GENES,
            'RNA_normalization':'log1p(10000*original_GEO_gene_UMI/original_GEO_total_UMI)',
            'RNA_zero_coding_MT_targets':int((values[:,len(donors):].sum(0)==0).sum()),
            'parameters':{'Neff':'LR e=.001 frac_tran=FALSE','RNA':'Human LOO epsilon=.1','cutoff':'fixed 50%',
                          'maxima':'published original matched-N table; no new null'},'seconds':round(time.monotonic()-started,3)}
    (a.output/'ready.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
