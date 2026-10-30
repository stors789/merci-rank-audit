"""Rebuild the fixed P2 discovery/evaluation-half source qualification."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from scipy import sparse
from half_counts import summarize, evaluation

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'gbml'))
from upstream import check as check_upstream

LENGTH=16569;BASES='ACGT'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('roles','consensus','alignment','reference','upstream','r-library','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--rscript',default='Rscript');a=p.parse_args();started=time.monotonic()
    check_upstream(a.upstream.resolve())
    meta=pd.read_csv(a.roles,sep='\t',dtype={'barcode':str}).set_index('barcode')
    cells=json.loads((a.consensus/'allele_cell_barcodes.json').read_text())
    consensus=json.loads((a.consensus/'consensus_receipt.json').read_text())
    alignment=json.loads((a.alignment/'alignment_receipt.json').read_text())
    reads=Counter(json.loads((a.alignment/'accepted_MT_reads_by_cell.json').read_text()))
    if consensus['status']!='COMPLETE' or consensus['matrix_shape']!=[len(cells),8*LENGTH] or not meta.index.is_unique or len(set(cells))!=len(cells) or not set(cells)<=set(meta.index):
        raise ValueError('New RNA and consensus cell axes differ')
    if sum(reads.values())!=alignment['accepted_MT_reads']:
        raise ValueError('Accepted alignment counts differ')
    meta['full_library_accepted_MT_reads']=[reads[b] for b in meta.index]
    meta['full_library_MT_reads_gt1000']=meta.full_library_accepted_MT_reads.gt(1000)
    reference=''.join(x.strip() for x in a.reference.read_text().splitlines() if not x.startswith('>')).upper()
    if len(reference)!=LENGTH:raise ValueError('Unexpected chrM reference')
    def groups():
        for number,start in enumerate(range(0,len(cells),16)):
            part=sparse.load_npz(a.consensus/'groups'/f'group_{number:04d}.npz').tocsr()
            if part.shape!=(len(cells[start:start+16]),8*LENGTH) or part.dtype.kind!='i':
                raise ValueError('New consensus group axis or dtype differs')
            yield start,part
    agg=np.zeros((LENGTH,8),dtype=np.int64);mutated_cells=np.zeros_like(agg)
    umi_positions=np.zeros(len(cells),dtype=np.int64);sites_depth5=np.zeros(len(cells),dtype=np.int64)
    for start,part in groups():
        total,mutated,umi,sites=summarize(part,LENGTH);end=start+part.shape[0]
        agg+=total;mutated_cells+=mutated;umi_positions[start:end]=umi;sites_depth5[start:end]=sites
    if int(umi_positions.sum())!=consensus['UMI_positions'] or int(agg.sum())!=consensus['UMI_positions']:
        raise ValueError('Consensus UMI-position totals differ')
    meta['sequence_candidate']=meta.index.isin(cells)
    meta['UMI_positions']=pd.Series(umi_positions,index=cells).reindex(meta.index,fill_value=0)
    meta['sites_UMI_depth_ge5']=pd.Series(sites_depth5,index=cells).reindex(meta.index,fill_value=0)
    meta['UMI_sequence_observed']=meta.UMI_positions.gt(0)
    candidates=[]
    for pos in range(LENGTH):
        depth=int(agg[pos,:4].sum())
        if depth<30:continue
        for b,alt in enumerate(BASES):
            count=int(agg[pos,b]);af=count/depth
            if alt==reference[pos] or count<10 or not .01<=af<=.99:continue
            ncell=int(mutated_cells[pos,b])
            if ncell>=3:candidates.append({'ID':f'MT_{pos+1}_{reference[pos]}-{alt}','position_1based':pos+1,
                                           'ref':reference[pos],'alt':alt,'discovery_depth':depth,
                                           'discovery_alt_UMI':count,'discovery_AF':af,'discovery_alt_cells':ncell})
    cand=pd.DataFrame(candidates,columns=['ID','position_1based','ref','alt','discovery_depth','discovery_alt_UMI','discovery_AF','discovery_alt_cells'])
    af=np.full((len(cand),len(meta)),np.nan);depths=np.zeros_like(af,dtype=np.int32);alts=np.zeros_like(depths)
    positions=cand.position_1based.to_numpy(dtype=int)-1;base_indices=np.array([BASES.index(x) for x in cand.alt],dtype=int)
    for start,part in groups():
        indices=meta.index.get_indexer(cells[start:start+part.shape[0]])
        d,alt=evaluation(part,positions,base_indices,LENGTH)
        depths[:,indices]=d;alts[:,indices]=alt
        values=np.full(d.shape,np.nan);np.divide(alt,d,out=values,where=d>=5);af[:,indices]=values
    mask=meta.full_library_MT_reads_gt1000.to_numpy();donors=np.flatnonzero(mask&meta.working_type.eq('myeloid').to_numpy())
    targets=np.flatnonzero(mask&meta.working_type.eq('hepatic').to_numpy());pair=np.r_[donors,targets];keep=[];stages=[]
    for j,v in enumerate(cand.itertuples()):
        da,ta=af[j,donors],af[j,targets];do=int(np.isfinite(da).sum());to=int(np.isfinite(ta).sum());dm=int((da>0).sum());tm=int((ta>0).sum())
        observed=bool(len(donors) and len(targets) and 10*do>=len(donors) and 10*to>=len(targets) and 10*(do+to)>=len(pair))
        eligible=observed and dm+tm>=5
        stages.append({'ID':v.ID,'donor_observed':do,'target_observed':to,'donor_mutated':dm,'target_mutated':tm,
                       'coverage_pass':observed,'mutated_cell_minimum_pass':dm+tm>=5,'test_eligible':eligible})
        if eligible:keep.append(j)
    a.output.mkdir(parents=True,mode=0o700,exist_ok=False)
    pd.DataFrame(af[keep][:,pair],index=cand.iloc[keep].ID,columns=meta.index[pair]).to_csv(a.output/'VAF.tsv',sep='\t',na_rep='NA')
    (a.output/'donors.txt').write_text(''.join(x+'\n' for x in meta.index[donors]))
    (a.output/'targets.txt').write_text(''.join(x+'\n' for x in meta.index[targets]))
    cand.to_csv(a.output/'discovery_candidates.csv',index=False);pd.DataFrame(stages).to_csv(a.output/'qualification_stages.csv',index=False)
    meta.to_csv(a.output/'full_metadata.tsv.gz',sep='\t')
    np.savez_compressed(a.output/'qualification_inputs.npz',AF=af,depth=depths,ALT=alts,
                        IDs=cand.ID.to_numpy(dtype=str),barcodes=meta.index.to_numpy(dtype=str))
    if not keep:raise ValueError('No test-eligible features; first generated inputs retained')
    output=a.output/'feature_statistics.csv';script=Path(__file__).resolve().parent/'source_features.R'
    with (a.output/'R.log').open('x') as log:
        subprocess.run([a.rscript,'--vanilla',str(script),str(a.output.resolve()),str(output.resolve()),
                        str(a.upstream.resolve()),str(a.r_library.resolve())],stdout=log,stderr=subprocess.STDOUT,check=True)
    stats=pd.read_csv(output);features=stats.loc[stats.catogary.eq('donor_enrich'),'ID'].tolist()
    result={'status':'COMPLETE','donor_cells':len(donors),'target_cells':len(targets),'discovered_candidates':len(cand),
            'test_eligible_variants':len(keep),'donor_enriched':len(features),'donor_depleted':int(stats.catogary.eq('donor_depleted').sum()),
            'feature_ids':features,'strict_MT_reads_gt':1000,'no_exploratory_rescue':True,
            'seconds':round(time.monotonic()-started,3)}
    (a.output/'qualification_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    if not features:raise ValueError('No qualified source features; downstream scoring is not run')


if __name__=='__main__':main()
