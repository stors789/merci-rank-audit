"""Create joint cell/UMI consensus from all newly aligned P2 candidate chunks."""
import argparse
from collections import Counter
import json
from pathlib import Path
import resource
import sys
import time

import pysam
from scipy import sparse

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'gbml'))
from pileup import consensus_matrix, LENGTH
from joint_concat import stack_on_disk


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('extraction','alignment','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();started=time.monotonic()
    alignment=json.loads((a.alignment/'alignment_receipt.json').read_text())
    if alignment['status']!='COMPLETE':raise ValueError('New full alignment is incomplete')
    cells=json.loads((a.extraction/'allele_cell_barcodes.json').read_text())
    if len(cells)!=len(set(cells)) or cells!=sorted(cells):raise ValueError('Candidate cell axis differs')
    a.output.mkdir(parents=True,mode=0o700,exist_ok=False);folder=a.output/'groups';folder.mkdir()
    chunks=[];shapes=[];total=Counter();molecules=positions=0
    with pysam.AlignmentFile(a.alignment/'all.accepted.CB.bam','rb') as source:
        iterator=source.fetch(until_eof=True);record=next(iterator,None)
        for number,start in enumerate(range(0,len(cells),16)):
            group=cells[start:start+16];allowed=set(group);sam=folder/f'group_{number:04d}.sam';npz=folder/f'group_{number:04d}.npz'
            with sam.open('x') as out:
                while record is not None and record.get_tag('CB')<=group[-1]:
                    if record.get_tag('CB') not in allowed:raise ValueError('Unexpected sorted cell')
                    out.write(record.to_string()+'\n');record=next(iterator,None)
            if sam.stat().st_size:matrix,counts,n,pos=consensus_matrix(sam,group)
            else:matrix=sparse.csr_matrix((len(group),8*LENGTH),dtype='int32');counts=Counter();n=pos=0
            sparse.save_npz(npz,matrix,compressed=False);chunks.append(npz);shapes.append({'rows':len(group),'nnz':matrix.nnz})
            total.update(counts);molecules+=n;positions+=pos
            (folder/f'group_{number:04d}.json').write_text(json.dumps({'cells':len(group),'counts':counts,'molecules':n,'positions':pos},indent=2)+'\n')
            if number%16==0:print(json.dumps({'completed_groups':number+1,'total_groups':(len(cells)+15)//16}),flush=True)
        if record is not None:raise ValueError('Unconsumed accepted reads')
    if total['accepted_chrM_alignments']!=alignment['accepted_MT_reads']:raise ValueError('Accepted read total differs')
    full=stack_on_disk(chunks,shapes,8*LENGTH,a.output,prefix='consensus')
    sparse.save_npz(a.output/'mt_UMI_alleles.npz',full)
    (a.output/'allele_cell_barcodes.json').write_text(json.dumps(cells)+'\n')
    result={'status':'COMPLETE','counts':total,'matrix_shape':list(full.shape),'UMI_molecules':molecules,
            'UMI_positions':positions,'cell_batch_size':16,'same_cell_UMI_joined_across_all_chunks':True,
            'independent_DNA_not_claimed':True,'seconds':round(time.monotonic()-started,3),
            'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (a.output/'consensus_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
