"""Align newly generated P2 candidate chunks against the full hg38 index."""
import argparse
from collections import Counter
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path
import re
import resource
import subprocess
import time

import pandas as pd
import pysam


def accept(line,counts):
    if line.startswith('@'):return None
    f=line.rstrip('\r\n').split('\t');flag=int(f[1])
    if flag & (4|256|2048):counts['unmapped_or_secondary_supplementary']+=1;return None
    counts['primary_'+f[2]]+=1
    if f[2]!='chrM' or int(f[4])<30:return None
    cigar=[(int(n),op) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',f[5])]
    if any(op=='N' for _,op in cigar):counts['chrM_spliced_alignment_skipped']+=1;return None
    if sum(n for n,op in cigar if op in 'M=X')<60:return None
    counts['accepted_chrM_alignments']+=1
    counts['accepted_strand_reverse' if flag&16 else 'accepted_strand_forward']+=1
    return f[:11]


def child_limit():
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('extraction','roles','index','minimap','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--resume',action='store_true');a=p.parse_args();a.output=a.output.resolve();start=time.monotonic()
    receipt=json.loads((a.extraction/'extraction_receipt.json').read_text())
    if receipt['status']!='COMPLETE':raise ValueError('New full extraction is incomplete')
    chunks=json.loads((a.extraction/'candidate_chunks.json').read_text())
    if sum(x['candidate_records'] for x in chunks)!=receipt['counts']['mitochondrial_kmer_candidate']:
        raise ValueError('Candidate chunk totals differ')
    if a.index.stat().st_size!=15420201719:raise ValueError('Expected the original hg38 split index')
    version=subprocess.check_output([str(a.minimap),'--version'],text=True).strip()
    if version!='2.30-r1287':raise ValueError('Expected minimap2 2.30-r1287')
    roles=pd.read_csv(a.roles,sep='\t',dtype={'barcode':str});allowed=set(roles.barcode)
    def identity(path):
        stat=path.stat();return {'path':str(path.resolve()),'bytes':stat.st_size,'mtime_ns':stat.st_mtime_ns}
    inputs={'driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'extraction':str(a.extraction.resolve()),
            'chunks_sha256':hashlib.sha256((a.extraction/'candidate_chunks.json').read_bytes()).hexdigest(),
            'roles_sha256':hashlib.sha256(a.roles.read_bytes()).hexdigest(),
            'candidate_files':[identity(Path(chunk['path'])) for chunk in chunks],
            'index':identity(a.index),'minimap':identity(a.minimap),'minimap_version':version}
    folder=a.output/'aligned';input_receipt=a.output/'alignment_inputs.json'
    if a.resume:
        if (a.output/'alignment_receipt.json').exists():raise FileExistsError('Complete alignment is already retained')
        if json.loads(input_receipt.read_text())!=inputs:raise ValueError('Alignment inputs changed since the first attempt')
        if not folder.is_dir():raise FileNotFoundError(folder)
    else:
        a.output.mkdir(parents=True,mode=0o700,exist_ok=False);folder.mkdir()
        input_receipt.write_text(json.dumps(inputs,indent=2)+'\n')
    def preserve(path):
        if path.exists():path.rename(path.with_name(path.name+f'.interrupted_{time.time_ns()}'))
    header=pysam.AlignmentHeader.from_dict({'HD':{'VN':'1.6','SO':'unsorted'},'SQ':[{'SN':'chrM','LN':16569}],
           'CO':['Accepted chrM records only; full nuclear competition counts are retained separately.']})
    completed=[]
    for number,chunk in enumerate(chunks):
        bam=folder/f'chunk_{number:04d}.bam';log=folder/f'chunk_{number:04d}.log';done=folder/f'chunk_{number:04d}.json'
        if a.resume and done.exists():
            item=json.loads(done.read_text())
            if (item['chunk']!=number or item['candidate_records']!=chunk['candidate_records'] or
                item['BAM']!=str(bam) or not bam.is_file() or bam.stat().st_size!=item['BAM_bytes'] or
                sum(item['accepted_reads_by_cell'].values())!=item['counts'].get('accepted_chrM_alignments',0)):
                raise ValueError('Completed alignment chunk does not match its receipt')
            completed.append(item)
            print(json.dumps({'reused_completed_chunk':number+1,'chunks':len(chunks)}),flush=True)
            continue
        preserve(bam);preserve(log)
        counts=Counter();bycell=Counter();seen=bytearray((chunk['raw_pair_end']-chunk['raw_pair_start'])//8+1);primary=0
        attempt=1+len(list(folder.glob(f'chunk_{number:04d}.log.interrupted_*')))
        command=[str(a.minimap.resolve()),'-ax','splice:sr','-t','4','-K','20M','--split-prefix',
                 str(folder/f'chunk_{number:04d}.attempt_{attempt}.split'),'--sam-hit-only','-N','5',str(a.index.resolve()),chunk['path']]
        with ExitStack() as stack:
            output=stack.enter_context(pysam.AlignmentFile(bam,'wb',header=header,threads=2))
            error=stack.enter_context(log.open('x'))
            process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=error,text=True,bufsize=1024**2,preexec_fn=child_limit)
            try:
                for line in process.stdout:
                    if line.startswith('@'):continue
                    raw=line.rstrip('\n').split('\t');flag=int(raw[1]);cell,umi,read_index=raw[0].split('_');read_index=int(read_index)
                    if cell not in allowed or not chunk['raw_pair_start']<=read_index<=chunk['raw_pair_end']:
                        raise ValueError('Invalid candidate identity or read index')
                    if not flag&(4|256|2048):
                        slot,bit=divmod(read_index-chunk['raw_pair_start'],8)
                        if seen[slot]&(1<<bit):raise ValueError('Duplicate primary query')
                        seen[slot]|=1<<bit;primary+=1
                    fields=accept(line,counts)
                    if fields is not None:
                        rec=pysam.AlignedSegment.fromstring('\t'.join(fields),output.header)
                        rec.set_tag('CB',cell);rec.set_tag('UB',umi);output.write(rec);bycell[cell]+=1
                process.stdout.close();code=process.wait()
            finally:
                if process.poll() is None:process.terminate();process.wait()
        if code:raise RuntimeError(f'Alignment chunk {number} failed; partial BAM and log retained')
        if sum(bycell.values())!=counts['accepted_chrM_alignments'] or primary>chunk['candidate_records']:
            raise ValueError('Alignment counts disagree')
        item={'chunk':number,'candidate_records':chunk['candidate_records'],'primary_queries':primary,
              'counts':counts,'accepted_reads_by_cell':bycell,'BAM':str(bam),'BAM_bytes':bam.stat().st_size,
              'command':command,'attempt':attempt}
        writing=folder/f'chunk_{number:04d}.json.writing';preserve(writing)
        writing.write_text(json.dumps(item,indent=2)+'\n');writing.replace(done);completed.append(item)
        print(json.dumps({'aligned_chunks':number+1,'chunks':len(chunks),'accepted_chrM':counts['accepted_chrM_alignments']}),flush=True)
    counts=Counter();bycell=Counter()
    for item in completed:counts.update(item['counts']);bycell.update(item['accepted_reads_by_cell'])
    per_cell=a.output/'accepted_MT_reads_by_cell.json';preserve(per_cell)
    per_cell.write_text(json.dumps(bycell)+'\n')
    unsorted=a.output/'all.accepted.unsorted.bam';sorted_bam=a.output/'all.accepted.CB.bam'
    preserve(unsorted);preserve(sorted_bam)
    pysam.cat('--no-PG','-o',str(unsorted),*[x['BAM'] for x in completed])
    pysam.sort('-t','CB','-@','4','-m','128M','--no-PG','-T',str(a.output/'CBsort.tmp'),'-o',str(sorted_bam),str(unsorted))
    observed=Counter();previous=''
    with pysam.AlignmentFile(sorted_bam,'rb') as source:
        for record in source.fetch(until_eof=True):
            cell=record.get_tag('CB')
            if previous>cell:raise ValueError('CB sort order differs')
            previous=cell;observed[cell]+=1
    if observed!=bycell:raise ValueError('CB-sort counts differ')
    result={'status':'COMPLETE','chunks':len(completed),'candidate_records':sum(x['candidate_records'] for x in completed),
            'alignment_counts':counts,'accepted_MT_reads':sum(bycell.values()),'accepted_MT_read_cells':len(bycell),
            'CB_sort_count_and_order_exact':True,'full_nuclear_competition':True,'full_genome_BAM_written':False,
            'minimap_version':version,'resumed_completed_chunks':bool(a.resume),'seconds':round(time.monotonic()-start,3),
            'peak_child_RSS_KiB':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}
    (a.output/'alignment_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
