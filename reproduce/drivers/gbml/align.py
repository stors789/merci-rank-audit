"""Competitive hg38 alignment of the freshly extracted GBML181 candidates."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import time

import pysam

MAX_SPOT = 98159746 + 102126028


def accept(line, counts):
    f = line.rstrip('\r\n').split('\t'); flag = int(f[1])
    if flag & (4 | 256 | 2048):
        counts['unmapped_or_secondary_supplementary'] += 1
        return None
    counts['primary_' + f[2]] += 1
    if f[2] != 'chrM' or int(f[4]) < 30:
        return None
    cigar = [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', f[5])]
    if any(op == 'N' for _, op in cigar):
        counts['chrM_spliced_alignment_skipped'] += 1
        return None
    if sum(n for n, op in cigar if op in 'M=X') < 60:
        return None
    counts['accepted_chrM_alignments'] += 1
    counts['accepted_strand_reverse' if flag & 16 else 'accepted_strand_forward'] += 1
    return f[:11]


def child_limit():
    resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('extraction', 'roles', 'index', 'minimap', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args(); started = time.monotonic()
    receipt = json.loads((a.extraction / 'extraction_receipt.json').read_text())
    candidate = a.extraction / 'all_candidates.fastq.gz'
    h = hashlib.sha256()
    with candidate.open('rb') as f:
        while chunk := f.read(8*1024*1024):
            h.update(chunk)
    if receipt['status'] != 'COMPLETE' or h.hexdigest() != receipt['combined_sha256']:
        raise ValueError('Newly generated candidate bytes differ from their receipt')
    if a.index.stat().st_size != 15420201719:
        raise ValueError('Expected the original UCSC hg38 minimap2 2.30 split index')
    version = subprocess.check_output([str(a.minimap), '--version'], text=True).strip()
    if version != '2.30-r1287':
        raise ValueError('Expected minimap2 2.30-r1287')
    allowed = {line.split('\t')[0].removesuffix('-1') for line in a.roles.read_text().splitlines()[1:]}
    a.output.mkdir(parents=True, mode=0o700, exist_ok=False)
    split = a.output / 'minimap_split'; split.mkdir()
    final = a.output / 'accepted_chrM.unsorted.bam'; writing = a.output / 'accepted_chrM.unsorted.bam.writing'
    command = [str(a.minimap.resolve()), '-ax', 'splice:sr', '-t', '1', '-K', '20M',
               '--split-prefix', str(split / 'all'), '--sam-hit-only', '-N', '5',
               str(a.index.resolve()), str(candidate.resolve())]
    counts = Counter(); bycell = Counter(); seen = bytearray((MAX_SPOT+7)//8); primary = 0
    header = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6', 'SO': 'unsorted'},
             'SQ': [{'SN': 'chrM', 'LN': 16569}],
             'CO': ['Accepted chrM records only; full nuclear competition counts are retained separately.']})
    with (a.output / 'minimap.log').open('x') as err, pysam.AlignmentFile(writing, 'wb', header=header, threads=1) as bam:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=err, text=True, bufsize=1024**2, preexec_fn=child_limit)
        print(json.dumps({'status': 'started', 'aligner_pid': process.pid, 'candidate_reads': receipt['candidate_records']}), flush=True)
        try:
            for line in process.stdout:
                if line.startswith('@'):
                    continue
                fields = line.rstrip('\n').split('\t'); flag = int(fields[1]); cell, umi, number = fields[0].split('_'); number = int(number)
                if cell not in allowed or len(umi) != 12 or not 1 <= number <= MAX_SPOT:
                    raise ValueError('Invalid alignment query identity')
                if not flag & (4 | 256 | 2048):
                    slot, bit = divmod(number-1, 8)
                    if seen[slot] & (1 << bit):
                        raise ValueError('Duplicate primary query')
                    seen[slot] |= 1 << bit; primary += 1
                fields = accept(line, counts)
                if fields is not None:
                    record = pysam.AlignedSegment.fromstring('\t'.join(fields), bam.header)
                    record.set_tag('CB', cell); record.set_tag('UB', umi); bam.write(record); bycell[cell] += 1
            process.stdout.close(); code = process.wait()
        finally:
            if process.poll() is None:
                process.terminate(); process.wait()
    if code != 0:
        raise RuntimeError('Competitive alignment failed; partial BAM and log retained')
    if primary > receipt['candidate_records'] or sum(bycell.values()) != counts['accepted_chrM_alignments']:
        raise ValueError('Alignment totals disagree')
    writing.replace(final); (a.output / 'accepted_reads_by_cell.json').write_text(json.dumps(bycell) + '\n')
    result = {'status': 'COMPLETE', 'command': command, 'minimap_version': version,
              'candidate_reads': receipt['candidate_records'], 'primary_queries': primary,
              'alignment_counts': counts, 'accepted_chrM_reads': sum(bycell.values()), 'accepted_cells': len(bycell),
              'primary_query_unique': True, 'full_genome_BAM_written': False, 'accepted_chrM_bam_bytes': final.stat().st_size,
              'seconds': round(time.monotonic()-started, 3), 'peak_self_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'peak_child_RSS_KiB': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}
    (a.output / 'alignment_receipt.json').write_text(json.dumps(result, indent=2) + '\n'); print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
