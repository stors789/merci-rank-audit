"""Stream P2 paired FASTQs and generate the original GEXSCOPE candidate chunks."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import time

import pandas as pd


def read_record(handle, first):
    seq, plus, quality = [handle.readline().rstrip(b'\r\n') for _ in range(3)]
    if not first.startswith(b'@') or not plus.startswith(b'+') or not seq or len(seq) != len(quality):
        raise ValueError('Malformed FASTQ record')
    return [first.rstrip(b'\r\n'), seq, plus, quality]


def selected(a, b, whitelist, kmers, counts):
    barcode = a[1][0:8] + a[1][24:32] + a[1][48:56]
    if barcode not in whitelist:
        return None
    counts['exact_called_barcode'] += 1; umi = a[1][57:69]
    if not (len(umi) == 12 and all(x in b'ACGT' for x in umi) and min(a[3][57:69]) >= 53):
        return None
    counts['called_and_UMI_Q20'] += 1
    if a[1][69:87].count(b'T') < 10:
        return None
    counts['called_UMI_polyT'] += 1; seq = b[1].upper()
    if not (len(seq) >= 50 and (any(seq[j:j+21] in kmers for j in range(0,len(seq)-20,10)) or seq[-21:] in kmers)):
        return None
    counts['mitochondrial_kmer_candidate'] += 1
    return barcode, umi, seq, b[3]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('r1', 'r2', 'roles', 'reference', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args(); started = time.monotonic()
    for path in (a.r1, a.r2, a.roles, a.reference):
        if not path.is_file():
            raise FileNotFoundError(path)
    roles = pd.read_csv(a.roles, sep='\t', dtype={'barcode': str})
    if len(roles) != 4359 or roles.barcode.duplicated().any() or not roles.cell_id.eq('GSM7774400:'+roles.barcode).all():
        raise ValueError('New P2 roles must retain the exact 24-base barcode namespace')
    whitelist = {x.encode() for x in roles.barcode}
    seq = b''.join(x.strip() for x in a.reference.read_bytes().splitlines() if not x.startswith(b'>')).upper()
    if len(seq) != 16569:
        raise ValueError('Expected chrM extracted from the original full hg38 reference')
    reverse = seq.translate(bytes.maketrans(b'ACGT',b'TGCA'))[::-1]
    kmers = {s[i:i+21] for s in (seq,reverse) for i in range(len(s)-20)}
    a.output.mkdir(parents=True, mode=0o700, exist_ok=False); folder = a.output/'candidates'; folder.mkdir()
    counts = Counter(); bycell = Counter(); chunks = []; out = None; chunk_reads = 0; chunk_begin = 1; actual = 0

    def finish(end):
        nonlocal out, chunk_reads, chunk_begin
        if out is not None:
            out.close(); out = None
            chunks.append({'path': str(destination), 'candidate_records': chunk_reads,
                           'raw_pair_start': chunk_begin, 'raw_pair_end': end})
            chunk_reads = 0

    try:
        with gzip.open(a.r1, 'rb') as left, gzip.open(a.r2, 'rb') as right:
            while True:
                h1, h2 = left.readline(), right.readline()
                if not h1 or not h2:
                    if h1 or h2:
                        raise ValueError('Unequal mate EOF')
                    break
                r1, r2 = read_record(left, h1), read_record(right, h2); actual += 1
                if r1[0].split()[0] != r2[0].split()[0]:
                    raise ValueError('Mismatched pair IDs at index ' + str(actual))
                counts['paired_records_read'] += 1
                pick = selected(r1,r2,whitelist,kmers,counts)
                if pick:
                    barcode,umi,rna,quality = pick; bycell[barcode.decode()] += 1
                    if out is None:
                        destination = (folder/f'chunk_{len(chunks):04d}.fastq.gz').resolve()
                        out = gzip.open(destination,'wb',compresslevel=2)
                        chunk_begin = 1 if actual <= 20000000 and not chunks else actual
                    name = b'@'+barcode+b'_'+umi+b'_'+str(actual).encode()
                    out.write(b'\n'.join([name,rna,b'+',quality])+b'\n'); chunk_reads += 1
                    if actual > 20000000 and chunk_reads == 4000000:
                        finish(actual)
                if actual == 20000000:
                    finish(actual)
                if actual % 10000000 == 0:
                    print(json.dumps({'paired_records_read':actual,'candidates':counts['mitochondrial_kmer_candidate'],
                                      'seconds':round(time.monotonic()-started,3)}),flush=True)
        finish(actual)
    finally:
        if out is not None:
            out.close()
    if actual != 417883844 or sum(bycell.values()) != counts['mitochondrial_kmer_candidate'] or sum(x['candidate_records'] for x in chunks) != counts['mitochondrial_kmer_candidate']:
        raise ValueError('Full-run pair or candidate counts differ; first output retained')
    (a.output/'candidate_chunks.json').write_text(json.dumps(chunks,indent=2)+'\n')
    (a.output/'candidate_reads_by_cell.json').write_text(json.dumps(bycell)+'\n')
    (a.output/'allele_cell_barcodes.json').write_text(json.dumps(sorted(bycell))+'\n')
    receipt={'status':'COMPLETE','pairs':actual,'counts':counts,'candidate_cells':len(bycell),'chunks':len(chunks),
             'complete_pair_ID_and_gzip_EOF_checks':True,'new_prefix_and_tail_both_generated':True,
             'seconds':round(time.monotonic()-started,3)}
    (a.output/'extraction_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)


if __name__ == '__main__':
    main()
