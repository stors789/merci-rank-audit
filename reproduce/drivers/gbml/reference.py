"""Extract chrM from the supplied hg38 FASTA and encode the original 21-mer set."""
import argparse
import gzip
import hashlib
import json
import struct
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, mode=0o700, exist_ok=False)
    opener = gzip.open if a.reference.suffix == '.gz' else open
    chunks = []; found = False
    with opener(a.reference, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if found:
                    break
                found = line[1:].split()[0] == 'chrM'
            elif found:
                chunks.append(line.strip())
    seq = ''.join(chunks).upper().encode()
    if len(seq) != 16569 or set(seq) - set(b'ACGTN'):
        raise ValueError('Expected the UCSC hg38 chrM sequence')
    fasta = a.output / 'chrM.fa'
    fasta.write_bytes(b'>chrM\n' + b'\n'.join(seq[i:i+50] for i in range(0, len(seq), 50)) + b'\n')
    rc = seq.translate(bytes.maketrans(b'ACGT', b'TGCA'))[::-1]
    kmers = {s[i:i+21] for s in (seq, rc) for i in range(len(s)-20)}
    codes = []
    for k in kmers:
        value = 0
        for base in k:
            value = (value << 3) | b'ACGTN'.index(base)
        codes.append(value)
    if len(codes) != len(set(codes)):
        raise ValueError('3-bit kmer collision')
    out = a.output / 'chrM_21mers.u64'
    out.write_bytes(struct.pack('<Q', len(codes)) + b''.join(struct.pack('<Q', x) for x in sorted(codes)))
    result = {'status': 'COMPLETE', 'reference': 'UCSC hg38 chrM', 'length': len(seq), 'n_count': seq.count(b'N'),
              'sequence_sha256': hashlib.sha256(seq).hexdigest(), 'kmer_count': len(codes),
              'kmer_sha256': hashlib.sha256(out.read_bytes()).hexdigest(),
              'criterion': 'every 21-mer in both reference orientations; read offsets every 10 plus final window'}
    (a.output / 'reference_receipt.json').write_text(json.dumps(result, indent=2) + '\n'); print(json.dumps(result))


if __name__ == '__main__':
    main()
