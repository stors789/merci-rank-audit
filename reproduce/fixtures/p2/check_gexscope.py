"""Check the original 24-base barcode, UMI quality, poly-T and 21-mer boundaries."""
from collections import Counter
import io
from pathlib import Path
import sys

repro = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(repro/'drivers/p2'))
from extract_fastq import read_record, selected

barcode = b'ACGTACGT' * 3; umi = b'ACGTACGTACGT'; pattern = b'ACGTACGTACGTACGTACGTA'
seq = bytearray(b'C'*87);seq[0:8]=barcode[:8];seq[24:32]=barcode[8:16];seq[48:56]=barcode[16:24]
seq[57:69]=umi;seq[69:87]=b'T'*18
quality = b'I'*87; rna = b'G'*10 + pattern + b'C'*39


def record(s, q):
    return [b'@test',bytes(s),b'+',bytes(q)]


def run(a,b):
    return selected(a,b,{barcode},{pattern},Counter())


assert run(record(seq,quality),record(rna,b'I'*len(rna))) == (barcode,umi,rna,b'I'*len(rna))
bad = bytearray(seq);bad[0] = ord('T');assert run(record(bad,quality),record(rna,b'I'*len(rna))) is None
bad = bytearray(seq);bad[57] = ord('N');assert run(record(bad,quality),record(rna,b'I'*len(rna))) is None
low = bytearray(quality);low[57] = 52;assert run(record(seq,low),record(rna,b'I'*len(rna))) is None
bad = bytearray(seq);bad[69:87] = b'T'*9+b'C'*9;assert run(record(bad,quality),record(rna,b'I'*len(rna))) is None
bad_rna = b'G'*70;assert run(record(seq,quality),record(bad_rna,b'I'*len(bad_rna))) is None
end = b'G'*49 + pattern;assert run(record(seq,quality),record(end,b'I'*len(end))) is not None
assert read_record(io.BytesIO(b'ACGT\n+identifier\nIIII\n'),b'@test\n')[2] == b'+identifier'
try:
    read_record(io.BytesIO(b'ACGT\n-identifier\nIIII\n'),b'@test\n')
except ValueError:
    pass
else:
    raise AssertionError('Invalid FASTQ plus line was accepted')
print('PASS: GEXSCOPE barcode, UMI Q20, poly-T, kmer boundaries and original FASTQ plus-prefix rule')
