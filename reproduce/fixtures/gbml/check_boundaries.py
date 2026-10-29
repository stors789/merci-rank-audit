"""Synthetic checks for candidate kmers, alignment filters and joint UMI consensus."""
import argparse
from collections import Counter
import importlib.util
import os
from pathlib import Path
import struct
import subprocess
import sys

repro = Path(__file__).resolve().parents[2]
driver = repro / 'drivers/gbml'; sys.path.insert(0, str(driver))
from pileup import consensus_matrix
from align import accept


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--extractor', type=Path, required=True); p.add_argument('--vdb-lib', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=False)
    k = b'ACGTACGTACNTACGTACGTA'; value = 0
    for ch in k:
        value = (value << 3) | b'ACGTN'.index(ch)
    binary = a.output / 'kmers.u64'; binary.write_bytes(struct.pack('<QQ', 1, value))
    reads = [b'G'*10 + k + b'T'*39, b'T'*49 + k, b'A'*70, b'A'*20, b'G'*9 + k + b'T'*40]
    expected = [1, 1, 0, 0, 0]
    env = os.environ.copy(); env['LD_LIBRARY_PATH'] = str(a.vdb_lib.resolve())
    output = subprocess.check_output([str(a.extractor.resolve()), '--kmer-test', str(binary.resolve())],
                                     input=b'\n'.join(reads)+b'\n', env=env)
    assert [int(x) for x in output.splitlines()] == expected
    counts = Counter()
    def line(flag=0, chrom='chrM', mapq=30, cigar='60M'):
        return f'ACGTACGTACGTACGT_AAAAAAAAAAAA_1\t{flag}\t{chrom}\t1\t{mapq}\t{cigar}\t*\t0\t0\t' + 'A'*60 + '\t' + 'I'*60 + '\n'
    assert accept(line(), counts) is not None
    assert accept(line(mapq=29), counts) is None and accept(line(chrom='chr1'), counts) is None
    assert accept(line(cigar='59M'), counts) is None and accept(line(cigar='30M1N30M'), counts) is None
    assert accept(line(flag=256), counts) is None
    assert counts['accepted_chrM_alignments'] == 1 and counts['chrM_spliced_alignment_skipped'] == 1
    seq1 = 'A'*70; seq2 = seq1[:6] + 'C' + seq1[7:]
    sam = a.output / 'two_runs.sam'
    def record(number, seq):
        return f'ACGTACGTACGTACGT_AAAAAAAAAAAA_{number}\t0\tchrM\t1\t60\t70M\t*\t0\t0\t{seq}\t' + 'I'*70 + '\n'
    sam.write_text(record(1, seq1) + record(98159747, seq2))
    matrix, observed, molecules, positions = consensus_matrix(sam, ['ACGTACGTACGTACGT'])
    assert molecules == 1 and positions == 59 and int(matrix.sum()) == 59
    assert observed['accepted_chrM_alignments'] == 2 and observed['conflicting_UMI_position_discarded'] == 1
    print('PASS: 3-bit N-containing kmers, candidate offsets, competitive filters and cross-run UMI conflict')


if __name__ == '__main__':
    main()
