"""Check original alignment filters and UMI consensus boundary behavior."""
import argparse
from collections import Counter
from pathlib import Path
import sys

repro = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(repro / 'drivers/p2'))
sys.path.insert(0, str(repro / 'drivers/gbml'))
from align_candidates import accept
from pileup import consensus_matrix, LENGTH


def sam(flag=0, contig='chrM', mapq=30, cigar='70M', name='read', sequence='A' * 70, quality='I' * 70):
    return '\t'.join([name, str(flag), contig, '1', str(mapq), cigar, '*', '0', '0', sequence, quality]) + '\n'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, mode=0o700, exist_ok=False)
    counts = Counter()
    assert accept(sam(mapq=30, cigar='60M'), counts) is not None
    assert accept(sam(mapq=29), counts) is None
    assert accept(sam(cigar='59M'), counts) is None
    assert accept(sam(cigar='30M1N30M'), counts) is None
    assert accept(sam(flag=256), counts) is None
    assert accept(sam(contig='chr1'), counts) is None
    assert accept(sam(flag=16, cigar='30M1D30M'), counts) is not None
    assert counts['accepted_chrM_alignments'] == 2 and counts['accepted_strand_reverse'] == 1
    cell, umi = 'A' * 24, 'ACGT' * 3
    sequence_a = list('A' * 70); sequence_a[0] = 'C'
    sequence_b = list('A' * 70); sequence_b[10] = 'C'
    quality = list('I' * 70); quality[5] = '!'
    path = a.output / 'synthetic.sam'
    path.write_text(sam(name=cell + '_' + umi + '_1', sequence=''.join(sequence_a), quality=''.join(quality)) +
                    sam(name=cell + '_' + umi + '_2', sequence=''.join(sequence_b), quality=''.join(quality)))
    matrix, result, molecules, positions = consensus_matrix(path, [cell])
    assert molecules == 1 and positions == 58 and int(matrix.sum()) == 58
    assert result['accepted_chrM_alignments'] == 2 and result['conflicting_UMI_position_discarded'] == 1
    columns = set(matrix.indices.tolist())
    assert all((column % LENGTH) not in (0, 5, 10, 69) for column in columns)
    assert any((column % LENGTH) == 6 for column in columns)
    print('PASS: primary/MAPQ/CIGAR boundaries, base quality, end trimming and conflict discard')


if __name__ == '__main__':
    main()
