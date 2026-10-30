"""Check incomplete and mismatched synthetic FASTQs retain failure evidence."""
import argparse
import gzip
from pathlib import Path
import subprocess
import sys

import pandas as pd


def barcode(number):
    bases = []
    for _ in range(24):
        number, digit = divmod(number, 4)
        bases.append('ACGT'[digit])
    return ''.join(reversed(bases))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, mode=0o700, exist_ok=False)
    barcodes = [barcode(i) for i in range(4359)]
    pd.DataFrame({'barcode': barcodes, 'cell_id': ['GSM7774400:' + x for x in barcodes]}).to_csv(
        args.output / 'roles.tsv', sep='\t', index=False)
    (args.output / 'chrM.fa').write_text('>chrM\n' + 'A' * 16569 + '\n')
    r1 = list('C' * 87)
    cell = barcodes[0]
    r1[0:8], r1[24:32], r1[48:56] = cell[:8], cell[8:16], cell[16:24]
    r1[57:69] = 'ACGTACGTACGT'
    r1[69:87] = 'T' * 18
    sequence = ''.join(r1)
    script = Path(__file__).resolve().parents[2] / 'drivers/p2/extract_fastq.py'
    for name, mismatch in (('incomplete', False), ('mismatched', True)):
        left, right = args.output / (name + '_R1.gz'), args.output / (name + '_R2.gz')
        with gzip.open(left, 'wt') as first, gzip.open(right, 'wt') as second:
            for i in range(3):
                first.write(f'@pair{i}\n{sequence}\n+\n{"I" * 87}\n')
                header = 'wrong' if mismatch and i == 1 else f'pair{i}'
                second.write(f'@{header}\n{"A" * 70}\n+\n{"I" * 70}\n')
        out = args.output / name
        result = subprocess.run([sys.executable, str(script), '--r1', str(left), '--r2', str(right),
                                 '--roles', str(args.output / 'roles.tsv'), '--reference', str(args.output / 'chrM.fa'),
                                 '--output', str(out)], capture_output=True, text=True)
        (args.output / (name + '.log')).write_text(result.stdout + result.stderr)
        expected = 'Mismatched pair IDs' if mismatch else 'Full-run pair or candidate counts differ'
        if result.returncode == 0 or expected not in result.stderr or (out / 'extraction_receipt.json').exists():
            raise AssertionError(name + ' was not retained as a failed extraction')
        chunk = out / 'candidates/chunk_0000.fastq.gz'
        with gzip.open(chunk, 'rt') as f:
            lines = f.read().splitlines()
        expected_records = 1 if mismatch else 3
        if len(lines) != 4 * expected_records or not lines[0].startswith('@' + cell + '_'):
            raise AssertionError(name + ' did not retain its readable partial candidate chunk')
    print('PASS: incomplete and mismatched FASTQs retain readable partial chunks without a completion receipt')


if __name__ == '__main__':
    main()
