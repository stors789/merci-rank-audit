"""Test alignment recovery with two synthetic chunks and a fake minimap process."""
import argparse
import gzip
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

FAKE_MINIMAP = '''#!/usr/bin/env python3
import gzip
from pathlib import Path
import sys
if sys.argv[1:] == ['--version']:
    print('2.30-r1287')
    raise SystemExit(0)
query=Path(sys.argv[-1]); marker=Path(__file__).with_suffix('.failed')
with gzip.open(query,'rt') as source:
    while True:
        name=source.readline().strip()
        if not name: break
        sequence=source.readline().strip(); source.readline(); quality=source.readline().strip()
        print('\\t'.join([name[1:],'0','chrM','1','60',str(len(sequence))+'M','*','0','0',sequence,quality]))
        if query.name == 'chunk_0001.fastq.gz' and not marker.exists():
            marker.write_text('first failure intentionally retained')
            raise SystemExit(7)
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, mode=0o700, exist_ok=False)
    extraction = root / 'extraction'
    extraction.mkdir()
    cell, umi = 'A' * 24, 'ACGT' * 3
    chunks = []
    for number, indices in enumerate(((1, 2), (3,))):
        path = extraction / f'chunk_{number:04d}.fastq.gz'
        with gzip.open(path, 'wt') as f:
            for index in indices:
                f.write(f'@{cell}_{umi}_{index}\n{"A" * 70}\n+\n{"I" * 70}\n')
        chunks.append({'path': str(path), 'candidate_records': len(indices),
                       'raw_pair_start': min(indices), 'raw_pair_end': max(indices)})
    (extraction / 'candidate_chunks.json').write_text(json.dumps(chunks) + '\n')
    (extraction / 'allele_cell_barcodes.json').write_text(json.dumps([cell]) + '\n')
    (extraction / 'extraction_receipt.json').write_text(json.dumps(
        {'status': 'COMPLETE', 'counts': {'mitochondrial_kmer_candidate': 3}}) + '\n')
    roles = root / 'roles.tsv'
    pd.DataFrame({'barcode': [cell]}).to_csv(roles, sep='\t', index=False)
    index = root / 'fake_sparse_index.mmi'
    with index.open('wb') as f:
        f.truncate(15420201719)
    minimap = root / 'fake_minimap'
    minimap.write_text(FAKE_MINIMAP)
    minimap.chmod(0o755)
    script = Path(__file__).resolve().parents[2] / 'drivers/p2/align_candidates.py'
    output = root / 'alignment'
    command = [sys.executable, str(script), '--extraction', str(extraction), '--roles', str(roles),
               '--index', str(index), '--minimap', str(minimap), '--output', str(output)]
    first = subprocess.run(command, capture_output=True, text=True)
    (root / 'first.log').write_text(first.stdout + first.stderr)
    done = output / 'aligned/chunk_0000.json'
    bam = output / 'aligned/chunk_0000.bam'
    if first.returncode == 0 or not done.is_file() or not bam.is_file() or (output / 'alignment_receipt.json').exists():
        raise AssertionError('The first synthetic failure was not retained at a chunk boundary')
    old_receipt, old_bam = done.read_bytes(), bam.read_bytes()
    original_roles = roles.read_bytes()
    roles.write_bytes(original_roles + b'changed\n')
    changed = subprocess.run(command + ['--resume'], capture_output=True, text=True)
    (root / 'changed_input.log').write_text(changed.stdout + changed.stderr)
    if changed.returncode == 0 or 'inputs changed' not in changed.stderr:
        raise AssertionError('Changed input was not rejected before recovery')
    roles.write_bytes(original_roles)
    resumed = subprocess.run(command + ['--resume'], capture_output=True, text=True)
    (root / 'resumed.log').write_text(resumed.stdout + resumed.stderr)
    if resumed.returncode:
        raise AssertionError('Synthetic chunk recovery failed; inspect resumed.log')
    receipt = json.loads((output / 'alignment_receipt.json').read_text())
    if (receipt['status'] != 'COMPLETE' or receipt['accepted_MT_reads'] != 3 or
            done.read_bytes() != old_receipt or bam.read_bytes() != old_bam or
            not list((output / 'aligned').glob('chunk_0001.bam.interrupted_*'))):
        raise AssertionError('Recovered alignment or preserved prior artifacts differ')
    consensus = root / 'consensus'
    consensus_driver = script.parent / 'consensus.py'
    completed = subprocess.run([sys.executable, str(consensus_driver), '--extraction', str(extraction),
                                '--alignment', str(output), '--output', str(consensus)], capture_output=True, text=True)
    (root / 'consensus.log').write_text(completed.stdout + completed.stderr)
    if completed.returncode:
        raise AssertionError('Synthetic joint consensus failed; inspect consensus.log')
    consensus_receipt = json.loads((consensus / 'consensus_receipt.json').read_text())
    if (consensus_receipt['UMI_molecules'] != 1 or consensus_receipt['UMI_positions'] != 60 or
            consensus_receipt['counts']['accepted_chrM_alignments'] != 3):
        raise AssertionError('The repeated cross-chunk molecule was not joined once')
    print('PASS: alignment recovery and cross-chunk same-cell/UMI consensus')


if __name__ == '__main__':
    main()
