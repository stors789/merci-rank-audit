"""Combine accepted reads from both runs before cell/UMI/position consensus."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import resource
import time

import pysam
from scipy import sparse
from pileup import consensus_matrix, LENGTH


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--alignment', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); started = time.monotonic()
    align = json.loads((a.alignment / 'alignment_receipt.json').read_text())
    if align['status'] != 'COMPLETE':
        raise ValueError('New alignment is not complete')
    bycell = json.loads((a.alignment / 'accepted_reads_by_cell.json').read_text()); cells = sorted(bycell)
    a.output.mkdir(parents=True, mode=0o700, exist_ok=False)
    raw = a.alignment / 'accepted_chrM.unsorted.bam'; sorted_bam = a.output / 'accepted_chrM.CB.bam'
    writing = a.output / 'accepted_chrM.CB.bam.writing'
    pysam.sort('-t', 'CB', '-@', '1', '-m', '256M', '--no-PG', '-O', 'BAM',
               '-T', str(a.output / 'sort_tmp'), '-o', str(writing), str(raw))
    writing.replace(sorted_bam)
    actual = Counter(); previous = ''
    with pysam.AlignmentFile(sorted_bam, 'rb') as f:
        for record in f.fetch(until_eof=True):
            cell = record.get_tag('CB')
            if cell < previous:
                raise ValueError('CB order invalid')
            previous = cell; actual[cell] += 1
    if dict(actual) != {k: int(v) for k, v in bycell.items()}:
        raise ValueError('Accepted reads changed during BAM sorting')
    folder = a.output / 'consensus_batches'; folder.mkdir(); matrices = []; counts = Counter(); molecules = positions = 0; groups = []
    with pysam.AlignmentFile(sorted_bam, 'rb') as f:
        iterator = f.fetch(until_eof=True); record = next(iterator, None)
        for start in range(0, len(cells), 16):
            group = cells[start:start+16]; allowed = set(group); idx = start//16
            sam = folder / f'batch_{idx:04d}.sam'; npz = folder / f'batch_{idx:04d}.npz'
            with sam.open('x') as handle:
                while record is not None and record.get_tag('CB') <= group[-1]:
                    if record.get_tag('CB') not in allowed:
                        raise ValueError('Unexpected sorted cell')
                    handle.write(record.to_string() + '\n'); record = next(iterator, None)
            if sam.stat().st_size:
                matrix, batch_counts, n, pos = consensus_matrix(sam, group)
            else:
                matrix = sparse.csr_matrix((len(group), 8*LENGTH), dtype='int32'); batch_counts = Counter(); n = pos = 0
            sparse.save_npz(npz, matrix); matrices.append(matrix); counts.update(batch_counts); molecules += n; positions += pos
            groups.append({'batch': idx, 'cells': len(group), 'molecules': n, 'UMI_positions': pos,
                           'accepted_chrM_reads': batch_counts['accepted_chrM_alignments']})
        if record is not None:
            raise ValueError('Unconsumed accepted records')
    if counts['accepted_chrM_alignments'] != align['accepted_chrM_reads']:
        raise ValueError('Consensus read count mismatch')
    full = sparse.vstack(matrices, format='csr') if matrices else sparse.csr_matrix((0, 8*LENGTH), dtype='int32')
    if full.shape != (len(cells), 8*LENGTH):
        raise ValueError('Consensus matrix dimensions')
    sparse.save_npz(a.output / 'mt_UMI_alleles.npz', full)
    (a.output / 'allele_cell_barcodes.json').write_text(json.dumps(cells) + '\n')
    with (a.output / 'consensus_batches.csv').open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['batch','cells','molecules','UMI_positions','accepted_chrM_reads'], lineterminator='\n')
        writer.writeheader(); writer.writerows(groups)
    result = {'status': 'COMPLETE', 'accepted_chrM_reads': align['accepted_chrM_reads'], 'accepted_cells': len(cells),
              'CB_sorted_counts_match': True, 'same_cell_UMI_combined_across_runs_before_consensus': True,
              'consensus_counts': counts, 'molecules': molecules, 'UMI_positions': positions, 'matrix_shape': list(full.shape),
              'split_note': 'two deterministic complementary columns are combined, not independent samples',
              'seconds': round(time.monotonic()-started, 3), 'max_rss_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (a.output / 'consensus_receipt.json').write_text(json.dumps(result, indent=2) + '\n'); print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
