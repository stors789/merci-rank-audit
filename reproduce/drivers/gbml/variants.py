"""Adapt full-run UMI consensuses to the original MERCI 10x variant format."""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import os
import resource
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from upstream import check as check_upstream

LENGTH = 16569
BASES = np.array(list('ACGT'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('consensus', 'alignment', 'roles', 'reference', 'upstream', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    output = args.output.resolve(); receipt = output / 'variant_input_receipt.json'
    if output.exists():
        raise FileExistsError('Protect first MERCI input preparation')
    if json.loads((args.consensus / 'consensus_receipt.json').read_text())['status'] != 'COMPLETE':
        raise ValueError('New UMI consensus has not completed')
    roles_path = args.roles.resolve()
    fresh = json.loads((roles_path.parent / 'roles_receipt.json').read_text())
    if hashlib.sha256(roles_path.read_bytes()).hexdigest() != fresh['role_sha256']:
        raise ValueError('Newly generated roles changed')
    check_upstream(args.upstream.resolve())
    upstream = args.upstream.resolve() / 'MERCI-mtSNP.py'
    reference_path = args.reference.resolve()
    roles = pd.read_csv(roles_path, sep='\t', dtype=str)
    roles['cell'] = roles.barcode.str.removesuffix('-1')
    if roles.cell.duplicated().any() or roles.role.value_counts().to_dict() != {'source': 108, 'target': 2001}:
        raise ValueError('Unexpected frozen roles')
    role = roles.set_index('cell').role.to_dict()
    cells = json.loads((args.consensus / 'allele_cell_barcodes.json').read_text())
    reads = json.loads((args.alignment / 'accepted_reads_by_cell.json').read_text())
    if len(set(cells)) != len(cells) or set(cells) != set(reads) or not set(cells) <= set(role):
        raise ValueError('Consensus cell axis does not match accepted alignments and frozen roles')
    matrix = sparse.load_npz(args.consensus / 'mt_UMI_alleles.npz').tocsr()
    if matrix.shape != (len(cells), 8 * LENGTH):
        raise ValueError('Unexpected consensus shape')
    bases = matrix[:, :4 * LENGTH] + matrix[:, 4 * LENGTH:]
    depth = sum((bases[:, i * LENGTH:(i + 1) * LENGTH] for i in range(4)),
                sparse.csr_matrix((len(cells), LENGTH), dtype=np.int32)).tocsr()
    reference = ''.join(x.strip() for x in reference_path.read_text().splitlines() if not x.startswith('>'))
    if len(reference) != LENGTH or set(reference) - set('ACGTN'):
        raise ValueError('Unexpected chrM reference')
    reference = np.array(list(reference))
    mutations = {}
    for i, cell in enumerate(cells):
        base_row = bases.getrow(i)
        position = base_row.indices % LENGTH
        alternate = BASES[base_row.indices // LENGTH]
        selected = alternate != reference[position]
        pos, alt, supports = position[selected], alternate[selected], base_row.data[selected]
        coverage = np.asarray(depth.getrow(i)[:, pos].toarray()).ravel()
        if np.any(supports <= 0) or np.any(coverage < supports):
            raise ValueError('ALT counts exceed coverage')
        records = {}
        for p, a, n, d in zip(pos, alt, supports, coverage):
            p, n, d = int(p), int(n), int(d)
            identifier = f'MT_{p}_{reference[p]}-{a}'
            if identifier in records:
                raise ValueError('Duplicate cell-variant observation')
            records[identifier] = [p, n, d, n / d]
        mutations[cell] = records
    before = sum(map(len, mutations.values()))
    output.mkdir()
    os.chdir(output)
    spec = importlib.util.spec_from_file_location('original_merci_mtsnp', upstream)
    caller = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(caller)
    mutations = caller.MTmutations_filters(mutations, 'GBML181_UMI', ln=5)
    after = sum(map(len, mutations.values()))
    variants_path = output / 'MT_variants.tsv'
    positions = set()
    variant_cells = set()
    retained_records = 0
    with variants_path.open('x', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['Cell', 'ID', 'pos', 's_reads', 't_reads', 'AF', 'Cell_reads'])
        for cell in cells:
            for identifier, (pos, alt, cov, af) in mutations[cell].items():
                writer.writerow([cell, identifier, pos, alt, cov, format(af, '.17g'), reads[cell]])
                if reads[cell] > 10:
                    positions.add(pos)
                    variant_cells.add(cell)
                    retained_records += 1
    positions = sorted(positions)
    if positions:
        cov = pd.DataFrame(depth[:, positions].toarray().T, index=positions, columns=cells)
        cov.to_csv(output / 'Coverage.tsv.gz', sep='\t', index_label='position0', compression='gzip')
    for name, group in [('donors.txt', 'source'), ('targets.txt', 'target')]:
        (output / name).write_text(''.join(x + '\n' for x in sorted(roles.loc[roles.role == group, 'cell'])))
    counts = {}
    for group in ['source', 'target']:
        members = roles.loc[roles.role == group, 'cell'].tolist()
        counts[group] = dict(frozen=len(members), accepted_alignment=int(sum(x in reads for x in members)),
                             accepted_reads_gt10=int(sum(reads.get(x, 0) > 10 for x in members)),
                             variant_record_after_minReads=int(sum(x in variant_cells for x in members)))
    summary = dict(complete=True, roles=counts, alt_records_before_original_nearby_filter=before,
                   alt_records_after_original_nearby_filter=after,
                   alt_records_after_strict_minReads=retained_records,
                   covered_positions_used=len(positions), reference_coordinates='zero-based original MERCI-mtSNP convention',
                   support_and_coverage='unanimous cell-UMI-position counts, not original Cell Ranger reads',
                   cell_reads='accepted competitive chrM alignments before UMI consensus; duplicate flags unavailable',
                   original_MERCI_mtsnp_nearby_filter='ln=5; singleton support near another variant removed',
                   additional_discovery_depth_alt_or_AF_thresholds=False,
                   seconds=round(time.monotonic() - started, 3),
                   max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    receipt.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
