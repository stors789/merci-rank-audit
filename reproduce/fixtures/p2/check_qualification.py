"""Exercise discovery, evaluation and the strict read gate on synthetic UMI counts."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('upstream', 'r-library', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--rscript', default='Rscript')
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, mode=0o700, exist_ok=False)
    consensus = root / 'consensus'
    (consensus / 'groups').mkdir(parents=True)
    alignment = root / 'alignment'
    alignment.mkdir()
    donors = [f'D{i:02d}' for i in range(20)]
    targets = [f'R{i:02d}' for i in range(40)]
    cells = donors + targets + ['Z_excluded']
    roles = pd.DataFrame({'barcode': cells, 'cell_id': ['GSM7774400:' + x for x in cells], 'UMI': [91] * len(cells),
                          'working_type': ['myeloid'] * 20 + ['hepatic'] * 40 + ['myeloid']})
    roles.to_csv(root / 'roles.tsv', sep='\t', index=False)
    reads = {cell: 1001 for cell in cells}
    reads['Z_excluded'] = 1000
    (alignment / 'accepted_MT_reads_by_cell.json').write_text(json.dumps(reads))
    (alignment / 'alignment_receipt.json').write_text(json.dumps({'accepted_MT_reads': sum(reads.values())}))
    (consensus / 'allele_cell_barcodes.json').write_text(json.dumps(cells))
    (root / 'chrM.fa').write_text('>chrM\n' + 'A' * 16569 + '\n')
    mutated = set(donors[:19] + targets[:1])
    total = 0
    for number, start in enumerate(range(0, len(cells), 16)):
        group = cells[start:start + 16]
        row, col, data = [], [], []
        for i, cell in enumerate(group):
            values = [(0, 2)] + ([] if cell == 'R39' else [(4 * 16569, 10)])
            for column, value in values:
                row.append(i); col.append(column); data.append(value)
            if cell in mutated:
                for column, value in [(16569, 1), (5 * 16569, 5)]:
                    row.append(i); col.append(column); data.append(value)
        matrix = sparse.coo_matrix((np.array(data, dtype=np.int32), (row, col)),
                                   shape=(len(group), 8 * 16569)).tocsr()
        total += int(matrix.sum())
        sparse.save_npz(consensus / 'groups' / f'group_{number:04d}.npz', matrix)
    (consensus / 'consensus_receipt.json').write_text(json.dumps(
        {'status': 'COMPLETE', 'matrix_shape': [len(cells), 8 * 16569], 'UMI_positions': total}))
    output = root / 'qualification'
    driver = Path(__file__).resolve().parents[2] / 'drivers/p2/qualify.py'
    command = [sys.executable, str(driver), '--roles', str(root / 'roles.tsv'), '--consensus', str(consensus),
               '--alignment', str(alignment), '--reference', str(root / 'chrM.fa'),
               '--upstream', str(args.upstream.resolve()), '--r-library', str(args.r_library.resolve()),
               '--rscript', args.rscript, '--output', str(output)]
    result = subprocess.run(command, capture_output=True, text=True)
    (root / 'qualification.log').write_text(result.stdout + result.stderr)
    if result.returncode:
        raise AssertionError('Synthetic qualification failed; inspect qualification.log')
    receipt = json.loads((output / 'qualification_receipt.json').read_text())
    if (receipt['donor_cells'] != 20 or receipt['target_cells'] != 40 or receipt['discovered_candidates'] != 1 or
            receipt['test_eligible_variants'] != 1 or receipt['feature_ids'] != ['MT_1_A-C']):
        raise AssertionError('Synthetic discovery or source qualification differs')
    with np.load(output / 'qualification_inputs.npz', allow_pickle=False) as data:
        if (data['ALT'].shape != (1, 61) or data['depth'][0, 0] != 15 or data['ALT'][0, 0] != 5 or
                data['depth'][0, 21] != 10 or data['ALT'][0, 21] != 0):
            raise AssertionError('Synthetic evaluation-half counts differ')
    if 'Z_excluded' in (output / 'donors.txt').read_text().splitlines():
        raise AssertionError('The strict >1000 read boundary was not enforced')
    genes = sorted(['MT-ATP6', 'MT-ATP8', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-CYB', 'MT-ND1',
                    'MT-ND2', 'MT-ND3', 'MT-ND4', 'MT-ND4L', 'MT-ND5', 'MT-ND6'])
    matrix = sparse.csr_matrix(np.tile(np.arange(1, 14, dtype=np.int32), (len(cells), 1)))
    counts = ad.AnnData(X=matrix, obs=pd.DataFrame(index=roles.cell_id),
                        var=pd.DataFrame({'gene_symbol': genes}, index=[f'gene{i:02d}' for i in range(13)]))
    counts.write_h5ad(root / 'counts.h5ad')
    prepared = root / 'prepared'
    prepare = driver.parent / 'prepare_merci.py'
    result = subprocess.run([sys.executable, str(prepare), '--qualification', str(output),
                             '--counts', str(root / 'counts.h5ad'), '--output', str(prepared)],
                            capture_output=True, text=True)
    (root / 'prepare.log').write_text(result.stdout + result.stderr)
    if result.returncode:
        raise AssertionError('Synthetic MERCI input preparation failed; inspect prepare.log')
    vaf = pd.read_csv(prepared / 'VAF.tsv', sep='\t', index_col=0)
    expression = pd.read_csv(prepared / 'MT_expression.tsv', sep='\t', index_col=0)
    if (vaf.shape != (1, 60) or not np.isnan(vaf.loc['MT_1_A-C', 'R39']) or
            not np.isclose(vaf.loc['MT_1_A-C', 'D00'], 1/3) or
            not np.isclose(expression.loc[genes[0], 'D00'], np.log1p(10000/91))):
        raise AssertionError('Prepared missingness or original-GEO normalization differs')
    print('PASS: discovery, strict >1000 gate, evaluation missingness and GEO-normalized MERCI inputs')


if __name__ == '__main__':
    main()
