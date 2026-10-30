"""Check 52×2 nuclear scores against analytically known synthetic counts."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

CONTACT = 'RAC1 CDC42 RHOB MYH9 VIM FLNA'.split()
HUMAN = 'MKI67 PCNA MCM2 MCM5 TYMS TOP2A CCNB1 CDK1 UBE2C STMN1'.split()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, mode=0o700, exist_ok=False)
    n = 3428
    barcodes = [f'cell{i:04d}' for i in range(n)]
    ids = ['GSM7774400:' + cell for cell in barcodes]
    symbols = ['GENE1', 'GENE2'] + CONTACT + HUMAN
    values = np.ones((n, len(symbols)), dtype=np.int32)
    values[:, 0] = np.arange(1, n + 1)
    values[:, 1] = 10000 - values[:, 0] - len(CONTACT) - len(HUMAN)
    counts = ad.AnnData(X=sparse.csr_matrix(values), obs=pd.DataFrame(index=ids),
                        var=pd.DataFrame({'gene_symbol': symbols}, index=[f'gene{i:02d}' for i in range(len(symbols))]))
    counts.write_h5ad(root / 'counts.h5ad')
    meta = pd.DataFrame({'barcode': barcodes, 'cell_id': ids, 'working_type': 'hepatic', 'identity_cluster': 24,
                         'nuclear_UMI': 10000, 'UMI': 10000,
                         'full_library_accepted_MT_reads': np.arange(n) + 2000,
                         'contact_transport': np.log1p(6)})
    meta.to_csv(root / 'metadata.tsv', sep='\t', index=False)
    score = pd.DataFrame({'barcode': barcodes[:3265], 'source_marker_observable': np.arange(3265) < 2741,
                          'observed_source_feature_sites': 1 + np.arange(3265) % 2,
                          'raw_prediction': np.where(np.arange(3265) % 2 == 0, 'Receiver', 'non-Receiver')})
    score.to_csv(root / 'cells.tsv', sep='\t', index=False)
    (root / 'synthetic.gmt').write_text(''.join(f'HALLMARK_FAKE_{i:02d}\tfixture\tGENE1\n' for i in range(50)))
    driver = Path(__file__).resolve().parents[2] / 'drivers/p2/nuclear.py'
    output = root / 'nuclear'
    command = [sys.executable, str(driver), '--counts', str(root / 'counts.h5ad'), '--metadata', str(root / 'metadata.tsv'),
               '--cells', str(root / 'cells.tsv'), '--gmt', str(root / 'synthetic.gmt'), '--output', str(output)]
    result = subprocess.run(command, capture_output=True, text=True)
    (root / 'nuclear.log').write_text(result.stdout + result.stderr)
    if result.returncode:
        raise AssertionError('Synthetic nuclear calculation failed; inspect nuclear.log')
    receipt = json.loads((output / 'nuclear_receipt.json').read_text())
    if receipt['nuclear_features'] != 52 or receipt['representations'] != 2 or receipt['observable_domain'] != 2741:
        raise AssertionError('Nuclear axes differ')
    with np.load(output / 'nuclear.npz', allow_pickle=False) as data:
        expected = np.arange(2741)
        if not data['original'].flags.f_contiguous or not data['alternative'].flags.f_contiguous:
            raise AssertionError('Original column-major nuclear layout was not preserved')
        if not (np.allclose(data['original'][:, 0], np.log1p(expected + 1), rtol=0, atol=1e-14) and
                np.allclose(data['alternative'][:, 0], expected / (n - 1), rtol=0, atol=1e-14) and
                np.allclose(data['original'][:, 50], np.log1p(6), rtol=0, atol=1e-14) and
                np.allclose(data['alternative'][:, 50], .5, rtol=0, atol=1e-14) and
                np.allclose(data['original'][:, 51], .5, rtol=0, atol=1e-14) and
                np.allclose(data['alternative'][:, 51], np.log1p(1), rtol=0, atol=1e-14)):
            raise AssertionError('Analytical nuclear readout or full-domain tied percentile differs')
        constant = np.full(len(expected), np.log1p(10000))
        retained = bool(constant.std() > 0)
        expected_x = np.column_stack([np.ones(len(expected)), np.log1p(expected + 2000)] +
                                     ([constant] if retained else []) + [1 + expected % 2])
        if data['X'].shape != expected_x.shape or not np.array_equal(data['X'], expected_x):
            raise AssertionError('Synthetic covariate ordering differs from the original std>0 rule')
    print('PASS: synthetic 52×2 readouts, full-domain ranks, tied scores and original covariate rule; '
          f'floating constant retained={retained}')


if __name__ == '__main__':
    main()
