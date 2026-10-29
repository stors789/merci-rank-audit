"""Prepare MT13 expression using only the original gene-expression library denominator."""
import argparse
import gzip
import json
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

GENES = sorted(['MT-ATP6','MT-ATP8','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ND1',
                'MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('features', 'barcodes', 'matrix', 'roles', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    folder = args.output.resolve(); receipt = folder / 'RNA_input_receipt.json'
    if folder.exists():
        raise FileExistsError('Protect first RNA input')
    with gzip.open(args.features, 'rt') as handle:
        features = [x.rstrip('\n').split('\t') for x in handle]
    with gzip.open(args.barcodes, 'rt') as handle:
        cells = [x.strip() for x in handle]
    roles = pd.read_csv(args.roles, sep='\t', dtype=str)
    indices = pd.Index(cells).get_indexer(roles.barcode)
    if (indices < 0).any() or len(set(cells)) != len(cells):
        raise ValueError('Frozen roles not on the original expression axis')
    gene_rows = [i for i, feature in enumerate(features) if feature[2] == 'Gene Expression']
    symbols = np.array([features[i][1] for i in gene_rows])
    if len(gene_rows) != 36601 or any(np.sum(symbols == x) != 1 for x in GENES):
        raise ValueError('Gene Expression axis or MT13 genes changed')
    with gzip.open(args.matrix, 'rb') as handle:
        matrix = mmread(handle).tocsr()
    matrix.sum_duplicates()
    if matrix.shape != (93702, 4379) or matrix.dtype.kind not in 'iu':
        raise ValueError('Original multimodal matrix shape or counts changed')
    gene_matrix = matrix[gene_rows, :][:, indices]
    total = np.asarray(gene_matrix.sum(axis=0)).ravel()
    counts = gene_matrix[[int(np.flatnonzero(symbols == x)[0]) for x in GENES], :].toarray()
    if np.any(total <= 0) or np.any(counts < 0):
        raise ValueError('Invalid RNA library denominator or MT counts')
    normalized = np.log1p(10000 * counts / total[np.newaxis, :])
    names = roles.barcode.str.removesuffix('-1').tolist()
    if len(set(names)) != len(names) or not np.isfinite(normalized).all():
        raise ValueError('RNA names or normalized values invalid')
    folder.mkdir()
    pd.DataFrame(normalized,index=GENES,columns=names).to_csv(folder / 'MT_expression.tsv',sep='\t',index_label='gene')
    pd.DataFrame(counts,index=GENES,columns=names).to_csv(folder / 'MT_counts.tsv',sep='\t',index_label='gene')
    pd.DataFrame({'cell':names,'RNA_gene_total':total}).to_csv(folder / 'RNA_totals.tsv',sep='\t',index=False)
    group_counts = {}
    for role in ['source','target']:
        selected = roles.role.to_numpy() == role
        group_counts[role] = {'frozen':int(selected.sum()),
                              'zero_MT13_expression':int(np.sum(counts[:, selected].sum(axis=0)==0))}
    result = dict(complete=True, counts_are_provided_matrix=True, gene_expression_rows=len(gene_rows),
                  excluded_ATAC_peak_rows=sum(x[2]=='Peaks' for x in features), MT_genes=GENES,
                  normalization='log1p(10000 * MT gene count / all Gene Expression UMI counts)',
                  roles=group_counts, zero_expression_kept_in_frozen_roles=True,
                  seconds=round(time.monotonic()-started,3),
                  max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    receipt.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main()
