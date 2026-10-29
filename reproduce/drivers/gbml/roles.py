"""Reconstruct the fixed GBML181 RNA-only QC and cluster roles from GEO files."""
import argparse
import csv
import gc
import gzip
import hashlib
import json
import resource
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from threadpoolctl import threadpool_limits


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('matrix', 'features', 'barcodes', 'output'):
        p.add_argument('--' + name, required=True, type=Path)
    a = p.parse_args(); started = time.monotonic()
    a.output.mkdir(parents=True, mode=0o700, exist_ok=False)
    with gzip.open(a.features, 'rt') as f:
        features = [x.rstrip('\n').split('\t') for x in f]
    with gzip.open(a.barcodes, 'rt') as f:
        cells = [x.rstrip('\n') for x in f]
    with threadpool_limits(4), gzip.open(a.matrix, 'rb') as f:
        raw = mmread(f).tocsr()
    if raw.shape != (len(features), len(cells)) or raw.dtype.kind not in 'iu' or len(set(cells)) != len(cells):
        raise ValueError('Original feature/barcode matrix axes or count types differ')
    raw.sum_duplicates()
    ids = np.flatnonzero(np.array([v[2] == 'Gene Expression' for v in features]))
    symbols = np.array([features[i][1] for i in ids], dtype=str)
    mt = np.char.startswith(symbols, 'MT-')
    peaks = sum(v[2] == 'Peaks' for v in features)
    if len(cells) != 4379 or len(ids) != 36601 or peaks != 57101 or len(features) != len(ids) + peaks or mt.sum() != 13:
        raise ValueError('Expected GBML181 gene-expression and peak dimensions')
    x = raw[ids, :].T.tocsr(); del raw; gc.collect()
    total = np.asarray(x.sum(axis=1)).ravel(); ngenes = np.asarray(x.getnnz(axis=1)).ravel()
    mt_counts = np.asarray(x[:, mt].sum(axis=1)).ravel()
    pct = np.divide(100 * mt_counts, total, out=np.full(total.shape, np.nan, dtype=float), where=total > 0)
    eligible = (ngenes >= 500) & (ngenes <= 5000) & (pct < 15)
    names = np.array(cells)[eligible]
    if len(names) < 51:
        raise ValueError('Insufficient cells for fixed 50-component analysis')
    genes = symbols[~mt]; X = x[eligible, :][:, ~mt].tocsr(); del x; gc.collect()
    obj = ad.AnnData(X); obj.obs_names = names; obj.var_names = genes; obj.var_names_make_unique()
    with threadpool_limits(4):
        sc.pp.normalize_total(obj, target_sum=10000); sc.pp.log1p(obj)
        sc.pp.pca(obj, n_comps=50, zero_center=False, svd_solver='randomized', random_state=0)
        sc.pp.neighbors(obj, n_neighbors=15, n_pcs=50, random_state=0)
        sc.tl.leiden(obj, resolution=1, random_state=0, flavor='igraph', n_iterations=2, directed=False)
    labels = obj.obs['leiden'].astype(str).to_numpy()
    membership = a.output / 'role_clusters.tsv.gz'
    pd.DataFrame({'barcode': names, 'cluster': labels}).to_csv(membership, sep='\t', index=False, compression='gzip')
    rows = list(zip(names, labels))
    source = sorted((r for r in rows if r[1] == '13'), key=lambda x: x[0])
    target = sorted((r for r in rows if r[1] in {'0', '1', '2', '3', '4', '6', '8'}), key=lambda x: x[0])
    if len(source) != 108 or len(target) != 2001 or set(x[0] for x in source) & set(x[0] for x in target):
        raise ValueError('Predeclared RNA cluster role counts or disjointness differ')
    roles = a.output / 'fixed_roles.tsv'
    with roles.open('x', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n'); writer.writerow(['barcode', 'cluster', 'role'])
        for role, members in (('source', source), ('target', target)):
            for barcode, cluster in members:
                writer.writerow([barcode, cluster, role])
    receipt = {'status': 'COMPLETE', 'cells': len(cells), 'gene_expression_features': len(ids),
               'atac_peak_features': peaks, 'mt_gene_features': int(mt.sum()),
               'RNA_QC_n': int(eligible.sum()), 'source_n': len(source), 'target_n': len(target),
               'clusters': {str(k): int((labels == k).sum()) for k in sorted(set(labels), key=int)},
               'role_sha256': hashlib.sha256(roles.read_bytes()).hexdigest(),
               'roles': 'myeloid-compatible source and glioma-compatible target, not malignant truth',
               'numeric_threads': 4, 'seconds': round(time.monotonic() - started, 3),
               'max_rss_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (a.output / 'roles_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
