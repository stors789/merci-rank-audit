"""Joint nuclear identity representation of the five paired HCC specimens."""
import argparse
import gc
import json
import tarfile
import re
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse
import scanpy as sc
from joint_concat import align_sparse, stack_on_disk

CACHE = OUT = GMT = None
CONTACT = ['RAC1', 'CDC42', 'RHOB', 'MYH9', 'VIM', 'FLNA']
MARKERS = {
    'hepatic': ['ALB', 'APOA1', 'APOA2', 'ASGR1', 'CPS1', 'TF', 'KRT18'],
    'epithelial': ['EPCAM', 'KRT19', 'KRT7', 'KRT8', 'KRT18', 'MUC1'],
    'fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'PDGFRA', 'COL3A1'],
    'mural': ['RGS5', 'MCAM', 'NOTCH3', 'CSPG4', 'PDGFRB', 'DES', 'MYH11'],
    'endothelial': ['PECAM1', 'VWF', 'KDR', 'EMCN', 'ENG', 'RAMP2', 'PLVAP'],
    'myeloid': ['LST1', 'TYROBP', 'FCER1G', 'AIF1', 'CTSS', 'CD68'],
    'T_NK': ['CD3D', 'CD3E', 'TRAC', 'NKG7', 'KLRD1', 'CD247'],
    'B': ['MS4A1', 'CD79A', 'CD79B', 'CD22', 'BANK1'],
    'plasma': ['MZB1', 'JCHAIN', 'SDC1', 'DERL3', 'IGHG1'],
    'mast': ['KIT', 'TPSAB1', 'TPSB2', 'CPA3'],
}


def resolve_symbols(samples):
    features = []
    for sample in samples:
        f = pd.read_csv(next((CACHE / 'counts' / sample['GSM']).rglob('genes.tsv')),
                        sep='\t', header=None, dtype=str, keep_default_na=False)
        features.append(f.rename(columns={0: 'gene_ID', 1: 'raw_symbol'}))
    table = pd.concat(features, ignore_index=True)
    symbols, resolutions = {}, []
    for gene_id, sub in table.groupby('gene_ID', sort=False):
        names = sorted(sub['raw_symbol'].unique(), key=lambda name: (len(name), name))
        base = names[0]
        assert all(n == base or re.fullmatch(re.escape(base) + r'_\d+', n) for n in names), (gene_id, names)
        symbols[gene_id] = base
        if len(names) > 1:
            resolutions.append(dict(gene_ID=gene_id, raw_symbols=';'.join(names), resolved_symbol=base,
                                    rule='same stable ID; other names equal base plus underscore numeric suffix'))
    pd.DataFrame(resolutions).to_csv(OUT / 'annotation_symbol_resolution.csv', index=False)
    return symbols


def read_counts(sample, symbol_by_id):
    folder = CACHE / 'counts' / sample['GSM']
    matrix = next(folder.rglob('matrix.mtx'))
    feature_path = next(folder.rglob('genes.tsv'))
    features = pd.read_csv(feature_path, sep='\t', header=None, dtype=str, keep_default_na=False)
    features[1] = features[0].map(symbol_by_id)
    assert features[1].notna().all()
    barcodes = pd.read_csv(next(folder.rglob('barcodes.tsv')), sep='\t', header=None, dtype=str)[0]
    assert features[0].is_unique and barcodes.is_unique
    X = io.mmread(matrix).T.tocsr().astype(np.int32)
    assert X.shape == (len(barcodes), len(features))
    genes = features[1].to_numpy()
    mt = np.char.startswith(genes.astype(str), 'MT-')
    total = np.asarray(X.sum(axis=1)).ravel()
    mitochondrial = np.asarray(X[:, mt].sum(axis=1)).ravel()
    detected = X.getnnz(axis=1)
    keep = (detected >= 300) & (total >= 1000) & (mitochondrial / np.maximum(total, 1) <= .2)
    obs = pd.DataFrame({'barcode': barcodes, 'GSM': sample['GSM'], 'pair': sample['pair'],
                        'tissue': sample['tissue'], 'UMI': total, 'nuclear_UMI': total - mitochondrial,
                        'MT_UMI': mitochondrial, 'genes': detected})
    obs.index = sample['GSM'] + ':' + obs['barcode']
    obs.index.name = 'cell_id'
    a = ad.AnnData(X=X[keep], obs=obs.loc[keep],
                   var=pd.DataFrame({'gene_symbol': genes}, index=features[0].to_numpy()))
    return a, dict(GSM=sample['GSM'], pair=sample['pair'], tissue=sample['tissue'],
                   raw_cells=X.shape[0], raw_genes=X.shape[1], passing_cells=int(keep.sum()),
                   nonzero_entries=int(X.nnz), raw_UMI=int(total.sum()))


def main():
    global CACHE, OUT, GMT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archives', type=Path, required=True)
    parser.add_argument('--gmt', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    CACHE = args.output.resolve(); OUT = CACHE / 'results'; GMT = args.gmt.resolve()
    if CACHE.exists():
        raise FileExistsError('Protect existing output')
    samples = json.loads((Path(__file__).resolve().parent / 'sample_manifest.json').read_text())
    if len(samples) != 10 or not GMT.is_file():
        raise ValueError('Expected all ten public samples and the explicit Hallmark GMT')
    for sample in samples:
        if not (args.archives / sample['file']).is_file():
            raise FileNotFoundError(args.archives / sample['file'])
    OUT.mkdir(parents=True, exist_ok=False)
    for sample in samples:
        target = CACHE / 'counts' / sample['GSM']; target.mkdir(parents=True)
        with tarfile.open(args.archives / sample['file'], mode='r:*') as archive:
            archive.extractall(target, filter='data')
        for name in ('genes.tsv', 'matrix.mtx', 'barcodes.tsv'):
            if len(list(target.rglob(name))) != 1:
                raise ValueError('Expected one ' + name + ' in ' + sample['GSM'])

    symbol_by_id = resolve_symbols(samples)
    union_ids = pd.Index(sorted(symbol_by_id))
    chunks, shapes, observations, summaries = [], [], [], []
    chunk_dir = CACHE / 'aligned_chunks'; chunk_dir.mkdir()
    for sample in samples:
        a, summary = read_counts(sample, symbol_by_id)
        matrix = align_sparse(a.X, a.var_names, union_ids)
        path = chunk_dir / (sample['GSM'] + '.npz'); sparse.save_npz(path, matrix, compressed=False)
        chunks.append(path); shapes.append({'rows': matrix.shape[0], 'nnz': matrix.nnz})
        observations.append(a.obs.copy())
        summaries.append(summary)
        print(summary, flush=True)
        del a, matrix; gc.collect()
    X = stack_on_disk(chunks, shapes, len(union_ids), CACHE)
    a = ad.AnnData(X=X, obs=pd.concat(observations),
                   var=pd.DataFrame({'gene_symbol': [symbol_by_id[x] for x in union_ids]}, index=union_ids))
    del observations, X
    assert a.obs_names.is_unique and a.var_names.is_unique
    genes = a.var['gene_symbol'].to_numpy().astype(str)
    mt = np.char.startswith(genes, 'MT-')
    hallmark = {line.split('\t')[0]: line.rstrip().split('\t')[2:] for line in GMT.read_text().splitlines()}
    ox = np.isin(genes, hallmark['HALLMARK_OXIDATIVE_PHOSPHORYLATION']) & ~mt
    contact = np.isin(genes, CONTACT)
    identity = ~mt & ~ox & ~contact & (genes != 'ACTA2')
    X = a.X
    for name, mask in [('nuclear_OXPHOS', ox), ('contact_transport', contact)]:
        counts = np.asarray(X[:, mask].sum(axis=1)).ravel()
        a.obs[name + '_UMI'] = counts
        a.obs[name] = np.log1p(1e4 * counts / a.obs['nuclear_UMI'].to_numpy())
    pd.DataFrame(summaries).to_csv(OUT / 'input_QC_summary.csv', index=False)
    a.write_h5ad(CACHE / 'raw_HCC_counts.h5ad', compression='gzip')
    scope = dict(status='COUNTS_READY', identity_seed=0, HVG=2000, PCA=30, neighbors=15, leiden_resolution=.5,
                 identity_excludes=['MT','nuclear_OXPHOS','six_contact_genes','ACTA2_observation'],
                 contact_symbols=CONTACT, nuclear_OXPHOS_symbols=genes[ox].tolist(),
                 passing_cells=a.n_obs, unique_gene_IDs=a.n_vars,
                 duplicate_symbols=int(pd.Series(genes).duplicated().sum()), patient_pairs=5,
                 symbol_resolution='same stable ID; numeric suffix variants reconciled; no rows merged',
                 clinical_mapping='not used in identity or scoring')
    (OUT / 'preparation_scope.json').write_text(json.dumps(scope,indent=2)+'\n')
    print('COUNTS_READY', a.shape, flush=True)



if __name__ == '__main__':
    main()
