"""Reconstruct the original joint identity clusters from newly generated counts."""
import argparse
import gc
import json
from pathlib import Path
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc

from joint_concat import stack_on_disk
from prepare_counts import CONTACT, MARKERS


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--counts', type=Path, required=True); p.add_argument('--gmt', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); started = time.monotonic()
    a.output.mkdir(parents=True, mode=0o700, exist_ok=False); folder = a.output/'normalized_chunks'; folder.mkdir()
    data = ad.read_h5ad(a.counts, backed='r'); obs = data.obs.copy(); var = data.var.copy()
    genes = var.gene_symbol.to_numpy().astype(str); mt = np.char.startswith(genes, 'MT-')
    hallmark = {line.split('\t')[0]: line.rstrip().split('\t')[2:] for line in a.gmt.read_text().splitlines()}
    ox = np.isin(genes, hallmark['HALLMARK_OXIDATIVE_PHOSPHORYLATION']) & ~mt
    identity = ~mt & ~ox & ~np.isin(genes, CONTACT) & (genes != 'ACTA2')
    nuc = obs.nuclear_UMI.to_numpy(); paths = []; shapes = []
    for start in range(0, len(obs), 512):
        end = min(start+512, len(obs)); x = data.X[start:end,:].tocsr()
        normalized = (sparse.diags(1e4/nuc[start:end]) @ x.astype(np.float32)).tocsr()
        normalized.data = np.log1p(normalized.data)
        part = normalized[:,identity].copy(); path = folder/f'{start:06d}.npz'
        sparse.save_npz(path, part, compressed=False); paths.append(path)
        shapes.append({'rows': part.shape[0], 'nnz': part.nnz})
        del x, normalized, part
    data.file.close(); del data; gc.collect()
    matrix = stack_on_disk(paths, shapes, int(identity.sum()), a.output, data_dtype=np.float64, prefix='identity')
    b = ad.AnnData(X=matrix, obs=obs.copy(), var=var.loc[identity].copy())
    identity_genes = genes[identity]
    all_markers = {g for members in MARKERS.values() for g in members}
    marker_mask = np.isin(identity_genes, list(all_markers)); marker_genes = identity_genes[marker_mask]
    markers = b.X[:,marker_mask].copy()
    sc.pp.highly_variable_genes(b, n_top_genes=2000, flavor='seurat')
    reduced = b[:,b.var.highly_variable].copy(); del b, matrix; gc.collect(); b = reduced; del reduced
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=30, svd_solver='arpack', random_state=0)
    sc.pp.neighbors(b, n_neighbors=15, n_pcs=30, random_state=0)
    sc.tl.leiden(b, resolution=.5, random_state=0, flavor='igraph', n_iterations=2)
    obs['identity_cluster'] = b.obs['leiden'].to_numpy()
    clusters = sorted(obs.identity_cluster.unique(), key=int)
    means = np.stack([np.asarray(markers[obs.identity_cluster.eq(c).to_numpy()].mean(axis=0)).ravel() for c in clusters])
    std = means.std(axis=0); z = (means-means.mean(axis=0))/np.where(std>0,std,1)
    scores = np.stack([z[:,np.isin(marker_genes, geneset)].mean(axis=1) for geneset in MARKERS.values()], axis=1)
    order = np.argsort(-scores, axis=1); labels = list(MARKERS)
    mapping = {c: labels[order[i,0]] for i,c in enumerate(clusters)}
    obs['working_type'] = obs.identity_cluster.map(mapping)
    obs.to_csv(a.output/'global_roles.tsv.gz', sep='\t')
    counts = obs.identity_cluster.value_counts().to_dict()
    result = {'status':'COMPLETE','cells':len(obs),'genes':len(genes),'identity_genes':int(identity.sum()),
              'clusters':counts,'mapping':mapping,'HVG':2000,'PCA':30,'neighbors':15,'leiden_resolution':.5,'seed':0,
              'normalization':'original per-cell nuclear UMI; chunked CSR with preserved row/column order',
              'seconds':round(time.monotonic()-started,3)}
    (a.output/'clustering_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main()
