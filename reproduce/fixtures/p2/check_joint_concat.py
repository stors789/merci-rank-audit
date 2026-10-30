"""Compare streaming gene-axis reindexing against the original AnnData outer join."""
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

repro = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(repro/'drivers/p2'))
from joint_concat import align_sparse, stack_on_disk

samples = [
    ad.AnnData(X=sparse.csr_matrix([[1,2,0],[0,3,4]],dtype=np.int32),
               obs=pd.DataFrame(index=['a','b']), var=pd.DataFrame(index=['z','a','c'])),
    ad.AnnData(X=sparse.csr_matrix([[5,6,7]],dtype=np.int32),
               obs=pd.DataFrame(index=['c']), var=pd.DataFrame(index=['d','c','a'])),
]
original = ad.concat(samples, join='outer', merge=None)
union = pd.Index(sorted(set().union(*(set(x.var_names) for x in samples))))
actual = sparse.vstack([align_sparse(x.X, x.var_names, union) for x in samples], format='csr')
assert list(original.var_names) == list(union)
assert np.array_equal(original.X.data, actual.data)
assert np.array_equal(original.X.indices, actual.indices)
assert np.array_equal(original.X.indptr, actual.indptr)
assert original.obs_names.tolist() == ['a','b','c']
with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp); paths = []; shapes = []
    for i, sample in enumerate(samples):
        part = align_sparse(sample.X, sample.var_names, union)
        path = folder/f'{i}.npz'; sparse.save_npz(path, part, compressed=False)
        paths.append(path); shapes.append({'rows': part.shape[0], 'nnz': part.nnz})
    disk = stack_on_disk(paths, shapes, len(union), folder)
    assert np.array_equal(original.X.data, disk.data)
    assert np.array_equal(original.X.indices, disk.indices)
    assert np.array_equal(original.X.indptr, disk.indptr)
    nuclear = np.array([13,17,19]); identity = np.array([True,False,True,True])
    baseline = (sparse.diags(1e4/nuclear) @ original.X.astype(np.float32)).tocsr()
    baseline.data = np.log1p(baseline.data); baseline = baseline[:,identity].copy()
    paths=[]; shapes=[]
    for start, stop in ((0,2),(2,3)):
        part = (sparse.diags(1e4/nuclear[start:stop]) @ original.X[start:stop].astype(np.float32)).tocsr()
        part.data = np.log1p(part.data);part = part[:,identity].copy()
        path = folder/f'norm_{start}.npz';sparse.save_npz(path,part,compressed=False)
        paths.append(path);shapes.append({'rows':part.shape[0],'nnz':part.nnz})
    actual = stack_on_disk(paths,shapes,int(identity.sum()),folder,data_dtype=np.float64,prefix='norm')
    assert np.array_equal(baseline.data,actual.data) and np.array_equal(baseline.indices,actual.indices)
    assert np.array_equal(baseline.indptr,actual.indptr)
    assert np.array_equal(np.asarray(baseline[[0,2]].mean(axis=0)),np.asarray(actual[[0,2]].mean(axis=0)))
print('PASS: original outer CSR join, chunked normalized identity storage, means and cell order')
