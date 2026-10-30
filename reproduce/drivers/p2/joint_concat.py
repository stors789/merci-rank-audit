"""Sequentially apply the same sparse outer-axis reindex used by AnnData."""
import numpy as np
import pandas as pd
from scipy import sparse


def align_sparse(matrix, old_ids, union_ids):
    old = pd.Index(old_ids); union = pd.Index(union_ids)
    if not old.is_unique or not union.is_unique:
        raise ValueError('Stable gene IDs must be unique')
    if old.equals(union):
        return matrix
    new_pos = union.get_indexer(old)
    if (new_pos < 0).any():
        raise ValueError('A sample gene ID is absent from the union')
    mapping = sparse.coo_matrix((np.ones(len(old), dtype=bool), (np.arange(len(old)), new_pos)),
                                shape=(len(old), len(union)))
    return matrix @ mapping


def stack_on_disk(chunks, rows, ncols, folder, data_dtype=np.int32, prefix='joint'):
    """Stack CSR chunks without holding a second full copy."""
    total_rows = sum(x['rows'] for x in rows); total_nnz = sum(x['nnz'] for x in rows)
    if total_nnz == 0 or total_nnz >= 2**31:
        raise ValueError('This int32 reconstruction requires 0 < nnz < 2^31')
    data = np.memmap(folder/(prefix+'_data.bin'), mode='w+', dtype=data_dtype, shape=(total_nnz,))
    indices = np.memmap(folder/(prefix+'_indices.bin'), mode='w+', dtype=np.int32, shape=(total_nnz,))
    indptr = np.memmap(folder/(prefix+'_indptr.bin'), mode='w+', dtype=np.int32, shape=(total_rows+1,))
    row_offset = offset = 0
    for path, spec in zip(chunks, rows):
        part = sparse.load_npz(path).tocsr()
        if part.shape != (spec['rows'], ncols) or part.nnz != spec['nnz'] or part.dtype != np.dtype(data_dtype):
            raise ValueError('Aligned sparse chunk differs from its recorded shape')
        stop = offset + part.nnz; end_row = row_offset + part.shape[0]
        data[offset:stop] = part.data; indices[offset:stop] = part.indices
        indptr[row_offset:end_row] = part.indptr[:-1] + offset
        offset = stop; row_offset = end_row
        del part
    indptr[total_rows] = offset
    data.flush(); indices.flush(); indptr.flush()
    return sparse.csr_matrix((data, indices, indptr), shape=(total_rows, ncols), copy=False)
