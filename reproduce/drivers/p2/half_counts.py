"""Small-block summaries of the original discovery/evaluation UMI matrix."""
import numpy as np


def summarize(part, length):
    aggregate = np.column_stack([part[:,b*length:(b+1)*length].sum(0).A1 for b in range(8)])
    mutated = np.column_stack([(part[:,b*length:(b+1)*length]>0).sum(0).A1 for b in range(8)])
    depth = sum(part[:,b*length:(b+1)*length] for b in range(8)).tocsr()
    return aggregate, mutated, np.asarray(part.sum(1)).ravel(), np.asarray((depth>=5).sum(1)).ravel()


def evaluation(part, positions, base_indices, length):
    depth = sum(part[:,(4+b)*length:(5+b)*length] for b in range(4)).tocsr()
    columns = (4+np.asarray(base_indices))*length+np.asarray(positions)
    return depth[:,positions].toarray().T, part[:,columns].toarray().T
