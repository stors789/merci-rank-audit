"""Functions copied from scripts/rank_identifiability_20261024/bounds.py at 9084e149."""
import numpy as np


def top_k(n, cutoff):
    # R quantile type 7: 1 + (n-1)*p. Verified against the actual R function.
    return (100 + (n - 1) * cutoff + 99) // 100

def top_mask(rank, k):
    return np.asarray(rank) > len(rank) - k

def legal_rank(neff, rank):
    n = len(neff)
    return (np.array_equal(np.sort(rank), np.arange(1, n + 1))
            and np.all(np.diff(neff[np.argsort(rank)]) >= 0))


def global_witness(neff, rna_rank, upper):
    order = np.lexsort((-rna_rank if upper else rna_rank, -neff))
    rank = np.empty(len(neff), dtype=int)
    rank[order] = np.arange(len(neff), 0, -1)
    assert legal_rank(neff, rank)
    return rank

def intersection_bound(neff, rna_rank, k):
    value = np.partition(neff, len(neff)-k)[len(neff)-k]
    block = neff == value
    b = int(block.sum())
    strict = neff > value
    before = int(strict.sum())
    s = k - before
    rna = top_mask(rna_rank, k)
    m = int((block & rna).sum())
    c = int((strict & rna).sum())
    must = rna & (strict | (block & (s == b)))
    possible = rna & (strict | (block & (s > 0)))
    return dict(boundary_Neff=int(value), boundary_b=b, boundary_s=s,
                boundary_RNA_m=m, strict_intersection_c=c,
                lower=c + max(0, s - (b - m)), upper=c + min(s, m),
                must=must, possible=possible, boundary=block, strict=strict,
                rna=rna)
