"""Original deterministic cell/UMI/position consensus for accepted chrM reads."""
from collections import Counter
import hashlib
import re
import numpy as np
from scipy import sparse

LENGTH = 16569
BASE = {'A': 1, 'C': 2, 'G': 4, 'T': 8}
BASE_INDEX = {1: 0, 2: 1, 4: 2, 8: 3}


def consensus_matrix(sam, cells):
    cell_index = {cell: i for i, cell in enumerate(cells)}
    masks = {}
    molecule_split = {}
    counts = Counter()
    with sam.open() as source:
        for line in source:
            if line.startswith('@'):
                continue
            fields = line.rstrip('\n').split('\t')
            flag, mapq = int(fields[1]), int(fields[4])
            if flag & (4 | 256 | 2048):
                counts['unmapped_or_secondary_supplementary'] += 1
                continue
            counts['primary_'+fields[2]] += 1
            if fields[2] != 'chrM' or mapq < 30:
                continue
            cigar = [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', fields[5])]
            if any(op == 'N' for _, op in cigar):
                counts['chrM_spliced_alignment_skipped'] += 1
                continue
            if sum(n for n, op in cigar if op in 'M=X') < 60:
                continue
            cell, umi, read_number = fields[0].split('_')
            ci = cell_index[cell]
            ui = int(umi.translate(str.maketrans('ACGT', '0123')), 4)
            molecule = (ci << 24) | ui
            if molecule not in molecule_split:
                molecule_split[molecule] = hashlib.sha256((cell+':'+umi).encode()).digest()[0] & 1
            sequence, quality = fields[9], fields[10]
            q, r = 0, int(fields[3])-1
            counts['accepted_chrM_alignments'] += 1
            for n, op in cigar:
                if op in 'M=X':
                    for offset in range(n):
                        qi, ri = q+offset, r+offset
                        if 5 <= qi < len(sequence)-5 and 0 <= ri < LENGTH and ord(quality[qi])-33 >= 30:
                            bit = BASE.get(sequence[qi])
                            if bit is not None:
                                key = (molecule << 15) | ri
                                masks[key] = masks.get(key, 0) | bit
                    q += n
                    r += n
                elif op in 'IS':
                    q += n
                elif op in 'DN':
                    r += n
            assert q == len(sequence), fields[5]
    row, col = [], []
    for key, mask in masks.items():
        if mask not in BASE_INDEX:
            counts['conflicting_UMI_position_discarded'] += 1
            continue
        molecule, position = key >> 15, key & 32767
        ci = molecule >> 24
        split = molecule_split[molecule]
        row.append(ci)
        col.append((split*4+BASE_INDEX[mask])*LENGTH+position)
    data = np.ones(len(row), dtype=np.int32)
    matrix = sparse.coo_matrix((data, (row, col)), shape=(len(cells), 8*LENGTH)).tocsr()
    return matrix, counts, len(molecule_split), len(row)
