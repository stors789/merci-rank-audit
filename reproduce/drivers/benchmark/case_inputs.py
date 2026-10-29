"""Load explicit prepared NPZ cases and verify saved array identities."""
import hashlib
from pathlib import Path
import numpy as np

def sha_file(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def digest_arrays(ids,neff,rna,maxima):
    h=hashlib.sha256()
    for value in ids:
        raw=str(value).encode();h.update(len(raw).to_bytes(4,"little"));h.update(raw)
    for values in (neff,rna,maxima):
        a=np.asarray(values,dtype="<i8");h.update(len(a).to_bytes(8,"little"));h.update(a.tobytes())
    return h.hexdigest()
def resolve_case(case,base):
    case=dict(case)
    for name in ("prepared_input","score","thresholds"):
        if name in case:case[name]=str((base/case[name]).resolve())
    return case
def load_prepared(case):
    path=Path(case["prepared_input"])
    if sha_file(path)!=case["prepared_sha256"]:raise ValueError("Prepared input identity changed: "+case["case"])
    if case["scope"]=="real":
        if sha_file(case["score"])!=case["score_sha256"] or sha_file(case["thresholds"])!=case["thresholds_sha256"]:
            raise ValueError("Original score or thresholds identity changed: "+case["case"])
    with np.load(path,allow_pickle=False) as d:
        ids=d["ids"];neff=d["neff"];rna=d["rna"];maxima=d["maxima"]
    if digest_arrays(ids,neff,rna,maxima)!=case["canonical_array_sha256"]:raise ValueError("Canonical array identity changed")
    return ids,neff,rna,maxima
