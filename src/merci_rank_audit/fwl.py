"""Functions copied from scripts/merci_nuclear_robustness_20261026/regression.py at 9084e149."""
import numpy as np


def projection(x,y):
    u,s,_=np.linalg.svd(x,full_matrices=False)
    rank=int((s>s[0]*max(x.shape)*np.finfo(float).eps).sum())
    q=u[:,:rank]
    return q,y-q@(q.T@y)

def beta_batch(y,x,labels,block=256):
    q,yr=projection(x,y);out=np.full((len(labels),y.shape[1]),np.nan)
    for start in range(0,len(labels),block):
        z=labels[start:start+block].astype(float);zr=z-(z@q)@q.T;den=(zr*zr).sum(axis=1)
        valid=(z.sum(axis=1)>=5)&((len(x)-z.sum(axis=1))>=5)&(den>1e-10)
        values=np.full((len(z),y.shape[1]),np.nan);values[valid]=(zr[valid]@yr)/den[valid,None]
        out[start:start+len(z)]=values
    return out
