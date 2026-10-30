"""Compare chunked half-count calculations with direct whole-matrix formulas."""
from pathlib import Path
import sys

import numpy as np
from scipy import sparse

repro=Path(__file__).resolve().parents[2];sys.path.insert(0,str(repro/'drivers/p2'))
from half_counts import summarize,evaluation

length=5;dense=np.zeros((4,8*length),dtype=np.int32)
for row,col,value in [(0,0,3),(0,5,2),(0,20,4),(0,25,1),(1,2,5),(1,22,6),(1,27,2),
                      (2,4,1),(2,39,5),(3,0,7),(3,20,2),(3,25,3),(3,39,1)]:dense[row,col]=value
full=sparse.csr_matrix(dense);parts=[full[:2],full[2:]]
direct=summarize(full,length);piece=[summarize(x,length) for x in parts]
assert np.array_equal(direct[0],sum(x[0] for x in piece))
assert np.array_equal(direct[1],sum(x[1] for x in piece))
assert np.array_equal(direct[2],np.concatenate([x[2] for x in piece]))
assert np.array_equal(direct[3],np.concatenate([x[3] for x in piece]))
positions=np.array([0,2,4]);bases=np.array([1,0,3])
whole=evaluation(full,positions,bases,length);small=[evaluation(x,positions,bases,length) for x in parts]
assert np.array_equal(whole[0],np.concatenate([x[0] for x in small],axis=1))
assert np.array_equal(whole[1],np.concatenate([x[1] for x in small],axis=1))
af=np.full(whole[0].shape,np.nan);np.divide(whole[1],whole[0],out=af,where=whole[0]>=5)
assert np.isnan(af[0,1]) and af[0,0]==.2 and af[0,3]==.6
print('PASS: chunked discovery counts, observed depth, evaluation ALT/depth and missingness')
