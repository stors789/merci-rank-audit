"""Independent finite hypergeometric reference for zero, certain, or single-gate cases."""
import math
from fractions import Fraction
import numpy as np


def k_for(n,percent):
    return (100+(n-1)*int(percent)+99)//100


def event(n,neff,rna,maximum,percent,sorted_scores):
    k=k_for(n,percent);value=sorted_scores[n-k];boundary=neff==value;strict=neff>value
    b=int(boundary.sum());s=k-int(strict.sum());rna_top=rna>n-k
    m=int((boundary&rna_top).sum());c=int((strict&rna_top).sum())
    lo=max(0,s-(b-m));hi=min(s,m);minimum=max(lo,int(maximum)-c+1)
    result=dict(cutoff=int(percent),k=k,value=int(value),b=b,s=s,m=m,c=c,lo=lo,hi=hi,minimum=minimum)
    if minimum>hi:return dict(result,status='impossible',probability=Fraction())
    if minimum<=lo:return dict(result,status='certain',probability=Fraction(1))
    weight=math.comb(m,minimum)*math.comb(b-m,s-minimum);total=weighted=0
    for x in range(minimum,hi+1):
        total+=weight;weighted+=x*weight
        if x<hi:weight=weight*(m-x)*(s-x)//((x+1)*(b-m-s+x+1))
    return dict(result,status='variable',probability=Fraction(total,math.comb(b,s)),conditional_successes=Fraction(weighted,total))


def analytic(neff,rna,maxima,cutoffs=tuple(range(100,0,-10))):
    """Return raw-possible mask and independent exact reference when tractable.

    For one nontrivial gate, the event is a hypergeometric tail in one score
    block.  Conditional within each gate-RNA class, positions are exchangeable.
    This also handles the gate and 50% cutoff sharing only part of that block.
    """
    neff=np.asarray(neff,dtype=np.int64);rna=np.asarray(rna,dtype=np.int64);n=len(neff)
    if n<1 or len(rna)!=n or len(maxima)!=len(cutoffs):raise ValueError('Invalid dimensions')
    sorted_scores=np.sort(neff);k50=k_for(n,50);value50=sorted_scores[n-k50]
    block50=neff==value50;strict50=neff>value50;b50=int(block50.sum());s50=k50-int(strict50.sum())
    rna50=rna>n-k50;possible=rna50&(strict50|block50);q=np.zeros(n,dtype=float)
    m50=int((block50&rna50).sum());c50=int((strict50&rna50).sum())
    raw_range=(c50+max(0,s50-(b50-m50)),c50+min(s50,m50))
    q[rna50&strict50]=1.;q[rna50&block50]=s50/b50
    events=[event(n,neff,rna,maximum,pct,sorted_scores) for pct,maximum in zip(cutoffs,maxima)]
    if any(e['status']=='certain' for e in events):
        return dict(kind='exact_certain',gate_probability=Fraction(1),members=q,raw_possible=possible,raw_count_range=raw_range,events=events)
    active=[e for e in events if e['status']=='variable']
    if not active:return dict(kind='exact_impossible',gate_probability=Fraction(),members=None,raw_possible=possible,raw_count_range=raw_range,events=events)
    if len(active)>1:return dict(kind='multi_gate_reference_required',gate_probability=None,members=None,raw_possible=possible,raw_count_range=raw_range,events=events)
    e=active[0]
    if e['value']==int(value50):
        expected=e['conditional_successes'];b=e['b'];m=e['m'];sg=e['s']
        for gate_class in (False,True):
            group=rna50&block50&((rna>n-e['k'])==gate_class)
            if not group.any():continue
            gate_probability=expected/m if gate_class else (sg-expected)/(b-m)
            if s50<=sg:conditional=gate_probability*Fraction(s50,sg)
            else:conditional=gate_probability+(1-gate_probability)*Fraction(s50-sg,b-sg)
            q[group]=float(conditional)
    return dict(kind='exact_single_gate',gate_probability=e['probability'],members=q,raw_possible=possible,raw_count_range=raw_range,events=events)
