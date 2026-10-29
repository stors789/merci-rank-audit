"""Uniform Neff tie orders conditional on a union of frozen cutoff gates."""
from fractions import Fraction
from bisect import bisect_right
from decimal import Decimal, InvalidOperation
from numbers import Integral
import math
import numpy as np
import pandas as pd

from .bounds import top_k, intersection_bound, legal_rank

CUTOFFS = np.arange(100, 0, -10)
MAX_EXACT_INTEGER = 2**53 - 1


def exact_integers(values, name):
    """Validate original integer values without an intermediate float conversion."""
    if isinstance(values, np.ndarray) and values.ndim != 1:
        raise ValueError(f'{name} must be one-dimensional')
    result = []
    for value in values:
        if isinstance(value, (bool, np.bool_, float, np.floating)):
            raise ValueError(f'{name} requires exact integers, not floating-point values')
        try:
            number = Decimal(int(value)) if isinstance(value, Integral) else Decimal(value)
            if not number.is_finite() or number != number.to_integral_value() or not 0 <= number <= MAX_EXACT_INTEGER:
                raise ValueError
            result.append(int(number))
        except (InvalidOperation, TypeError, ValueError, OverflowError):
            raise ValueError(f'{name} must contain exact integers in 0..{MAX_EXACT_INTEGER}') from None
    return np.asarray(result, dtype=np.int64)


def exact_rank(values, n):
    """Accept an exact integer rank array, never a previously rounded float."""
    rank = values if isinstance(values, np.ndarray) and values.ndim == 1 and values.dtype.kind in 'iu' else exact_integers(values, 'DNA_rank')
    if rank.shape != (n,) or np.any(rank < 1) or np.any(rank > n):
        raise ValueError('DNA_rank must contain integers in 1..N')
    return rank


class GateSampler:
    def __init__(self, ids, neff, rna_rank, maxima, include_events=True):
        ids=np.asarray(ids,dtype=object)
        if ids.ndim!=1:raise ValueError('IDs must be a one-dimensional array')
        if np.any(pd.isna(ids)):raise ValueError('IDs must not contain missing values')
        ids=ids.astype(str)
        neff=exact_integers(neff,'Neff');rna_rank=exact_integers(rna_rank,'RNA_rank');maxima=exact_integers(maxima,'maxima')
        n=len(ids)
        if n<1 or np.any(ids=='') or len(np.unique(ids))!=n or len(neff)!=n or len(rna_rank)!=n:
            raise ValueError('IDs must be nonempty, unique and have the same length as scores/ranks')
        if not np.array_equal(np.sort(rna_rank),np.arange(1,n+1)):
            raise ValueError('RNA ranks must be a unique permutation 1..N')
        if maxima.shape!=(len(CUTOFFS),):
            raise ValueError('All ten matched-N integer maxima are required')
        self.order=np.argsort(ids,kind='stable');self.ids=ids[self.order];self.neff=neff[self.order]
        self.rna=rna_rank[self.order];self.maxima=maxima;self.n=n
        self.ks=np.array([top_k(n,int(c)) for c in CUTOFFS]);self.rna_masks=self.rna[:,None]>n-self.ks
        if np.any(self.maxima<0) or np.any(self.maxima>self.ks):raise ValueError('Maxima must lie within 0..K')
        by_score=np.argsort(self.neff,kind='stable');ordered_scores=self.neff[by_score]
        self.blocks=np.split(by_score,np.flatnonzero(ordered_scores[1:]!=ordered_scores[:-1])+1)
        self.offsets=np.r_[0,np.cumsum([len(b) for b in self.blocks])]
        if not include_events:
            self.events=None;self.active=None;self.S=None;self.logS=None;self.certain=None
            return
        self.events=[]
        for j,k in enumerate(self.ks):
            bound=intersection_bound(self.neff,self.rna,int(k))
            b,s,m,c=(bound[x] for x in ['boundary_b','boundary_s','boundary_RNA_m','strict_intersection_c'])
            lo=max(0,s-(b-m));hi=min(s,m);minimum=max(lo,int(self.maxima[j])-c+1)
            if minimum>hi:
                self.events.append(dict(j=j,k=int(k),bound=bound,logp=-np.inf))
                continue
            support=range(minimum,hi+1)
            weight=math.comb(m,minimum)*math.comb(b-m,s-minimum)
            cumulative=[];total=0
            for x in support:
                total+=weight;cumulative.append(total)
                if x<hi:
                    weight=weight*(m-x)*(s-x)//((x+1)*(b-m-s+x+1))
            probability=Fraction(total,math.comb(b,s))
            logp=math.log(probability.numerator)-math.log(probability.denominator)
            self.events.append(dict(j=j,k=int(k),bound=bound,logp=logp,probability=probability,
                                    minimum=minimum,cumulative=cumulative,total=total,certain=minimum==lo))
        self.active=[e for e in self.events if np.isfinite(e['logp'])]
        self.certain=any(e.get('certain',False) for e in self.active)
        self.S=sum((e['probability'] for e in self.active),Fraction())
        self.logS=math.log(self.S.numerator)-math.log(self.S.denominator) if self.active else -np.inf
        if self.active:
            common=math.lcm(*(e['probability'].denominator for e in self.active))
            weights=[e['probability'].numerator*(common//e['probability'].denominator) for e in self.active]
            self.mixture_cumulative=[];total=0
            for weight in weights:total+=weight;self.mixture_cumulative.append(total)
            self.mixture_total=total

    @staticmethod
    def randbelow(rng,upper):
        bits=upper.bit_length();mask=(1<<bits)-1
        while True:
            value=int.from_bytes(rng.bytes((bits+7)//8),'little')&mask
            if value<upper:return value

    def counts(self,rank):
        rank=exact_rank(rank,self.n)
        return np.sum((rank[:,None]>self.n-self.ks)&self.rna_masks,axis=0)

    def uniform_order(self,rng):
        rank=np.empty(self.n,dtype=int)
        for block,lo,hi in zip(self.blocks,self.offsets[:-1],self.offsets[1:]):
            rank[rng.permutation(block)]=np.arange(lo+1,hi+1)
        return rank

    @staticmethod
    def choose(rng,group,n):
        if n==0:return np.empty(0,dtype=int)
        if n==len(group):return group.copy()
        return rng.choice(group,size=n,replace=False)

    def single_event_order(self,event,rng):
        bnd=event['bound'];s=bnd['boundary_s']
        x=event['minimum']+bisect_right(event['cumulative'],self.randbelow(rng,event['total']))
        rank=np.empty(self.n,dtype=int)
        for block,lo,hi in zip(self.blocks,self.offsets[:-1],self.offsets[1:]):
            if not bnd['boundary'][block[0]]:
                rank[rng.permutation(block)]=np.arange(lo+1,hi+1)
                continue
            inside=block[bnd['rna'][block]];outside=block[~bnd['rna'][block]]
            chosen=np.r_[self.choose(rng,inside,x),self.choose(rng,outside,s-x)]
            unchosen=np.setdiff1d(block,chosen,assume_unique=True)
            rank[rng.permutation(unchosen)]=np.arange(lo+1,hi-s+1)
            rank[rng.permutation(chosen)]=np.arange(hi-s+1,hi+1)
        return rank

    def proposal(self,rng):
        if self.active is None:raise ValueError('Conditional events were not prepared')
        if not self.active:raise ValueError('G is impossible; conditional distribution is undefined')
        if self.certain:
            rank=self.uniform_order(rng);return rank,self.counts(rank),True
        event=self.active[bisect_right(self.mixture_cumulative,self.randbelow(rng,self.mixture_total))]
        rank=self.single_event_order(event,rng);counts=self.counts(rank);multiplicity=int((counts>self.maxima).sum())
        if multiplicity<1 or counts[event['j']]<=self.maxima[event['j']]:raise AssertionError('conditional event failed')
        return rank,counts,bool(rng.integers(multiplicity)==0)

    def validate_order(self,rank):
        return legal_rank(self.neff,exact_rank(rank,self.n))
