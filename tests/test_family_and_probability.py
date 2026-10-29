import itertools
import math
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from merci_rank_audit.analytic import analytic
from merci_rank_audit.conditioned import evaluate, solve
from merci_rank_audit.core import exact_number
from merci_rank_audit.sampling import sample_case
from merci_rank_audit.sampler import GateSampler

FIXTURES=Path(__file__).resolve().parents[1]/'src/merci_rank_audit/examples'


@pytest.fixture(scope='module')
def all_eight_orders():
    ranks=np.asarray(list(itertools.permutations(range(1,9))),dtype=int)
    rna=np.arange(1,9);counts=[]
    for cutoff in range(100,0,-10):
        k=math.ceil(1+7*cutoff/100)
        counts.append(((ranks>8-k)&(rna>8-k)).sum(axis=1))
    k=math.ceil(1+7*.5)
    return ranks,np.column_stack(counts),(ranks>8-k)&(rna>8-k)


@pytest.mark.parametrize('case,kind',[
    ('all_tied_certain','exact_certain'),
    ('all_tied_single','exact_single_gate'),
    ('all_tied_overlapping','multi_gate_reference_required')])
def test_exact_families_probabilities_and_witnesses_against_enumeration(all_eight_orders,case,kind):
    scores=pd.read_csv(FIXTURES/(case+'_score.tsv'),sep='\t')
    thresholds=pd.read_csv(FIXTURES/(case+'_thresholds.csv'))
    s=GateSampler(scores.cell.to_numpy(str),scores.Neff.to_numpy(int),scores.RNA_rank.to_numpy(int),thresholds.maximum.to_numpy(int))
    ranks,counts,labels=all_eight_orders
    passed=(counts>s.maxima).any(axis=1)
    gated=labels&passed[:,None]
    ctx,records=solve(s);answer=analytic(s.neff,s.rna,s.maxima)
    assert answer['kind']==kind
    assert np.array_equal([r['raw_core'] for r in records],labels.all(axis=0))
    assert np.array_equal([r['raw_union'] for r in records],labels.any(axis=0))
    assert np.array_equal([r['conditioned_core'] for r in records],labels[passed].all(axis=0))
    assert np.array_equal([r['conditioned_union'] for r in records],labels[passed].any(axis=0))
    assert np.array_equal(gated.all(axis=0),labels.all(axis=0) if passed.all() else np.zeros(8,dtype=bool))
    if answer['members'] is not None:
        assert answer['gate_probability']==Fraction(int(passed.sum()),len(passed))
        for i,probability in enumerate(answer['members']):
            assert probability==Fraction(int(labels[passed,i].sum()),int(passed.sum()))
    if case=='all_tied_single':assert answer['gate_probability']==Fraction(1,28)
    for i in range(8):
        for selected in [False,True]:
            expected=bool(np.any(labels[passed,i]==selected))
            w=evaluate(ctx,i,selected)
            assert w['reachable']==expected
            if expected:
                rank=w['rank'];assert s.validate_order(rank)
                actual_counts=np.array([int(np.sum((rank>8-math.ceil(1+7*c/100))&
                    (s.rna>8-math.ceil(1+7*c/100)))) for c in range(100,0,-10)])
                assert (actual_counts>s.maxima).any()
                assert bool((rank[i]>8-5) and (s.rna[i]>8-5))==selected


def test_impossible_family_is_undefined_and_original_rank_optional():
    n=5;k=[math.ceil(1+(n-1)*c/100) for c in range(100,0,-10)]
    s=GateSampler([f'id{i}' for i in range(n)],[0]*n,list(range(1,n+1)),k)
    context,rows=solve(s)
    assert not context['possible']
    assert all(x['conditioned_core'] is None and x['conditioned_union'] is None for x in rows)
    assert evaluate(context,0,True)['status']=='undefined_G_impossible'
    assert analytic(s.neff,s.rna,s.maxima)['members'] is None


@pytest.mark.parametrize('scores,rna,cutoff',[
    ([0]*8,list(range(1,9)),70),
    ([0,0,0,1,1,1],[1,2,3,4,5,6],10)])
def test_single_gate_across_cutoff_boundary_geometries(scores,rna,cutoff):
    n=len(scores);ks=[math.ceil(1+(n-1)*c/100) for c in range(100,0,-10)]
    maxima=ks.copy();maxima[list(range(100,0,-10)).index(cutoff)]-=1
    sampler=GateSampler([f'id{i}' for i in range(n)],scores,rna,maxima)
    blocks=[np.flatnonzero(np.asarray(scores)==value).tolist() for value in sorted(set(scores))]
    orders=[]
    for parts in itertools.product(*(itertools.permutations(block) for block in blocks)):
        position=[i for part in parts for i in part]
        rank=np.empty(n,dtype=int);rank[position]=np.arange(1,n+1);orders.append(rank)
    ranks=np.asarray(orders);rna=np.asarray(rna)
    passed=np.zeros(len(ranks),dtype=bool)
    for k,maximum in zip(ks,maxima):passed|=(((ranks>n-k)&(rna>n-k)).sum(axis=1)>maximum)
    labels=(ranks>n-ks[5])&(rna>n-ks[5])
    result=analytic(sampler.neff,sampler.rna,sampler.maxima)
    assert result['kind']=='exact_single_gate' and passed.any()
    for i in range(n):
        assert result['members'][i]==Fraction(int(labels[passed,i].sum()),int(passed.sum()))


def test_extreme_positive_gate_never_becomes_zero():
    n=2000;ks=[math.ceil(1+(n-1)*c/100) for c in range(100,0,-10)]
    maximum=ks.copy();maximum[5]-=1
    s=GateSampler([f'id{i:04d}' for i in range(n)],[0]*n,list(range(1,n+1)),maximum)
    answer=analytic(s.neff,s.rna,s.maxima)
    assert answer['kind']=='exact_single_gate'
    assert answer['gate_probability']==Fraction(1,math.comb(n,ks[5]))
    number=exact_number(answer['gate_probability'])
    assert float(answer['gate_probability'])==0.
    assert Decimal(number['decimal'])>0 and number['log'] < -1000
    assert set(answer['members'])=={Fraction(0),Fraction(1)}
    huge=Fraction(1,math.comb(100000,50000));representation=exact_number(huge)
    assert representation['integer_base']==16 and int(representation['numerator'],16)==huge.numerator
    assert int(representation['denominator'],16)==huge.denominator and Decimal(representation['decimal'])>0
    assert math.isfinite(representation['log']) and representation['log'] < -60000


def test_overlapping_conditional_sampling_against_small_exact_marginals(all_eight_orders):
    scores=pd.read_csv(FIXTURES/'all_tied_overlapping_score.tsv',sep='\t')
    thresholds=pd.read_csv(FIXTURES/'all_tied_overlapping_thresholds.csv')
    s=GateSampler(scores.cell.to_numpy(str),scores.Neff.to_numpy(int),scores.RNA_rank.to_numpy(int),thresholds.maximum.to_numpy(int))
    _,counts,labels=all_eight_orders;passed=(counts>s.maxima).any(axis=1)
    truth=labels[passed].mean(axis=0)
    context,records=solve(s);union=np.array([r['conditioned_union'] for r in records],dtype=bool)
    meta,selected,prefixes,_=sample_case({'sampler':s,'nuclear':None},'conditional',2000,31,[100,2000],union)
    assert meta['accepted']>0 and meta['gate_probability']['method']=='overlap_corrected_MC'
    assert meta['gate_probability']['MC_SE_decimal'] is not None
    assert np.max(np.abs(selected/meta['accepted']-truth))<.1
    assert all(x['union_unseen']>=0 for x in prefixes)


def test_integer_contract_and_id_reordering():
    ids=['b','a'];rna=[2,1];maxima=[2,2,2,2,2,2,2,2,2,1]
    first=GateSampler(ids,[0,2**53-1],rna,maxima)
    second=GateSampler(ids[::-1],[2**53-1,0],rna[::-1],maxima)
    assert np.array_equal(first.ids,second.ids)
    assert np.array_equal(first.neff,second.neff)
    assert np.array_equal(first.rna,second.rna)
    for bad in ([0.,1.],[0,2**53],[0,float('nan')],[0,-1]):
        with pytest.raises(ValueError):GateSampler(ids,bad,rna,maxima)
    with pytest.raises(ValueError):GateSampler(['','valid'],[0,1],rna,maxima)
    for missing in (None,float('nan'),pd.NA):
        with pytest.raises(ValueError):GateSampler([missing,'valid'],[0,1],rna,maxima)
    with pytest.raises(ValueError):analytic([0.,1.],rna,maxima)
    with pytest.raises(ValueError):analytic([0,1],rna,maxima,[100.,90.,80.,70.,60.,50.,40.,30.,20.,10.])
    assert GateSampler(ids,['0','1.0'],rna,maxima).n==2
