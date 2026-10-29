"""Exact 50%-member reachability conditional on a nested OR gate."""
import numpy as np
from .bounds import global_witness, intersection_bound, top_k


def setup(sampler):
    s=sampler;k=top_k(s.n,50);bound=intersection_bound(s.neff,s.rna,k)
    upper=global_witness(s.neff,s.rna,True);upper_counts=s.counts(upper)
    possible=bool((upper_counts>s.maxima).any())
    assert possible==bool(s.active)
    block=np.flatnonzero(bound['boundary']);ordered=block[np.argsort(-s.rna[block],kind='stable')]
    positions=np.full(s.n,-1,dtype=int);positions[ordered]=np.arange(len(ordered))
    prefix=np.clip(s.ks-int(bound['strict'].sum()),0,len(block))
    return dict(s=s,k=k,bound=bound,upper=upper,upper_counts=upper_counts,possible=possible,
                block=block,ordered=ordered,positions=positions,prefix=prefix)


def constrained_best(context,index,selected):
    """Return a complete order maximizing all gates under the member constraint."""
    s=context['s'];b=context['bound'];i=int(index)
    if selected and not b['possible'][i]:return None,'structurally_cannot_select'
    if not selected and b['must'][i]:return None,'structurally_cannot_omit'
    if (selected and b['must'][i]) or (not selected and not b['possible'][i]):
        return context['upper'].copy(),'global_upper_satisfies_constraint'
    block=context['block'];size=int(b['boundary_s'])
    assert i in block and 0<size<len(block)
    ordered=context['ordered']
    others=ordered[ordered!=i];chosen=np.zeros(s.n,bool)
    if selected:
        chosen[i]=True;chosen[others[:size-1]]=True
    else:
        chosen[others[:size]]=True
    high=ordered[chosen[ordered]];low=ordered[~chosen[ordered]]
    rank=context['upper'].copy();hi=int(rank[block].max());lo=hi-len(block)
    rank[high]=np.arange(hi,hi-size,-1);rank[low]=np.arange(hi-size,lo,-1)
    assert s.validate_order(rank)
    actual=bool(rank[i]>s.n-context['k'] and s.rna[i]>s.n-context['k'])
    assert actual==selected
    return rank,'constrained_upper'


def evaluate_summary(context,index,selected):
    """Exact constrained gate counts without materializing a full witness."""
    if not context['possible']:
        return dict(status='undefined_G_impossible',reachable=None,counts=None,reason='condition_has_no_orders')
    s=context['s'];b=context['bound'];i=int(index)
    if selected and not b['possible'][i]:
        return dict(status='not_reachable',reachable=False,counts=None,reason='structurally_cannot_select')
    if not selected and b['must'][i]:
        return dict(status='not_reachable',reachable=False,counts=None,reason='structurally_cannot_omit')
    counts=context['upper_counts'].copy()
    fixed=(selected and b['must'][i]) or (not selected and not b['possible'][i])
    if not fixed:
        pos=int(context['positions'][i]);size=int(b['boundary_s']);ordered=context['ordered']
        assert pos>=0 and 0<size<len(ordered)
        for j,t in enumerate(context['prefix']):
            if selected and pos>=size and size<=t<=pos:
                other=ordered[t-1]
                counts[j]+=int(s.rna[i]>s.n-s.ks[j])-int(s.rna[other]>s.n-s.ks[j])
            elif not selected and pos<size and pos<t<=size:
                other=ordered[t]
                counts[j]+=int(s.rna[other]>s.n-s.ks[j])-int(s.rna[i]>s.n-s.ks[j])
    passed=bool((counts>s.maxima).any())
    return dict(status='reachable' if passed else 'not_reachable',reachable=passed,counts=counts,
                reason='passing_constrained_witness' if passed else 'all_constrained_gate_maxima_fail')


def evaluate(context,index,selected):
    if not context['possible']:
        return dict(status='undefined_G_impossible',reachable=None,rank=None,counts=None,reason='condition_has_no_orders')
    rank,reason=constrained_best(context,index,selected)
    if rank is None:
        return dict(status='not_reachable',reachable=False,rank=None,counts=None,reason=reason)
    s=context['s'];counts=s.counts(rank);passed=bool((counts>s.maxima).any())
    return dict(status='reachable' if passed else 'not_reachable',reachable=passed,rank=rank if passed else None,
                counts=counts,reason='passing_constrained_witness' if passed else 'all_constrained_gate_maxima_fail')


def solve(sampler):
    context=setup(sampler);s=sampler;b=context['bound'];records=[]
    if not context['possible']:
        for i in range(s.n):
            records.append(dict(cell=s.ids[i],raw_core=bool(b['must'][i]),raw_union=bool(b['possible'][i]),
                can_select=None,can_omit=None,conditioned_core=None,conditioned_union=None,status='undefined_G_impossible'))
        return context,records
    for i in range(s.n):
        inside=evaluate_summary(context,i,True);outside=evaluate_summary(context,i,False)
        assert inside['reachable'] or outside['reachable']
        records.append(dict(cell=s.ids[i],raw_core=bool(b['must'][i]),raw_union=bool(b['possible'][i]),
            can_select=inside['reachable'],can_omit=outside['reachable'],
            conditioned_core=not outside['reachable'],conditioned_union=inside['reachable'],status='exact',
            selected_reason=inside['reason'],omitted_reason=outside['reason'],
            selected_max_gate_margin=int(np.max(inside['counts']-s.maxima)) if inside['counts'] is not None else None,
            omitted_max_gate_margin=int(np.max(outside['counts']-s.maxima)) if outside['counts'] is not None else None))
    return context,records
