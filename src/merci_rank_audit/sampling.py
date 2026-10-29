"""Streaming uniform or overlap-corrected conditional samples."""
import math
import time
from decimal import Decimal, localcontext

import numpy as np

from .fwl import beta_batch
from .inputs import ordinary_interval, percentile_with_MC_interval


def sample_case(case,method,budget,seed,prefixes=None,passing_union=None,batch_size=256,on_prefix=None):
    if method not in ('uniform','conditional') or budget<1 or batch_size<1:
        raise ValueError('Invalid sampling arguments')
    sampler=case['sampler']; n=sampler.n; nuclear=case['nuclear']
    prefixes=sorted(set(prefixes or [budget]))
    if not prefixes or prefixes[-1]>budget or prefixes[0]<1:
        raise ValueError('Invalid prefixes')
    rng=np.random.default_rng(seed); k=int(sampler.ks[5]); selected_count=np.zeros(n,dtype=np.int64)
    accepted=0; weights=[]; prefix_rows=[]; group=[]; beta_parts={'original':[],'alternative':[]}
    start=time.monotonic()
    impossible=method=='conditional' and not sampler.active
    if impossible:
        return {'method':method,'status':'undefined_G_impossible','budget':budget,'proposals':0,'accepted':0,
                'seed':int(seed),'seconds':time.monotonic()-start,'gate_probability':None},None,[],[]

    def process_nuclear():
        if not group:return
        labels=np.asarray(group,dtype=bool);group.clear()
        for representation in ('original','alternative'):
            y=nuclear[representation];result=beta_batch(y,nuclear['X'],labels)
            valid=np.isfinite(y).all(axis=0)&(y.std(axis=0,ddof=1)>1e-14)
            result[:,~valid]=np.nan;beta_parts[representation].append(result)

    for proposal in range(1,budget+1):
        if method=='conditional':
            rank,counts,yes=sampler.proposal(rng)
            multiplicity=int((counts>sampler.maxima).sum())
            weights.append(1/multiplicity)
        else:
            rank=sampler.uniform_order(rng)
            counts=sampler.counts(rank)
            yes=bool((counts>sampler.maxima).any())
        if yes:
            accepted+=1
            selected=(rank>n-k)&(sampler.rna>n-k)
            selected_count+=selected
            if nuclear is not None:
                group.append(selected[nuclear['index']][nuclear['observable']])
                if len(group)>=batch_size:process_nuclear()
        if proposal in prefixes:
            seen=selected_count>0
            if passing_union is not None and np.any(seen & ~passing_union):
                raise ValueError('A sampled member is outside the exact passing-only union')
            if on_prefix is not None:on_prefix(proposal,seen.copy())
            prefix_rows.append({'method':method,'budget':proposal,'accepted':accepted,'seen_members':int(seen.sum()),
                                'exact_passing_union':int(passing_union.sum()) if passing_union is not None else None,
                                'union_unseen':int((passing_union&~seen).sum()) if passing_union is not None else None,
                                'frequency_status':'defined' if accepted else 'NA_no_accepted_order',
                                'seconds':time.monotonic()-start})
    if nuclear is not None:process_nuclear()
    rows=[]
    if nuclear is not None:
        for representation in ('original','alternative'):
            values=np.concatenate(beta_parts[representation],axis=0) if beta_parts[representation] else np.empty((0,len(nuclear['names'])))
            for j,name in enumerate(nuclear['names']):
                v=values[:,j];valid=v[np.isfinite(v)]
                q=[percentile_with_MC_interval(valid,p) if len(valid)>=500 else (None,None,None) for p in (.025,.5,.975)]
                rows.append({'method':method,'representation':representation,'feature':str(name),'accepted':accepted,
                             'valid_betas':len(valid),'invalid_betas':len(v)-len(valid),
                             'quantile_status':'defined' if len(valid)>=500 else 'not_reported_under_500_valid',
                             'q025':q[0][0],'q025_MC_low':q[0][1],'q025_MC_high':q[0][2],
                             'median':q[1][0],'median_MC_low':q[1][1],'median_MC_high':q[1][2],
                             'q975':q[2][0],'q975_MC_low':q[2][1],'q975_MC_high':q[2][2]})
    if method=='uniform':
        low,high=ordinary_interval(accepted,budget)
        gate={'estimate':accepted/budget,'interval_low':low,'interval_high':high,'method':'uniform_binomial'}
    elif sampler.certain:
        gate={'estimate':1.0,'method':'exact_certain'}
    else:
        mean=float(np.mean(weights))
        relative_se=float(np.std(weights,ddof=1)/math.sqrt(budget)/mean) if budget>1 else None
        with localcontext() as ctx:
            ctx.prec=40
            probability=(Decimal(sampler.S.numerator)/Decimal(sampler.S.denominator))*Decimal(str(mean))
            se=probability*Decimal(str(relative_se)) if relative_se is not None else None
        gate={'estimate_decimal':str(probability),'log_estimate':sampler.logS+math.log(mean),
              'MC_SE_decimal':str(se) if se is not None else None,'MC_relative_SE':relative_se,
              'method':'overlap_corrected_MC','mean_inverse_multiplicity':mean}
    meta={'method':method,'status':'complete' if accepted else 'NA_no_accepted_order','budget':budget,
          'proposals':budget,'accepted':accepted,'seed':int(seed),'seconds':time.monotonic()-start,'gate_probability':gate}
    return meta,selected_count if accepted else None,prefix_rows,rows
