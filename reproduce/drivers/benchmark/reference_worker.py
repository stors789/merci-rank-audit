"""Create one independent analytic or bounded streaming multi-gate reference."""
import argparse
import json
import math
import os
import resource
import sys
import time
from decimal import Decimal,localcontext
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.stats import binom

from merci_rank_audit.sampler import GateSampler
from merci_rank_audit.conditioned import solve
from analytic_reference import analytic
from case_inputs import load_prepared,sha_file,resolve_case


def fraction_text(value):
    with localcontext() as context:
        context.prec=40
        return str(Decimal(value.numerator)/Decimal(value.denominator))


def wilson_half(counts,n):
    z=1.959963984540054;p=counts/n
    return z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)


def count_quantiles(hist,n,raw_range):
    cumulative=np.cumsum(hist);rows=[];limit=max(1.,.05*(raw_range[1]-raw_range[0]))
    def order_stat(rank):return int(np.searchsorted(cumulative,rank,side='left'))
    for q in (.025,.5,.975):
        point=order_stat(max(1,math.ceil(q*n)));low=int(binom.ppf(.025,n,q));high=int(binom.ppf(.975,n,q))+1
        lo=order_stat(low) if low>=1 else None;hi=order_stat(high) if high<=n else None
        rows.append(dict(q=q,estimate=point,CI95_low=lo,CI95_high=hi,width=hi-lo if lo is not None and hi is not None else None))
    return rows,limit


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--case',required=True);p.add_argument('--output-root',type=Path,required=True);p.add_argument('--spawn-time',type=float,required=True)
    args=p.parse_args();ready=time.perf_counter();config=json.loads(args.config.read_text());case=next(x for x in config['cases'] if x['case']==args.case)
    case=resolve_case(case,args.config.parent)
    output=args.output_root/args.case
    if output.exists():raise FileExistsError('Protect first reference output')
    output.mkdir(parents=True,mode=0o700);started=time.perf_counter()
    meta=dict(case=args.case,N=case['N'],scope=case['scope'],config_sha256=sha_file(args.config),code_commit=config['code_commit'],
        prepared_sha256=case['prepared_sha256'],canonical_array_sha256=case['canonical_array_sha256'],
        cold_start_seconds=round(ready-args.spawn_time,6),status='STARTED')
    try:
        if len(os.listdir('/proc/self/task'))!=1:raise RuntimeError('Reference worker is not single-threaded')
        t=time.perf_counter();ids,neff,rna,maxima=load_prepared(case);meta['input_read_verify_seconds']=round(time.perf_counter()-t,6)
        t=time.perf_counter();s=GateSampler(ids,neff,rna,maxima)
        meta['sampler_preprocess_seconds']=round(time.perf_counter()-t,6);t=time.perf_counter()
        context,records=solve(s);meta['exact_member_summary_seconds']=round(time.perf_counter()-t,6)
        meta.update(gate_state='impossible' if not s.active else 'certain' if s.certain else 'possible_not_certain',
            active_events=len(s.active),raw_core=int(context['bound']['must'].sum()),raw_union=int(context['bound']['possible'].sum()),
            passing_core=sum(r['conditioned_core'] for r in records) if context['possible'] else None,
            passing_union=sum(r['conditioned_union'] for r in records) if context['possible'] else None)
        del records;t=time.perf_counter();ind=analytic(s.neff,s.rna,s.maxima);domain=ind['raw_possible']
        assert np.array_equal(domain,context['bound']['possible'])
        assert tuple(ind['raw_count_range'])==(context['bound']['lower'],context['bound']['upper'])
        for observed,expected in zip(s.events,ind['events']):
            assert observed.get('probability',Fraction())==expected['probability']
        meta['analytic_reference_seconds']=round(time.perf_counter()-t,6);meta['reference_kind']=ind['kind'];meta['raw_possible_N']=int(domain.sum())
        n=s.n;q=np.full(n,np.nan);q_se=np.full(n,np.nan);reference_p=None;p_se=None;p_half=None;max_half=None;hist=np.empty(0,dtype=np.int64)
        if ind['kind']!='multi_gate_reference_required':
            p_exact=ind['gate_probability'];p_value=float(p_exact);reference_p=p_value if p_value>0 or not p_exact else None
            meta['P_G_exact_decimal']=fraction_text(p_exact);meta['P_G_exact_log']=math.log(p_exact.numerator)-math.log(p_exact.denominator) if p_exact else None
            if ind['members'] is not None:q=ind['members'];q_se=np.zeros(n)
            meta.update(reference_status='EXACT',gate_qualified=True,member_qualified=bool(ind['members'] is not None and domain.any()),
                        reference_proposals=0,reference_accepted=None,reference_weight_ESS=None,gate_MC_SE=0.,gate_95_halfwidth=0.,
                        member_max_Wilson_95_halfwidth=0. if ind['members'] is not None and domain.any() else None,
                        member_reference_rms_MC_SE=0. if ind['members'] is not None and domain.any() else None)
            if ind['members'] is None:
                meta.update(count_quantile_status='NA_G_IMPOSSIBLE',count_quantiles=None)
            elif ind['raw_count_range'][0]==ind['raw_count_range'][1]:
                meta.update(count_quantile_status='EXACT_CONSTANT',count_quantiles=[dict(q=q,estimate=ind['raw_count_range'][0],CI95_low=ind['raw_count_range'][0],CI95_high=ind['raw_count_range'][0],width=0) for q in (.025,.5,.975)])
            else:meta.update(count_quantile_status='NOT_NEEDED_FOR_GATE_MEMBER_ANALYTIC_REFERENCE',count_quantiles=None)
        else:
            rng=np.random.default_rng(case['reference_seed']);target=config['reference']['multi_proposals'];limit=config['reference']['multi_wall_limit_seconds']
            hit=0;nprop=0;sum_w=sum_w2=0.;selected=np.zeros(int(domain.sum()),dtype=np.int64);hist=np.zeros(n+1,dtype=np.int64)
            k=int(s.ks[5]);rna_top=s.rna>s.n-k;start=time.perf_counter()
            for _ in range(target):
                if time.perf_counter()-start>=limit:break
                rank,counts,accepted=s.proposal(rng);multiplicity=int((counts>s.maxima).sum());assert multiplicity>=1
                weight=1/multiplicity;sum_w+=weight;sum_w2+=weight*weight;nprop+=1
                if nprop==1 and len(os.listdir('/proc/self/task'))!=1:raise RuntimeError('Reference worker started extra threads')
                if accepted:
                    hit+=1;labels=(rank>s.n-k)&rna_top;selected+=labels[domain];hist[int(labels.sum())]+=1
                if nprop%10000==0:print(json.dumps({'case':args.case,'phase':'reference','proposals':nprop,'seconds':round(time.perf_counter()-start,2)}),flush=True)
            elapsed=time.perf_counter()-start;meta['reference_sampling_seconds']=round(elapsed,6)
            if nprop:
                scale=float(s.S);estimate=scale*sum_w/nprop
                reference_p=estimate if estimate>0 else None
                variance=max(0.,(sum_w2-sum_w*sum_w/nprop)/(nprop-1)) if nprop>1 else None
                p_se=scale*math.sqrt(variance/nprop) if variance is not None and reference_p is not None else None
                p_half=1.959963984540054*p_se if p_se is not None else None
            if hit:
                q[domain]=selected/hit;q_se[domain]=np.sqrt(q[domain]*(1-q[domain])/hit)
                if domain.any():max_half=float(np.max(wilson_half(selected,hit)))
            full=nprop==target;gate_ok=full and reference_p is not None and p_half is not None and p_half<=max(1e-5,.1*reference_p)
            member_ok=full and hit>=config['reference']['minimum_accepted'] and domain.any() and max_half<=config['reference']['max_marginal_wilson_95_halfwidth']
            qtiles,qlimit=count_quantiles(hist,hit,ind['raw_count_range']) if hit else (None,None)
            count_ok=full and hit>=config['reference']['minimum_accepted'] and qtiles is not None and all(x['width'] is not None and x['width']<=qlimit for x in qtiles)
            meta.update(reference_status='MC_COMPLETE' if full else 'MC_TIMEOUT_INCOMPLETE',gate_qualified=bool(gate_ok),
                member_qualified=bool(member_ok and gate_ok),reference_proposals=nprop,reference_accepted=hit,
                reference_weight_ESS=sum_w*sum_w/sum_w2 if sum_w2 else None,gate_MC_SE=p_se,gate_95_halfwidth=p_half,
                member_max_Wilson_95_halfwidth=max_half,
                member_reference_rms_MC_SE=float(np.sqrt(np.mean(q_se[domain]**2))) if hit and domain.any() else None,
                count_quantile_status='QUALIFIED' if count_ok else 'NOT_QUALIFIED',count_quantile_width_limit=qlimit,count_quantiles=qtiles,
                P_G_exact_decimal=None,P_G_exact_log=s.logS+math.log(sum_w/nprop) if nprop else None)
        if len(os.listdir('/proc/self/task'))!=1:raise RuntimeError('Reference worker started extra threads')
        meta['reference_P_G']=reference_p;meta['input_and_event_probabilities_checked']=True
        t=time.perf_counter();np.savez_compressed(output/'reference.npz',member_probability=q,member_MC_SE=q_se,raw_possible=domain,accepted_count_histogram=hist)
        meta['reference_output_seconds']=round(time.perf_counter()-t,6)
        meta['reference_npz_sha256']=sha_file(output/'reference.npz');meta['status']='COMPLETE'
    except Exception as error:
        meta.update(status='ERROR',error=repr(error));(output/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
        print(json.dumps({'case':args.case,'status':'ERROR','error':repr(error)}),flush=True);raise
    meta['process_elapsed_seconds']=round(time.perf_counter()-args.spawn_time,6);meta['work_seconds']=round(time.perf_counter()-started,6)
    meta['threads']=len(os.listdir('/proc/self/task'));meta['actual_core']=sorted(os.sched_getaffinity(0))
    meta['nice']=os.getpriority(os.PRIO_PROCESS,0);meta['AS_limit_bytes']=resource.getrlimit(resource.RLIMIT_AS)[0]
    meta['peak_RSS_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (output/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps({'case':args.case,'status':meta['status'],'reference_status':meta['reference_status'],
                      'kind':meta['reference_kind'],'qualified':meta['member_qualified']}),flush=True)


if __name__=='__main__':main()
