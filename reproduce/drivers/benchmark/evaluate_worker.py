"""Stream twenty fixed replicates for one case, retaining every checkpoint."""
import argparse
import hashlib
import json
import math
import os
import resource
import sys
import time
from datetime import datetime
from pathlib import Path
import numpy as np

from merci_rank_audit.sampler import GateSampler
from merci_rank_audit.conditioned import setup
from case_inputs import load_prepared,sha_file,resolve_case


def atomic_json(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    if path.exists() or temporary.exists():raise FileExistsError(path)
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temporary.rename(path)


def rmse(estimate,truth):
    return float(np.sqrt(np.mean((estimate-truth)**2))) if len(estimate) else None


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--case',required=True);parser.add_argument('--reference-root',type=Path,required=True)
    parser.add_argument('--output-root',type=Path,required=True);parser.add_argument('--spawn-time',type=float,required=True)
    args=parser.parse_args();ready=time.perf_counter();config=json.loads(args.config.read_text());case=next(c for c in config['cases'] if c['case']==args.case)
    case=resolve_case(case,args.config.parent)
    if len(os.listdir('/proc/self/task'))!=1:raise RuntimeError('Evaluation worker is not single-threaded')
    config_hash=sha_file(args.config);folder=args.output_root/args.case;folder.mkdir(parents=True,exist_ok=True,mode=0o700)
    identity=dict(case=args.case,config_sha256=config_hash,code_commit=config['code_commit'],prepared_sha256=case['prepared_sha256'],
                  canonical_array_sha256=case['canonical_array_sha256'])
    identity_path=folder/'identity.json'
    if identity_path.exists():
        if json.loads(identity_path.read_text())!=identity:raise ValueError('Existing case output identity differs')
    else:atomic_json(identity_path,identity)
    start=time.perf_counter();ids,neff,rna,maxima=load_prepared(case);input_seconds=time.perf_counter()-start
    start=time.perf_counter();s=GateSampler(ids,neff,rna,maxima);sampler_seconds=time.perf_counter()-start
    start=time.perf_counter();context=setup(s);domain=context['bound']['possible'];setup_seconds=time.perf_counter()-start
    ref_path=args.reference_root/args.case;reference=None;truth=None
    if (ref_path/'metadata.json').exists():
        reference=json.loads((ref_path/'metadata.json').read_text())
        if reference.get('config_sha256')!=config_hash or reference.get('prepared_sha256')!=case['prepared_sha256']:
            raise ValueError('Reference input or configuration differs')
        if reference.get('status')=='COMPLETE':
            npz=ref_path/'reference.npz'
            if sha_file(npz)!=reference['reference_npz_sha256']:raise ValueError('Reference array integrity mismatch')
            with np.load(npz,allow_pickle=False) as data:
                if not np.array_equal(data['raw_possible'],domain):raise ValueError('Reference member domain differs')
                truth=data['member_probability'][domain]
    if reference is None:reference=dict(status='REFERENCE_MISSING',gate_qualified=False,member_qualified=False)
    k=int(s.ks[5]);rna_top=s.rna>s.n-k;domain_n=int(domain.sum());gate_ref=reference.get('reference_P_G')
    gate_ok=bool(reference.get('gate_qualified')) and gate_ref is not None
    member_ok=bool(reference.get('member_qualified') and truth is not None and domain_n>0 and np.isfinite(truth).all())
    ref_noise=reference.get('member_reference_rms_MC_SE') if member_ok else None
    deadline=datetime.fromisoformat(config['deadline']['compute_stop_CST']).timestamp()
    case_meta=dict(**identity,N=s.n,scope=case['scope'],gate_state='impossible' if not s.active else 'certain' if s.certain else 'possible_not_certain',
        raw_possible_N=domain_n,reference_status=reference.get('reference_status',reference.get('status')),
        gate_reference_qualified=gate_ok,member_reference_qualified=member_ok,reference_P_G=gate_ref,
        reference_gate_95_halfwidth=reference.get('gate_95_halfwidth'),reference_member_rms_MC_SE=ref_noise,
        cold_start_seconds=round(ready-args.spawn_time,6),input_read_verify_seconds=round(input_seconds,6),
        sampler_preprocess_seconds=round(sampler_seconds,6),member_setup_seconds=round(setup_seconds,6),
        initial_threads=len(os.listdir('/proc/self/task')),actual_core=sorted(os.sched_getaffinity(0)),
        nice=os.getpriority(os.PRIO_PROCESS,0),AS_limit_bytes=resource.getrlimit(resource.RLIMIT_AS)[0])
    meta_path=folder/'case_metadata.json'
    if not meta_path.exists():atomic_json(meta_path,case_meta)
    max_seconds=float(config['per_flow_max_seconds']);prefixes=list(config['proposal_prefixes']);times=list(config['wall_checkpoints_seconds'])
    completed=0;partial=0;deadline_stopped=False
    for rep in range(config['replicates']):
        done=folder/f'rep_{rep+1:02d}.done.json';journal=folder/f'rep_{rep+1:02d}.jsonl'
        if done.exists():
            marker=json.loads(done.read_text())
            if marker['config_sha256']!=config_hash or not journal.exists() or sha_file(journal)!=marker['journal_sha256']:
                raise ValueError('Existing replicate completion identity mismatch')
            completed+=1;continue
        if journal.exists():partial+=1;continue
        if time.time()+2*max_seconds+5>=deadline:deadline_stopped=True;break
        previous={};stream=journal.open('x',buffering=1);rep_start=time.perf_counter();rows=0;output_seconds=0.
        def append(row):
            nonlocal rows,output_seconds
            before=time.perf_counter();stream.write(json.dumps(row,allow_nan=False)+'\n');output_seconds+=time.perf_counter()-before;rows+=1
        def snapshot(flow,axis,target,n,elapsed,gate,member_a,member_b,accepted,ess,status,seed,gate_log=None):
            if len(os.listdir('/proc/self/task'))!=1:raise RuntimeError('Evaluation worker started extra threads')
            total=time.perf_counter()-args.spawn_time;refstatus='QUALIFIED' if member_ok else 'REFERENCE_NOT_QUALIFIED'
            base=dict(case=args.case,N=s.n,scope=case['scope'],replicate=rep+1,flow=flow,axis=axis,target=target,
                proposals=n,sampling_seconds=round(elapsed,6),process_elapsed_seconds=round(total,6),
                status=status,seed=seed,gate_estimate=gate,
                gate_estimator_id='uniform' if flow=='uniform' else 'shared_conditional',
                gate_log_estimate=gate_log if gate_log is not None else math.log(gate) if gate is not None and gate>0 else None,
                gate_abs_error=abs(gate-gate_ref) if gate_ok and gate is not None else None,
                gate_reference_qualified=gate_ok,raw_possible_N=domain_n,accepted=accepted,weight_ESS=ess,
                reference_status=reference.get('reference_status',reference.get('status')),
                reference_member_rms_MC_SE=ref_noise,reference_gate_95_halfwidth=reference.get('gate_95_halfwidth'),
                reference_log_P_G=reference.get('P_G_exact_log'))
            def member_state(value):
                if status=='NOT_REACHED_60S':return 'NOT_REACHED_60S'
                if status=='EXACT_IMPOSSIBLE_SHORTCUT':return 'NA_G_IMPOSSIBLE'
                if not domain_n:return 'NA_EMPTY_RAW_POSSIBLE_DOMAIN'
                if value is None:return 'NA_ZERO_GATE_HITS' if flow=='uniform' else 'NA_NO_ACCEPTED_PROPOSALS'
                return 'EXACT_SHORTCUT' if status=='EXACT_NO_TIES_SHORTCUT' else 'ESTIMATED'
            if flow=='uniform':
                previous[('uniform',axis,target)]=member_a.copy() if member_a is not None else None
                state=member_state(member_a)
                row=dict(base,method='uniform',member_mean=float(member_a.mean()) if member_a is not None and domain_n else None,
                    member_RMSE=rmse(member_a,truth) if member_ok and member_a is not None else None,
                    member_status=state,RMSE_status=refstatus if state in ('ESTIMATED','EXACT_SHORTCUT') else state,
                    agreement_with_uniform_RMSE=None,weighted_accepted_RMSE=None)
                append(row);return
            uniform=previous.get(('uniform',axis,target));paired=rmse(member_a,member_b) if member_a is not None and member_b is not None else None
            for method,estimate in (('conditional_weighted',member_a),('conditional_accepted',member_b)):
                previous[(method,axis,target)]=estimate.copy() if estimate is not None else None
                state=member_state(estimate)
                row=dict(base,method=method,member_mean=float(estimate.mean()) if estimate is not None and domain_n else None,
                    member_RMSE=rmse(estimate,truth) if member_ok and estimate is not None else None,
                    member_status=state,RMSE_status=refstatus if state in ('ESTIMATED','EXACT_SHORTCUT') else state,
                    agreement_with_uniform_RMSE=rmse(estimate,uniform) if estimate is not None and uniform is not None else None,
                    weighted_accepted_RMSE=paired)
                append(row)
        def shortcut(flow,seed):
            if not s.active:
                gate=0.;q=None;status='EXACT_IMPOSSIBLE_SHORTCUT'
            else:
                rank=np.empty(s.n,dtype=int);rank[np.argsort(s.neff,kind='stable')]=np.arange(1,s.n+1)
                passed=bool((s.counts(rank)>s.maxima).any());assert passed
                gate=1.;q=((rank>s.n-k)&rna_top)[domain].astype(float);status='EXACT_NO_TIES_SHORTCUT'
            for axis,targets in (('proposals',prefixes),('seconds',times)):
                for target in targets:snapshot(flow,axis,target,0,0.,gate,q,q if flow=='conditional' else None,
                    None if not s.active else 0,None,status,seed)
        def run_flow(flow,seed):
            if not s.active or len(s.blocks)==s.n:shortcut(flow,seed);return
            rng=np.random.default_rng(seed);n=hits=0;sum_w=sum_w2=0.;u=np.zeros(domain_n,dtype=np.int64)
            weighted=np.zeros(domain_n,float);accepted=np.zeros(domain_n,dtype=np.int64)
            seen_p=set();seen_t=set();t0=time.perf_counter()
            while True:
                elapsed=time.perf_counter()-t0
                if elapsed>=max_seconds:break
                if flow=='uniform':
                    rank=s.uniform_order(rng);counts=s.counts(rank);yes=bool((counts>s.maxima).any())
                    if yes:hits+=1;u+=((rank>s.n-k)&rna_top)[domain]
                else:
                    rank,counts,yes=s.proposal(rng);multiplicity=1 if s.certain else int((counts>s.maxima).sum())
                    if multiplicity<1:raise AssertionError('Conditional proposal did not pass')
                    weight=1/multiplicity;sum_w+=weight;sum_w2+=weight*weight
                    selected=((rank>s.n-k)&rna_top)[domain];weighted+=weight*selected
                    if yes:hits+=1;accepted+=selected
                n+=1;elapsed=time.perf_counter()-t0
                reached_p=[x for x in prefixes if x not in seen_p and n>=x]
                reached_t=[x for x in times if x not in seen_t and elapsed>=x]
                if reached_p or reached_t:
                    if flow=='uniform':gate=hits/n;gate_log=None;a=u/hits if hits else None;b=None;ess=None
                    else:
                        gate_value=1. if s.certain else float(s.S)*sum_w/n
                        gate=gate_value if gate_value>0 else None
                        gate_log=0. if s.certain else s.logS+math.log(sum_w/n)
                        a=weighted/sum_w;b=accepted/hits if hits else None;ess=sum_w*sum_w/sum_w2
                    for target in reached_p:
                        snapshot(flow,'proposals',target,n,elapsed,gate,a,b,hits,ess,'MEASURED',seed,gate_log);seen_p.add(target)
                    for target in reached_t:
                        snapshot(flow,'seconds',target,n,elapsed,gate,a,b,hits,ess,'MEASURED',seed,gate_log);seen_t.add(target)
                if n>=prefixes[-1] and elapsed>=times[-1]:break
            elapsed=time.perf_counter()-t0
            for axis,targets,seen in (('proposals',prefixes,seen_p),('seconds',times,seen_t)):
                for target in targets:
                    if target not in seen:snapshot(flow,axis,target,n,elapsed,None,None,None,hits,None,'NOT_REACHED_60S',seed)
        try:
            order=('uniform','conditional') if rep%2==0 else ('conditional','uniform')
            for flow in order:
                seed=case['uniform_replicate_seeds'][rep] if flow=='uniform' else case['conditional_replicate_seeds'][rep]
                run_flow(flow,seed)
            for axis,targets in (('proposals',prefixes),('seconds',times)):
                for target in targets:
                    u=previous.get(('uniform',axis,target))
                    for method in ('conditional_weighted','conditional_accepted'):
                        estimate=previous.get((method,axis,target))
                        append(dict(case=args.case,N=s.n,scope=case['scope'],replicate=rep+1,flow='comparison',axis=axis,target=target,
                            method='uniform_vs_'+method,agreement_RMSE=rmse(u,estimate) if u is not None and estimate is not None else None,
                            status='MEASURED' if u is not None and estimate is not None else 'NA_ONE_ESTIMATE_UNDEFINED',
                            raw_possible_N=domain_n))
        except Exception as error:
            stream.close();atomic_json(folder/f'rep_{rep+1:02d}.error.json',dict(**identity,error=repr(error),rows_saved=rows));raise
        stream.close();atomic_json(done,dict(**identity,replicate=rep+1,journal_sha256=sha_file(journal),rows=rows,
            flow_order=list(order),replicate_seconds=round(time.perf_counter()-rep_start,6),
            output_write_seconds=round(output_seconds,6),actual_core=sorted(os.sched_getaffinity(0)),
            threads=len(os.listdir('/proc/self/task')),peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,status='COMPLETE'))
        completed+=1
        print(json.dumps({'case':args.case,'replicate':rep+1,'of':config['replicates'],'seconds':round(time.perf_counter()-rep_start,2)}),flush=True)
    final=dict(**identity,status='COMPLETE' if completed==config['replicates'] else 'PARTIAL',completed_replicates=completed,
        partial_existing=partial,deadline_stopped=deadline_stopped,process_elapsed_seconds=round(time.perf_counter()-args.spawn_time,6),
        final_threads=len(os.listdir('/proc/self/task')),peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    target=folder/'case_status.json'
    if target.exists():
        previous_status=folder/('case_status_previous_'+str(int(time.time()))+'.json');target.rename(previous_status)
    atomic_json(target,final);print(json.dumps(final),flush=True)


if __name__=='__main__':main()
