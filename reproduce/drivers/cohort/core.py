"""Compute and retain exact family states with explicit undefined conditions."""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from merci_rank_audit.bounds import global_witness, top_k
from merci_rank_audit.sampler import GateSampler
from merci_rank_audit.conditioned import solve, evaluate

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def digest_ids(ids):
    h=hashlib.sha256()
    for item in ids:
        raw=str(item).encode('utf-8');h.update(len(raw).to_bytes(4,'little'));h.update(raw)
    return h.hexdigest()


def ratio(a,b):return float(a/b) if b else None


def states(core,union):
    if np.any(core&~union):raise ValueError('Core is not contained in union')
    return np.where(core,0,np.where(union,1,2)).astype(np.int8)


def model_case(arrays,maxima):
    s=GateSampler(arrays['ids'],arrays['neff'],arrays['rna_rank'],maxima);order=s.order;n=s.n;k=top_k(n,50)
    dna=np.asarray(arrays['dna_rank'])[order];truth=np.asarray(arrays['truth'],bool)[order]
    observed=np.asarray(arrays['observed_source_features'],float)[order]
    if not s.validate_order(dna):raise ValueError('Projected DNA order is not legal')
    context,records=solve(s);bound=context['bound'];can=context['possible']
    lower=global_witness(s.neff,s.rna,False);lower_counts=s.counts(lower);upper_counts=context['upper_counts']
    original_counts=s.counts(dna);must=bool((lower_counts>s.maxima).any());gate=bool((original_counts>s.maxima).any())
    if must!=s.certain or can!=bool(s.active):raise ValueError('Gate state differs from exact bounds')
    raw_core=bound['must'];raw_union=bound['possible'];raw_call=(dna>n-k)&(s.rna>n-k)
    default_core=raw_core if must else np.zeros(n,bool);default_union=raw_union if can else np.zeros(n,bool)
    raw_state=states(raw_core,raw_union);default_state=states(default_core,default_union);default_call=raw_call&gate
    if np.any(raw_core&~raw_call) or np.any(raw_call&~raw_union) or np.any(default_core&~default_call) or np.any(default_call&~default_union):
        raise ValueError('Original call is outside its raw/default family')
    conditional_state=None
    if can:
        conditional_core=np.array([r['conditioned_core'] for r in records],bool)
        conditional_union=np.array([r['conditioned_union'] for r in records],bool)
        if np.any(raw_core&~conditional_core) or np.any(conditional_union&~raw_union):raise ValueError('Conditional family does not restrict raw family')
        conditional_state=states(conditional_core,conditional_union)
    dna_core=bound['strict']|(bound['boundary']&(bound['boundary_s']==bound['boundary_b']))
    dna_union=bound['strict']|(bound['boundary']&(bound['boundary_s']>0));dna_state=states(dna_core,dna_union)
    rna_top=s.rna>n-k
    if not np.array_equal(raw_state,np.where(rna_top,dna_state,2)):raise ValueError('DNA/RNA state composition differs from raw family')
    return dict(sampler=s,context=context,ids=s.ids,neff=s.neff,rna_rank=s.rna,dna_rank=dna,truth=truth,observed=observed,
        raw_state=raw_state,default_state=default_state,conditional_state=conditional_state,dna_state=dna_state,rna_top=rna_top,
        raw_call=raw_call,default_call=default_call,gate=gate,can=can,must=must,N=n,K=k,maxima=s.maxima,
        original_counts=original_counts,lower_counts=lower_counts,upper_counts=upper_counts)


def confusion(selected,truth):
    if selected is None:return dict(selected=None,TP=None,FP=None,TN=None,FN=None,precision=None,sensitivity=None,FPR=None)
    selected=np.asarray(selected,bool);truth=np.asarray(truth,bool);tp=int((selected&truth).sum());fp=int((selected&~truth).sum())
    tn=int((~selected&~truth).sum());fn=int((~selected&truth).sum())
    return dict(selected=int(selected.sum()),TP=tp,FP=fp,TN=tn,FN=fn,precision=ratio(tp,tp+fp),sensitivity=ratio(tp,tp+fn),FPR=ratio(fp,fp+tn))


def write_members(model,path):
    if path.exists():raise FileExistsError(path)
    m=model;n=m['N'];cond=m['conditional_state'] if m['conditional_state'] is not None else np.full(n,-1,dtype=np.int8)
    np.savez_compressed(path,ids=m['ids'],neff=m['neff'],dna_rank=m['dna_rank'],rna_rank=m['rna_rank'],truth=m['truth'],
        observed_source_features=m['observed'],raw_state=m['raw_state'],default_state=m['default_state'],conditional_state=cond,
        dna_state=m['dna_state'],rna_top=m['rna_top'],raw_call=m['raw_call'],default_call=m['default_call'],maxima=m['maxima'])
    return sha(path)


def witness(model,index,selected,path):
    value=evaluate(model['context'],int(index),bool(selected));rank=value['rank'];s=model['sampler']
    result=dict(cell_ID_sha256=digest_ids([s.ids[index]]),selected=bool(selected),status=value['status'],reason=value['reason'],
                gate_counts=value['counts'].tolist() if value['counts'] is not None else None,rank_sha256=None)
    if rank is not None:
        if not s.validate_order(rank) or not bool((s.counts(rank)>s.maxima).any()):raise ValueError('Invalid passing witness')
        if bool(rank[index]>s.n-model['K'] and s.rna[index]>s.n-model['K'])!=bool(selected):raise ValueError('Witness membership differs')
        if path.exists():raise FileExistsError(path)
        np.savez_compressed(path,ids=s.ids,rank=rank,selected_index=int(index),maxima=s.maxima);result['rank_sha256']=sha(path)
    return result
