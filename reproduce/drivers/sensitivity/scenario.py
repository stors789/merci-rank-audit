"""Portable original BIO-06 scenario and transition functions."""
import hashlib
import numpy as np
from merci_rank_audit.sampler import GateSampler
from merci_rank_audit.conditioned import solve, evaluate
from merci_rank_audit.bounds import global_witness

def array_hash(x):
    return hashlib.sha256(np.asarray(x, dtype="<i8").tobytes()).hexdigest()

FAMILIES=("raw","default","passing")
STATES=("outside","variable","core")
def family_array(core,union):
    if core is None:return None
    core=np.asarray(core,dtype=bool);union=np.asarray(union,dtype=bool)
    assert np.all(~core|union)
    return np.where(core,2,np.where(union,1,0)).astype(np.int8)

def evaluate_case(ids,neff,rna,dna,maxima):
    s=GateSampler(ids,neff,rna,maxima)
    assert np.array_equal(s.ids,ids) and s.validate_order(dna)
    counts=s.counts(dna)
    independent=((dna[:,None]>s.n-s.ks)&(rna[:,None]>s.n-s.ks)).sum(axis=0)
    assert np.array_equal(counts,independent)
    ctx,recs=solve(s)
    low=s.counts(global_witness(s.neff,s.rna,False))
    must=bool((low>s.maxima).any());can=ctx["possible"]
    fam={}
    for f,c,u in [("raw","raw_core","raw_union"),("default","default_core","default_union"),
                  ("passing","conditioned_core","conditioned_union")]:
        if f=="default":
            co=[bool(r["raw_core"] and must) for r in recs]
            un=[bool(r["raw_union"] and can) for r in recs]
        elif f=="passing" and not can:co=un=None
        else:
            co=[r[c] for r in recs];un=[r[u] for r in recs]
        fam[f]=family_array(co,un)
    raw=(dna>s.n-s.ks[5])&(rna>s.n-s.ks[5]);original_gate=bool((counts>s.maxima).any())
    return dict(s=s,ctx=ctx,family=fam,raw_call=raw,default_call=raw&original_gate,
                original_gate=original_gate,original_counts=counts,upper=ctx["upper_counts"],lower=low,
                gate_state="impossible" if not can else "certain" if must else "possible_not_certain")

def transition(base,new,family):
    if base is None or new is None:
        d=dict(status="NA_G_IMPOSSIBLE",base_core=None,base_outside=None,
               changed=None,core_lost=None,core_gained=None,union_lost=None,union_gained=None)
        d.update({f"{a}_to_{b}":None for a in STATES for b in STATES})
        return d
    d=dict(status="DEFINED",base_core=int((base==2).sum()),base_outside=int((base==0).sum()),
       changed=int((base!=new).sum()),core_lost=int(((base==2)&(new!=2)).sum()),
       core_gained=int(((base!=2)&(new==2)).sum()),
       union_lost=int(((base!=0)&(new==0)).sum()),union_gained=int(((base==0)&(new!=0)).sum()),
       core_to_outside=int(((base==2)&(new==0)).sum()),
       outside_to_core=int(((base==0)&(new==2)).sum()))
    for i,a in enumerate(STATES):
        for j,b in enumerate(STATES):d[f"{a}_to_{b}"]=int(((base==i)&(new==j)).sum())
    assert sum(d[f"{a}_to_{b}"] for a in STATES for b in STATES)==len(base)
    return d

def witness(base,current,ids,private):
    s=current["s"];ctx=current["ctx"];dna_counts=current["original_counts"]
    out=dict(original_gate=current["original_gate"],original_counts=dna_counts.tolist(),
             maxima=s.maxima.tolist(),upper_counts=current["upper"].tolist(),lower_counts=current["lower"].tolist(),
             global_upper_legal=bool(s.validate_order(ctx["upper"])),first_member=None)
    changed=np.zeros(s.n,dtype=bool)
    for f in FAMILIES:
        a=base["family"][f];b=current["family"][f]
        if a is not None and b is not None:changed|=(a!=b)
    changed|=(base["default_call"]!=current["default_call"])
    if changed.any() and ctx["possible"]:
        index=int(np.flatnonzero(changed)[0])
        out["first_member"]=dict(ID_sha256=hashlib.sha256(ids[index].encode()).hexdigest(),
                                 family_states={f:None if current["family"][f] is None else int(current["family"][f][index]) for f in FAMILIES},
                                 baseline_states={f:None if base["family"][f] is None else int(base["family"][f][index]) for f in FAMILIES},
                                 checks=[])
        saved={}
        for label,case in [("baseline",base),("scenario",current)]:
            for selected in (True,False):
                ans=evaluate(case["ctx"],index,selected);rec=dict(context=label,selected=selected,
                    status=ans["status"],reachable=ans["reachable"],reason=ans["reason"],
                    counts=None if ans["counts"] is None else ans["counts"].tolist())
                if ans["rank"] is not None:
                    rank=ans["rank"];sampler=case["s"]
                    assert sampler.validate_order(rank) and np.array_equal(sampler.counts(rank),ans["counts"])
                    assert ((rank[index]>sampler.n-sampler.ks[5]) and (sampler.rna[index]>sampler.n-sampler.ks[5]))==selected
                    assert (ans["counts"]>sampler.maxima).any()
                    saved[label+"_"+str(selected)]=rank
                    rec["rank_sha256"]=array_hash(rank)
                out["first_member"]["checks"].append(rec)
        if saved:np.savez_compressed(private,**saved)
    return out
