"""Portable extraction of the original composition-comparison functions."""
import csv
import json
from pathlib import Path
import numpy as np
from core import confusion, model_case, ratio, write_members

LEVELS=(0,25,50,75,100);FAMILIES=('raw','default','conditional');STATES=('core','variable','outside')
REPRESENTATIVE_SPLITS={('bench1',1),('bench1',21),('bench2',1)}


def csv_write(path,rows):
    if not rows:return
    with path.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),extrasaction='ignore',lineterminator='\n');writer.writeheader();writer.writerows(rows)


def projected_rank(rank,indices):
    subset=np.asarray(rank)[indices];order=np.argsort(subset,kind='stable');result=np.empty(len(indices),dtype=int);result[order]=np.arange(1,len(indices)+1)
    return result


def indices_for(ids,wanted):
    idx=np.searchsorted(ids,wanted)
    if np.any(idx>=len(ids)) or not np.array_equal(ids[idx],wanted):raise ValueError('Reference IDs are not constant across conditions')
    return idx


def make_condition(data,split,mode,level,maxima,members):
    m=split['M'];cc=data['cc_order'];anchor=data['anchor'];background=data['background'];b=len(background)
    if mode=='nested':
        take=len(cc)*level//100;indices=np.concatenate((anchor,background,cc[:take]));ref=np.concatenate((anchor,background));background_replaced=None
    else:
        take=b*level//100;indices=np.concatenate((anchor,background[:b-take],cc[:take]));ref=anchor;background_replaced=take
    if len(np.unique(indices))!=len(indices) or (mode=='fixedN' and len(indices)!=m):raise ValueError('Cohort selection violates fixed membership or N')
    arrays=dict(ids=data['ids'][indices],neff=data['neff'][indices],rna_rank=projected_rank(data['rna_rank'],indices),
        dna_rank=projected_rank(data['dna_rank'],indices),truth=data['truth'][indices],observed_source_features=data['observed_source_features'][indices])
    model=model_case(arrays,maxima);model['rna_score']=data['rna_score'][indices][model['sampler'].order]
    ref_ids=np.sort(data['ids'][ref]);untied_ids=np.sort(data['ids'][ref[data['untied_endpoints'][ref]]])
    model['reference_ids']=ref_ids;model['untied_ids']=untied_ids
    model['member_sha256']=write_members(model,members/(f"{split['experiment']}_{split['split']:03d}_{mode}_{level:03d}.npz"))
    model['CC_added']=take;model['background_replaced']=background_replaced
    return model


def condition_rows(split,mode,level,model):
    n=model['N'];truth=model['truth'];ref=np.zeros(n,bool);ref[indices_for(model['ids'],model['reference_ids'])]=True
    untied=np.zeros(n,bool);untied[indices_for(model['ids'],model['untied_ids'])]=True
    base=dict(experiment=split['experiment'],stage=split['stage'],split=split['split'],design=mode,level=level,
        CC_added=model['CC_added'],MC_n=int((~truth).sum()),N=n,K=model['K'],background_slots_replaced=model['background_replaced'],
        background_replaced_fraction=ratio(model['background_replaced'],split['background_slots']) if model['background_replaced'] is not None else None,
        cohort_CC_fraction=ratio(int(truth.sum()),n),ordering_origin='parent_order_projection',projected_gate=model['gate'],gate_can_pass=model['can'],gate_must_pass=model['must'],
        projected_counts=json.dumps(model['original_counts'].tolist()),maxima=json.dumps(model['maxima'].tolist()),member_output_sha256=model['member_sha256'])
    rows=[]
    for group,mask in (('cohort',np.ones(n,bool)),('reference',ref),('RNA_untied_endpoints',untied)):
        group_truth=truth[mask];raw=confusion(model['raw_call'][mask],group_truth);default=confusion(model['default_call'][mask],group_truth)
        for family in FAMILIES:
            state=model[family+'_state'];defined=state is not None
            state_group=state[mask] if defined else None;core=(state_group==0) if defined else None
            core_metrics=confusion(core,group_truth)
            row=dict(**base,group=group,family=family,status='DEFINED' if defined else 'NA_G_IMPOSSIBLE',group_N=int(mask.sum()),
                group_CC=int(group_truth.sum()),group_MC=int((~group_truth).sum()),
                core_N=int((state_group==0).sum()) if defined else None,variable_N=int((state_group==1).sum()) if defined else None,
                outside_N=int((state_group==2).sum()) if defined else None,
                core_CC=int(((state_group==0)&group_truth).sum()) if defined else None,core_MC=int(((state_group==0)&~group_truth).sum()) if defined else None,
                union_CC=int(((state_group!=2)&group_truth).sum()) if defined else None,union_MC=int(((state_group!=2)&~group_truth).sum()) if defined else None)
            row.update({'core_'+k:v for k,v in core_metrics.items()});row.update({'projected_raw_'+k:v for k,v in raw.items()});row.update({'projected_default_'+k:v for k,v in default.items()})
            rows.append(row)
    return rows


def compare_levels(split,mode,left,right,a,b):
    if not np.array_equal(a['reference_ids'],b['reference_ids']) or not np.array_equal(a['untied_ids'],b['untied_ids']):raise ValueError('Reference set changed')
    ids=a['reference_ids'];ia=indices_for(a['ids'],ids);ib=indices_for(b['ids'],ids);untied=np.isin(ids,a['untied_ids'])
    base=dict(experiment=split['experiment'],stage=split['stage'],split=split['split'],design=mode,from_level=left,to_level=right,
        from_N=a['N'],to_N=b['N'],from_K=a['K'],to_K=b['K'],from_projected_gate=a['gate'],to_projected_gate=b['gate'],
        from_can=a['can'],to_can=b['can'],from_must=a['must'],to_must=b['must'])
    transitions=[];channels=[];joint=[];calls=[]
    for group,mask in (('reference',np.ones(len(ids),bool)),('RNA_untied_endpoints',untied)):
        group_n=int(mask.sum());da=a['dna_state'][ia];db=b['dna_state'][ib];ra=a['rna_top'][ia];rb=b['rna_top'][ib]
        dna_changed=da!=db;rna_changed=ra!=rb
        categories={'DNA_only':dna_changed&~rna_changed,'RNA_only':~dna_changed&rna_changed,'both':dna_changed&rna_changed,'unchanged':~dna_changed&~rna_changed}
        raw_a=a['raw_state'][ia];raw_b=b['raw_state'][ib]
        for category,cmask in categories.items():
            q=mask&cmask;channels.append(dict(**base,group=group,group_N=group_n,channel=category,count=int(q.sum()),
                raw_core_to_outside=int((q&(raw_a==0)&(raw_b==2)).sum()),raw_outside_to_core=int((q&(raw_a==2)&(raw_b==0)).sum())))
        if sum(int((mask&v).sum()) for v in categories.values())!=group_n:raise ValueError('Channel partition not exhaustive')
        index=((da*3+db)*4+ra.astype(int)*2+rb.astype(int));counts=np.bincount(index[mask],minlength=36)
        for x in range(3):
            for y in range(3):
                for u in range(2):
                    for v in range(2):joint.append(dict(**base,group=group,group_N=group_n,DNA_from=STATES[x],DNA_to=STATES[y],
                        RNA_from=bool(u),RNA_to=bool(v),count=int(counts[(x*3+y)*4+u*2+v])))
        if int(counts.sum())!=group_n:raise ValueError('DNA/RNA cross-tab not exhaustive')
        for family in FAMILIES:
            sa=a[family+'_state'];sb=b[family+'_state'];defined=sa is not None and sb is not None
            table=np.bincount((sa[ia][mask]*3+sb[ib][mask]),minlength=9) if defined else None
            for x in range(3):
                for y in range(3):transitions.append(dict(**base,group=group,group_N=group_n,family=family,status='DEFINED' if defined else 'NA_G_IMPOSSIBLE',
                    from_state=STATES[x],to_state=STATES[y],count=int(table[x*3+y]) if defined else None,
                    fraction=ratio(int(table[x*3+y]),group_n) if defined else None))
            if defined and int(table.sum())!=group_n:raise ValueError('Family cross-tab not exhaustive')
        for policy in ('raw','default'):
            ca=a[policy+'_call'][ia];cb=b[policy+'_call'][ib];counts=np.bincount(ca[mask].astype(int)*2+cb[mask].astype(int),minlength=4)
            calls.append(dict(**base,group=group,group_N=group_n,policy=policy,ordering_origin='parent_order_projection',both_unselected=int(counts[0]),
                from_only=int(counts[2]),to_only=int(counts[1]),both_selected=int(counts[3]),flipped=int(counts[1]+counts[2])))
    return transitions,channels,joint,calls
