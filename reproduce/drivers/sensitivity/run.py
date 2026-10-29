"""Run frozen RNA-order/finite-threshold scenarios and gate rarity streams."""
import argparse
import hashlib
import itertools
import json
import math
from decimal import Decimal
from pathlib import Path

import numpy as np
import pandas as pd
from merci_rank_audit.sampler import GateSampler
from scenario import FAMILIES, array_hash, evaluate_case, transition, witness


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def load_case(spec, base):
    score = (base/spec['score']).resolve(); threshold = (base/spec['thresholds']).resolve()
    if spec.get('score_sha256') and sha(score) != spec['score_sha256']: raise ValueError('score digest differs')
    if spec.get('threshold_sha256') and sha(threshold) != spec['threshold_sha256']: raise ValueError('threshold digest differs')
    data = pd.read_csv(score, sep='\t', dtype={'cell':str,'Neff':str,'RNA_rank':str,'DNA_rank':str,'RNA_score':str},keep_default_na=False)
    data = data.sort_values('cell',kind='stable').reset_index(drop=True)
    required={'cell','Neff','RNA_rank','DNA_rank'}
    if spec['rna_orders'] > 1: required.add('RNA_score')
    if not required.issubset(data.columns): raise ValueError('missing score fields: '+str(required-set(data.columns)))
    ids=data.cell.to_numpy(str);neff=data.Neff.to_numpy(str);rna=data.RNA_rank.to_numpy(str);dna=data.DNA_rank.to_numpy(str)
    if len(ids)!=spec['N'] or len(set(ids))!=len(ids): raise ValueError('IDs or N differ')
    return data,ids,neff,np.asarray(rna,dtype=np.int64),np.asarray(dna,dtype=np.int64),threshold


def rna_orders(data, rank, count):
    if count==1: return [rank.copy()]
    groups={}
    for i,value in enumerate(data.RNA_score): groups.setdefault(Decimal(value),[]).append(i)
    blocks=sorted((v for v in groups.values() if len(v)>1),key=lambda z:min(rank[z]))
    sizes=[list(itertools.permutations(sorted(rank[v].tolist()))) for v in blocks]
    orders=[]
    for assignment in itertools.product(*sizes):
        order=rank.copy()
        for idx,values in zip(blocks,assignment):order[idx]=values
        orders.append(order)
    if len(orders)!=count or len({x.tobytes() for x in orders})!=count:raise ValueError('RNA tie orders differ from frozen config')
    original=next(i for i,x in enumerate(orders) if np.array_equal(x,rank))
    return [orders[original]]+[x for i,x in enumerate(orders) if i!=original]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--stage',choices=['both','scenarios','rarity'],default='both')
    a=p.parse_args();c=json.loads(a.config.read_text());base=a.config.resolve().parent
    thresholds=pd.read_csv((base/c['all_thresholds']).resolve());a.output.mkdir(parents=True,exist_ok=False)
    private=a.output/'private_members';private.mkdir(mode=0o700);scenario_rows=[];transition_rows=[];rarity_rows=[];first=None
    for spec in c['cases']:
        data,ids,neff,rna,dna,original_path=load_case(spec,base); orders=rna_orders(data,rna,spec['rna_orders'])
        null_ids=['original']+[x['id'] for x in spec['nulls']]; maxima={}
        for null_id in null_ids:
            z=thresholds[(thresholds.case==spec['name'])&(thresholds.null_id==null_id)].sort_values('cutoff',ascending=False)
            if len(z)!=10 or z.cutoff.tolist()!=list(range(100,0,-10)) or not z.N.eq(len(ids)).all():raise ValueError('threshold table mismatch')
            maxima[null_id]=z.maximum.to_numpy(dtype=np.int64)
        original=pd.read_csv(original_path).sort_values('cutoff',ascending=False)
        if not np.array_equal(original.maximum.to_numpy(np.int64),maxima['original']):raise ValueError('original thresholds disagree')
        if a.stage in ('both','scenarios'):
            baseline=evaluate_case(ids,neff,rna,dna,maxima['original'])
            if spec.get('baseline_expected'):
                e=spec['baseline_expected']
                if (baseline['original_gate']!=e['original_gate'] or baseline['original_counts'].tolist()!=e['original_counts']
                    or [int((baseline['family']['raw']==2).sum()),int((baseline['family']['raw']!=0).sum())]!=e['raw_core_union']):
                    raise ValueError('historical baseline differs')
            for i,order in enumerate(orders):
                states=[];axis=[]
                for null_id in null_ids:
                    cur=evaluate_case(ids,neff,order,dna,maxima[null_id]);sid=f"{spec['name']}_rna{i:02d}_{null_id}"
                    mode='baseline' if i==0 and null_id=='original' else 'RNA_only' if null_id=='original' else 'maxima_only' if i==0 else 'combined'
                    row=dict(scenario=sid,case=spec['name'],N=len(ids),rna_index=i,rna_sha256=array_hash(order),null_id=null_id,
                             comparison=mode,gate_state=cur['gate_state'],original_gate=cur['original_gate'],
                             original_counts=json.dumps(cur['original_counts'].tolist()),maxima=json.dumps(cur['s'].maxima.tolist()),
                             lower_counts=json.dumps(cur['lower'].tolist()),upper_counts=json.dumps(cur['upper'].tolist()),
                             raw_call=int(cur['raw_call'].sum()),default_call=int(cur['default_call'].sum()),
                             raw_call_changed=int((cur['raw_call']!=baseline['raw_call']).sum()),
                             default_call_changed=int((cur['default_call']!=baseline['default_call']).sum()))
                    changed=(cur['gate_state']!=baseline['gate_state'] or cur['original_gate']!=baseline['original_gate']
                             or row['raw_call_changed'] or row['default_call_changed'])
                    for family in FAMILIES:
                        v=cur['family'][family]; old=baseline['family'][family]
                        row[family+'_core']=None if v is None else int((v==2).sum())
                        row[family+'_variable']=None if v is None else int((v==1).sum())
                        row[family+'_outside']=None if v is None else int((v==0).sum())
                        row[family+'_union']=None if v is None else int((v!=0).sum())
                        t=transition(old,v,family);transition_rows.append(dict(scenario=sid,case=spec['name'],N=len(ids),rna_index=i,
                                                                             null_id=null_id,comparison=mode,family=family,**t))
                        changed |= (v is None)!=(old is None) or (v is not None and old is not None and not np.array_equal(v,old))
                    if changed and first is None:
                        first={'scenario':sid,'case':spec['name'],'comparison':mode,
                               'details':witness(baseline,cur,ids,private/'first_witness_ranks.npz')}
                    scenario_rows.append(row);axis.append(sid)
                    states.append(np.stack([cur['family'][f] if cur['family'][f] is not None else np.full(len(ids),-1,dtype=np.int8) for f in FAMILIES]))
                np.savez_compressed(private/(spec['name']+f'_rna{i:02d}.npz'),ids=ids,scenario=np.asarray(axis),states=np.asarray(states,dtype=np.int8))
        if a.stage in ('both','rarity'):
            for item in spec['sampling']:
                null_id=item['id'];sampler=GateSampler(ids,neff,rna,maxima[null_id]);record=dict(case=spec['name'],N=len(ids),null_id=null_id,seed=item['seed'],active_events=len(sampler.active))
                if not sampler.active:record.update(status='EXACT_IMPOSSIBLE',actual_proposals=0,accepted=None,P_G=0.,MC_SE=0.,MC_95_halfwidth=0.)
                elif sampler.certain:record.update(status='EXACT_CERTAIN',actual_proposals=0,accepted=None,P_G=1.,MC_SE=0.,MC_95_halfwidth=0.)
                else:
                    n=c['conditional_proposals'];rng=np.random.default_rng(item['seed']);multiplicity=np.empty(n,dtype=np.int16)
                    accepted=np.empty(n,dtype=bool);counts=np.empty((n,10),dtype=np.int32)
                    for j in range(n):
                        rank,counts[j],accepted[j]=sampler.proposal(rng);multiplicity[j]=int((counts[j]>sampler.maxima).sum())
                        if multiplicity[j]<1:raise ValueError('conditional proposal did not pass')
                    weights=1./multiplicity.astype(float);scale=float(sampler.S);se=scale*float(weights.std(ddof=1))/math.sqrt(n) if n>1 else None
                    np.savez_compressed(private/(spec['name']+'_'+null_id+'_draws.npz'),counts=counts,accepted=accepted,multiplicity=multiplicity)
                    record.update(status='SAMPLED',actual_proposals=n,accepted=int(accepted.sum()),P_G=scale*float(weights.mean()),
                                  MC_SE=se,MC_95_halfwidth=1.96*se if se is not None else None)
                rarity_rows.append(record)
    if scenario_rows:pd.DataFrame(scenario_rows).to_csv(a.output/'scenarios.csv',index=False)
    if transition_rows:pd.DataFrame(transition_rows).to_csv(a.output/'transitions.csv',index=False)
    if rarity_rows:
        frame=pd.DataFrame(rarity_rows)
        frame['delta_from_original']=np.nan;frame['combined_MC_95_halfwidth']=np.nan
        frame['difference_exceeds_combined_MC_95']=pd.Series([pd.NA]*len(frame),dtype='boolean')
        for case_name,group in frame.groupby('case'):
            original=group[group.null_id=='original'].iloc[0]
            for idx,row in group[group.null_id!='original'].iterrows():
                delta=float(row.P_G-original.P_G);frame.loc[idx,'delta_from_original']=delta
                if pd.notna(row.MC_95_halfwidth) and pd.notna(original.MC_95_halfwidth):
                    error=math.hypot(float(row.MC_95_halfwidth),float(original.MC_95_halfwidth))
                    frame.loc[idx,'combined_MC_95_halfwidth']=error
                    frame.loc[idx,'difference_exceeds_combined_MC_95']=abs(delta)>error
        frame.to_csv(a.output/'rarity.csv',index=False)
    (a.output/'first_change.json').write_text(json.dumps(first if first else {'status':'NO_SUBSTANTIVE_CHANGE'},indent=2)+'\n')
    receipt={'status':'COMPLETE','scenarios':len(scenario_rows),'family_transitions':len(transition_rows),'rarity_streams':len(rarity_rows),
             'scope':'member files contain private IDs; finite input scenarios only'}
    (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))

if __name__=='__main__':main()
