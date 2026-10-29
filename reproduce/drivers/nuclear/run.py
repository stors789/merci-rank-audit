"""Run a frozen P2 conditional-proposal stream with exact input identity."""
import argparse
import csv
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from merci_rank_audit.fwl import projection
from merci_rank_audit.inputs import percentile_with_MC_interval, validate_manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--stream', choices=['reference', 'repeat'], required=True)
    p.add_argument('--repeat', type=int, default=1)
    p.add_argument('--all-outcomes', action='store_true', help='Compute both representations for every supplied outcome')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); config = json.loads(a.config.read_text())
    _, cases = validate_manifest(a.manifest)
    if len(cases) != 1 or cases[0]['nuclear'] is None:
        raise ValueError('one score/threshold/nuclear case is required')
    case = cases[0]; sampler = case['sampler']; data = case['nuclear']
    if config.get('input_identity') is not None and config['input_identity'] != case['input_identity']:
        raise ValueError('input bytes differ from the frozen P2 identities')
    if not sampler.active:
        raise ValueError('G impossible: conditional nuclear distribution undefined')
    if a.stream == 'reference':
        seed = config['reference_seed_entropy']; budgets = [config['reference_proposals']]
        limit = config['reference_seconds_limit']
    else:
        if not 1 <= a.repeat <= len(config['repeat_seed_entropy']):
            raise ValueError('repeat number outside the frozen range')
        seed = config['repeat_seed_entropy'][a.repeat - 1]; budgets = config['prefixes']
        limit = config['repeat_seconds_limit']
    series = ([(rep, str(name)) for rep in ('original','alternative') for name in data['names']]
              if a.all_outcomes else config['series'])
    columns = [(rep, int(np.flatnonzero(data['names'] == name)[0])) for rep, name in series]
    y = np.column_stack([data[rep][:, index] for rep, index in columns]); q, yr = projection(data['X'], y)
    valid_y = np.isfinite(y).all(axis=0) & (y.std(axis=0, ddof=1) > 1e-14)
    rng = np.random.default_rng(np.random.SeedSequence(seed)); k = int(sampler.ks[5])
    a.output.mkdir(parents=True, exist_ok=False)
    numbers=[]; betas=[]; denominators=[]; total=[]; observed=[]; first_rank=[]; first_counts=[]
    proposals=accepted=0; start=time.monotonic(); prefixes=[]; stop=None
    for budget in budgets:
        while proposals < budget:
            if time.monotonic()-start >= limit:
                stop = 'TIME_LIMIT'; break
            labels=[]; batch_numbers=[]; batch_total=[]; batch_obs=[]
            for _ in range(min(config['batch_proposals'],budget-proposals)):
                rank, counts, yes = sampler.proposal(rng); proposals += 1
                if yes:
                    full = (rank > sampler.n-k) & (sampler.rna > sampler.n-k)
                    z = full[data['index']][data['observable']]; accepted += 1
                    labels.append(z); batch_numbers.append(proposals); batch_total.append(int(full.sum())); batch_obs.append(int(z.sum()))
                    if len(first_rank) < 5: first_rank.append(rank.copy()); first_counts.append(counts.copy())
            if labels:
                z=np.asarray(labels,dtype=float); zr=z-(z@q)@q.T; den=np.sum(zr*zr,axis=1); selected=z.sum(axis=1)
                valid=(selected>=5)&((z.shape[1]-selected)>=5)&(den>1e-10)
                values=np.full((len(z),len(series)),np.nan); values[valid]=(zr[valid]@yr)/den[valid,None]; values[:,~valid_y]=np.nan
                betas.extend(values); denominators.extend(den); numbers.extend(batch_numbers); total.extend(batch_total); observed.extend(batch_obs)
        values=np.asarray(betas,dtype=float).reshape(-1,len(series))
        prefixes.append({'budget':budget,'status':'COMPLETE' if proposals==budget else stop,
                         'proposals':proposals,'accepted':accepted,'valid_betas':np.isfinite(values).sum(axis=0).tolist()})
        if stop: break
    values=np.asarray(betas,dtype=float).reshape(-1,len(series)); indices=np.asarray(numbers,dtype=np.int32)
    for budget in budgets:
        if not any(x['budget']==budget for x in prefixes):
            prefixes.append({'budget':budget,'status':stop or 'NOT_REACHED','proposals':proposals,
                             'accepted':accepted,'valid_betas':np.isfinite(values).sum(axis=0).tolist()})
    sample=a.output/'samples.npz'
    np.savez_compressed(sample,proposal_index=indices,betas=values,denominator=np.asarray(denominators),
                        selected_total=np.asarray(total,dtype=np.int32),selected_observed=np.asarray(observed,dtype=np.int32),
                        first_ranks=np.asarray(first_rank,dtype=np.int32).reshape(-1,sampler.n),
                        first_counts=np.asarray(first_counts,dtype=np.int32).reshape(-1,10))
    rows=[]
    for prefix in prefixes:
        for j,(rep,name) in enumerate(series):
            x=values[indices<=prefix['budget'],j]; x=x[np.isfinite(x)]
            for probability in config['quantiles']:
                defined=prefix['status']=='COMPLETE' and len(x)>=config['min_valid_beta']
                v=percentile_with_MC_interval(x,probability) if defined else (None,None,None)
                rows.append(dict(budget=prefix['budget'],representation=rep,feature=name,quantile=probability,
                                 valid_betas=len(x),status='DEFINED' if defined else ('LOW_VALID_N' if prefix['status']=='COMPLETE' else prefix['status']),
                                 estimate=v[0],MC95_low=v[1],MC95_high=v[2]))
    with (a.output/'quantiles.csv').open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0],lineterminator='\n');w.writeheader();w.writerows(rows)
    receipt={'status':'COMPLETE' if not stop else stop,'stream':a.stream,'seed_entropy':seed,'proposals':proposals,
             'accepted':accepted,'series':series,'prefixes':prefixes,'input_identity':case['input_identity'],
             'samples_sha256':hashlib.sha256(sample.read_bytes()).hexdigest(),'seconds':time.monotonic()-start}
    (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'proposals':proposals,'accepted':accepted}))

if __name__=='__main__': main()
