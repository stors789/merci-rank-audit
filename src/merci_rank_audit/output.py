"""Run the explicit-input audit and write its small, inspectable outputs."""
import hashlib
import json
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from . import __version__
from .core import analyze, exact_number, member_witness
from .inputs import require_identity, validate_manifest
from .sampling import sample_case


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def methods_for(mode):
    return {'none':[],'uniform':['uniform'],'conditional':['conditional'],
            'both':['uniform','conditional']}[mode]


def sampling_contract():
    package=Path(__file__).resolve().parent
    components=['bounds.py','sampler.py','inputs.py','fwl.py','sampling.py','output.py']
    return {'software_version':__version__,
            'component_sha256':{name:sha(package/name) for name in components},
            'runtime_versions':{'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__}}


def required_files(methods,nuclear):
    names=[]
    for method in methods:
        names.extend([f'{method}_members.tsv',f'{method}_prefixes.tsv'])
        if nuclear:names.append(f'{method}_nuclear.tsv')
    return names


def load_reuse(root,case,budget,methods):
    source=root/'cases'/case['name']/'samples'
    metadata_path=source/'metadata.json'
    if not metadata_path.is_file():raise ValueError('Reusable sample metadata is missing')
    meta=json.loads(metadata_path.read_text())
    require_identity(case,meta,require_nuclear=case['nuclear'] is not None)
    if meta.get('schema')!='merci-rank-audit-sample-v1' or not isinstance(meta.get('results'),dict):
        raise ValueError('Reusable sample schema is unsupported')
    if meta.get('sampling_contract')!=sampling_contract():
        raise ValueError('Reusable sample software or runtime identity differs or is missing')
    seeds={m:case['seed_'+m] for m in methods}
    if meta.get('budget')!=budget or meta.get('methods')!=methods or meta.get('seeds')!=seeds:
        raise ValueError('Reusable sample budget, methods, or seeds differ')
    names=required_files(methods,case['nuclear'] is not None)
    if set(meta.get('file_sha256',{}))!=set(names):
        raise ValueError('Reusable sample file list differs')
    for name in names:
        p=source/name
        if not p.is_file() or sha(p)!=meta['file_sha256'][name]:
            raise ValueError('Reusable sample output missing or changed: '+name)
    values={}
    for method in methods:
        frame=pd.read_csv(source/f'{method}_members.tsv',sep='\t',dtype={'cell':str})
        if frame.cell.duplicated().any() or set(frame.cell)!=set(case['sampler'].ids):
            raise ValueError('Reusable sample IDs differ')
        frame=frame.set_index('cell').loc[case['sampler'].ids]
        accepted=meta['results'][method]['accepted']
        counts=frame.selected_count.to_numpy(float)
        if accepted:
            if not np.isfinite(counts).all() or np.any(counts!=np.floor(counts)) or np.any((counts<0)|(counts>accepted)):
                raise ValueError('Reusable sample counts invalid')
            counts=counts.astype(np.int64)
        elif not np.isnan(counts).all():
            raise ValueError('Zero accepted orders require undefined counts')
        else:counts=None
        values[method]=(meta['results'][method],counts)
    return source,meta,values


def write_sample(case,out,budget,methods,union):
    out.mkdir()
    results={};values={}
    prefixes=sorted(set([x for x in (100,1000,5000,20000,budget) if x<=budget]))
    for method in methods:
        seed=case['seed_'+method]
        meta,counts,rows,nuclear=sample_case(case,method,budget,seed,prefixes,union)
        results[method]=meta;values[method]=(meta,counts)
        members=pd.DataFrame({'cell':case['sampler'].ids,'accepted_orders':meta['accepted'],
                              'selected_count':counts if counts is not None else [None]*case['sampler'].n})
        members.to_csv(out/f'{method}_members.tsv',sep='\t',index=False)
        pd.DataFrame(rows,columns=['method','budget','accepted','seen_members','exact_passing_union',
                                  'union_unseen','frequency_status','seconds']).to_csv(out/f'{method}_prefixes.tsv',sep='\t',index=False)
        if case['nuclear'] is not None:
            columns=['method','representation','feature','accepted','valid_betas','invalid_betas','quantile_status',
                     'q025','q025_MC_low','q025_MC_high','median','median_MC_low','median_MC_high',
                     'q975','q975_MC_low','q975_MC_high']
            pd.DataFrame(nuclear,columns=columns).to_csv(out/f'{method}_nuclear.tsv',sep='\t',index=False)
    names=required_files(methods,case['nuclear'] is not None)
    meta={'schema':'merci-rank-audit-sample-v1','input_identity':case['input_identity'],
          'sampling_contract':sampling_contract(),'budget':budget,
          'methods':methods,'seeds':{m:case['seed_'+m] for m in methods},'results':results,
          'file_sha256':{name:sha(out/name) for name in names}}
    (out/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta,values


def copy_reuse(source,meta,out):
    out.mkdir()
    for name in meta['file_sha256']:shutil.copyfile(source/name,out/name)
    shutil.copyfile(source/'metadata.json',out/'metadata.json')


def run_audit(manifest,output,sample_mode='none',case_name=None,witness=None,reuse_samples=None):
    start=time.monotonic();manifest=Path(manifest).resolve();out=Path(output).resolve()
    budget,cases=validate_manifest(manifest)
    if case_name is not None:
        cases=[case for case in cases if case['name']==case_name]
        if len(cases)!=1:raise ValueError('Requested case does not exist')
    if witness is not None and len(cases)!=1:raise ValueError('A member witness requires exactly one case')
    if witness is not None and witness[0] not in cases[0]['sampler'].ids:
        raise ValueError('Requested member ID is not in the case')
    methods=methods_for(sample_mode)
    if reuse_samples is not None and not methods:raise ValueError('Reusing samples requires a sampling mode')
    reuse={}
    if reuse_samples is not None:
        source=Path(reuse_samples).resolve()
        for case in cases:reuse[case['name']]=load_reuse(source,case,budget,methods)
    out.mkdir(parents=True,exist_ok=False)
    (out/'cases').mkdir()
    summaries=[]
    for case in cases:
        case_out=out/'cases'/case['name'];case_out.mkdir()
        sampler=case['sampler'];context,records,analytic,summary=analyze(case)
        union=np.array([r['conditioned_union'] for r in records],dtype=bool) if context['possible'] else None
        samples={}
        if methods:
            if reuse_samples is not None:
                src,meta,samples=reuse[case['name']]
                copy_reuse(src,meta,case_out/'samples')
            else:
                meta,samples=write_sample(case,case_out/'samples',budget,methods,union)
        if analytic['kind']=='multi_gate_reference_required' and samples:
            method=next((m for m in ('conditional','uniform') if m in samples),None)
            summary['probability_method']='MC_'+method
            summary['gate_probability']=samples[method][0]['gate_probability']
        rows=[];probability_values={};probability_keys={}
        for i,row in enumerate(records):
            item=dict(row)
            if analytic['members'] is not None:
                value=analytic['members'][i];identity=id(value)
                if identity not in probability_keys:
                    key=f'exact_{len(probability_values):03d}'
                    probability_keys[identity]=key;probability_values[key]=exact_number(value)
                key=probability_keys[identity];number=probability_values[key]
                short=(len(number['numerator'])+len(number['denominator']))<=256
                item.update(probability_method=analytic['kind'],probability_status='defined',
                            probability_decimal=number['decimal'],probability_numerator=number['numerator'],
                            probability_denominator=number['denominator'],probability_integer_base=number['integer_base'],
                            probability_exact_id=key,sample_accepted=None)
                if not short:item['probability_numerator']=item['probability_denominator']=None
            elif not context['possible']:
                item.update(probability_method='undefined_G_impossible',probability_status='undefined',
                            probability_decimal=None,probability_numerator=None,probability_denominator=None,
                            probability_integer_base=None,probability_exact_id=None,sample_accepted=None)
            else:
                method=next((m for m in ('conditional','uniform') if m in samples),None)
                m,count=samples[method] if method else (None,None)
                item.update(probability_method='MC_'+method if method else 'sampling_required',
                            probability_status='sample_estimate' if count is not None else ('NA_no_accepted_order' if method else 'not_requested'),
                            probability_decimal=str(count[i]/m['accepted']) if count is not None else None,
                            probability_numerator=str(count[i]) if count is not None else None,
                            probability_denominator=str(m['accepted']) if count is not None else None,
                            probability_integer_base=10 if count is not None else None,probability_exact_id=None,
                            sample_accepted=m['accepted'] if m else None)
            rows.append(item)
        pd.DataFrame(rows).to_csv(case_out/'members.tsv',sep='\t',index=False)
        if probability_values:
            (case_out/'probability_values.json').write_text(json.dumps(probability_values,indent=2)+'\n')
        if witness is not None:
            certificate=member_witness(context,witness[0],witness[1])
            rank=certificate.pop('rank')
            if rank is not None:
                pd.DataFrame({'cell':sampler.ids,'DNA_rank':rank}).to_csv(case_out/'witness.tsv',sep='\t',index=False)
            (case_out/'witness.json').write_text(json.dumps(certificate,indent=2)+'\n')
        nuclear_status=('not_supplied' if case['nuclear'] is None else 'not_requested' if not methods else
                        'computed' if any(value[0]['accepted'] for value in samples.values()) else
                        'not_evaluable_condition_impossible' if not context['possible'] else 'not_evaluable_no_accepted')
        summary.update({'case':case['name'],'N':sampler.n,'input_identity':case['input_identity'],
                        'nuclear_status':nuclear_status,
                        'samples':{key:value[0] for key,value in samples.items()}})
        (case_out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        summaries.append(summary)
    run={'version':__version__,'manifest_sha256':sha(manifest),'budget':budget,'sample_mode':sample_mode,
         'sample_source':str(Path(reuse_samples).resolve()) if reuse_samples is not None else 'new_or_not_requested',
         'cases':summaries,'elapsed_seconds':time.monotonic()-start}
    (out/'run.json').write_text(json.dumps(run,indent=2)+'\n')
    lines=['# Fixed-score audit','']
    for row in summaries:
        original=row['original_gate_status'] if row['original_gate_pass'] is None else str(row['original_gate_pass'])
        passing=('undefined (G impossible)' if row['passing_core'] is None else
                 f"{row['passing_core']}/{row['passing_union']} (exact)")
        lines.extend([f"## {row['case']}",f"- N: {row['N']}; original gate: {original}; family gate: {row['gate_state']}.",
                      f"- Raw core/union: {row['raw_core']}/{row['raw_union']}; failed-as-empty: {row['default_core']}/{row['default_union']}.",
                      f"- Passing-only core/union: {passing}.",
                      f"- Member probability method: {row['probability_method']}; nuclear: {row['nuclear_status']}.",''])
    lines.append('Conditional results describe the specified ranking mechanism, not physical transfer or a patient-level probability.')
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return run
