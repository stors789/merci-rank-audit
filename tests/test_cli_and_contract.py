import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from merci_rank_audit.cli import main
from merci_rank_audit.output import run_audit

EXAMPLES=Path(__file__).resolve().parents[1]/'src/merci_rank_audit/examples'


def prepare(tmp_path,budget=20,nuclear=True,case='synthetic_certain'):
    inputs=tmp_path/'inputs';shutil.copytree(EXAMPLES,inputs)
    manifest={'budget':budget,'cases':[{'name':case,'score':'score.tsv',
               'thresholds':'certain_thresholds.csv','seed_uniform':21,'seed_conditional':22}]}
    if nuclear:manifest['cases'][0]['nuclear']='nuclear.npz'
    (inputs/'test_manifest.json').write_text(json.dumps(manifest)+'\n')
    return inputs/'test_manifest.json'


def test_exact_only_no_nuclear_and_missing_original_rank(tmp_path):
    manifest=prepare(tmp_path,nuclear=False)
    score=manifest.parent/'score.tsv'
    d=pd.read_csv(score,sep='\t').drop(columns=['DNA_rank']);d.to_csv(score,sep='\t',index=False)
    out=tmp_path/'out';run=run_audit(manifest,out)
    case=run['cases'][0]
    assert case['original_gate_pass'] is None and case['original_gate_status']=='not_evaluated_missing_DNA_rank'
    assert case['nuclear_status']=='not_supplied' and not (out/'cases/synthetic_certain/samples').exists()
    members=pd.read_csv(out/'cases/synthetic_certain/members.tsv',sep='\t')
    values=json.loads((out/'cases/synthetic_certain/probability_values.json').read_text())
    assert members.probability_exact_id.notna().all()
    assert all(key in values for key in members.probability_exact_id)
    with pytest.raises(FileExistsError):run_audit(manifest,out)


def test_witness_and_low_valid_nuclear_quantile(tmp_path):
    manifest=prepare(tmp_path,budget=20)
    out=tmp_path/'out'
    run=run_audit(manifest,out,'conditional',witness=('synthetic_19',True))
    assert run['cases'][0]['samples']['conditional']['gate_probability']['method']=='exact_certain'
    assert run['cases'][0]['samples']['conditional']['gate_probability']['estimate']==1.0
    w=json.loads((out/'cases/synthetic_certain/witness.json').read_text())
    assert w['reachable'] is True and (out/'cases/synthetic_certain/witness.tsv').is_file()
    nuclear=pd.read_csv(out/'cases/synthetic_certain/samples/conditional_nuclear.tsv',sep='\t')
    assert len(nuclear)==4 and nuclear.quantile_status.eq('not_reported_under_500_valid').all()
    assert nuclear[['q025','median','q975']].isna().all().all()


def test_reuse_rejects_missing_changed_and_unbound_outputs(tmp_path):
    manifest=prepare(tmp_path,budget=20,nuclear=False)
    first=tmp_path/'first';run_audit(manifest,first,'conditional')
    reuse=tmp_path/'reuse';run_audit(manifest,reuse,'conditional',reuse_samples=first)
    metadata=first/'cases/synthetic_certain/samples/metadata.json'
    saved_metadata=metadata.read_bytes();contract=json.loads(saved_metadata)
    contract['sampling_contract']['software_version']='different'
    metadata.write_text(json.dumps(contract))
    with pytest.raises(ValueError,match='software or runtime identity'):run_audit(manifest,tmp_path/'changed_code','conditional',reuse_samples=first)
    assert not (tmp_path/'changed_code').exists()
    metadata.write_bytes(saved_metadata)
    score=manifest.parent/'score.tsv';saved_score=score.read_bytes();score.write_bytes(saved_score+b'\n')
    with pytest.raises(ValueError,match='saved sample input identity'):run_audit(manifest,tmp_path/'changed_input','conditional',reuse_samples=first)
    assert not (tmp_path/'changed_input').exists()
    score.write_bytes(saved_score)
    source=first/'cases/synthetic_certain/samples/conditional_members.tsv'
    original=source.read_bytes();source.write_bytes(original+b'\n')
    with pytest.raises(ValueError,match='missing or changed'):run_audit(manifest,tmp_path/'changed','conditional',reuse_samples=first)
    assert not (tmp_path/'changed').exists()
    source.write_bytes(original);source.unlink()
    with pytest.raises(ValueError,match='missing or changed'):run_audit(manifest,tmp_path/'missing','conditional',reuse_samples=first)
    meta=first/'cases/synthetic_certain/samples/metadata.json'
    data=json.loads(meta.read_text());data.pop('input_identity');meta.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='no supported input identity'):run_audit(manifest,tmp_path/'unbound','conditional',reuse_samples=first)


def test_impossible_gate_and_zero_sample_are_not_zero_probability(tmp_path):
    manifest=prepare(tmp_path,budget=20,nuclear=False)
    d=json.loads(manifest.read_text());d['cases'][0]['thresholds']='impossible_thresholds.csv'
    manifest.write_text(json.dumps(d))
    out=tmp_path/'out';result=run_audit(manifest,out,'both')
    row=result['cases'][0]
    assert row['gate_state']=='impossible' and row['passing_core'] is None and row['passing_union'] is None
    members=pd.read_csv(out/'cases/synthetic_certain/samples/uniform_members.tsv',sep='\t')
    assert members.accepted_orders.eq(0).all() and members.selected_count.isna().all()
    table=pd.read_csv(out/'cases/synthetic_certain/members.tsv',sep='\t')
    assert table.probability_decimal.isna().all() and table.probability_status.eq('undefined').all()


def test_possible_multigate_zero_uniform_passes_have_undefined_conditional_estimate(tmp_path):
    inputs=tmp_path/'inputs';shutil.copytree(EXAMPLES,inputs)
    manifest={'budget':1,'cases':[{'name':'multi','score':'all_tied_overlapping_score.tsv',
              'thresholds':'all_tied_overlapping_thresholds.csv','seed_uniform':1,'seed_conditional':2}]}
    path=inputs/'manifest_test.json';path.write_text(json.dumps(manifest))
    out=tmp_path/'out';run=run_audit(path,out,'uniform')
    assert run['cases'][0]['gate_state']=='possible_not_certain' and run['cases'][0]['passing_union']>0
    members=pd.read_csv(out/'cases/multi/members.tsv',sep='\t')
    assert members.probability_status.eq('NA_no_accepted_order').all()
    assert members.probability_decimal.isna().all()
    prefixes=pd.read_csv(out/'cases/multi/samples/uniform_prefixes.tsv',sep='\t')
    assert int(prefixes.iloc[0].seen_members)==0 and int(prefixes.iloc[0].union_unseen)>0


def test_constant_outcome_and_degenerate_label_are_not_valid_beta(tmp_path):
    manifest=prepare(tmp_path,budget=20)
    path=manifest.parent/'nuclear.npz'
    with np.load(path,allow_pickle=False) as a:data={k:a[k] for k in a.files}
    data['original'][:,0]=1.;np.savez_compressed(path,**data)
    out=tmp_path/'out';run_audit(manifest,out,'conditional')
    nuclear=pd.read_csv(out/'cases/synthetic_certain/samples/conditional_nuclear.tsv',sep='\t')
    assert int(nuclear[(nuclear.representation=='original')&(nuclear.feature==data['names'][0])].iloc[0].valid_betas)==0
    from merci_rank_audit.fwl import beta_batch
    x=np.ones((20,1));y=np.arange(20.)[:,None];z=np.zeros((1,20))
    assert np.isnan(beta_batch(y,x,z)).all()


@pytest.mark.parametrize('problem',['duplicate_id','duplicate_rna','wrong_N','nuclear_id'])
def test_invalid_domain_is_rejected_before_output(tmp_path,problem):
    manifest=prepare(tmp_path,budget=10)
    score=manifest.parent/'score.tsv'
    if problem in ('duplicate_id','duplicate_rna'):
        d=pd.read_csv(score,sep='\t')
        column='cell' if problem=='duplicate_id' else 'RNA_rank'
        d.loc[1,column]=d.loc[0,column];d.to_csv(score,sep='\t',index=False)
    elif problem=='wrong_N':
        path=manifest.parent/'certain_thresholds.csv'
        d=pd.read_csv(path);d.loc[0,'N']+=1;d.to_csv(path,index=False)
    else:
        path=manifest.parent/'nuclear.npz'
        with np.load(path,allow_pickle=False) as d:a={k:d[k] for k in d.files}
        a['ids'][0]='missing';np.savez_compressed(path,**a)
    with pytest.raises(ValueError):run_audit(manifest,tmp_path/'not_created')
    assert not (tmp_path/'not_created').exists()


def test_nuclear_quantiles_defined_when_valid_count_is_sufficient(tmp_path):
    manifest=prepare(tmp_path,budget=1000)
    out=tmp_path/'out';run_audit(manifest,out,'conditional')
    d=pd.read_csv(out/'cases/synthetic_certain/samples/conditional_nuclear.tsv',sep='\t')
    assert d.valid_betas.ge(500).all() and d.quantile_status.eq('defined').all()
    assert d[['q025','median','q975']].notna().all().all()


def test_long_exact_member_rationals_are_shared_not_repeated(tmp_path):
    n=2000;inputs=tmp_path/'inputs';inputs.mkdir()
    score=pd.DataFrame({'cell':[f'id{i:04d}' for i in range(n)],'Neff':np.zeros(n,dtype=int),'RNA_rank':np.arange(1,n+1)})
    score.to_csv(inputs/'score.tsv',sep='\t',index=False)
    cutoffs=list(range(100,0,-10));k=[(100+(n-1)*c+99)//100 for c in cutoffs]
    maxima=k.copy();maxima[5]=700
    pd.DataFrame({'N':n,'cutoff':cutoffs,'maximum':maxima}).to_csv(inputs/'thresholds.csv',index=False)
    manifest={'budget':1,'cases':[{'name':'long','score':'score.tsv','thresholds':'thresholds.csv',
                                   'seed_uniform':1,'seed_conditional':2}]}
    path=inputs/'manifest.json';path.write_text(json.dumps(manifest))
    out=tmp_path/'out';run_audit(path,out)
    members=pd.read_csv(out/'cases/long/members.tsv',sep='\t')
    values=json.loads((out/'cases/long/probability_values.json').read_text())
    long={key:v for key,v in values.items() if len(v['numerator'])+len(v['denominator'])>256}
    assert long and len(values)<10 and all(int(v['denominator'],v['integer_base'])>0 for v in long.values())
    assert members.probability_exact_id.notna().all()
    assert members.loc[members.probability_exact_id.isin(long),'probability_denominator'].isna().all()
    assert (out/'cases/long/members.tsv').stat().st_size<2_000_000


def test_cli_rejects_partial_witness_and_missing_inputs(tmp_path):
    manifest=prepare(tmp_path,nuclear=False)
    assert main(['audit','--manifest',str(manifest),'--output',str(tmp_path/'x'),'--member','synthetic_00'])==2
    assert not (tmp_path/'x').exists()
    assert main(['bench2','--info','no','--expression','no','--variants','no','--coverage','no',
                 '--upstream','no','--r-library','no','--output',str(tmp_path/'y')])==2
    assert not (tmp_path/'y').exists()
