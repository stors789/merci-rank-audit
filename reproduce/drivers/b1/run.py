"""Rebuild a frozen B1 strict split from four explicit processed source files."""
import argparse
import csv
import hashlib
from importlib.resources import as_file, files
import json
import os
from pathlib import Path
import subprocess
import sys

UPSTREAM_COMMIT='79e105aaf13159bee325c3ceceb6ee679e430341'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        while chunk:=stream.read(1024*1024):h.update(chunk)
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('info','expression','variants','coverage','upstream','r-library','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--rscript',default='Rscript')
    parser.add_argument('--split',type=int,required=True)
    parser.add_argument('--config',type=Path,default=Path(__file__).resolve().parents[2]/'configs/cohort_formal_portable.json')
    args=parser.parse_args()
    if not 21<=args.split<=100:raise ValueError('B1 frozen split must be in 21..100')
    for path in (args.info,args.expression,args.variants,args.coverage):
        if not path.is_file():raise FileNotFoundError(path)
    upstream=args.upstream.resolve();library=args.r_library.resolve()
    if not library.is_dir():raise FileNotFoundError(library)
    revision=subprocess.check_output(['git','-C',str(upstream),'rev-parse','HEAD'],text=True).strip()
    if revision!=UPSTREAM_COMMIT:raise ValueError('Upstream checkout does not match the pinned MERCI commit')
    subprocess.run(['git','-C',str(upstream),'diff','--quiet','HEAD','--',
                    'MERCI/R/MERCI_MT_est.r','MERCI/R/MERCI_ReceiverPre.r'],check=True)
    config=json.loads(args.config.read_text());original=config['cases']
    cases=[dict(x) for x in original if x['experiment']=='bench1' and x['split']==args.split
           and x['case'] in ('strict_regular_mixed','strict_regular_negative')]
    splits=[dict(x) for x in config['splits'] if x['experiment']=='bench1' and x['split']==args.split]
    if len(cases)!=2 or len(splits)!=1:raise ValueError('Frozen split or two strict cases are absent from the config')
    output=args.output.resolve();output.mkdir(parents=True,mode=0o700,exist_ok=False)
    generated=output/'generated';generated.mkdir(mode=0o700)
    maxima=output/'original_maxima.csv'
    with maxima.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['case','N','cutoff','maximum'],lineterminator='\n');writer.writeheader()
        for case in cases:
            if len(case['maxima'])!=10:raise ValueError('Expected ten original maxima')
            for cutoff,value in zip(range(100,0,-10),case['maxima']):
                writer.writerow({'case':case['case'],'N':case['N'],'cutoff':cutoff,'maximum':value})
    env=os.environ.copy()
    env.update({'R_LIBS_USER':str(library),'R_LIBS_SITE':'','R_LIBS':'',
                'R_PROFILE_USER':os.devnull,'R_ENVIRON_USER':os.devnull,
                'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1',
                'NUMEXPR_NUM_THREADS':'1','NUMBA_NUM_THREADS':'1'})
    here=Path(__file__).resolve().parent
    with as_file(files('merci_rank_audit').joinpath('r_helpers','engine.R')) as engine:
        command=[args.rscript,'--vanilla',str(here/'score_strict.R'),
                 *(str(x.resolve()) for x in (args.info,args.expression,args.variants,args.coverage)),
                 str(upstream),str(engine),str(here/'read_coverage.R'),str(generated),str(args.split),str(maxima),str(library)]
        subprocess.run(command,check=True,env=env)
    # Bind newly produced files to new byte hashes. Historical hashes remain separate provenance.
    for case in cases:
        case['historical_score_sha256']=case['score_sha256']
        case['source']=str(generated/(case['case']+'_cells.tsv.gz'))
        case['score_sha256']=sha(Path(case['source']))
    split=splits[0];split['historical_prepared_sha256']=split['prepared_sha256']
    split['prepared']=str(output/'prepared.npz')
    null=(args.config.resolve().parent/config['null_maxima']).resolve()
    binding={'status':'NEW_FILES_FROM_EXPLICIT_PUBLIC_INPUTS','cases':cases,'splits':splits,'null_maxima':str(null),
             'upstream_commit':revision,'original_config':str(args.config.resolve())}
    binding_path=output/'binding_before_prepare.json';binding_path.write_text(json.dumps(binding,indent=2)+'\n')
    cohort=here.parent/'cohort';manifest=output/'cohort_manifest.json'
    subprocess.run([sys.executable,str(cohort/'prepare_split.py'),'--manifest',str(binding_path),
                    '--experiment','bench1','--split',str(args.split),
                    '--mixed',str(generated/'strict_regular_mixed_cells.tsv.gz'),
                    '--negative',str(generated/'strict_regular_negative_cells.tsv.gz'),
                    '--output',str(output/'prepared.npz'),'--output-manifest',str(manifest)],check=True,env=env)
    subprocess.run([sys.executable,str(cohort/'run.py'),'--manifest',str(manifest),
                    '--output',str(output/'cohort')],check=True,env=env)
    prepared_sha=sha(output/'prepared.npz')
    score_match={x['case']:x['score_sha256']==x['historical_score_sha256'] for x in cases}
    prepared_match=prepared_sha==split['historical_prepared_sha256']
    result={'status':'COMPLETE','split':args.split,'upstream_commit':revision,
            'generated_score_sha256':{x['case']:x['score_sha256'] for x in cases},
            'historical_score_sha256':{x['case']:x['historical_score_sha256'] for x in cases},
            'prepared_sha256':prepared_sha,'historical_prepared_sha256':split['historical_prepared_sha256'],
            'historical_byte_match':{'scores':score_match,'prepared':prepared_match},
            'historical_identity_status':('BYTE_IDENTICAL' if all(score_match.values()) and prepared_match
                                          else 'DIFFERENT_BYTES_REQUIRES_COMPARISON'),
            'cohort_manifest_sha256':sha(manifest),
            'scope':'new files only; cell and member outputs remain private'}
    (output/'run_complete.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)

if __name__=='__main__':main()
