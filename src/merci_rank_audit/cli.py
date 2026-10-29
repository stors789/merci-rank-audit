"""Command-line interface for installed fixed-score audits."""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from importlib.resources import as_file, files
from pathlib import Path

from .output import run_audit

UPSTREAM_COMMIT='79e105aaf13159bee325c3ceceb6ee679e430341'


def build_parser():
    parser=argparse.ArgumentParser(prog='merci-rank-audit',description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    example=sub.add_parser('init-example',help='Copy the small synthetic quickstart inputs')
    example.add_argument('--output',required=True,type=Path)
    audit=sub.add_parser('audit',help='Exact families; optional sampling and member witness')
    audit.add_argument('--manifest',required=True,type=Path)
    audit.add_argument('--output',required=True,type=Path)
    audit.add_argument('--case')
    audit.add_argument('--sample',choices=['none','uniform','conditional','both'],default='none')
    audit.add_argument('--reuse-samples',type=Path)
    audit.add_argument('--member')
    audit.add_argument('--decision',choices=['selected','omitted'])
    bench=sub.add_parser('bench2',help='Run the explicit external R preparation, then audit its output')
    for name in ('info','expression','variants','coverage','upstream','r-library','output'):
        bench.add_argument('--'+name,required=True,type=Path)
    bench.add_argument('--rscript',default='Rscript')
    bench.add_argument('--sample',choices=['none','uniform','conditional','both'],default='none')
    return parser


def run_bench2(args):
    inputs=[args.info,args.expression,args.variants,args.coverage]
    if not all(x.is_file() for x in inputs):raise FileNotFoundError('All four processed input files are required')
    upstream=args.upstream.resolve()
    revision=subprocess.check_output(['git','-C',str(upstream),'rev-parse','HEAD'],text=True).strip()
    if revision!=UPSTREAM_COMMIT:raise ValueError('Upstream checkout must match the pinned MERCI commit')
    r_library=args.r_library.resolve()
    if not r_library.is_dir():raise ValueError('An explicit isolated R user library is required')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    run={'status':'started','case':'bench2_split1_original_regular_negative','upstream_commit':revision,
         'input_paths':[str(p.resolve()) for p in inputs],'Rscript':args.rscript,'R_user_library':str(r_library)}
    (out/'run_started.json').write_text(json.dumps(run,indent=2)+'\n')
    env=os.environ.copy()
    env.update({'R_LIBS_USER':str(r_library),'R_LIBS_SITE':'','R_LIBS':'','R_PROFILE_USER':'/dev/null',
                'R_ENVIRON_USER':'/dev/null','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    try:
        with as_file(files('merci_rank_audit').joinpath('r_helpers','score_one.R')) as script:
            with as_file(files('merci_rank_audit').joinpath('r_helpers','engine.R')) as engine:
                command=[args.rscript,'--vanilla',str(script),*(str(x.resolve()) for x in inputs),
                         str(upstream),str(engine),str(out)]
                subprocess.run(command,check=True,env=env)
        manifest={'budget':5000,'cases':[{'name':run['case'],'score':'score.tsv','thresholds':'thresholds.csv',
                  'seed_uniform':20261027103,'seed_conditional':20261027003}]}
        (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        run_audit(out/'manifest.json',out/'audit',sample_mode=args.sample)
    except Exception as exc:
        (out/'run_failure.json').write_text(json.dumps({'status':'failed','error':str(exc)},indent=2)+'\n')
        raise
    result={'status':'generated_and_audited','case':run['case'],'upstream_commit':revision,
            'seconds':time.monotonic()-started,'historical_outputs_used_for_generation':False}
    (out/'run_complete.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main(argv=None):
    parser=build_parser();args=parser.parse_args(argv)
    try:
        if args.command=='init-example':
            with as_file(files('merci_rank_audit').joinpath('examples')) as source:
                shutil.copytree(source,args.output)
            result={'status':'created','output':str(args.output)}
        elif args.command=='bench2':
            result=run_bench2(args)
        else:
            if (args.member is None)!=(args.decision is None):
                raise ValueError('--member and --decision must be supplied together')
            witness=(args.member,args.decision=='selected') if args.member is not None else None
            result=run_audit(args.manifest,args.output,args.sample,args.case,witness,args.reuse_samples)
        print(json.dumps(result,indent=2))
        return 0
    except (ValueError,FileNotFoundError,FileExistsError,subprocess.CalledProcessError) as exc:
        print(f'error: {exc}',file=sys.stderr)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
