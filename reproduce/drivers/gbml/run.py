"""Rebuild the fixed GBML181 score and exact audit from public input files."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

from upstream import check as check_upstream

GIB = 1024**3


def check_thresholds(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != '7d75d473b3861ec69ed5b3368305730fdcb8ee190d517329014b2381dbdeb0d6':
        raise ValueError('The supplied GBML181 threshold table differs from the published frozen table')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('matrix', 'features', 'barcodes', 'sra1', 'sra2', 'reference', 'index', 'minimap',
                 'extractor', 'vdb-lib', 'upstream', 'r-library', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--rscript', default='Rscript')
    p.add_argument('--audit-python', default=sys.executable)
    p.add_argument('--thresholds', type=Path, default=Path(__file__).resolve().parents[2] / 'data/bio06/all_thresholds.csv')
    a = p.parse_args(); here = Path(__file__).resolve().parent
    for path in (a.matrix, a.features, a.barcodes, a.sra1, a.sra2, a.reference, a.index, a.minimap, a.extractor, a.thresholds):
        if not path.is_file():
            raise FileNotFoundError(path)
    if not a.vdb_lib.is_dir() or not a.r_library.is_dir():
        raise FileNotFoundError('Explicit VDB and R library directories are required')
    check_thresholds(a.thresholds)
    revision = check_upstream(a.upstream.resolve())
    cpus = sorted(os.sched_getaffinity(0))
    if len(cpus) < 4:
        raise ValueError('RNA role stage requires access to four CPUs for the recorded numerical setup')
    hard = resource.getrlimit(resource.RLIMIT_AS)[1]
    if hard != resource.RLIM_INFINITY and hard < 8*GIB:
        raise ValueError('The minimap child needs an 8 GiB address-space hard limit')
    out = a.output.resolve(); out.mkdir(parents=True, mode=0o700, exist_ok=False)
    base_env = os.environ.copy(); base_env.update({'PYTHONNOUSERSITE': '1', 'PYTHONDONTWRITEBYTECODE': '1', 'PYTHONPATH': '',
        'R_LIBS_USER': str(a.r_library.resolve()), 'R_LIBS_SITE': '', 'R_LIBS': '',
        'R_PROFILE_USER': os.devnull, 'R_ENVIRON_USER': os.devnull})
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):
        base_env[name] = '1'
    elapsed = {}

    def stage(name, args, rna=False, align=False):
        env = base_env.copy()
        if rna:
            for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
                env[key] = '4'
        def limits():
            os.sched_setaffinity(0, cpus[:4] if rna else cpus[:1])
            resource.setrlimit(resource.RLIMIT_AS, (4*GIB, 8*GIB if align else 4*GIB))
        clock = time.monotonic()
        with (out / (name + '.log')).open('x') as log:
            subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, env=env, check=True, preexec_fn=limits)
        elapsed[name] = round(time.monotonic()-clock, 3)
        print(json.dumps({'stage': name, 'status': 'COMPLETE', 'seconds': elapsed[name]}), flush=True)

    py = sys.executable; roles = out / 'roles'; reference = out / 'reference'; extraction = out / 'extraction'
    alignment = out / 'alignment'; consensus = out / 'consensus'; variants = out / 'variants'; rna = out / 'rna'
    qualification = out / 'qualification'; scores = out / 'scores'; audit = out / 'audit_input'
    stage('roles', [py, str(here/'roles.py'), '--matrix', str(a.matrix.resolve()), '--features', str(a.features.resolve()),
          '--barcodes', str(a.barcodes.resolve()), '--output', str(roles)], rna=True)
    stage('reference', [py, str(here/'reference.py'), '--reference', str(a.reference.resolve()), '--output', str(reference)])
    stage('extraction', [py, str(here/'extract.py'), '--sra1', str(a.sra1.resolve()), '--sra2', str(a.sra2.resolve()),
          '--barcodes', str(a.barcodes.resolve()), '--roles', str(roles/'fixed_roles.tsv'),
          '--kmers', str(reference/'chrM_21mers.u64'), '--extractor', str(a.extractor.resolve()),
          '--vdb-lib', str(a.vdb_lib.resolve()), '--output', str(extraction)])
    stage('alignment', [py, str(here/'align.py'), '--extraction', str(extraction), '--roles', str(roles/'fixed_roles.tsv'),
          '--index', str(a.index.resolve()), '--minimap', str(a.minimap.resolve()), '--output', str(alignment)], align=True)
    stage('consensus', [py, str(here/'consensus.py'), '--alignment', str(alignment), '--output', str(consensus)])
    stage('variants', [py, str(here/'variants.py'), '--consensus', str(consensus), '--alignment', str(alignment),
          '--roles', str(roles/'fixed_roles.tsv'), '--reference', str(reference/'chrM.fa'),
          '--upstream', str(a.upstream.resolve()), '--output', str(variants)])
    stage('expression', [py, str(here/'expression.py'), '--matrix', str(a.matrix.resolve()), '--features', str(a.features.resolve()),
          '--barcodes', str(a.barcodes.resolve()), '--roles', str(roles/'fixed_roles.tsv'), '--output', str(rna)])
    stage('qualification', [a.rscript, '--vanilla', str(here/'qualify.R'), str(a.upstream.resolve()),
          str(variants), str(qualification), str(a.r_library.resolve())])
    stage('scores', [a.rscript, '--vanilla', str(here/'score.R'), str(a.upstream.resolve()), str(variants),
          str(qualification), str(rna), str(a.thresholds.resolve()), str(scores), str(a.r_library.resolve())])
    stage('summary', [py, str(here/'summarize.py'), '--root', str(out), '--output', str(out/'source_and_score_summary.json')])
    stage('audit', [py, str(here/'audit.py'), '--scores', str(scores), '--audit-python', a.audit_python, '--output', str(audit)])
    result = {'status': 'COMPLETE', 'upstream_commit': revision, 'seconds_by_stage': elapsed,
              'scope': 'new private outputs; no historical score, role, variant or sample inputs',
              'exact_audit_receipt': str(audit/'audit_receipt.json')}
    (out/'run_complete.json').write_text(json.dumps(result, indent=2)+'\n'); print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
