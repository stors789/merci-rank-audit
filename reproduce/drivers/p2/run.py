"""Regenerate the fixed P2 input chain and its exact audit from public inputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'gbml'))
from upstream import check as check_upstream

GIB = 1024**3
TABLE_SHA = '7d75d473b3861ec69ed5b3368305730fdcb8ee190d517329014b2381dbdeb0d6'
GMT_SHA = 'ee2463540042078bfa3f67828e1e223bb354446d9fbb4d22845866835ba5c772'
CHRM_SHA = 'f156ff3f65bbcc80c7ebb9936dceb96b1477b4f8f535c4e1dbe7baea225cbc66'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('archives', 'r1', 'r2', 'reference', 'index', 'minimap', 'gmt',
                 'upstream', 'r-library', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--rscript', default='Rscript')
    parser.add_argument('--audit-python', default=sys.executable)
    parser.add_argument('--thresholds', type=Path,
                        default=Path(__file__).resolve().parents[2] / 'data/bio06/all_thresholds.csv')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    gbml = here.parent / 'gbml'
    for path in (args.r1, args.r2, args.reference, args.index, args.minimap, args.gmt, args.thresholds):
        if not path.is_file():
            raise FileNotFoundError(path)
    if not args.archives.is_dir() or not args.r_library.is_dir():
        raise FileNotFoundError('Explicit GEO archive and R library directories are required')
    for sample in json.loads((here / 'sample_manifest.json').read_text()):
        if not (args.archives / sample['file']).is_file():
            raise FileNotFoundError(args.archives / sample['file'])
    if hashlib.sha256(args.thresholds.read_bytes()).hexdigest() != TABLE_SHA:
        raise ValueError('The original fixed threshold table differs')
    if hashlib.sha256(args.gmt.read_bytes()).hexdigest() != GMT_SHA:
        raise ValueError('The fixed Hallmark GMT differs from the original local input')
    if args.index.stat().st_size != 15420201719:
        raise ValueError('The original hg38 index has a different size')
    if subprocess.check_output([str(args.minimap.resolve()), '--version'], text=True).strip() != '2.30-r1287':
        raise ValueError('Expected minimap2 2.30-r1287')
    revision = check_upstream(args.upstream.resolve())
    cpus = sorted(os.sched_getaffinity(0))
    if len(cpus) < 4:
        raise ValueError('The recorded identity and alignment stages require four CPUs')
    hard = resource.getrlimit(resource.RLIMIT_AS)[1]
    if hard != resource.RLIM_INFINITY and hard < 8 * GIB:
        raise ValueError('The alignment child requires an 8 GiB address-space hard limit')
    out = args.output.resolve()
    out.mkdir(parents=True, mode=0o700, exist_ok=False)
    env = os.environ.copy()
    env.update({'PYTHONNOUSERSITE': '1', 'PYTHONDONTWRITEBYTECODE': '1', 'PYTHONPATH': '',
                'R_LIBS_USER': str(args.r_library.resolve()), 'R_LIBS_SITE': '', 'R_LIBS': '',
                'R_PROFILE_USER': os.devnull, 'R_ENVIRON_USER': os.devnull})
    elapsed = {}

    def stage(name, command, threads=4, align=False):
        child_env = env.copy()
        for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                    'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
            child_env[key] = str(threads)

        def limits():
            os.sched_setaffinity(0, cpus[:4] if threads == 4 else cpus[:1])
            resource.setrlimit(resource.RLIMIT_AS, (4 * GIB, (8 if align else 4) * GIB))

        started = time.monotonic()
        with (out / (name + '.log')).open('x') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=child_env,
                           check=True, preexec_fn=limits)
        elapsed[name] = round(time.monotonic() - started, 3)
        print(json.dumps({'stage': name, 'status': 'COMPLETE', 'seconds': elapsed[name]}), flush=True)

    py = sys.executable
    counts = out / 'counts'
    clusters = out / 'clusters'
    refined = out / 'refined'
    roles = out / 'roles'
    reference = out / 'reference'
    extraction = out / 'extraction'
    alignment = out / 'alignment'
    consensus = out / 'consensus'
    qualification = out / 'qualification'
    inputs = out / 'merci_inputs'
    scores = out / 'scores'
    nuclear = out / 'nuclear'
    binding = out / 'audit_input'
    stage('counts', [py, str(here / 'prepare_counts.py'), '--archives', str(args.archives.resolve()),
                     '--gmt', str(args.gmt.resolve()), '--output', str(counts)])
    matrix = counts / 'raw_HCC_counts.h5ad'
    stage('clusters', [py, str(here / 'cluster_joint.py'), '--counts', str(matrix),
                       '--gmt', str(args.gmt.resolve()), '--output', str(clusters)])
    stage('refined', [py, str(here / 'refine_stream.py'), '--counts', str(matrix),
                      '--global-roles', str(clusters / 'global_roles.tsv.gz'), '--gmt', str(args.gmt.resolve()),
                      '--output', str(refined)])
    stage('roles', [py, str(here / 'finalize_roles.py'), '--roles', str(refined / 'refined_roles.tsv.gz'),
                    '--output', str(roles)])
    role_file = roles / 'P2_roles.tsv.gz'
    stage('reference', [py, str(gbml / 'reference.py'), '--reference', str(args.reference.resolve()),
                        '--output', str(reference)], threads=1)
    if json.loads((reference / 'reference_receipt.json').read_text())['sequence_sha256'] != CHRM_SHA:
        raise ValueError('The newly extracted chrM differs from the original reference')
    chr_m = reference / 'chrM.fa'
    stage('extraction', [py, str(here / 'extract_fastq.py'), '--r1', str(args.r1.resolve()),
                         '--r2', str(args.r2.resolve()), '--roles', str(role_file), '--reference', str(chr_m),
                         '--output', str(extraction)])
    stage('alignment', [py, str(here / 'align_candidates.py'), '--extraction', str(extraction),
                        '--roles', str(role_file), '--index', str(args.index.resolve()),
                        '--minimap', str(args.minimap.resolve()), '--output', str(alignment)], align=True)
    stage('consensus', [py, str(here / 'consensus.py'), '--extraction', str(extraction),
                        '--alignment', str(alignment), '--output', str(consensus)])
    stage('qualification', [py, str(here / 'qualify.py'), '--roles', str(role_file), '--consensus', str(consensus),
                            '--alignment', str(alignment), '--reference', str(chr_m),
                            '--upstream', str(args.upstream.resolve()), '--r-library', str(args.r_library.resolve()),
                            '--rscript', args.rscript, '--output', str(qualification)])
    stage('merci_inputs', [py, str(here / 'prepare_merci.py'), '--qualification', str(qualification),
                           '--counts', str(matrix), '--output', str(inputs)])
    stage('scores', [args.rscript, '--vanilla', str(here / 'score_merci.R'), str(inputs), str(scores),
                     str(args.upstream.resolve()), str(args.r_library.resolve()), str(args.thresholds.resolve())])
    stage('nuclear', [py, str(here / 'nuclear.py'), '--counts', str(matrix),
                      '--metadata', str(qualification / 'full_metadata.tsv.gz'), '--cells', str(scores / 'cells.tsv'),
                      '--gmt', str(args.gmt.resolve()), '--output', str(nuclear)])
    stage('audit_input', [py, str(here / 'prepare_audit.py'), '--cells', str(scores / 'cells.tsv'),
                          '--threshold-table', str(args.thresholds.resolve()), '--nuclear', str(nuclear / 'nuclear.npz'),
                          '--output', str(binding)], threads=1)
    manifest = binding / 'manifest.json'
    stage('exact_audit', [args.audit_python, '-m', 'merci_rank_audit', 'audit', '--manifest', str(manifest),
                          '--output', str(out / 'exact_audit'), '--sample', 'none'], threads=1)
    stage('reference_prefix', [py, str(here / 'short_prefix.py'), '--manifest', str(manifest),
                               '--audit-python', args.audit_python, '--output', str(out / 'reference_prefix')], threads=1)
    result = {'status': 'COMPLETE', 'upstream_commit': revision, 'seconds_by_stage': elapsed,
              'scope': 'new private outputs; no historical roles, scores, consensus or samples used as input',
              'binding_receipt': str(binding / 'binding_receipt.json'),
              'exact_audit': str(out / 'exact_audit' / 'run.json'),
              'reference_prefix': str(out / 'reference_prefix' / 'prefix_receipt.json')}
    (out / 'run_complete.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
