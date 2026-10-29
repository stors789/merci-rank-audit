"""Bind freshly generated score and fixed thresholds to a new exact audit."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--audit-python', default=sys.executable)
    a = p.parse_args(); score = (a.scores / 'score.tsv').resolve(); thresholds = (a.scores / 'thresholds.csv').resolve()
    if not score.is_file() or not thresholds.is_file():
        raise FileNotFoundError('New score and threshold files are required')
    a.output.mkdir(parents=True, mode=0o700, exist_ok=False)
    manifest = {'budget': 5000, 'cases': [{'name': 'GBML181', 'score': str(score), 'thresholds': str(thresholds),
                'seed_uniform': 20261028101, 'seed_conditional': 20261028102}]}
    path = a.output / 'manifest.json'; path.write_text(json.dumps(manifest, indent=2) + '\n')
    with (a.output / 'audit.log').open('x') as log:
        subprocess.run([a.audit_python, '-m', 'merci_rank_audit', 'audit', '--manifest', str(path.resolve()),
                        '--output', str((a.output / 'audit').resolve()), '--sample', 'none'], stdout=log, stderr=subprocess.STDOUT, check=True)
    result = {'status': 'COMPLETE', 'score_sha256': sha(score), 'threshold_sha256': sha(thresholds),
              'manifest_sha256': sha(path), 'random_sampling': 'not_requested',
              'audit': json.loads((a.output / 'audit/run.json').read_text())}
    (a.output / 'audit_receipt.json').write_text(json.dumps(result, indent=2) + '\n'); print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
