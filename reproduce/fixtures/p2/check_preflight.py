"""Ensure changed fixed inputs are rejected before a new output is created."""
import argparse
import json
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, mode=0o700, exist_ok=False)
    driver = Path(__file__).resolve().parents[2] / 'drivers/p2'
    archives = root / 'archives'
    archives.mkdir()
    for sample in json.loads((driver / 'sample_manifest.json').read_text()):
        (archives / sample['file']).touch()
    for name in ('r1', 'r2', 'reference', 'index', 'minimap', 'gmt', 'thresholds'):
        (root / name).write_text('artificial input\n')
    (root / 'r-library').mkdir()
    base = [sys.executable, str(driver / 'run.py'), '--archives', str(archives),
            '--r1', str(root / 'r1'), '--r2', str(root / 'r2'), '--reference', str(root / 'reference'),
            '--index', str(root / 'index'), '--minimap', str(root / 'minimap'), '--gmt', str(root / 'gmt'),
            '--upstream', str(root / 'missing_upstream'), '--r-library', str(root / 'r-library')]
    for name, extra, expected in (('bad_table', ['--thresholds', str(root / 'thresholds')], 'threshold table differs'),
                                  ('bad_gmt', [], 'Hallmark GMT differs')):
        destination = root / name
        result = subprocess.run(base + ['--output', str(destination)] + extra, capture_output=True, text=True)
        (root / (name + '.log')).write_text(result.stdout + result.stderr)
        if result.returncode == 0 or expected not in result.stderr or destination.exists():
            raise AssertionError(name + ' was not rejected before output creation')
    print('PASS: changed threshold table and GMT rejected before starting the workflow')


if __name__ == '__main__':
    main()
