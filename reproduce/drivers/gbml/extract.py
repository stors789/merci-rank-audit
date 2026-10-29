"""Extract and reconcile both native SRA streams with explicit generated roles."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import time

RUNS = (('SRR24283915', 98159746, 0), ('SRR24283916', 102126028, 98159746))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('sra1', 'sra2', 'barcodes', 'roles', 'kmers', 'extractor', 'vdb-lib', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args(); started = time.monotonic()
    for path in (a.sra1, a.sra2, a.barcodes, a.roles, a.kmers, a.extractor):
        if not path.is_file():
            raise FileNotFoundError(path)
    if not a.vdb_lib.is_dir():
        raise FileNotFoundError(a.vdb_lib)
    a.output.mkdir(parents=True, mode=0o700, exist_ok=False); batches = a.output / 'candidate_batches'; batches.mkdir()
    env = os.environ.copy(); env['LD_LIBRARY_PATH'] = str(a.vdb_lib.resolve())
    # A user-supplied local SRA archive is read directly; no accession lookup is needed.
    env['NCBI_SETTINGS'] = str((a.output / 'local_settings.mkfg').resolve())
    (a.output / 'local_settings.mkfg').write_text('/repository/remote/disabled = "true"\n')
    times = []
    for (run, spots, offset), archive in zip(RUNS, (a.sra1, a.sra2)):
        command = [str(a.extractor.resolve()), str(archive.resolve()), run, str(spots), str(offset),
                   str(a.barcodes.resolve()), str(a.roles.resolve()), str(a.kmers.resolve()), str(batches.resolve()), '10000000']
        clock = time.monotonic()
        with (a.output / (run + '.log')).open('x') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env)
        times.append({'run': run, 'returncode': result.returncode, 'seconds': round(time.monotonic()-clock, 3)})
        if result.returncode:
            raise RuntimeError(f'{run} extraction failed; first partial files and log retained')
        print(json.dumps(times[-1]), flush=True)
    roles = {line.split('\t')[0].removesuffix('-1') for line in a.roles.read_text().splitlines()[1:]}
    totals = {}; candidate_files = []
    for run, expected, _ in RUNS:
        with (batches / (run + '_batches.csv')).open() as f:
            rows = list(csv.DictReader(f))
        end = 0
        for row in rows:
            first, last, n = int(row['first_spot']), int(row['last_spot']), int(row['spots'])
            if first != end+1 or last-first+1 != n or int(row['valid_structure'])+int(row['unexpected_structure']) != n:
                raise ValueError('Noncontiguous spots or unexpected ledger dimensions')
            if int(row['source'])+int(row['target']) != int(row['role_match']):
                raise ValueError('Role totals differ')
            count, name = int(row['candidates']), row['candidate_file']
            if bool(count) != bool(name):
                raise ValueError('Candidate batch filename/count mismatch')
            if name:
                path = batches / name
                if not path.is_file():
                    raise FileNotFoundError(path)
                candidate_files.append(path)
            end = last
        if end != expected:
            raise ValueError('Full run EOF not reached')
        cols = ('spots','valid_structure','unexpected_structure','barcode_acgt','axis_match','barcode_q20',
                'role_match','source','target','umi_acgt','umi_q20','candidates')
        totals[run] = {k: sum(int(row[k]) for row in rows) for k in cols}
        if totals[run]['unexpected_structure']:
            raise ValueError('Unexpected full-run read structures')
    combined = a.output / 'all_candidates.fastq.gz'
    with combined.open('xb') as out:
        for path in candidate_files:
            with path.open('rb') as f:
                shutil.copyfileobj(f, out, 8*1024*1024)
    count = 0; last_id = 0; byrun = {run: 0 for run, _, _ in RUNS}
    with gzip.open(combined, 'rt') as f:
        while header := f.readline():
            seq, plus, quality = (f.readline().rstrip('\n') for _ in range(3))
            if not header.startswith('@') or plus != '+' or len(seq) != 91 or len(quality) != 91:
                raise ValueError('Malformed candidate FASTQ')
            cell, umi, number = header.strip()[1:].split('_'); number = int(number)
            if cell not in roles or len(umi) != 12 or set(umi)-set('ACGT') or number <= last_id:
                raise ValueError('Invalid role, UMI, or duplicate/out-of-order spot ID')
            last_id = number; count += 1
            byrun['SRR24283915' if number <= RUNS[0][1] else 'SRR24283916'] += 1
    if count != sum(x['candidates'] for x in totals.values()) or any(byrun[k] != v['candidates'] for k, v in totals.items()):
        raise ValueError('Candidate counts differ from per-run ledgers')
    h = hashlib.sha256()
    with combined.open('rb') as f:
        while chunk := f.read(8*1024*1024):
            h.update(chunk)
    result = {'status': 'COMPLETE', 'jobs': times, 'full_run_spot_counts': totals, 'candidate_records': count,
              'candidate_by_run': byrun, 'combined_sha256': h.hexdigest(), 'combined_bytes': combined.stat().st_size,
              'candidate_gzip_crc_passed': True, 'candidate_ID_order_and_uniqueness_passed': True,
              'no_full_fastq': True, 'seconds': round(time.monotonic()-started, 3),
              'max_child_RSS_KiB': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}
    (a.output / 'extraction_receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
