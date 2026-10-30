"""Apply the previously fixed C5 decisions to newly reconstructed joint RNA roles."""
import argparse
import hashlib
import json
from pathlib import Path
import re

import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--roles', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, mode=0o700, exist_ok=False)
    obs = pd.read_csv(a.roles,sep='\t',dtype={'cell_id':str,'GSM':str,'barcode':str,'identity_cluster':str})
    p2 = obs[obs.GSM.eq('GSM7774400')].copy()
    if len(p2) != 4359 or not p2.cell_id.is_unique or not p2.barcode.is_unique:
        raise ValueError('P2 joint-QC cell axis differs')
    if not p2.cell_id.eq('GSM7774400:' + p2.barcode).all() or not p2.barcode.map(lambda x: bool(re.fullmatch('[ACGT]{24}', x))).all():
        raise ValueError('P2 requires the exact 24-base barcode and GSM namespace')
    roles = p2.working_type.value_counts().to_dict()
    path = a.output/'P2_roles.tsv.gz'; p2.to_csv(path, sep='\t', index=False, compression='gzip')
    receipt = {'status': 'COMPLETE', 'GSM': 'GSM7774400', 'RNA_QC_n': len(p2), 'working_roles': roles,
               'hepatic_clusters': p2[p2.working_type.eq('hepatic')].identity_cluster.value_counts().to_dict(),
               'roles_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
               'role_status': 'working RNA identity; no patient mixing or independent malignant labels'}
    (a.output/'roles_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n');print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
