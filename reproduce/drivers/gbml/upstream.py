"""Check the separately obtained MERCI checkout before invoking its functions."""
import subprocess

COMMIT = '79e105aaf13159bee325c3ceceb6ee679e430341'
FILES = ('MERCI-mtSNP.py', 'MERCI/R/readMTvar_10x.r', 'MERCI/R/MTmutMatrix_refined.r',
         'MERCI/R/MTmutMatrix_transform.r', 'MERCI/R/Denrich_mtMut_extr.r', 'MERCI/R/filter_MTmuts.r',
         'MERCI/R/Donor_enriched_vars.r', 'MERCI/R/Cell_Neff_cal.r', 'MERCI/R/w_cal2.r',
         'MERCI/R/mut_inf_extr.r', 'MERCI/R/MERCI_LOO_MT_est.r', 'MERCI/R/MERCI_ReceiverPre.r')


def check(path):
    if any(not (path / name).is_file() for name in FILES):
        raise FileNotFoundError('A required MERCI source file is missing from the pinned checkout')
    revision = subprocess.check_output(['git', '-C', str(path), 'rev-parse', 'HEAD'], text=True).strip()
    if revision != COMMIT:
        raise ValueError('The separate MERCI checkout is not at the pinned commit')
    subprocess.run(['git', '-C', str(path), 'diff', '--quiet', 'HEAD', '--', *FILES], check=True)
    return revision
