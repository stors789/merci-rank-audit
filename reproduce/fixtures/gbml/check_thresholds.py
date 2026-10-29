"""Reject a changed GBML181 maximum before any full-input work starts."""
import importlib.util
from pathlib import Path
import sys
import tempfile

repro = Path(__file__).resolve().parents[2]; driver = repro/'drivers/gbml'; sys.path.insert(0, str(driver))
spec = importlib.util.spec_from_file_location('gbml_run', driver/'run.py')
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
original = repro/'data/bio06/all_thresholds.csv'; module.check_thresholds(original)
source = original.read_bytes()
changed = source.replace(b'GBML181,original,1343,50,672,372,', b'GBML181,original,1343,50,672,373,')
assert changed != source
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder)/'changed.csv'; path.write_bytes(changed)
    try:
        module.check_thresholds(path)
    except ValueError:
        pass
    else:
        raise AssertionError('A changed frozen maximum was accepted')
print('PASS: original GBML181 table accepted and altered maximum rejected')
