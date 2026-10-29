"""Run the adapted formal reference and evaluation workers for explicit cases."""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument("--config",type=Path,required=True)
p.add_argument("--case",required=True,help="Case name or 'all'")
p.add_argument("--output",type=Path,required=True)
p.add_argument("--stage",choices=["reference","evaluate","both"],default="both")
p.add_argument("--reference-root",type=Path,help="Existing references for evaluation-only runs")
a=p.parse_args();config=json.loads(a.config.read_text())
names=[x["case"] for x in config["cases"]] if a.case=="all" else [a.case]
if a.output.exists():raise FileExistsError(a.output)
a.output.mkdir(parents=True)
here=Path(__file__).resolve().parent
for name in names:
 for stage in (["reference","evaluate"] if a.stage=="both" else [a.stage]):
    script=here/("reference_worker.py" if stage=="reference" else "evaluate_worker.py")
    command=[sys.executable,str(script),"--config",str(a.config.resolve()),"--case",name,
             "--output-root",str((a.output/("references" if stage=="reference" else "results")).resolve()),
             "--spawn-time",str(time.perf_counter())]
    if stage=="evaluate":command+=["--reference-root",str((a.reference_root or (a.output/"references")).resolve())]
    subprocess.run(command,check=True)
