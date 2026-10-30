"""Run the first 100 proposals of the previously fixed P2 reference stream."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--audit-python',default=sys.executable)
    p.add_argument('--config',type=Path,default=Path(__file__).resolve().parents[2]/'configs/nuclear_formal_portable.json')
    a=p.parse_args();manifest=json.loads(a.manifest.read_text());config=json.loads(a.config.read_text())
    if len(manifest['cases'])!=1 or manifest['cases'][0]['name']!='P2' or config['reference_seed_entropy']!=[1573,5005,0]:
        raise ValueError('Expected the original P2 reference stream')
    case=manifest['cases'][0];base=a.manifest.resolve().parent
    identity={'scheme':'sha256-file-bytes-v1','integer_contract':'exact-integer-v2',
              'score':sha((base/case['score']).resolve()),'thresholds':sha((base/case['thresholds']).resolve()),
              'nuclear':sha((base/case['nuclear']).resolve())}
    a.output.mkdir(parents=True,mode=0o700,exist_ok=False)
    config['historical_input_identity']=config['input_identity'];config['input_identity']=identity
    config['reference_proposals']=100;config['reference_seconds_limit']=120
    local=a.output/'config.json';local.write_text(json.dumps(config,indent=2)+'\n')
    driver=Path(__file__).resolve().parents[1]/'nuclear/run.py'
    with (a.output/'stream.log').open('x') as log:
        subprocess.run([a.audit_python,str(driver),'--manifest',str(a.manifest.resolve()),'--config',str(local.resolve()),
                        '--stream','reference','--output',str((a.output/'stream').resolve())],stdout=log,stderr=subprocess.STDOUT,check=True)
    stream=json.loads((a.output/'stream/receipt.json').read_text())
    if stream['status']!='COMPLETE' or stream['proposals']!=100 or stream['seed_entropy']!=[1573,5005,0]:
        raise ValueError('The original reference prefix did not finish; its actual receipt is retained')
    result={'status':'COMPLETE','new_input_identity':identity,'original_reference_seed_entropy':[1573,5005,0],
            'proposals':stream['proposals'],'accepted':stream['accepted'],'historical_samples_used':False}
    (a.output/'prefix_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
