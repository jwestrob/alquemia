#!/usr/bin/env python3
"""Two explicit isolated installation phases; no numerical/model calls."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import time


def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def save(p,x):
    Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')


def run(mode,config_path):
    c=json.loads(Path(config_path).read_text());stage=Path(c['stage']);prefix=Path(c['prefix'])
    assert c['implementation']==pin(__file__) and c['mamba']==pin(c['mamba']['path'])
    out=stage/(mode+'_v1');out.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy();env.update(c['environment'])
    if mode=='resolve':
        assert not prefix.exists()
        argv=[c['mamba']['path'],'create','--dry-run','--json','--yes','--no-default-packages','--strict-channel-priority','--override-channels','-c','conda-forge','-p',str(prefix),*c['requested_packages']]
    else:
        previous=json.loads((stage/'resolve_v1/RECEIPT.json').read_text());assert previous['exit_code']==0
        plan=json.loads((stage/'resolve_v1/stdout.json').read_text());assert plan['success'] is True
        assert c['exact_plan']==pin(c['exact_plan']['path'])
        exact=json.loads(Path(c['exact_plan']['path']).read_text())
        assert exact['source_plan']==pin(stage/'resolve_v1/stdout.json')
        links=exact['packages'];assert links and not prefix.exists()
        assert [(r['name'],r['version'],r['build']) for r in links]==[(r['name'],r['version'],r['build_string']) for r in plan['actions']['LINK']]
        urls=[]
        for p in links:
            assert p['url'].startswith('https://conda.anaconda.org/conda-forge/') and len(p['md5'])==32
            assert len(p['sha256'])==64
            urls.append(p['url']+'#'+p['md5'])
        explicit=stage/'EXPLICIT.txt';explicit.write_text('@EXPLICIT\n'+'\n'.join(urls)+'\n')
        save(stage/'SELECTED_PACKAGES.json',links)
        argv=[c['mamba']['path'],'create','--json','--yes','--no-default-packages','--copy','-p',str(prefix),'--file',str(explicit)]
    before=resource.getrusage(resource.RUSAGE_CHILDREN);tick=time.monotonic();cpu=time.process_time()
    rec={'mode':mode,'implementation':pin(__file__),'configuration':pin(config_path),'argv':argv,'environment_overrides':c['environment'],'status':'started','native_calls':0,'factorizations':0}
    save(out/'RECEIPT.json',rec)
    try:
        with open(out/'stdout.json','w') as stdout,open(out/'stderr.txt','w') as stderr:
            process=subprocess.run(argv,env=env,stdout=stdout,stderr=stderr,timeout=c['timeouts_seconds'][mode])
        rec.update(exit_code=process.returncode,status='complete' if process.returncode==0 else 'failed')
    except Exception as exc:rec.update(status='failed',exit_code=None,error=repr(exc))
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    rec['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,'child_user_CPU_seconds':after.ru_utime-before.ru_utime,'child_system_CPU_seconds':after.ru_stime-before.ru_stime,'child_max_RSS_KiB':after.ru_maxrss,'numerical_threads':1}
    rec['stdout']=pin(out/'stdout.json');rec['stderr']=pin(out/'stderr.txt')
    if mode=='install':rec.update(explicit=pin(stage/'EXPLICIT.txt'),selected_packages=pin(stage/'SELECTED_PACKAGES.json'))
    save(out/'RECEIPT.json',rec);print(json.dumps({'receipt':pin(out/'RECEIPT.json'),'status':rec['status'],'exit_code':rec['exit_code'],'cost':rec['cost']}))
    return rec


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('resolve','install'));p.add_argument('--config',type=Path,required=True);a=p.parse_args()
    result=run(a.mode,a.config);raise SystemExit(result['status']!='complete')
