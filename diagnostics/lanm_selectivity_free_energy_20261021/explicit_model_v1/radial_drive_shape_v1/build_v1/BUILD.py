#!/usr/bin/env python3
"""Pin/release adapter for the unchanged qualified compiler/targets recipe."""
import hashlib,json,os,runpy,sys
from pathlib import Path

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def load(p):
    assert pin(p['path'])==p,p;return json.loads(Path(p['path']).read_text())
def fingerprint(cfg):
    body={k:v for k,v in cfg.items() if k!='root_release'}
    return hashlib.sha256(json.dumps(body,sort_keys=True,separators=(',',':')).encode()).hexdigest()
config_path,out=sys.argv[1:];cfg=json.loads(Path(config_path).read_text())
assert cfg['implementation']==pin(__file__)
assert cfg['compiler_workers']==12 and os.environ.get('SLURM_CPUS_PER_TASK')=='12'
assert os.environ.get('SLURM_JOB_ID','').isdigit() and os.environ.get('SLURM_NTASKS','1')=='1'
assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
ready=load(cfg['source_ready']);review=load(cfg['source_review']);release=load(cfg['root_release'])
assert cfg['source_manifest']==ready['source_manifest']
assert all(review[k]==v for k,v in ready['source_review_contract'].items())
assert release['status']=='released_radial_drive_isolated_build' and release['configuration_fingerprint']==fingerprint(cfg)
assert release['source_manifest']==cfg['source_manifest'] and release['source_review']==cfg['source_review']
assert release['output']==str(Path(out).resolve()) and release['maximum_build_attempts']==1
assert pin(cfg['parent_builder']['path'])==cfg['parent_builder']
# No import of native source or scientific execution. Parent owns immutable output
# creation, exact compiler commands, logs, build receipts and actual child costs.
sys.argv=[cfg['parent_builder']['path'],cfg['source_manifest']['path'],cfg['source_manifest']['sha256'],out]
runpy.run_path(cfg['parent_builder']['path'],run_name='__main__')
