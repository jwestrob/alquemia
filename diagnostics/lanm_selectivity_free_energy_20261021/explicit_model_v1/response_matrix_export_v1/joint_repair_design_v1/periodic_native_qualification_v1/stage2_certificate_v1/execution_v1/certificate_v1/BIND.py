"""Bind completed state initialization and scout; no assembly/factor/native call."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import CERTIFY as C

def run(a):
    sources=C.read(C.pin(a.sources));m=C.read(sources['execution_manifest'])
    for p in sources['certificate_sources'].values():C.verify(p)
    C.require(sources['certificate_sources']['BIND.py']==C.pin(__file__),'exact binder source')
    C.require(a.state_id in m['rows'] and a.state_id!='central' and m['rows'][a.state_id]['role']=='ACTIVE','one original32 changed-state slot')
    cfg={'status':'bound_stage2_certificate_configuration','implementation':sources['certificate_sources']['CERTIFY.py'],'manifest':sources['execution_manifest'],'certificate_sources':sources['certificate_sources'],'authorities':sources['authorities'],'state_id':a.state_id,'build':C.pin(a.build),'central_scout_review':C.pin(a.scout_review),'initializer':C.pin(a.initializer),'caps':C.CAPS,'nominal_contract':m['nominal_contract']}
    build=C.read(cfg['build']);scout=C.read(cfg['central_scout_review'])
    C.require(build['status']=='stage2_build_complete' and build['manifest']==cfg['manifest'],'exact successful build')
    C.require(scout['status']=='admitted_stage2_central_and_bypass_scout' and scout['manifest']==cfg['manifest'] and scout['build']==cfg['build'],'actual scout identity')
    C.require(scout['all_three_bypass_tokens_equal'] is True and scout['central_dual_energy_passed'] is True,'actual scout gates')
    C.initialization(cfg,m)
    a.out.mkdir(parents=True,exist_ok=False);C.save(a.out/'CONFIG.json',cfg)
    C.save(a.out/'ROOT_RELEASE_UNBOUND.json',{'status':None,'implementation':cfg['implementation'],'configuration':C.pin(a.out/'CONFIG.json'),'output':None,'source_review':None,'state_id':a.state_id,'maximum_assemblies':1,'maximum_numeric_factors':1,'maximum_solve_probes':10})
    print(C.json.dumps({'configuration':C.pin(a.out/'CONFIG.json'),'native_calls':0,'assemblies':0,'factors':0,'root_release_required':True}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--sources',type=Path,required=True);p.add_argument('--state-id',required=True);p.add_argument('--build',type=Path,required=True);p.add_argument('--scout-review',type=Path,required=True);p.add_argument('--initializer',type=Path,required=True);p.add_argument('--out',type=Path,required=True);run(p.parse_args())
