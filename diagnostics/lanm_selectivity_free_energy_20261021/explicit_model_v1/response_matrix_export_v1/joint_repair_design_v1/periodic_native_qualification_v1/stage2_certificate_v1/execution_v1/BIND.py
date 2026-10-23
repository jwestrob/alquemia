"""Bind exact completed dependencies for one Stage2 slot; no launch or release."""
import argparse
from COMMON import *

def run(a):
    mp=pin(a.manifest);m=read(mp);cfg={'manifest':mp}
    if a.kind=='build':
        cfg.update(status='bound_stage2_build_configuration',implementation=m['execution_sources']['BUILD.py'])
    else:
        require(a.state_id in m['rows'] and a.operation in ('initialize','gradient'),'declared slot')
        require(a.build is not None,'actual build required');bp=pin(a.build);build=read(bp)
        require(build['status']=='stage2_build_complete' and build['manifest']==mp,'exact completed build')
        cfg.update(status='bound_stage2_native_configuration',implementation=m['execution_sources']['NATIVE.py'],state_id=a.state_id,operation=a.operation,build=bp,central_scout_review=None,certificate=None,certificate_admission=None)
        if a.state_id not in ('central','BYPASS_OLD_ABSENT','BYPASS_CLONE_ABSENT','BYPASS_CLONE_ZERO'):
            require(a.scout_review is not None,'actual central scout review required');cfg['central_scout_review']=pin(a.scout_review)
            sr=read(cfg['central_scout_review']);require(sr['status']=='admitted_stage2_central_and_bypass_scout' and sr['manifest']==mp and sr['build']==bp,'scout review identity')
        if a.operation=='initialize':require(a.state_id!='central' and m['rows'][a.state_id]['role']=='ACTIVE','only32 declared initializer slots')
        elif m['rows'][a.state_id]['role']=='ACTIVE':
            require(a.certificate is not None,'state certificate required');cfg['certificate']=pin(a.certificate)
            if a.state_id!='central':
                require(a.certificate_admission is not None,'root state certificate admission required');cfg['certificate_admission']=pin(a.certificate_admission)
            certificate(cfg,m)
    a.out.mkdir(parents=True,exist_ok=False);save(a.out/'CONFIG.json',cfg)
    save(a.out/'ROOT_RELEASE_UNBOUND.json',{'status':None,'implementation':cfg['implementation'],'configuration':pin(a.out/'CONFIG.json'),'output':None,'source_review':None,**({'operation':a.operation,'state_id':a.state_id,'maximum_native_programs':1} if a.kind=='native' else {})})
    print(json.dumps({'configuration':pin(a.out/'CONFIG.json'),'release_required':True,'native_calls':0}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--kind',choices=['build','native'],required=True)
    p.add_argument('--state-id');p.add_argument('--operation',choices=['initialize','gradient']);p.add_argument('--build',type=Path)
    p.add_argument('--scout-review',type=Path);p.add_argument('--certificate',type=Path);p.add_argument('--certificate-admission',type=Path)
    p.add_argument('--out',type=Path,required=True);run(p.parse_args())
