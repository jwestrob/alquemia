"""Bind all nine inactive native configs after an actual separately admitted build."""
import argparse
from BINDING import *

def bind(manifest,build,review,out):
    mp=pin(manifest);m=read(mp);validate_manifest(m);bp=pin(build);rp=pin(review);r=read(rp)
    require(r['status']=='admitted_cutoff12_baseline_build' and r['result']==bp and r['manifest']==mp,'Root actual build admission')
    native_build({'build':bp,'manifest':mp},m)
    out.mkdir(parents=True,exist_ok=False);configs={}
    for sid in IDS:
        folder=out/sid;folder.mkdir()
        cfg={'status':'bound_cutoff12_native_configuration','manifest':mp,'implementation':m['execution_sources']['RUN.py'],
             'state_id':sid,'operation':'gradient','build':bp,'build_root_review':rp,'certificate':None,'certificate_admission':None}
        save(folder/'CONFIG.json',cfg);configs[sid]=pin(folder/'CONFIG.json')
        save(folder/'ROOT_RELEASE_UNBOUND.json',{'status':None,'implementation':cfg['implementation'],'configuration':configs[sid],
             'output':m['rows'][sid]['output'],'source_review':None,'operation':'gradient','state_id':sid,'maximum_native_programs':1})
    result={'status':'bound_nine_cutoff12_baseline_configs_no_native','implementation':pin(__file__),'manifest':mp,'build':bp,'build_root_review':rp,
            'configurations':configs,'native_calls':0,'submissions':0,'maximum_native_programs':9}
    save(out/'RESULT.json',result);print(json.dumps({'result':pin(out/'RESULT.json'),'slots':9,'native_calls':0}))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for n in ('manifest','build','review','out'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();bind(a.manifest,a.build,a.review,a.out)
