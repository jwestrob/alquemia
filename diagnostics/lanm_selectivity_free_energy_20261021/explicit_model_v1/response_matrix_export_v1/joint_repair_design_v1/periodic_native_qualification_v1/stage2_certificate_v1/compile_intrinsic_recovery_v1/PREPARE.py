"""One declaration-only source recovery and later exact metadata rebinding."""
import argparse
import copy
import difflib
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'execution_v1'))
from COMMON import *

LINE='      intrinsic :: maxval\n'
ANCHOR='   subroutine after_direct(a,b)\n'

def source(out):
    tick=time.monotonic();cpu=time.process_time();out.mkdir(parents=True,exist_ok=False)
    # Derive workspace from the existing source authority, without changing it.
    parent_inputs=json.loads((HERE.parent/'SOURCE_INPUTS_v2.json').read_text())
    old_source=pin(Path(parent_inputs['source_output'])/'MANIFEST.json');old=read(old_source)
    wstage=Path(parent_inputs['source_output']).parent
    failure=pin(wstage/'execution_v1/build_run_v1/execution/RESULT.json');failed=read(failure)
    require(failed['status']=='stage2_build_failed' and failed['compile_attempts']==1 and failed['link_attempts']==0 and failed['native_programs']==0,'exact first-compile failure')
    stderr=failed['commands']['compile_GATE_f90']['stderr'];error=verify(stderr).read_text()
    require('108:30:' in error and 'Syntax error in argument list' in error and 'stage2_gate_state.mod' in error,'actual compiler symptom')
    parent_dir=Path(parent_inputs['induce']['path']).parent
    sizes=pin(parent_dir/'sizes.f');atoms=pin(parent_dir/'atoms.f')
    require('parameter (maxval=8)' in verify(sizes).read_text() and 'use sizes' in verify(atoms).read_text(),'host-associated intrinsic collision authority')
    files={};checks=[];difference=''
    for name,item in old['files'].items():
        original=verify(item).read_bytes();current=original
        if name=='GATE.f90':
            require(original.count(ANCHOR.encode())==1 and LINE.encode() not in original,'single unpatched routine')
            current=original.replace(ANCHOR.encode(),(ANCHOR+LINE).encode())
            require(current.replace((ANCHOR+LINE).encode(),ANCHOR.encode())==original,'inverse patch exactly original bytes')
            difference=''.join(difflib.unified_diff(original.decode().splitlines(True),current.decode().splitlines(True),fromfile=item['path'],tofile=str(out/name)))
        (out/name).write_bytes(current);files[name]=pin(out/name)
        checks.append({'file':name,'identical_to_original':current==original,'inverse_patch_identity':name=='GATE.f90' and current.replace((ANCHOR+LINE).encode(),ANCHOR.encode())==original})
    require(sum(not x['identical_to_original'] for x in checks)==1,'exactly one source changed')
    (out/'PATCH.diff').write_text(difference)
    manifest={'status':'prepared_stage2_intrinsic_declaration_recovery_unbuilt','producer':pin(__file__),'parent_source':old_source,'failed_build':failure,'failed_stderr':stderr,'name_collision_sources':{'atoms':atoms,'sizes':sizes},'files':files,'patch':pin(out/'PATCH.diff'),'checks':checks,'executable_statements_unchanged':True,'native_calls':0,'compile_calls':0,'actual':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    save(out/'MANIFEST.json',manifest);print(json.dumps({'source_manifest':pin(out/'MANIFEST.json'),'patch':manifest['patch'],'actual':manifest['actual']}))

def bindings(source_manifest,native_review,out):
    tick=time.monotonic();cpu=time.process_time();sp=pin(source_manifest);source=read(sp);nrp=pin(native_review);nr=read(nrp)
    require(source['status']=='prepared_stage2_intrinsic_declaration_recovery_unbuilt' and source['executable_statements_unchanged'] is True,'exact source-only recovery')
    require(nr['status']=='root_reviewed_stage2_native_source_design' and nr['source_manifest']==sp,'new root native-source review')
    for p in source['files'].values():verify(p)
    parent=json.loads((HERE.parent/'SOURCE_INPUTS_v2.json').read_text());wstage=Path(parent['source_output']).parent
    oldp=pin(wstage/'execution_v1/preparation_v1/MANIFEST.json');old=read(oldp);new=copy.deepcopy(old)
    new['native_sources']=sp;new['native_source_review']=nrp
    new['binding_recovery']={'implementation':pin(__file__),'parent_manifest':oldp}
    out.mkdir(parents=True,exist_ok=False);save(out/'MANIFEST.json',new);mp=pin(out/'MANIFEST.json')
    center=read(pin(wstage/'execution_v1/preparation_v1/CENTRAL_CERTIFICATE.json'));center['manifest']=mp;center['source_manifest']=sp
    center['parent_certificate']=pin(wstage/'execution_v1/preparation_v1/CENTRAL_CERTIFICATE.json');center['implementation']=pin(__file__)
    save(out/'CENTRAL_CERTIFICATE.json',center)
    config={'status':'bound_stage2_build_configuration','implementation':new['execution_sources']['BUILD.py'],'manifest':mp}
    save(out/'CONFIG_BUILD.json',config)
    save(out/'ROOT_RELEASE_BUILD_UNBOUND.json',{'status':None,'implementation':config['implementation'],'configuration':pin(out/'CONFIG_BUILD.json'),'output':None,'source_review':None})
    cert=read(pin(HERE.parent/'execution_v1/certificate_v1/SOURCE_INPUTS.json'));cert['execution_manifest']=mp;cert['authorities']['native_source_review']=nrp
    save(out/'CERTIFICATE_SOURCE_INPUTS.json',cert)
    result={'status':'prepared_stage2_intrinsic_recovery_bindings','implementation':pin(__file__),'old_manifest':oldp,'manifest':mp,'source_manifest':sp,'native_source_review':nrp,'central_certificate':pin(out/'CENTRAL_CERTIFICATE.json'),'configuration':pin(out/'CONFIG_BUILD.json'),'certificate_source_inputs':pin(out/'CERTIFICATE_SOURCE_INPUTS.json'),'checks':['native source/review fields changed; explicit metadata-recovery lineage added','all36 runtime/source input and execution/helper pins unchanged','central expected state PME and factor pins unchanged','certificate producer and physics helpers unchanged'],'actual':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'native_calls':0,'compile_calls':0,'assemblies':0,'factors':0,'solves':0}}
    changed=[key for key in old if old[key]!=new[key]];require(changed==['native_sources','native_source_review'] and set(new)-set(old)=={'binding_recovery'},'exact manifest field delta')
    save(out/'RESULT.json',result);print(json.dumps(result))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--phase',choices=['source','bindings'],required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--source-manifest',type=Path);p.add_argument('--native-review',type=Path);a=p.parse_args()
    if a.phase=='source':source(a.out.resolve())
    else:bindings(a.source_manifest,a.native_review,a.out.resolve())
