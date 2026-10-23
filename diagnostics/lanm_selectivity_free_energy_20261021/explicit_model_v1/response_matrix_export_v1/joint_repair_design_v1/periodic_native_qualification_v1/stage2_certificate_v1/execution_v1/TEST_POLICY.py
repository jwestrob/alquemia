"""Tiny metadata/mock controls only; no compiler, native or numerical backend."""
import argparse
import ast
import copy
import tempfile
from unittest.mock import patch
from COMMON import *

HERE=Path(__file__).resolve().parent

def run(out):
    tick=time.monotonic();cpu=time.process_time();checks=[]
    def check(name,fn,reject=False):
        try:fn()
        except (ValueError,KeyError,FileNotFoundError):
            if not reject:raise
        else:require(not reject,'expected rejection: '+name)
        checks.append(name)
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)
        def doc(name,x):save(p/name,x);return pin(p/name)
        s=p/'implementation.py';s.write_text('# immutable fixture\n');ip=pin(s)
        token=p/'state.dat';token.write_text('one exact state\n');tp=pin(token)
        env=doc('env.json',{'prefix':sys.prefix,'runtime_files':[]})
        stage1=doc('stage1.json',{'status':'admitted_finite_periodic_operator_parity','sampled_operator_parity_passed':True,'candidate':tp,'static_input':tp,'expected_central_state':tp,'initialized_PME':tp,'static_root_admission':tp})
        nr=doc('native_review.json',{'status':'root_reviewed_stage2_native_source_design'})
        m={'status':'prepared_stage2_execution_inputs_no_numerical_calls','execution_sources':{'implementation':ip},'helpers':{},'stage1_actual_admission':stage1,'candidate':tp,'central_static_input':tp,'native_source_review':nr,'environment':env,'native_sources':tp,'rows':{'central':{'inputs':{'xyz':tp}},'changed':{'inputs':{'xyz':tp}}}}
        mp=doc('manifest.json',m)
        review={'status':'passed_stage2_execution_source_review','manifest':mp,'execution_sources':m['execution_sources'],'implementations':{'build':ip,'native':ip}}
        rp=doc('review.json',review)
        cfg={'status':'bound_stage2_build_configuration','implementation':ip,'manifest':mp};cp=doc('config.json',cfg)
        rel={'status':'released_stage2_build','implementation':ip,'configuration':cp,'output':str(p/'out'),'source_review':rp};lp=doc('release.json',rel)
        def context(c=cfg,r=rel,kind='build'):
            doc('config.json',c);r=copy.deepcopy(r);r['configuration']=pin(p/'config.json');doc('release.json',r)
            with patch.dict(os.environ,{k:'1' for k in THREADS}),patch('COMMON.allocation',return_value={'CPUs':1}):
                return source_context(p/'config.json',p/'release.json',p/'out',s,kind,1)
        check('matching finite build binding accepted',context)
        check('unbound configuration rejected',lambda:context(dict(cfg,status=None)),True)
        check('unreleased source rejected',lambda:context(r=dict(rel,status=None)),True)
        check('different output rejected',lambda:context(r=dict(rel,output=str(p/'other'))),True)
        check('wrong implementation rejected',lambda:context(dict(cfg,implementation=tp)),True)
        bad=copy.deepcopy(review);bad['implementations']['build']=tp
        br=doc('bad_review.json',bad)
        check('future phase cannot inherit another implementation review',lambda:context(r=dict(rel,source_review=br)),True)
        b=doc('build.json',{'status':'stage2_build_complete','manifest':mp,'source_manifest':tp})
        nc=dict(cfg,status='bound_stage2_native_configuration',build=b,state_id='changed',central_scout_review=None)
        nl=dict(rel,status='released_stage2_native')
        check('changed state requires completed central scout',lambda:context(nc,nl,'native'),True)
        scout=doc('scout.json',{'status':'admitted_stage2_central_and_bypass_scout','build':b,'manifest':mp,'all_three_bypass_tokens_equal':True,'central_dual_energy_passed':True})
        nc['central_scout_review']=scout
        check('matching actual scout permits changed-state phase',lambda:context(nc,nl,'native'))
        c={'status':'stage2_state_nominal_certificate_passed','state_id':'central','inputs':m['rows']['central']['inputs'],'manifest':mp,'candidate':tp,'e_app':None,'application_error_bound_established':False,'sufficient_numerical_test_passed':True,'expected_state':tp,'pme_export':tp,'root_admission':stage1,'factor_root_admission':tp}
        cc={'certificate':doc('cert.json',c),'state_id':'central','manifest':mp}
        check('central certificate reuses exact admitted state',lambda:certificate(cc,m))
        def altered(**change):
            obj=dict(c,**change);return certificate(dict(cc,certificate=doc('changed_cert.json',obj)),m)
        check('wrong state certificate rejected',lambda:altered(state_id='changed'),True)
        check('failed nominal certificate rejected',lambda:altered(sufficient_numerical_test_passed=False),True)
        check('unbound application-error claim rejected',lambda:altered(application_error_bound_established=True),True)
        check('other initialized raw state rejected',lambda:altered(expected_state=ip),True)
        token.write_text('mutated\n')
        check('changed pinned input rejected',lambda:verify(tp),True)
        outmock=p/'mock';outmock.mkdir();r={'commands':{}}
        with patch('COMMON.subprocess.Popen',side_effect=OSError('manufactured pre-launch failure')):
            try:invoke(r,outmock,'native',['never-launched'],p,1,1024,1)
            except OSError:pass
        check('pre-launch failure persisted without PID',lambda:require(r['commands']['native']['status']=='failed' and r['commands']['native']['PID'] is None,'failure ledger'))
        check('same command cannot be retried',lambda:invoke(r,outmock,'native',['never-launched'],p,1,1024,1),True)
    build=(HERE/'BUILD.py').read_text();native=(HERE/'NATIVE.py').read_text();prep=(HERE/'PREPARE.py').read_text();wrapper=(HERE/'RUN.sbatch').read_text()
    for name in ['COMMON.py','BUILD.py','NATIVE.py','BIND.py','PREPARE.py']:ast.parse((HERE/name).read_text())
    check('fixed five compiles and three links',lambda:require("<=5" in build and "<=3" in build and "('gradient_old','GRADIENT.f90',False)" in build,'build counts'))
    check('single native invocation source',lambda:require(native.count("invoke(result,out,'native'")==1 and "automatic_retry':False" in native,'one native'))
    check('observer-only runtime key construction',lambda:require("key+b'\\nPME-NORM-EXPORT 1\\n'" in prep,'observer delta'))
    check('Slurm default-memory one-CPU no-requeue wrapper',lambda:require('--cpus-per-task=1' in wrapper and '--no-requeue' in wrapper and '--mem' not in wrapper and '--time' not in wrapper,'resource scope'))
    result={'status':'passed_stage2_execution_policy_controls','implementation':pin(__file__),'sources':{n:pin(HERE/n) for n in ('COMMON.py','BUILD.py','NATIVE.py','BIND.py','PREPARE.py','RUN.sbatch')},'checks':checks,'count':len(checks),'actual':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'native_calls':0,'compiler_calls':0,'assemblies':0,'factors':0,'solves':0}}
    out.mkdir(parents=True,exist_ok=False);save(out/'RESULT.json',result);print(json.dumps(result))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.out)
