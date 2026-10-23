"""Narrow scheduling/lineage controls; no native, arrays, probes or scheduler."""
import ast
import copy
import importlib.util
import json
from pathlib import Path
import time
import PREPARE as P
import BIND as B


def run(prepared,out):
    tick,cpu=time.monotonic(),time.process_time();out.mkdir(parents=True,exist_ok=False)
    r={'status':'started','implementation':P.pin(__file__),'checks':[],'error':None,
       'counts':{'native_compiler_solver_calls':0,'numerical_arrays':0,'scheduler_queries_or_mutations':0}}
    def check(label,ok):assert ok,label;r['checks'].append(label)
    def rejects(label,fn):
        try:fn()
        except (ValueError,AssertionError):check(label,True)
        else:raise AssertionError(label)
    try:
        a=json.loads((P.HERE/'AUTHORITIES.json').read_text());oldready=P.read(a['parent_ready']);old=P.read(oldready['inputs'])
        x=P.read(P.pin(prepared/'INPUTS.json'));m=P.read(x['central_manifest']);probe=P.read(x['probe_addendum'])
        before=P.verify(a['parent_binder']).read_text();after=(P.HERE/'BIND.py').read_text();restored=after
        for oldtext,newtext in reversed(P.REPLACEMENTS):
            assert restored.count(newtext)==1;restored=restored.replace(newtext,oldtext,1)
        check('inverse binder patch restores exact original bytes',restored==before)
        old_ast,new_ast=ast.parse(before),ast.parse(after)
        changed={'validate_build','bound_manifest','context'}
        check('all binder AST outside three declared metadata functions unchanged',len(old_ast.body)==len(new_ast.body) and all(ast.dump(l)==ast.dump(rr) for l,rr in zip(old_ast.body,new_ast.body) if not isinstance(l,ast.FunctionDef) or l.name not in changed))
        for key,cpus in [('slurm1',1),('slurm12',12)]:
            ot=P.verify(m['helpers'][key]).read_text();nt=P.verify(x['standard_helpers'][key]).read_text()
            expected=ot.replace(P.QOS_OLD,P.QOS_NEW,1)
            if cpus==1:expected=expected.replace("qos='standard_expand'","qos='standard'",1)
            check(key+' source only changes expected/reported QOS',nt==expected)
            spec=importlib.util.spec_from_file_location('mock_'+key,x['standard_helpers'][key]['path']);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
            env={'SLURM_JOB_ID':'123','SLURM_CPUS_PER_TASK':str(cpus),'SLURM_NTASKS':'1','SLURM_JOB_NUM_NODES':'1'}
            raw=f'JobId=123 JobState=RUNNING Account=banfield QOS=standard Partition=standard NumNodes=1 NumCPUs={cpus} NumTasks=1 CPUs/Task={cpus} Requeue=0 OverSubscribe=NO BatchHost=mock UserId=mock(601) AllocTRES=cpu={cpus},mem={cpus*4000}M,node=1 MinMemoryCPU=4000M'
            rec=mod.parse_allocation(raw,env,'mock',601)
            check(key+' accepts exact standard and reports truthful QOS',rec['qos']=='standard' and rec['allocated_CPUs']==cpus)
            rejects(key+' rejects expand',lambda:mod.parse_allocation(raw.replace('QOS=standard ','QOS=standard_expand '),env,'mock',601))
        for target,authority in [('RUN.sbatch',oldready['wrapper']),('CERTIFICATE.sbatch',a['old_certificate_wrapper'])]:
            text=P.verify(authority).read_text();new=(P.HERE/target).read_text()
            check(target+' exact QOS-only wrapper delta',new==text.replace('#SBATCH --qos=standard_expand\n','#SBATCH --qos=standard\n',1))
        for name,p in m['execution_sources'].items():P.verify(p)
        check('all PR execution/native/source identity pins preserved',x['native_implementation']==old['native_implementation'] and x['central_manifest']==old['central_manifest'] and x['runtime_ready']==old['runtime_ready'])
        req=json.loads((P.HERE/'REQUEST_TEMPLATE.json').read_text());fake={'path':'/manufactured/review.json','sha256':'a'*64};req['binder_source_review']=fake
        before_m=copy.deepcopy(m);bound=B.bound_manifest(m,req,probe,x)
        check('both helpers bound before initializer without changing rows/model/native or execution',m==before_m and bound['helpers']['slurm1']==x['standard_helpers']['slurm1'] and bound['helpers']['slurm12']==x['standard_helpers']['slurm12'] and all(bound[k]==m[k] for k in ('rows','model','native_sources','execution_sources','identity_sources')))
        check('explicit scheduling authority/review and fixed probes retained',bound['scheduling_recovery']['declaration']==x['scheduling_recovery_declaration'] and bound['scheduling_recovery']['source_review']==fake and bound['probe_data']==probe['probe_data'] and bound['probe_source']==probe['probe_source'])
        check('no initializer certificate parity or scout reused',all(bound[k] is None for k in ('runtime_release','initializer','certificate','central_parity_admission','central_scout_admission')))
        build=P.read(a['build']);admission=P.read(a['build_admission']);submission=P.read(a['build_submission'])
        B.validate_build(build,admission,req,submission,m)
        check('actual admitted20813 metadata joins exact release field and four roles',True)
        wrong=copy.deepcopy(build);wrong['allocation']['job_id']='20739'
        rejects('cancelled old20739 cannot be bound',lambda:B.validate_build(wrong,admission,req,submission,m))
        wrongsub=copy.deepcopy(submission);wrongsub['release']=fake
        rejects('different replacement release rejected',lambda:B.validate_build(build,admission,req,wrongsub,m))
        restored_inputs=copy.deepcopy(x)
        for k in ('implementation','build_submission','sources'):restored_inputs[k]=old[k]
        for k in ('scheduling_recovery_declaration','admitted_build','admitted_build_review','standard_helpers','scheduling_recovery_evidence'):del restored_inputs[k]
        check('all remaining original binder inputs exactly preserved',restored_inputs==old)
        check('request remains unbound pending review/output with exact admitted replacement',req['status']=='unbound_single_central_initializer_binding' and req['native_output'] is None and req['build']==a['build'] and req['build_admission']==a['build_admission'])
        check('original twenty controls reused without rerun',P.read(oldready['controls'])['check_count']==20)
        r.update(status='passed_standard_central_initializer_recovery_controls',inputs=P.pin(prepared/'INPUTS.json'),preparation=P.pin(prepared/'RESULT.json'),old_controls_reused=oldready['controls'])
    except BaseException as exc:r.update(status='failed_standard_central_initializer_recovery_controls',error=repr(exc))
    r.update(check_count=len(r['checks']),cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1})
    P.write(out/'RESULT.json',r);print(json.dumps({k:r[k] for k in ('status','check_count','error','cost')}));return 0 if r['error'] is None else 1


if __name__=='__main__':
    import sys
    raise SystemExit(run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve()))
