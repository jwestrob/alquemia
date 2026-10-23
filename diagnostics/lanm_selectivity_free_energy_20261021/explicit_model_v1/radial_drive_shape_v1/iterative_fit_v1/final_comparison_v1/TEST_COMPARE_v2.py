#!/usr/bin/env python3
"""Manufactured scalar/persistence checks only. Never open the running campaign."""
import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import time
import COMPARE_v2 as C


def main(out):
    require=lambda x,m: C.require(x,m)
    require(all(os.environ.get(k)=='1' for k in C.THREADS),'One test thread')
    out=Path(out);out.mkdir(parents=True,exist_ok=False);tick,cpu=time.monotonic(),time.process_time();passed=[]
    def check(name,fn):fn();passed.append(name)
    def rejected(fn):
        try:fn()
        except (ValueError,RuntimeError):return
        raise AssertionError('Expected rejection')
    fixed=dict(q_denominators=dict(J_P=1.,J_B=1.),electrical_actual_scores=dict(pair_mu=1.,C_mu=1.,chi=1.),signed_W9_reference_Debye=-.1,signed_W9_absolute_ceiling_Debye=.25)
    gates={**{'energy/'+str(i):dict(error_kcal_mol=.1,passed=True) for i in range(10)},**{'force/'+str(i):dict(rigid_RMS_kcal_mol_A=1.,maximum_body_component_error_kcal_mol_A=2.,passed=True) for i in range(15)}}
    metrics=dict(component_errors_kcal_mol=dict(P8=.6,B8=.8,P9=.6,B9=.8),J_P=.36,J_B=.64,q=.64,original25=gates,all25_original_physical_gates=True,physical_constraints_pass=True,
                 original_block_scores=dict(old_EF=.2,pair_mu=.4,C_mu=.4,chi=.4),fixed_electrical_safeguards={k:dict(score=.4,ceiling=1.,passed=True) for k in C.ELECTRICAL},
                 signed_W9_z_safeguard=dict(model_Dy_minus_La_z_Debye=.1,reference_Dy_minus_La_z_Debye=-.1,model_minus_reference_Debye=.2,absolute_ceiling_Debye=.25,passed=True),all_declared_safeguards_pass=True)
    check('q and exact25 versus three electrical plus signed safeguard',lambda:require(C.decisions(metrics,fixed)['all_declared_safeguards_pass'],'Combined'))
    bad=copy.deepcopy(metrics);bad['q']=.65
    check('wrong q arithmetic rejected',lambda:rejected(lambda:C.decisions(bad,fixed)))
    lastbit=copy.deepcopy(metrics);lastbit['q']=C.math.nextafter(.64,0.)
    check('producer literal q survives tolerated reconstruction difference',lambda:require(C.decisions(lastbit,fixed)['q']==lastbit['q'] and C.decisions(lastbit,fixed)['independent_reconstruction']['q']!=lastbit['q'],'Exact decision authority'))
    edge=copy.deepcopy(metrics);edge['original25']['energy/0']['error_kcal_mol']=1.
    check('exact physical boundary passes',lambda:require(C.decisions(edge,fixed)['physical25_pass'],'Boundary'))
    outside=copy.deepcopy(edge);outside['original25']['energy/0']['error_kcal_mol']=1.+1e-12
    check('no arithmetic closeness relaxation of physical gate',lambda:rejected(lambda:C.decisions(outside,fixed)))
    electrical=copy.deepcopy(metrics);electrical['original_block_scores']['chi']=2.;electrical['fixed_electrical_safeguards']['chi']=dict(score=2.,ceiling=1.,passed=False);electrical['all_declared_safeguards_pass']=False
    check('electrical failure does not relabel physical25',lambda:require(C.decisions(electrical,fixed)['physical25_pass'] and not C.decisions(electrical,fixed)['all_declared_safeguards_pass'],'Separate meanings'))
    signed=copy.deepcopy(metrics);signed['signed_W9_z_safeguard'].update(model_minus_reference_Debye=.3,passed=False);signed['all_declared_safeguards_pass']=False
    check('absolute signed residual ceiling enforced',lambda:require(not C.decisions(signed,fixed)['all_declared_safeguards_pass'],'Signed'))
    slots=[dict(fraction=f,eligible=True,execution_status='complete',actual=dict(q=q,all_declared_safeguards_pass=True)) for f,q in zip(C.FRACTIONS,(.8,.8,.9,.95))]
    check('earlier fraction wins ties',lambda:require(C.recorded_choice(1.,slots)[0]==0,'Earlier'))
    tie=copy.deepcopy(slots)
    for x in tie:x['actual']['q']=1.-5e-11
    check('current epsilon wins near tie',lambda:require(C.recorded_choice(1.,tie)[0] is None,'Current'))
    excluded=copy.deepcopy(slots);excluded[0]['actual']['all_declared_safeguards_pass']=False
    check('better unsafeguarded point excluded',lambda:require(C.recorded_choice(1.,excluded)[0]==1,'Safeguards'))
    incomplete=copy.deepcopy(slots);incomplete[-1]['execution_status']='not_attempted'
    check('incomplete eligible batch cannot select',lambda:rejected(lambda:C.recorded_choice(1.,incomplete)))
    header='JobIDRaw|State|ExitCode|ElapsedRaw|AllocCPUS|TotalCPU|MaxRSS\n'
    sacct=header+'20605|COMPLETED|0:0|100|13|00:12:00|\n20605.batch|COMPLETED|0:0|100|13|00:12:00|321K\n'
    check('parent cost alone no child double count',lambda:require(C.accounting(sacct,'20605')['parent_allocated_CPU_seconds']==1300,'Parent'))
    check('running scheduler evidence rejected',lambda:rejected(lambda:C.accounting(sacct.replace('COMPLETED','RUNNING'),'20605')))
    check('wrong CPU grant rejected',lambda:rejected(lambda:C.accounting(sacct.replace('|13|','|12|'),'20605')))
    domain=[dict(name='p'+str(i),lower=-1.,upper=1.) for i in range(14)];v=[0.]*14;v[7]=-1.
    check('parameter coordinate boundary is literal and descriptive',lambda:require(C.parameter_report(v,[0.]*14,domain)['box_coordinates'][7]['at_literal_boundary'],'Boundary'))
    dependencies=json.loads((Path(__file__).parent/'DEPENDENCIES.json').read_text());jp=dependencies['journal'];s=importlib.util.spec_from_file_location('test_saved_journal',C.verify(jp));J=importlib.util.module_from_spec(s);s.loader.exec_module(J)
    with tempfile.TemporaryDirectory() as tmp:
        tmp=Path(tmp);raw=tmp/'native.out';raw.write_text('manufactured failed output\n');inputs={k:C.pin(raw) for k in ('key','xyz','parameters')}
        cfg=dict(adapter=dict(path='manufactured',sha256='a'*64),root_release=dict(path='manufactured-release',sha256='b'*64),limits=dict(new_parameter_vectors=228,new_native_programs=7280))
        result=dict(configuration_fingerprint='c'*64,job_id='20605',native_programs=1)
        ledger=dict(native_program_invocations=0,native_returncodes_observed=0,vector_keys=[],requests=[],stopped=False,limits=cfg['limits'])
        values=[0.]*14;request=dict(values=values,vector_key=C.token(values),role=dict(parameters=values,role='trial',outer_step=1),status='started',required_new_programs=1,cells={'c'+str(i):dict(status='not_attempted') for i in range(46)})
        attempt=dict(number=1,cache_key='k',cell_id='c0',program='analyze',folder=str(tmp),status='invoke_boundary_entered')
        events=[dict(kind='initialize',data=dict(ledger=ledger,policy=J.POLICY,authority=dict(implementation=cfg['adapter'],configuration_fingerprint=result['configuration_fingerprint'],root_release=cfg['root_release']))),
                dict(kind='request_start',data=dict(index=0,request=request)),dict(kind='attempt_intent',data=dict(attempt=attempt,request_index=0,program_inputs=inputs)),
                dict(kind='attempt_return',data=dict(attempt=dict(attempt,status='failed_native_attempt',error='manufactured stop',returncode_observed=True),receipt_delta=[dict(returncode=1,output=C.pin(raw))],program_inputs=inputs)),
                dict(kind='request_update',data=dict(index=0,updates=dict(status='failed',error='manufactured stop'))),dict(kind='stop',data=dict(error='manufactured stop'))]
        path=tmp/'EVENTS.jsonl';previous='0'*64;lines=[]
        for i,event in enumerate(events,1):
            body=dict(sequence=i,previous_sha256=previous,event=event);previous=J.sha(body);lines.append(J.encode(dict(body,sha256=previous))+b'\n')
        path.write_bytes(b''.join(lines));_,summary=C.stream_journal(path,J,cfg,result,False)
        check('failed attempt remains charged exactly once',lambda:require(summary['charged_program_attempts']==1 and summary['observed_returns']==1 and summary['returned_success']==0,'Failure charged'))
        check('legacy673/24000 not needed for current journal caps',lambda:require(summary['unique_new_vectors']==1,'Current caps'))
        path.write_bytes(b''.join(lines)+b'{partial');_,partial=C.stream_journal(path,J,cfg,result,False)
        check('partial committed tail preserved never repaired',lambda:require(partial['partial_suffix_bytes']==8 and path.read_bytes().endswith(b'{partial'),'Partial'))
        path.write_bytes(b''.join(lines)+lines[-1]);check('duplicate/reordered journal event rejected',lambda:rejected(lambda:C.stream_journal(path,J,cfg,result,False)))
    answer=dict(status='passed_manufactured_terminal_comparison_controls',implementation=C.pin(C.__file__),test_implementation=C.pin(__file__),checks=passed,count=len(passed),
                actual_campaign_reads=0,actual_scientific_vectors=0,native_calls=0,solver_calls=0,Jacobian_replays=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
    C.save(out/'RESULT.json',answer);print(json.dumps(dict(result=C.pin(out/'RESULT.json'),count=len(passed))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args()
    try:main(a.out)
    except Exception as e:
        dest=Path(a.out);dest.mkdir(parents=True,exist_ok=True)
        failure=dest/'FAILURE.json'
        if not failure.exists():C.save(failure,dict(status='failed_manufactured_control',error=repr(e),implementation=C.pin(C.__file__),test_implementation=C.pin(__file__)))
        raise
