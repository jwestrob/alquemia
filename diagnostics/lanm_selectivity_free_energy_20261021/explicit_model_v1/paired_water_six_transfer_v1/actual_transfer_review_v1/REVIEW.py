#!/usr/bin/env python3
"""Saved-only actual17612 raw/control/score replay and changes versus17578."""
import importlib.util
import json
import os
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4]
BASE=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_six_transfer_v1'
RUN=BASE/'run_v1';OUT=BASE/'actual_transfer_review_v1'
spec=importlib.util.spec_from_file_location('actual_six_transfer',HERE.parent/'RUN.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
start,cpu=time.monotonic(),time.process_time();checks=[]
def check(ok,name):
    checks.append({'name':name,'passed':bool(ok)})
    if not ok:raise AssertionError(name)
check(all(os.environ.get(k)=='1'for k in m.d.e.THREADS),'one-thread saved review')
rp=m.pin(RUN/'execution/RESULT.json');r=m.read(rp)
check(r['status']=='consumed_six_parameter_transfer_complete' and r['error']is None and not r['global_stop'] and r['job_id']=='17612','actual terminal17612 admitted')
OUT.mkdir(parents=True,exist_ok=False)
bound=m.read(r['bound']);release=m.read(r['release']);data=m.read(bound['data']);previous=m.read(bound['previous_transfer_result'])
check(r['implementation']==m.pin(HERE.parent/'RUN.py')==bound['implementation']==release['implementation'],'exact executed source')
check(release['status']=='approved_consumed_six_parameter_water_transfer' and release['bound']==r['bound'],'exact root release')
for key in ('candidate_result','candidate_config','fit_parent_result','specification','renderer','pin_inventory'):m.verify(bound[key])
check(set(r['rows'])==set(m.d.CELLS) and all(x['status']=='complete'for x in r['rows'].values()),'complete32 denominator')
check(set(bound['reuse'])==set(m.WATER+m.d.DRY),'eight selected water and six unchanged dry reused')
check(r['stage_A_passed'] and all(r['candidate_flags'][k]for k in ('development_checks_passed','selected_score_development_passed')),'water qualification gate passed')
native=m.module(data['config']['native_helper'],'actual_transfer_native');body=m.module(data['config']['body_helper'],'actual_transfer_body');score=m.module(bound['score_source'],'actual_transfer_score')
native.invoke=lambda *a,**k: (_ for _ in ()).throw(AssertionError('No native calls in saved review'))
receipts=[];cartesian={};reused=[];fresh=[]
for gid,row in r['rows'].items():
    original=data['rows'][gid];inp=bound['inputs'][gid]['files']
    for p in inp.values():m.verify(p)
    m.replay_native(row['model'],original,native);check(True,gid+':raw E/F replay exact')
    if gid in bound['reuse']:
        reuse=bound['reuse'][gid];m.verify(reuse['authority'])
        check(row['model']==reuse['value'],gid+':exact reused model and provenance');reused.append(gid)
    else:fresh.append(gid)
    if gid in m.d.DRY:
        check(row['model']==original['baseline']==previous['rows'][gid]['model'],gid+':dry model unchanged from original and17578')
        check(all(inp[k]['sha256']==original['input'][k]['sha256']for k in ('xyz','key','parameters')),gid+':dry physical input unchanged')
    else:
        check(row['input']==bound['inputs'][gid],gid+':bound input exact')
        for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
            check(m.pin(RUN/'execution'/gid/n)['sha256']==inp[k]['sha256'],gid+':executed input '+k)
    if gid in m.d.FD:
        saved=row['cartesian_check'];actual=m.gradient(m.verify(saved['output']),original['geometry']['atoms'],True,native)
        check(actual==saved and actual['fd_passed'] and actual['fd_total_separation_A']==1e-5 and actual['fd_each_side_displacement_A']==5e-6,gid+':all original Cartesian components replay/pass')
        check(abs(actual['energy']-row['model']['energy'])<=.0001 and max(abs(x-y)for a,b in zip(actual['forces'],row['model']['forces'])for x,y in zip(a,b))<=.001,gid+':Cartesian and admitted E/F agree')
        cartesian[gid]={'maximum_component_error_kcal_mol_A':actual['max_component_error_kcal_mol_A'],'atoms':original['geometry']['atoms'],'output':actual['output']}
    for receipt in row['receipts']:
        program=Path(receipt['argv'][0]).name;fd=program=='testgrad' and gid in m.d.FD
        args=['system.xyz','E']if program=='analyze'else['system.xyz','Y','Y','0.00001']if fd else['system.xyz','Y','N']
        check(receipt['returncode']==0 and program in ('analyze','testgrad') and receipt['argv']==[data['config']['executables'][program]['path']]+args and set(receipt['threads'].values())=={'1'},gid+':native command/threads/return')
        m.program_admission(m.verify(receipt['output']),program,native);receipts.append(receipt)
check(len(fresh)==18 and len(reused)==14 and len(cartesian)==10 and len(receipts)==40==r['native_program_invocations'],'18 new/14 reused/10 FD/40 programs')
check(r['internal_FD_energies_requested']==r['internal_FD_energies_completed']==1380==sum(6*x['atoms']for x in cartesian.values()),'1380 internal displaced energies')
now=score.evaluate(data,{k:v['model']for k,v in r['rows'].items()},body)
old=score.evaluate(data,{k:v['model']for k,v in previous['rows'].items()},body)
check(now==r['comparison'],'all new comparison arithmetic exact replay')
check(old==previous['comparison']==r['previous_candidate_comparison'],'all17578 arithmetic exact replay')
groups={}
for group in ('fitted_centers','water_radial','hydrated_donor','dry_baseline'):
    groups[group]={}
    for kind in ('energy','rigid'):
        a=[x for x in old['development_checks']if x['group']==group and x['kind']==kind];b=[x for x in now['development_checks']if x['group']==group and x['kind']==kind]
        check([x['id']for x in a]==[x['id']for x in b],group+':'+kind+':unchanged denominator')
        groups[group][kind]={'count':len(b),'failed17578':[x['id']for x in a if not x['candidate_passed']],'failed17612':[x['id']for x in b if not x['candidate_passed']]}
energy=[]
for kind in ('radial_work','balanced_exchange','matched_water_addition'):
    old_index={x['id']:x for x in old[kind]}
    for row in now[kind]:
        a=old_index[row['id']]['errors']['candidate'];b=row['errors']['candidate']
        energy.append({'kind':kind,'id':row['id'],'group':row['group'],'error17578_kcal_mol':a['signed_kcal_mol'],'error17612_kcal_mol':b['signed_kcal_mol'],'absolute_error_change_kcal_mol':abs(b['signed_kcal_mol'])-abs(a['signed_kcal_mol']),'passed17578':a['passed'],'passed17612':b['passed']})
force=[];old_index={(x['geometry'],x['kind']):x for x in old['rigid_response']}
for row in now['rigid_response']:
    a=old_index[(row['geometry'],row['kind'])]['models']['candidate'];b=row['models']['candidate']
    force.append({'geometry':row['geometry'],'kind':row['kind'],'group':row['group'],'RMS17578_kcal_mol_A':a['rigid_RMS_kcal_mol_A'],'RMS17612_kcal_mol_A':b['rigid_RMS_kcal_mol_A'],'max17578_kcal_mol_A':a['maximum_body_component_error_kcal_mol_A'],'max17612_kcal_mol_A':b['maximum_body_component_error_kcal_mol_A'],'passed17578':a['passed'],'passed17612':b['passed']})
check(all(x['RMS17578_kcal_mol_A']==x['RMS17612_kcal_mol_A'] and x['max17578_kcal_mol_A']==x['max17612_kcal_mol_A']for x in force if x['group']=='dry_baseline'),'all dry force failures identically unchanged')
event=json.loads((RUN/'completion_event.json').read_text())['records']['17612']
check(event['state']=='COMPLETED' and event['exit_code']=='0:0' and int(event['allocated_cpus'])==8,'terminal scheduler admission')
for key in ('candidate_result','candidate_config','previous_transfer_result','data'):m.verify(bound[key])
m.verify(rp)
result={'status':'passed_actual_raw_numerical_and_comparison_review','implementation':m.pin(__file__),'actual_result':rp,'bound':r['bound'],'release':r['release'],'previous17578':bound['previous_transfer_result'],'check_count':len(checks),'checks':checks,'counts':{'complete_cells':32,'new_energy_force_cells':18,'reused_energy_force_cells':14,'Cartesian_checks':10,'native_programs':40,'internal_FD_energies':1380},'Cartesian_checks':cartesian,'grouped_failures':groups,'energy_changes_vs17578':energy,'force_changes_vs17578':force,'all_development_checks_passed':now['all_development_checks_passed'],'dry_failures_unchanged':True,'scientific_scope':'Consumed diagnostic transfer against unchanged provisional PBE0 targets; numerical consistency does not qualify physical accuracy or exchange free energy.','cost':{'scheduler_elapsed_seconds':int(event['elapsed_seconds']),'allocated_CPUs':8,'allocated_CPU_seconds':8*int(event['elapsed_seconds']),'runner_wall_seconds':r['elapsed_seconds'],'parent_CPU_seconds':r['parent_cpu_seconds'],'aggregate_native_user_CPU_seconds':r['aggregate_child_user_cpu_seconds'],'aggregate_native_system_CPU_seconds':r['aggregate_child_system_cpu_seconds'],'overlapping_per_invocation_CPU_counters_not_summed':True},'review_cost':{'wall_seconds':time.monotonic()-start,'parent_CPU_seconds':time.process_time()-cpu,'maximum_threads':1,'new_native_calls':0},'reserved_validation_opened':False,'parameters_fit_or_changed':False}
m.save(OUT/'REVIEW.json',result)
print(json.dumps({'review':m.pin(OUT/'REVIEW.json'),'status':result['status'],'checks':len(checks),'groups':{g:{k:[len(v['failed17578']),len(v['failed17612']),v['count']]for k,v in z.items()}for g,z in groups.items()},'cost':result['review_cost']}))
