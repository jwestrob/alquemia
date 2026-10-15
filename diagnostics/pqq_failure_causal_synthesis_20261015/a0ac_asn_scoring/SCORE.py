#!/usr/bin/env python3
"""Ten-source fixed Asn flip: unchanged warm worker, explicit204-atom admission."""
import argparse
import copy
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import signal
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_scoring'
INVENTORY=HERE.parent/'a0ac_asn_inventory/EVIDENCE.json'
INVENTORY_SHA='c91b152b1fb8b548395c41aac72fb2370b5ba73d015abbfce14767ac2787d4d8'
PRIOR=ROOT/'workspaces/exaa_metal_position_20261014/scoring_v1/CONFIG.json'
BASES=('origin','adaptive_Ca','adaptive_La');NEW='amide_flip';ORDER=BASES+(NEW,);METALS=('Ca','La');EV=23.06054783061903
PROTOCOL='A0ACD6B9F2_intact_Asn290_180degree_origin_flip_common4_v1'
FORCE_DEFINITIONS=('negative Cartesian total vacuum OMOL energy gradient, eV/angstrom','negative Cartesian energy gradient, eV/angstrom','negative_Cartesian_gradient_of_total_vacuum_OMOL_energy')

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':str(p),'sha256':h.hexdigest()}
def check(p):
    if pin(p['path'])!=p:raise ValueError('Changed artifact: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    s=importlib.util.spec_from_file_location(name,check(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def require(v,reason):
    if not v:raise ValueError(reason)

def prepare(preparation,output):
    """Metadata only. Preserve all ten rows, dispatch only paired admitted flips."""
    prep=read(preparation);inv=read(prep['inventory']);prior=read(PRIOR);template=read(prior['template'])
    require(prep['inventory']['sha256']==pin(INVENTORY)['sha256']==INVENTORY_SHA and len(prep['rows'])==10,'Exact ten-source preparation/inventory required')
    require({r['source_case_id'] for r in prep['rows']}=={r['source_case_id'] for r in inv['rows']},'Preparation lost or duplicated source')
    tasks=[]
    for row in prep['rows']:
        source,=[r for r in inv['rows'] if r['source_case_id']==row['source_case_id']]
        if row['status']!='prepared':continue
        for metal in METALS:
            v=row['matrix'][metal][NEW]
            require(v['status']=='prepared' and v['atom_count']==204,'Admitted source must contain both204-atom cells')
            check(v['xyz']);check(v['geometry_admission'])
            tasks.append({'task_id':row['source_case_id']+'__'+metal+'__'+NEW,'case_id':row['case_id'],
                'source_case_id':row['source_case_id'],'conditioning':source['conditioning'],
                'sample':int(row['source_case_id'].rsplit('sample-',1)[1]),'metal':metal,'candidate':NEW,
                **{k:v[k] for k in ('xyz','charge','multiplicity','atom_count','geometry_admission')},
                **{k:row[k] for k in ('maps','state_signature','source_preparation')}})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('SCORE.py','RUN.sbatch'):shutil.copyfile(HERE/name,out/name)
    worker=copy.deepcopy(template);worker.update(protocol_id=PROTOCOL,tasks=tasks,cases=[],fixtures=[],origin_task_count=len(tasks),preparation=pin(preparation),
        finite_scope={'declared_new_cells':20,'admitted_new_cells':len(tasks),'MACE_calls_maximum':20,'qualification_calls':0,'optimizer_starts':0,'SCC_calls':0,'DFT_calls':0},
        authorization='Root-authorized finite Asn290 rotamer test; actual Slurm only')
    write(out/'worker_manifest.json',worker)
    c={'protocol':PROTOCOL,'implementation':pin(out/'SCORE.py'),'wrapper':pin(out/'RUN.sbatch'),
        'plan':pin(HERE.parent/'ASN_ROTAMER_PLAN.md'),'preparation':pin(preparation),'preparation_manifest':prep['manifest'],
        'inventory':prep['inventory'],'inventory_reference':pin(INVENTORY),'rows':prep['rows'],'cells':tasks,'worker_manifest':pin(out/'worker_manifest.json'),
        'template':prior['template'],'prior_adapter':prior['implementation'],'qualification':prior['qualification'],
        'qualified_worker_summary':prior['qualified_worker_summary'],'python':template['optimizer_executable'],
        'declared_sources':10,'declared_new_cells':20,'admitted_new_cells':len(tasks),'reused_original_cells':60,'CPUs':7,
        'original_order':list(BASES),'new_order':list(ORDER),'energy_conversion_kcal_per_eV':EV,
        'operational_origin_tolerance_kcal_mol':0.1,'solvent_review_absolute_math_pool_change_kcal_mol':1.0,
        'new_SCC_calls':0,'new_DFT_calls':0,'new_optimizations':0,'output':str(out/'execution')}
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def physical_result(result,request,task,model,observed_model):
    """One204-atom scientific state predicate shared by real archives and new calls."""
    require(result['status'] in ('computed','complete') and observed_model==model,'Result/model unavailable')
    mult=request.get('multiplicity',request.get('spin_multiplicity'))
    require((request['charge'],mult)==(task['charge'],task['multiplicity']),'Request state differs')
    require(request['xyz']['sha256']==task['xyz']['sha256'],'Exact request coordinates differ')
    check(request['xyz']);check(result['forces'])
    state=result.get('state',result.get('input_state_check'))
    require(task['atom_count']==204 and (state['charge'],state['spin_multiplicity'],state['atoms'],state['graph_count'],state['head'],state['selected_metal_index'])==(task['charge'],1,204,1,'omol',0),'Actual204-atom OMOL state differs')
    require(result['force_definition'] in FORCE_DEFINITIONS and math.isfinite(result['energy_eV']),'Force convention/nonfinite energy')

def admitted_result(result,request,task,model,versions):
    physical_result(result,request,task,model,result['model'])
    require(result['status']=='complete' and (result['device'],result['threads'],result['actual_model_forwards'],result['parameter_versions_unchanged'])==('cpu',7,1,True),'New native CPU execution differs')
    require(result['versions']==versions,'New native runtime versions differ')

def archived_result(cell,task,model,versions,np):
    result=read(cell['MACE']);check(cell['xyz']);check(cell['native_frozen_ALPB'])
    if result.get('request'):
        require(result['request']==cell['saved_request'],'Archived request pin differs')
        request=read(result['request']);observed_model=result['model'];kind='request_receipt'
    else:
        require(result['status']=='computed' and cell['saved_request'] is None,'Unknown archived receipt schema')
        manifest=read(result['manifest']);request,=[t for t in manifest['tasks'] if t['task_id']==result['task_id']]
        observed_model=manifest['model'];kind='manifest_task_receipt'
    physical_result(result,request,task,model,observed_model)
    require(result['energy_eV']==cell['MACE_eV'] and result['forces']==cell['saved_forces'],'Archived energy/force identity differs')
    force=np.load(check(cell['saved_forces']),allow_pickle=False)
    require(force.shape==(204,3) and bool(np.isfinite(force).all()),'Real archived204-atom force array invalid')
    if result.get('device')=='cpu':admitted_result(result,request,task,model,versions)
    else:require(result['device']=='NVIDIA H200','Unexpected historical execution device')
    return {'status':'complete','receipt':cell['MACE'],'xyz':cell['xyz'],'energy_eV':result['energy_eV'],'forces':cell['saved_forces'],
        'request':cell['saved_request'],'historical_manifest':result.get('manifest'),'receipt_schema':kind,'historical_device':result['device'],
        'charge':task['charge'],'multiplicity':1,'atom_count':204,'reused':True}

def preflight(c,b):
    import numpy as np
    require(pin(sys.executable)==c['python'],'Qualified numerical parent interpreter differs')
    for k in ('plan','wrapper'):check(c[k])
    m=read(c['worker_manifest']);old=read(c['template']);prep=read(c['preparation']);inv=read(c['inventory']);pm=read(c['preparation_manifest'])
    check(c['inventory_reference']);require(c['inventory']['sha256']==c['inventory_reference']['sha256']==INVENTORY_SHA,'Actual and reference inventory content differs')
    require(c['protocol']==PROTOCOL and c['declared_sources']==10 and c['declared_new_cells']==20 and c['reused_original_cells']==60,'Finite protocol changed')
    require(c['original_order']==list(BASES) and c['new_order']==list(ORDER) and c['energy_conversion_kcal_per_eV']==EV and c['operational_origin_tolerance_kcal_mol']==0.1 and c['solvent_review_absolute_math_pool_change_kcal_mol']==1.0,'Pool/unit/gate changed')
    require(prep['manifest']==c['preparation_manifest'] and prep['inventory']==c['inventory'] and pm['plan']==c['plan'] and prep['rows']==c['rows'],'Actual geometry/protocol join differs')
    require(len(c['rows'])==10 and len({r['source_case_id'] for r in c['rows']})==10 and {r['source_case_id'] for r in c['rows']}=={r['source_case_id'] for r in inv['rows']},'Lost source denominator')
    expected={f'a0acd6b9f2-pqq-la_model__conditioned_{z}__seed-1_sample-{i}' for z in METALS for i in range(5)}
    require({r['source_case_id'] for r in c['rows']}==expected,'Exact five samples of both conditioning groups required')
    for k in ('implementation','model','settings','pool_settings','optimizer_software','source_modules','source_snapshot','calculator_source','model_python','optimizer_python','model_executable','optimizer_executable','gates'):require(m[k]==old[k],'Qualified worker changed: '+k)
    for p in m['source_modules'].values():check(p)
    for k in ('implementation','calculator_source','model_executable','optimizer_executable'):check(m[k])
    check(m['model']['checkpoint']);b.modules(m)
    require(m['tasks']==c['cells'] and len(c['cells'])==c['admitted_new_cells']<=20 and m['threads_per_model']==7 and m['lanes']==1 and not m['fixtures'] and not m['cases'],'Worker finite scope differs')
    qualification=read(c['qualification']);summary=read(c['qualified_worker_summary'])
    require(qualification['pass'] and summary['threads']==7 and summary['manifest']==c['template'],'Existing actual model qualification unavailable')
    for q in qualification['checks']:
        result=read(q['receipt']);check(result['forces']);read(result['request'])
        require(q['pass'] and result['model']==m['model'] and result['versions']==summary['versions'],'Existing model qualification differs')
    baseline={};checks=[];count=0
    for row in c['rows']:
        sid=row['source_case_id'];source,=[r for r in inv['rows'] if r['source_case_id']==sid]
        state=read(source['state_signature']);check(source['source_preparation']);check(source['protonated_source'])
        for p in source['maps'].values():check(p)
        require(source['atom_count']==204 and state['endpoint_charges']=={'Ca':-4,'La':-3} and state['multiplicities']=={'Ca':1,'La':1} and not state['water_inventory'],'Declared source chemical state differs')
        require(all(state[k]==source['state'][k] for k in source['state']),'Inventory state differs')
        baseline[sid]={z:{} for z in METALS};origins={}
        for z in METALS:
            for q in BASES:
                cell=source['matrix'][z][q];task={'xyz':cell['xyz'],'charge':state['endpoint_charges'][z],'multiplicity':1,'atom_count':204}
                coords=b.xyz(check(cell['xyz']));require(len(coords)==204 and coords[0][0]==z,'Archived coordinate count/metal differs')
                baseline[sid][z][q]=archived_result(cell,task,m['model'],summary['versions'],np);count+=1
                if q=='origin':origins[z]=check(cell['xyz']).read_text().splitlines(keepends=True)
        if row['status']!='prepared':
            require(not any(t['source_case_id']==sid for t in c['cells']),'Unsupported source dispatched')
            checks.append({'source_case_id':sid,'status':'unavailable','reason':row.get('reason'),'six_original_cells_validated':True});continue
        for k in ('maps','state_signature','source_preparation','protonated_source'):require(row[k]==source[k],'New geometry source provenance differs: '+k)
        gate=read(row['geometry_admission']);require(gate['pass'] and gate['status']=='prepared','Actual geometry gate failed')
        moved={x['qm_index'] for x in source['ASN_source_to_qm'] if x['source']['atom'] in ('OD1','ND2','HD21','HD22')}
        require(len(moved)==4,'Four exact amide atom identities required');paired={}
        source_tasks=[t for t in c['cells'] if t['source_case_id']==sid];require(len(source_tasks)==2 and {t['metal'] for t in source_tasks}==set(METALS),'Paired source must dispatch both metals')
        for t in source_tasks:
            z=t['metal'];v=row['matrix'][z][NEW];require(t['xyz']==v['xyz'] and t['geometry_admission']==v['geometry_admission']==row['geometry_admission'] and v['status']=='prepared','Prepared cell join differs')
            require(v['original_xyz']==source['matrix'][z]['origin']['xyz'] and v['original_MACE']==source['matrix'][z]['origin']['MACE'],'Wrong flip origin')
            lines=check(t['xyz']).read_text().splitlines(keepends=True);coords=b.xyz(check(t['xyz']))
            require(len(coords)==204 and coords[0][0]==z and [a[0] for a in coords]==[line.split()[0] for line in origins[z][2:]],'Flip atom count/order differs')
            require(all(lines[i+2]==origins[z][i+2] for i in range(204) if i not in moved),'Unchanged atom coordinate line changed')
            require(t['charge']==state['endpoint_charges'][z] and t['multiplicity']==1 and t['atom_count']==204 and all(math.isfinite(x) for atom in coords for x in atom[1:]),'Flip state/coordinates invalid')
            paired[z]=coords
        require(paired['Ca'][1:]==paired['La'][1:] and paired['Ca'][0][1:]==paired['La'][0][1:],'Paired physical coordinates differ')
        checks.append({'source_case_id':sid,'status':'passed','geometry_admission':row['geometry_admission'],'moved_indices_zero_based':sorted(moved),'six_original_cells_validated':True})
    require(count==60,'All sixty actual204-atom baseline receipts required before calls')
    return {'status':'passed','checks':checks,'baseline':baseline,'validated_real204_atom_receipts':count,'versions':summary['versions'],
        'qualification':c['qualification'],'new_model_calls':0,'historical_GPU_receipts':'Preserved actual schema, not relabeled as CPU executions'}

def choose(energies,order):
    if any(energies[z].get(q) is None for z in METALS for q in order):return None
    works={z:{q:(energies[z][q]-energies[z]['origin'])*EV for q in order} for z in METALS};out={'works_kcal_mol':works}
    for policy in ('mathematical','operational'):
        selected={}
        for z in METALS:
            best=min(order,key=lambda q:works[z][q]);selected[z]=best if policy=='mathematical' or works[z][best]<-0.1 else 'origin'
        w={z:works[z][selected[z]] for z in METALS}
        out[policy]={'selected_geometries':selected,'per_metal_work_kcal_mol':w,'Ca_minus_La_work_kcal_mol':w['Ca']-w['La'],
            'paired_vacuum_contrast_kcal_mol':(energies['Ca'][selected['Ca']]-energies['La'][selected['La']])*EV}
    return out

def span(values):
    return {'count':len(values),'minimum_kcal_mol':min(values) if values else None,'maximum_kcal_mol':max(values) if values else None,'range_kcal_mol':max(values)-min(values) if values else None}

def collect(c,pf,results):
    rows=[]
    for source in c['rows']:
        sid=source['source_case_id'];conditioning='Ca' if '__conditioned_Ca__' in sid else 'La';sample=int(sid.rsplit('sample-',1)[1]);matrix=copy.deepcopy(pf['baseline'][sid])
        for z in METALS:
            tid=sid+'__'+z+'__'+NEW
            matrix[z][NEW]=results.get(tid,{'status':'unavailable','energy_eV':None,'reason':source.get('reason') or 'not dispatched'})
        energy={z:{q:(v['energy_eV'] if v['status']=='complete' else None) for q,v in matrix[z].items()} for z in METALS}
        pools={'original3':choose(energy,BASES),'expanded4':choose(energy,ORDER)};changes={}
        for policy in ('mathematical','operational'):
            before,after=pools['original3'],pools['expanded4'];changes[policy]=None
            if before and after:
                a,b=before[policy],after[policy];d={z:b['per_metal_work_kcal_mol'][z]-a['per_metal_work_kcal_mol'][z] for z in METALS}
                changes[policy]={'per_metal_change_kcal_mol':d,'paired_contrast_change_kcal_mol':d['Ca']-d['La'],'before':a['selected_geometries'],'after':b['selected_geometries']}
        flip={z:(energy[z][NEW]-energy[z]['origin'])*EV if energy[z][NEW] is not None else None for z in METALS}
        rows.append({'source_case_id':sid,'case_id':source['case_id'],'conditioning':conditioning,'sample':sample,'preparation_status':source['status'],'preparation_reason':source.get('reason'),
          'matrix':matrix,'flip_work_kcal_mol':flip,'Ca_minus_La_flip_work_kcal_mol':flip['Ca']-flip['La'] if all(v is not None for v in flip.values()) else None,
          'pools':pools,'common3_to4':changes,'classification':None})
    grouped={}
    for condition in METALS:
        group=[r for r in rows if r['conditioning']==condition];matched=[r for r in group if r['pools']['expanded4'] is not None]
        grouped[condition]={'declared_sources':5,'complete_expanded_sources':len(matched),'matched_source_ids':[r['source_case_id'] for r in matched],
          'by_policy':{p:{'original3_all5':span([r['pools']['original3'][p]['paired_vacuum_contrast_kcal_mol'] for r in group]),
             'original3_matched':span([r['pools']['original3'][p]['paired_vacuum_contrast_kcal_mol'] for r in matched]),
             'expanded4_matched':span([r['pools']['expanded4'][p]['paired_vacuum_contrast_kcal_mol'] for r in matched])} for p in ('mathematical','operational')},
          'scope':'Any partial-coverage ranges describe the explicitly listed common sources only'}
    valid=[r for r in rows if r['common3_to4']['mathematical'] is not None]
    qualifies=any(abs(r['common3_to4']['mathematical']['paired_contrast_change_kcal_mol'])>=1.0 for r in valid)
    incomplete=any(r['preparation_status']=='prepared' and r['pools']['expanded4'] is None for r in rows)
    decision='qualifies_for_root_solvent_review' if qualifies else ('unresolved_incomplete_scoring' if incomplete else 'close_finite_flip_test')
    return {'protocol':PROTOCOL,'terminal':True,'preparation':c['preparation'],'inventory':c['inventory'],'declared_sources':10,
        'declared_new_cells':20,'admitted_new_cells':len(c['cells']),'complete_new_cells':sum(v['status']=='complete' for v in results.values()),'reused_original_cells':60,
        'rows':rows,'conditioning_spreads':grouped,'pool_orders':{'original3':list(BASES),'expanded4':list(ORDER)},
        'stage_gate':{'threshold_kcal_mol':1.0,'criterion':'Absolute mathematical common3_to4 paired contrast change >=1kcal/mol in any complete source; either sign',
           'decision':decision,'source_ids_meeting_gate':[r['source_case_id'] for r in valid if abs(r['common3_to4']['mathematical']['paired_contrast_change_kcal_mol'])>=1.0],
           'unavailable_preparations_remain_visible':True,'solvent_executed':False,'automatic_solvent_launch':False},
        'conversion_kcal_per_eV':EV,'new_SCC_calls':0,'new_DFT_calls':0,'new_optimizations':0,'production_changed':False,
        'limitations':'Finite discrete rotamer/common-pool sensitivity; raw vacuum contrast is not calibrated preference or affinity; no continuous path/minimum or physical-state validation'}

def run(config,action):
    require(bool(os.environ.get('SLURM_JOB_ID')),'Actual Slurm required before scientific imports')
    c=read(config);require(c['implementation']==pin(__file__),'Wrong scorer snapshot');prior=load('_asn_prior_placement',c['prior_adapter']);placement=prior.placement()
    out=Path(c['output'] if action=='run' else str(Path(c['output']).parent/'preflight_'+os.environ['SLURM_JOB_ID']));start=time.monotonic();cpu=time.process_time();transport=None;handlers={};results={};pf=None
    receipt={'status':'failed','config':pin(config),'placement':placement,'action':action}
    with (Path(config).parent/'RUN.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);out.mkdir(parents=True,exist_ok=False)
        try:
            m=read(c['worker_manifest']);b=load('_asn_unchanged_warm_worker',m['implementation']);pf=preflight(c,b);write(out/'PREFLIGHT.json',pf)
            if action=='run':
                if c['cells']:
                    attempt=out/'mace';(attempt/'lane0').mkdir(parents=True);handlers=b.handlers();transport=b.WarmCPU(check(c['worker_manifest']),m,attempt,0)
                    require(transport.ready['threads']==7 and transport.ready['versions']==pf['versions'],'Loaded worker runtime differs')
                for t in c['cells']:
                    request={'manifest':c['worker_manifest'],'task_id':t['task_id'],**{k:t[k] for k in ('xyz','charge','multiplicity')}}
                    value,rp=b.cell(attempt/'cells'/t['task_id'],request,transport,m['model']);r={'status':'unavailable','receipt':rp,'xyz':t['xyz'],'energy_eV':None,'reason':value.get('reason')}
                    try:
                        admitted_result(value,read(value['request']),t,m['model'],pf['versions']);r.update(status='complete',energy_eV=value['energy_eV'],forces=value['forces'],request=value['request'])
                    except Exception as e:r['reason']=str(e)
                    results[t['task_id']]=r
                if transport:transport.close();transport=None
                receipt['status']='complete' if len(results)==len(c['cells']) and all(r['status']=='complete' for r in results.values()) else 'partial'
            else:receipt['status']='complete'
        except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc());raise
        finally:
            try:write(out/'PARTIAL_RESULTS.json',{'config':pin(config),'rows':results,'declared_cells':20,'admitted_cells':len(c['cells'])})
            except BaseException as e:receipt.update(status='failed',partial_reporting_failure={'reason':str(e),'traceback':traceback.format_exc()})
            try:
                if transport:transport.abort()
            except BaseException as e:receipt.update(status='failed',cleanup_failure={'reason':str(e),'traceback':traceback.format_exc()})
            try:
                for sig,handler in handlers.items():signal.signal(sig,handler)
                if pf and action=='run':write(out/'FINAL_COLLECTION.json',{'config':pin(config),**collect(c,pf,results)})
            except BaseException as e:receipt.update(status='failed',collection_failure={'reason':str(e),'traceback':traceback.format_exc()});raise
            finally:
                elapsed=time.monotonic()-start;receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=7*elapsed,parent_process_CPU_seconds=time.process_time()-cpu)
                try:
                    receipt.update(actual_model_forwards=sum(read(p).get('actual_model_forwards',0) for p in out.glob('mace/cells/*/mace_result.json')),
                        worker_resource_receipts=[pin(p) for p in out.glob('mace/lane*/CPU_WORKER_SUMMARY.json')])
                except BaseException as e:receipt.update(status='failed',resource_reporting_failure={'reason':str(e),'traceback':traceback.format_exc()})
                write(out/'EXECUTION.json',receipt)
    if receipt['status']=='failed':raise RuntimeError('Execution or reporting failed; inspect preserved EXECUTION.json')

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    q=sub.add_parser('prepare');q.add_argument('--preparation',required=True);q.add_argument('--output',required=True)
    for name in ('preflight','run'):
        q=sub.add_parser(name);q.add_argument('--config',required=True);q.add_argument('--config-sha',required=True)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.preparation,a.output)
    else:run(check({'path':str(Path(a.config).resolve()),'sha256':a.config_sha}),a.action)
