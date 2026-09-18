"""Source-mapped established H normalization for a matched vacuum hybrid."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,cache_key,energy,paired,read_json,record,verify,write_new,xyz
from affordable_response import source_key
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL,accepted_attempt,check_atoms,write_xyz
from mace_hydrogen import repair as repair_hydrogen
from mace_omol_ablation import ADAPTER,COMPONENT
from mace_omol_ablation_run import descriptor_model,SEMANTICS
from mace_omol_intact import energy_task,EVALUATION
from mace_omol_panel import ADAPTER_TASK
from mace_omol_vacuum import METHOD,scientific_input,parse_endpoint

PROTOCOL='masked_omol_matched_normalized_H_vacuum_hybrid_v1'
STAGE='matched_H_core'
AUDIT_SHA='94b7b1e5bf6186dce8ed4c6b5a60e158aed80177d78e55798b0284c900919ca2'
OLD_SHA='1d720cb8703f1137c215709acede071db1f4c58d4eb207600c828962cde0d6f3'
PREP_SHA='361d2ac85b93360c8560a761e0d76d622b8f314a7c238013f8fb7faa184ae99c'
CASES=('GGR_extended','GGR_connected','ALPHA_1F6S','ALPHA_6IP9')
COUNTS={'GGR_extended':21,'GGR_connected':47,'ALPHA_1F6S':16,'ALPHA_6IP9':18}
TOL={'accounting_model_kcal':.01,'algebra_kcal_scale':1e-7,'partition_kcal_scale':2.,'ordering_kcal_scale':.02}
UNAVAILABLE={'aqueous_score':None,'calibrated_class':None,'predictive_improvement_claimed':False,
 'response_status':'response_model_not_validated','relaxation_correction_kcal_mol':None,'baseline_changed':False}


def normalized_physical(prep,original):
    """Replay the already frozen radial rule; no forcefield or geometry search."""
    for key in ('source','source_preparation','forcefield','water_forcefield'):verify(prep[key])
    if prep['preparation_details']['terminal_additions']:raise InvalidArtifact('unexpected terminal additions')
    rebuilt,moves=repair_hydrogen(original,prep['preparation_details']['bonds'])
    byid={a['id']:a for a in rebuilt};wanted={a['id']:a for a in prep['physical_atoms']}
    if len(byid)!=len(rebuilt) or set(byid)!=set(wanted):raise InvalidArtifact('physical identity inventory changed')
    if len(moves)!=len(prep['preparation_details']['protein_H_moves']):raise InvalidArtifact('protein H inventory differs')
    for move in prep['water_H_moves']:
        atom=byid[move['id']]
        if atom['element']!='H' or atom['xyz_A']!=move['before_A']:raise InvalidArtifact('water H source differs')
        oid=move['id'].rsplit('/',1)[0]+'/O';oxygen=byid[oid]
        vector=np.array(atom['xyz_A'])-oxygen['xyz_A']
        after=np.array(oxygen['xyz_A'])+vector*move['new_length_A']/np.linalg.norm(vector)
        if not np.allclose(after,move['after_A'],atol=1e-12,rtol=0):raise InvalidArtifact('recorded water H rule differs')
        atom['xyz_A']=move['after_A']
    errors=[]
    for old,new in zip(original,rebuilt):
        target=wanted[new['id']]
        if new['element']!=target['element']:raise InvalidArtifact('normalized element differs')
        if old['element']!='H' and old['xyz_A']!=target['xyz_A']:raise InvalidArtifact('heavy/metal coordinate changed')
        error=float(np.max(np.abs(np.array(new['xyz_A'])-target['xyz_A'])));errors.append(error)
        if error>1e-12:raise InvalidArtifact('normalized H arithmetic does not replay')
    return {'maximum_replay_error_A':max(errors),'protein_H_count':len(moves),'water_H_count':len(prep['water_H_moves']),
            'heavy_coordinates_unchanged':True,'chemical_state_unchanged':True}


@cached_file_checks
def source_data(prepared,audit,old_report):
    if record(prepared)['sha256']!=PREP_SHA or record(audit)['sha256']!=AUDIT_SHA or record(old_report)['sha256']!=OLD_SHA:
        raise InvalidArtifact('declared normalized-H source pins differ')
    p=read_json(prepared);a=read_json(audit);old=read_json(old_report)
    if old['status']!='complete' or not old['numerical_gate_pass']:raise InvalidArtifact('actual old-H comparison absent')
    collection=read_json(verify(a['collection']));mp=verify(collection['manifest']);parent=read_json(mp)
    for pin in parent['implementation'].values():verify(pin)
    whole={};physical={};checks={}
    for gid,row in a['rows'].items():
        prep=read_json(verify(row['preparation']));original=read_json(verify(p['physical_cases'][gid]['physical_atoms']))
        if row['preparation']!=p['physical_cases'][gid]['original_global_preparation']:raise InvalidArtifact('whole source differs')
        checks[gid]=normalized_physical(prep,original);physical[gid]=(prep,original);whole[gid]={}
        for metal,ref in row['exact_existing_masked_whole_endpoints'].items():
            e=prep['endpoints'][metal];coords=xyz(verify(e['xyz']));matches=[]
            if ref['normalized_input']!=e['xyz'] or e['state']!=check_atoms(coords,e['charge']):raise InvalidArtifact('whole normalized state differs')
            expected=[(metal if x['element']=='M' else x['element'],*x['xyz_A']) for x in prep['physical_atoms']]
            if coords!=expected:raise InvalidArtifact('whole physical coordinate inventory differs')
            for entry in ref['matches']:
                t=next(t for t in parent['tasks'] if t['task_id']==entry['task_id']);r=collection['rows'][t['task_id']]
                if (r!=entry['result'] or xyz(verify(t['xyz']))!=coords or t['charge']!=e['charge'] or t['spin_multiplicity']!=1
                    or t['energy_component']!=COMPONENT or t['charge_feature_adapter']!=ADAPTER):
                    raise InvalidArtifact('reused whole task/physical state differs')
                if not any(accepted_attempt(d,t,mp)==r for d in (mp.parent/'execution'/t['task_id']).glob('attempt_*')):
                    raise InvalidArtifact('whole reuse lacks accepted receipt')
                if abs(r['native_readout']['component_sum_error_kcal_mol'])>TOL['accounting_model_kcal']:
                    raise InvalidArtifact('whole reuse readout failed')
                matches.append(entry)
            if not matches:raise InvalidArtifact('required whole reuse missing')
            spread=(max(x['result']['energy_eV'] for x in matches)-min(x['result']['energy_eV'] for x in matches))*EV_TO_KCAL
            if spread>TOL['accounting_model_kcal']:raise InvalidArtifact('duplicate actual whole outputs disagree')
            chosen=sorted(matches,key=lambda x:x['task_id'])[0]
            whole[gid][metal]={'collection':a['collection'],'selected':chosen,'all_matching_task_ids':sorted(x['task_id'] for x in matches),
                              'duplicate_spread_model_kcal':spread,'normalized_preparation':row['preparation'],'state':e}
    cores={}
    for name in CASES:
        c=read_json(verify(p['cases'][name]));prep,original=physical[c['global_id']]
        bykey={x['source_key']:x for x in original};byid={x['id']:x for x in prep['physical_atoms']};endpoints={};mapping=[]
        for metal in ('Ca','La'):
            e=c['grids']['center']['endpoints'][metal];q=xyz(verify(e['xyz']));new=list(q);seen={0};moves=[]
            for item in c['source_graph']['source_to_qm']:
                i=item['qm_index']
                if i in seen:raise InvalidArtifact('duplicate QM mapping')
                seen.add(i)
                if item['kind']=='source':
                    orig=bykey[source_key(item['source'])];target=byid[orig['id']]
                    if q[i][0]!=target['element'] or np.max(np.abs(np.array(q[i][1:])-orig['xyz_A']))>1e-9:
                        raise InvalidArtifact('source atom/core mismatch')
                    if q[i][0]=='H':
                        new[i]=('H',*target['xyz_A'])
                        if np.max(np.abs(np.array(q[i][1:])-new[i][1:]))>1e-9:
                            moves.append({'qm_index':i,'physical_id':orig['id'],'before_A':list(q[i][1:]),'after_A':list(new[i][1:])})
                    elif np.max(np.abs(np.array(q[i][1:])-target['xyz_A']))>1e-9:raise InvalidArtifact('heavy atom changed')
                    if metal=='La':mapping.append({'qm_index':i,'kind':'source','physical_id':orig['id'],'source':item['source']})
                elif item['kind']=='sigma_link_H':
                    if q[i][0]!='H':raise InvalidArtifact('cap element differs')
                    x,y=(np.array(byid[bykey[source_key(item[k])]['id']]['xyz_A']) for k in ('retained','omitted'))
                    target=x+item['length_A']*(y-x)/np.linalg.norm(y-x)
                    if np.max(np.abs(target-q[i][1:]))>1e-9:raise InvalidArtifact('unchanged heavy cap anchors do not reproduce cap')
                    if metal=='La':mapping.append(copy.deepcopy(item))
                else:raise InvalidArtifact('unsupported source mapping')
            if len(seen)!=len(q) or len(moves)!=COUNTS[name] or check_atoms(new,e['charge'])!=e['state']:
                raise InvalidArtifact('normalized core atom/state/H inventory differs')
            if any(q[i]!=new[i] for i in range(len(q)) if q[i][0]!='H'):raise InvalidArtifact('heavy coordinates not exactly preserved')
            endpoints[metal]={'atoms':new,'charge':e['charge'],'spin_multiplicity':1,'state':e['state'],'original_xyz':e['xyz'],'H_moves':moves}
        if endpoints['Ca']['atoms'][1:]!=endpoints['La']['atoms'][1:] or endpoints['Ca']['atoms'][0][1:]!=endpoints['La']['atoms'][0][1:]:
            raise InvalidArtifact('normalized paired coordinates differ')
        cores[name]={'endpoints':endpoints,'mapping':mapping,'source_mapping':p['cases'][name],'global_id':c['global_id'],
            'evidence':c['evidence'],'explicit_waters':c['explicit_waters'],'physical_atoms':p['physical_cases'][c['global_id']]['physical_atoms'],
            'normalized_global_preparation':a['rows'][c['global_id']]['preparation'],'cap_coordinates_unchanged':True}
    return cores,whole,checks,parent


def prepare(prepared,audit,old_report,agreement,output):
    from mace_global_benchmark import snapshot
    start=time.monotonic();cpu=time.process_time();cores,whole,checks,parent=source_data(prepared,audit,old_report)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    names=(*parent['implementation'],'mace_omol_matched_h.py','mace_omol_vacuum.py','mace_omol_vacuum_hybrid.py',
           'ggr_sensitivity.py','ggr_preparation.py','carve_generic.py','coordination_policy.py','pqq_microstates.py',
           'site_mechanics.py','run_orca_task_manifest.py','render_orca_runtime_input.py','result_protocol.py',
           *(p.name for p in Path(__file__).parent.glob('affordable_*.py')))
    pins=snapshot(out,names);stored={}
    for name,c in cores.items():
        d=out/name;d.mkdir();saved={k:copy.deepcopy(v) for k,v in c.items() if k!='endpoints'};saved['endpoints']={}
        for metal,e in c['endpoints'].items():
            xp=d/(metal+'.xyz');write_xyz(xp,e['atoms']);saved['endpoints'][metal]={k:copy.deepcopy(v) for k,v in e.items() if k!='atoms'};saved['endpoints'][metal]['xyz']=record(xp)
        write_new(d/'mapping.json',saved);stored[name]=record(d/'mapping.json')
    m={'protocol_id':PROTOCOL,'agreement':record(agreement),'source_prepared':record(prepared),'reuse_audit':record(audit),
       'original_H_report':record(old_report),'implementation':pins,'cases':stored,'whole_reuses':whole,'geometry_checks':checks,
       'inventory':parent['inventory'],'software':parent['software'],'tolerances':TOL,'new_DFT_calls':8,'new_MACE_core_calls':8,
       'new_MACE_whole_calls':0,'reused_MACE_whole_calls':6,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,**UNAVAILABLE}
    write_new(out/'preparation.json',m);return validate_prepared(out/'preparation.json')


@cached_file_checks
def validate_prepared(preparation):
    p=read_json(preparation);cores,whole,checks,parent=source_data(verify(p['source_prepared']),verify(p['reuse_audit']),verify(p['original_H_report']))
    if p['protocol_id']!=PROTOCOL or p['tolerances']!=TOL or p['whole_reuses']!=whole or set(p['cases'])!=set(CASES):
        raise InvalidArtifact('normalized preparation/model/reuse differs')
    for gid,actual in checks.items():
        saved=p['geometry_checks'][gid]
        if ({k:v for k,v in saved.items() if k!='maximum_replay_error_A'}!={k:v for k,v in actual.items() if k!='maximum_replay_error_A'}
            or not 0<=saved['maximum_replay_error_A']<=1e-12 or not 0<=actual['maximum_replay_error_A']<=1e-12):
            raise InvalidArtifact('normalized arithmetic gate or physical inventory changed')
    for pin in [p['agreement'],*p['implementation'].values()]:verify(pin)
    for name,ref in p['cases'].items():
        c=read_json(verify(ref));wanted=cores[name]
        for key in set(wanted)-{'endpoints'}:
            if c[key]!=wanted[key]:raise InvalidArtifact('source mapping changed')
        for metal,e in c['endpoints'].items():
            w=wanted['endpoints'][metal]
            if {k:v for k,v in e.items() if k!='xyz'}!={k:v for k,v in w.items() if k!='atoms'} or xyz(verify(e['xyz']))!=w['atoms']:
                raise InvalidArtifact('normalized coordinates or electronic state changed')
        paired(verify(c['endpoints']['La']['xyz']),verify(c['endpoints']['Ca']['xyz']),0,-1)
    return {'status':'pass','cases':4,'new_DFT_calls':8,'new_MACE_core_calls':8,'reused_whole_endpoints':6,'preparation':record(preparation)}


def quantum_tasks(preparation,out):
    from ggr_sensitivity import ORCA
    p=read_json(preparation);tasks=[]
    for name,ref in p['cases'].items():
        c=read_json(verify(ref))
        for metal,e in c['endpoints'].items():
            tid=f'{name}_{metal}';d=out/tid;d.mkdir();xp=d/'core.xyz';ip=d/'endpoint.inp'
            shutil.copyfile(verify(e['xyz']),xp);ip.write_text(scientific_input(e['charge']))
            t={'task_id':tid,'case_id':name,'metal':metal,'charge':e['charge'],'multiplicity':1,'input':record(ip),'xyz':record(xp),
               'source_mapping':ref,'output_path':str(d/'endpoint.out'),'engrad_path':str(d/'endpoint.engrad'),'task_type':'analytic_gradient'}
            tasks.append(t)
    # Largest cores start first to keep the four existing 16-rank slots busy.
    tasks.sort(key=lambda t:(-len(xyz(verify(t['xyz']))),t['task_id']))
    return tasks,record(ORCA)


def prepare_quantum(preparation,output):
    from mace_global_benchmark import snapshot
    validate_prepared(preparation);p=read_json(preparation);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    pins=snapshot(out,p['implementation']);tasks,orca=quantum_tasks(preparation,out)
    for t in tasks:t['cache_key']=cache_key({'task':t,'method':METHOD,'protocol':PROTOCOL,'implementation':pins})
    m={'protocol_id':PROTOCOL,'stage':'normalized_vacuum_DFT','preparation':record(preparation),'agreement':p['agreement'],
       'implementation':pins,'tasks':tasks,'orca':orca,'method':METHOD,'energy_scope':'isolated_vacuum_endpoint',
       'execution_policy':{'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},
       'cost_tracking':{'compute_budget':None,'wall_time_limit':None},**UNAVAILABLE}
    write_new(out/'manifest.json',m);return validate_quantum(out/'manifest.json')


@cached_file_checks
def validate_quantum(manifest):
    if read_json(manifest).get('stage')=='coupled_donor_path':
        from mace_site_path_native import validate_quantum as coupled
        return coupled(manifest)
    if read_json(manifest).get('stage')=='matched_hybrid_GGR_transfer':
        from mace_omol_hybrid_transfer import validate_quantum as transfer
        return transfer(manifest)
    if read_json(manifest).get('stage')=='matched_hybrid_response_validation':
        from mace_omol_hybrid_minimum import validate_quantum as validate_response
        return validate_response(manifest)
    from affordable_workflow import dry_run
    m=read_json(manifest);validate_prepared(verify(m['preparation']));p=read_json(verify(m['preparation']))
    if m['protocol_id']!=PROTOCOL or m['stage']!='normalized_vacuum_DFT' or m['method']!=METHOD or len(m['tasks'])!=8:
        raise InvalidArtifact('quantum method/inventory differs')
    for pin in m['implementation'].values():verify(pin)
    if {t['task_id'] for t in m['tasks']}!={f'{c}_{metal}' for c in CASES for metal in ('Ca','La')}:raise InvalidArtifact('quantum tasks differ')
    for t in m['tasks']:
        c=read_json(verify(p['cases'][t['case_id']]));e=c['endpoints'][t['metal']]
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if (t['xyz']['sha256']!=e['xyz']['sha256'] or verify(t['input']).read_text()!=scientific_input(e['charge'])
            or t['charge']!=e['charge'] or t['multiplicity']!=1 or t['source_mapping']!=p['cases'][t['case_id']]
            or t['cache_key']!=cache_key({'task':payload,'method':METHOD,'protocol':PROTOCOL,'implementation':m['implementation']})):
            raise InvalidArtifact('quantum geometry/state/cache differs')
    return dry_run(manifest)


def collect_quantum(manifest):
    if read_json(manifest).get('stage')=='coupled_donor_path':
        from mace_site_path_native import collect_quantum as coupled
        return coupled(manifest)
    if read_json(manifest).get('stage')=='matched_hybrid_GGR_transfer':
        from mace_omol_hybrid_transfer import collect_quantum as transfer
        return transfer(manifest)
    if read_json(manifest).get('stage')=='matched_hybrid_response_validation':
        from mace_omol_hybrid_minimum import collect_quantum as collect_response
        return collect_response(manifest)
    from ggr_sensitivity import executed
    validate_quantum(manifest);m,rows=executed(manifest)
    for t in m['tasks']:
        r=rows[t['task_id']]
        if r['status']!='complete':continue
        try:r.update(parse_endpoint(t,verify(r['output']),t['engrad_path']))
        except (ValueError,OSError) as exc:r.update(status='invalid',reason=str(exc),energy_hartree=None)
    return {'manifest':record(manifest),'rows':rows,'status':'complete' if all(r['status']=='complete' for r in rows.values()) else 'incomplete',**UNAVAILABLE}


def core_tasks(preparation):
    p=read_json(preparation);tasks=[]
    for name,ref in p['cases'].items():
        c=read_json(verify(ref))
        for metal,e in c['endpoints'].items():
            tasks.append(energy_task({'task_id':f'{name}_{metal}','case_id':name,'metal':metal,'metal_index':0,'kind':'core','variant':'primary',
                'xyz':e['xyz'],'charge':e['charge'],'spin_multiplicity':1,'state':e['state'],'source_mapping':ref,
                'explicit_waters':c['explicit_waters'],'evidence':c['evidence'],'evidence_use':'consumed_method_development',
                'energy_component':COMPONENT,'charge_feature_adapter':ADAPTER,'output_semantics':SEMANTICS,'edge_adapter':copy.deepcopy(ADAPTER_TASK)}))
    return tasks


def prepare_mace(preparation,output):
    from mace_omol import common,seal
    validate_prepared(preparation);p=read_json(preparation)
    _,out,m=common(verify(p['inventory']),verify(p['software']),verify(p['agreement']),output,STAGE)
    m.update(protocol_id=PROTOCOL,model=descriptor_model(verify(p['software'])),tasks=core_tasks(preparation),preparation=record(preparation),
             energy_evaluation=EVALUATION,output_semantics=SEMANTICS,hybrid_tolerances=TOL,**UNAVAILABLE)
    return seal(out,m)


@cached_file_checks
def validate_mace(manifest):
    m=read_json(manifest);validate_prepared(verify(m['preparation']));p=read_json(verify(m['preparation']));expected=core_tasks(verify(m['preparation']))
    if (m['protocol_id']!=PROTOCOL or m['stage']!=STAGE or m['reused'] or m['model']!=descriptor_model(verify(p['software']))
        or m['software']!=p['software'] or m['inventory']!=p['inventory'] or m['hybrid_tolerances']!=TOL or len(m['tasks'])!=8):
        raise InvalidArtifact('normalized MACE method/inventory differs')
    for pin in [m['agreement'],*m['implementation'].values()]:verify(pin)
    parent=read_json(verify(read_json(verify(p['reuse_audit']))['collection']));pm=read_json(verify(parent['manifest']))
    for name in ('mace_omol_ablation.py','mace_omol_products.py','mace_omol_edges.py','mace_omol_readout.py'):
        if m['implementation'][name]['sha256']!=pm['implementation'][name]['sha256']:raise InvalidArtifact('descriptor implementation changed')
    for t,w in zip(m['tasks'],expected):
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if payload!=w or t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('normalized MACE state/cache differs')
    return {'status':'pass','tasks':8,'manifest':record(manifest)}


@cached_file_checks
def collect_mace(manifest):
    validate_mace(manifest);mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[]
    for t in m['tasks']:
        good=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp)
            if r is not None:good.append(r)
            attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None,'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
        rows[t['task_id']]=good[-1] if good else {'status':'unavailable','energy_eV':None}
    return {'manifest':record(manifest),'rows':rows,'attempts':attempts,'status':'complete' if all(v['status']=='computed' for v in rows.values()) else 'incomplete',**UNAVAILABLE}


def report(quantum,mace,output):
    q=collect_quantum(quantum);r=collect_mace(mace);qm=read_json(quantum);mm=read_json(mace)
    if qm['preparation']!=mm['preparation']:raise InvalidArtifact('quantum/MACE geometry preparations differ')
    p=read_json(verify(qm['preparation']));old=read_json(verify(p['original_H_report']));complete=q['status']==r['status']=='complete'
    checks=[];cases={};contrasts=[];partition=None
    if complete:
        for tid,v in r['rows'].items():
            error=v['native_readout']['component_sum_error_kcal_mol'];checks.append({'name':tid+'_readout','error':error,'pass':abs(error)<=TOL['accounting_model_kcal']})
        for name,ref in p['cases'].items():
            case=read_json(verify(ref));ends={}
            for metal in ('Ca','La'):
                tid=f'{name}_{metal}';d=q['rows'][tid]['energy_hartree'];c=r['rows'][tid]['energy_eV'];f=p['whole_reuses'][case['global_id']][metal]['selected']['result']['energy_eV']
                ends[metal]={'DFT_vacuum_hartree':d,'learned_core_eV':c,'learned_full_eV':f,'context_model_kcal':(f-c)*EV_TO_KCAL,
                             'hybrid_kcal_scale':d*HA_TO_KCAL+(f-c)*EV_TO_KCAL}
            ca,la=ends['Ca'],ends['La'];dr=(ca['DFT_vacuum_hartree']-la['DFT_vacuum_hartree'])*HA_TO_KCAL
            cr=(ca['learned_core_eV']-la['learned_core_eV'])*EV_TO_KCAL;fr=(ca['learned_full_eV']-la['learned_full_eV'])*EV_TO_KCAL;hr=dr+fr-cr
            error=hr-(ca['hybrid_kcal_scale']-la['hybrid_kcal_scale']);checks.append({'name':name+'_algebra','error':error,'pass':abs(error)<=TOL['algebra_kcal_scale']})
            cases[name]={'endpoints':ends,'DFT_vacuum_R_kcal_mol':dr,'learned_core_R_model_kcal':cr,'learned_full_R_model_kcal':fr,
                         'hybrid_R_kcal_scale':hr,'context_R_model_kcal':fr-cr,'original_H_hybrid_R_kcal_scale':old['cases'][name]['hybrid_R_kcal_scale'],
                         'normalization_change_kcal_scale':hr-old['cases'][name]['hybrid_R_kcal_scale'],'evidence':case['evidence'],'global_id':case['global_id']}
        a,b=cases['GGR_connected'],cases['GGR_extended'];partition=a['hybrid_R_kcal_scale']-b['hybrid_R_kcal_scale']
        error=partition-((a['DFT_vacuum_R_kcal_mol']-b['DFT_vacuum_R_kcal_mol'])-(a['learned_core_R_model_kcal']-b['learned_core_R_model_kcal']))
        checks.append({'name':'partition_full_term_cancellation','error':error,'pass':abs(error)<=TOL['algebra_kcal_scale']})
        for a in ('ALPHA_1F6S','ALPHA_6IP9'):
            for b in ('GGR_extended','GGR_connected'):
                delta=cases[a]['hybrid_R_kcal_scale']-cases[b]['hybrid_R_kcal_scale'];contrasts.append({'alpha':a,'GGR':b,'difference_kcal_scale':delta,'pass':delta>TOL['ordering_kcal_scale']})
    numerical=complete and bool(checks) and all(x['pass'] for x in checks)
    result={'protocol_id':PROTOCOL,'preparation':qm['preparation'],'quantum':q,'MACE':r,'cases':cases,'checks':checks,'contrasts':contrasts,
        'status':'complete' if complete else 'incomplete','numerical_gate_pass':numerical,'partition_shift_kcal_scale':partition,
        'partition_gate_pass':partition is not None and abs(partition)<=TOL['partition_kcal_scale'],
        'ordering_gate_pass':numerical and len(contrasts)==4 and all(x['pass'] for x in contrasts),'tolerances':TOL,
        'new_DFT_calls':8,'new_MACE_core_calls':8,'reused_MACE_whole_calls':6,'report_implementation':record(__file__),**UNAVAILABLE}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);write_new(out/'result.json',result)
    return {k:result[k] for k in ('status','numerical_gate_pass','partition_shift_kcal_scale','partition_gate_pass','contrasts','ordering_gate_pass')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare')
    for k in ('prepared','audit','old-report','agreement','output'):a.add_argument('--'+k,required=True)
    for op in ('prepare-quantum','prepare-mace'):
        a=s.add_parser(op);a.add_argument('--preparation',required=True);a.add_argument('--output',required=True)
    for op in ('dry-run-quantum','execute-quantum','collect-quantum'):
        a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
    a=s.add_parser('report')
    for k in ('quantum','mace','output'):a.add_argument('--'+k,required=True)
    a=p.parse_args()
    if a.op=='prepare':r=prepare(a.prepared,a.audit,a.old_report,a.agreement,a.output)
    elif a.op=='prepare-quantum':r=prepare_quantum(a.preparation,a.output)
    elif a.op=='prepare-mace':r=prepare_mace(a.preparation,a.output)
    elif a.op=='dry-run-quantum':r=validate_quantum(a.manifest)
    elif a.op=='execute-quantum':
        from affordable_workflow import execute
        validate_quantum(a.manifest);r=execute(a.manifest)
    elif a.op=='collect-quantum':r=collect_quantum(a.manifest)
    else:r=report(a.quantum,a.mace,a.output)
    if a.op in ('dry-run-quantum','execute-quantum','collect-quantum') and a.output:write_new(a.output,r)
    print(json.dumps(r if a.op!='collect-quantum' else {'status':r['status'],'rows':len(r['rows'])},indent=2))
