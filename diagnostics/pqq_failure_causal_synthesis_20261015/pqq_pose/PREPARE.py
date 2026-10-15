#!/usr/bin/env python3
"""Metadata freeze and Slurm-only finite intact-PQQ rigid-pose preparation."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import socket
import subprocess
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/pqq_pose'
PROTOCOL='ExaA_intact_PQQ_e3_plus_minus_one_degree_common9_preparation_v1'
BASES=('origin','adaptive_Ca','adaptive_La');METALS=('Ca','La');SIGNS=('minus','plus')
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','OPENMM_CPU_THREADS')

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def check(p):
    if pin(p['path'])!=p:raise ValueError('Changed input: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def require(value,reason):
    if not value:raise ValueError(reason)

def freeze(output):
    inventory=pin(ROOT/'workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json')
    require(inventory['sha256']=='757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343','Wrong exact source inventory')
    force=pin(HERE.parent/'pqq_force_inventory/INVENTORY.json');fi=read(force)
    origin=pin(ROOT/'workspaces/native_cavity_shielding_20261014/preparation_v1/execution/PREPARATION.json')
    donor=pin(ROOT/'workspaces/native_cavity_accommodation_20261014/preparation_v1/execution/PREPARATION.json')
    origins=read(origin);donors=read(donor)
    old=read(ROOT/'workspaces/exaa_metal_position_20261014/preparation_v1/manifest.json')
    sources=[]
    for sample in (0,3):
        src,=[r for r in read(inventory)['rows'] if r['case_id'].endswith('sample-'+str(sample))]
        ctx,=[r for r in fi['contexts'] if r['sample']==sample]
        proposals={z:read(src['native_proposals'][z]) for z in METALS}
        require(proposals['Ca']['manifest']==proposals['La']['manifest'],'Proposal manifest differs')
        union={}
        for z in METALS:
            union[z]={}
            for q in BASES:
                candidates=origins['rows'] if q=='origin' else donors['rows']
                row,=[r for r in candidates if r['sample']==sample and r['metal']==z and r.get('donor_candidate','origin')==q]
                require(row['status']=='prepared','Saved source union unavailable')
                union[z][q]=row
        sources.append({'sample':sample,'role':src['role'],'source_case_id':src['case_id'],'record':src,
          'PQQ_context':ctx,'proposal_manifest':proposals['Ca']['manifest'],'source_union':union})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir()
    code={}
    for name in ('PREPARE.py','PREPARE.sbatch'):
        shutil.copyfile(HERE/name,impl/name);code[name]=pin(impl/name)
    m={'protocol':PROTOCOL,'implementation':code,'inventory':inventory,'force_inventory':force,
      'force_result':pin(ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/pqq_forces_v1/execution/RESULT.json'),
      'origin_union_preparation':origin,'donor_union_preparation':donor,'sources':sources,
      'source_snapshot':old['source_snapshot'],'source_modules':old['source_modules'],'python':old['python'],
      'original_worker_settings':old['original_worker_settings'],'plan':pin(HERE.parent/'PQQ_POSE_PLAN.md'),
      'settings':{'angle_degrees':1.0,'frame_tolerance':1e-10,'PQQ_internal_distance_tolerance_A':1e-10,
        'original_kinematic_replay_A':1e-8,'source_bond_tolerance_A':1e-9,'maximum_heavy_displacement_A':0.8,
        'new_overlap_heavy_A':1.0,'new_overlap_H_A':0.55,'serialization':'.17g; non-PQQ coordinate lines unchanged'},
      'base_order':list(BASES),'sign_order':list(SIGNS),'declared_sources':2,'declared_new_cells':24,
      'reused_original_cells':12,'new_energy_calls':0,'new_optimizations':0,'output':str(out/'execution')}
    check(m['force_result']);write(out/'manifest.json',m);print(json.dumps(pin(out/'manifest.json')))

def placement():
    job=os.environ.get('SLURM_JOB_ID');require(job and job.isdigit(),'Actual Slurm required before scientific imports')
    text=subprocess.check_output(['scontrol','show','job','-o',job],text=True)
    fields=dict(x.split('=',1) for x in text.split() if '=' in x)
    require((fields.get('JobState'),fields.get('NumCPUs'),fields.get('NumNodes'))==('RUNNING','1','1'),'One actual CPU/node required')
    require(int(os.environ['SLURM_CPUS_ON_NODE'])==1 and all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread required')
    return {'job_id':job,'hostname':socket.gethostname(),'scontrol':text,'affinity':sorted(os.sched_getaffinity(0)),'CPUs':1}

def overlaps(np,before,after,symbols,ids,environment=None,environment_symbols=None):
    """Exact inherited new<limit and old>=limit predicate; no tuned contacts."""
    failures=[];minimum=None
    for i in ids:
        js=range(i) if environment is None else range(len(environment))
        for j in js:
            a=before[j] if environment is None else environment[j]
            b=after[j] if environment is None else environment[j]
            sym=symbols[j] if environment is None else environment_symbols[j]
            limit=0.55 if 'H' in (symbols[i],sym) else 1.0
            old=float(np.linalg.norm(before[i]-a));new=float(np.linalg.norm(after[i]-b))
            record={'local_index':i,'other_index':j if environment is None else j+192,'before_A':old,'after_A':new,'threshold_A':limit}
            if minimum is None or new<minimum['after_A']:minimum=record
            if new<limit and old>=limit:failures.append(record)
    return {'pass':not failures,'new_severe_overlaps':failures,'nearest_pair':minimum,
      'policy':'new<limit AND saved_base>=limit; inherited heavy1.0A/H0.55A'}

def execute(path,sha):
    place=placement();mp={'path':str(Path(path).resolve()),'sha256':sha};m=read(mp)
    require(m['protocol']==PROTOCOL and m['implementation']['PREPARE.py']==pin(__file__),'Wrong preparation snapshot')
    require(m['settings']['angle_degrees']==1.0 and m['base_order']==list(BASES) and m['sign_order']==list(SIGNS),'Frozen finite pose changed')
    out=Path(m['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    receipt={'status':'failed','manifest':mp,'placement':place,'new_energy_calls':0};rows=[]
    try:
        for key in ('inventory','force_inventory','force_result','origin_union_preparation','donor_union_preparation','python','plan'):check(m[key])
        require(pin(sys.executable)==m['python'],'Wrong numerical Python')
        for p in m['source_modules'].values():check(p)
        sys.path.insert(0,m['source_snapshot'])
        import numpy as np
        import adaptive_metal_proposals as geom
        import mace_site_kinematics as kinematic
        import mace_hybrid as atomic
        import affordable_common as common
        for module in (geom,kinematic,atomic,common):require(pin(module.__file__)==m['source_modules'][Path(module.__file__).name],'Wrong frozen helper')
        inv=read(m['inventory']);fi=read(m['force_inventory'])
        for s in m['sources']:
            src=s['record'];sample=s['sample'];ctx=s['PQQ_context'];directory=out/f'sample{sample}';directory.mkdir()
            row={'sample':sample,'role':s['role'],'case_id':s['source_case_id'],'source_case_id':s['source_case_id'],
              'status':'unavailable','reason':None,'matrix':{z:{f'{q}__PQQ_{sign}1deg':{'status':'unavailable'} for q in BASES for sign in SIGNS} for z in METALS},
              **{k:src[k] for k in ('selection_id','maps','state_signature','source_preparation','prepared_admission')}}
            rows.append(row);admission={'status':'unavailable','pass':False,'sample':sample,'variants':{},'declared_new_cells':12}
            try:
                actual,=[r for r in inv['rows'] if r['case_id']==s['source_case_id']];actualctx,=[r for r in fi['contexts'] if r['sample']==sample]
                require(src==actual and ctx==actualctx,'Changed source/context record')
                cp=read(src['source_preparation']);state=read(src['state_signature']);a=read(src['prepared_admission'])
                maps={z:read(src['maps'][z]) for z in METALS}
                require(maps['Ca']==maps['La'] and state==src['state'],'Source maps/state differs')
                require(a['status']=='prepared' and a['mapping_eligible'] and a['paired_maps_equal'] and (a['atom_count'],a['cap_count'])==(192,14),'Source admission differs')
                require(state['endpoint_charges']=={'Ca':-3,'La':-2} and state['multiplicities']=={'Ca':1,'La':1} and not state['water_inventory'],'Dry state changed')
                named={x['name']:x['qm_index'] for x in cp['mapping']['source_to_qm'] if x.get('fragment')==ctx['PQQ_fragment']}
                require(named==ctx['PQQ_indices_zero_based'] and len(named)==27 and ctx['microstate_id']=='pqq_ox_3minus_v1','Intact namedPQQ mapping changed')
                ids=list(named.values());identity={x['qm_index']:x for x in ctx['PQQ_atom_identity']};heavy=[i for i in ids if identity[i]['element']!='H']
                require(len(heavy)==24 and {identity[i]['name'] for i in ids if identity[i]['element']=='H'}=={'H_N1','H_C3','H_C8'},'PQQ physical H inventory changed')
                capids=[x['qm_index'] for x in cp['mapping']['source_to_qm'] if x['kind']=='sigma_link_H'];require(not set(ids)&set(capids),'PQQ includes synthetic cap')
                kin=kinematic.Kinematics(maps['Ca']['context']);pm=read(s['proposal_manifest']);tasks={}
                for z in METALS:
                    pr=read(src['native_proposals'][z])
                    task,=[t for t in pm['tasks'] if t['task_id']==pr['task_id'] and t['case_id']==pr['case_id'] and t['metal']==z]
                    require(task['mapping']==src['maps'][z] and task['source_preparation']==src['source_preparation'],'Original donor task provenance differs')
                    tasks[z]=task
                links={link[0]:link for link in kin.links};physical={i:links[i][2] for i in heavy if i in links and links[i][1]=='source'}
                require(len(physical)==24,'PQQ physical heavy mapping incomplete')
                physids=set(physical.values())
                require(not any((i in physids)!=(j in physids) for i,j in kin.data['bonds']),'Unsupported covalent PQQ-to-exterior boundary')
                require(not any(link[1]=='cap' and (link[2] in physids or link[3] in physids) for link in kin.links),'PQQ participates in cap anchor')
                origin_atoms=common.xyz(check(src['matrix']['Ca']['origin']['xyz']));origin=np.asarray([x[1:] for x in origin_atoms])
                center=np.mean(origin[heavy],axis=0);e1=origin[named['O1']]-origin[named['C5']];e1/=np.linalg.norm(e1)
                e2=origin[named['C6']]-origin[named['C5']];e2-=e1*np.dot(e1,e2);e2/=np.linalg.norm(e2);e3=np.cross(e1,e2);frame=np.stack((e1,e2,e3))
                require(np.isfinite(frame).all() and np.max(np.abs(frame@frame.T-np.eye(3)))<1e-10 and abs(np.linalg.det(frame)-1)<1e-10,'Invalid right-handed chemical frame')
                row.update(PQQ_indices_zero_based=named,PQQ_center_A=center.tolist(),PQQ_frame=frame.tolist())
                admission.update(PQQ_indices_zero_based=named,PQQ_center_A=center.tolist(),PQQ_frame=frame.tolist(),
                  PQQ_hydrogen_parents=ctx['PQQ_hydrogen_parents'],cap_indices=capids,source_unions={})
                for base in BASES:
                    owner='Ca' if base=='origin' else base.rsplit('_',1)[1];task=copy.deepcopy(tasks[owner])
                    if base=='origin':q=np.zeros(len(kin.modes))
                    else:
                        proposal=read(src['native_proposals'][owner]);require(proposal['status']=='proposal_available' and proposal['proposal']['geometry_checks']['pass'],'Saved proposal unavailable')
                        require(proposal['manifest']==s['proposal_manifest'],'Wrong proposal source');q=np.asarray(proposal['proposal']['full_q'])
                    metal_ids=[next(i for i,mode in enumerate(kin.modes) if mode['id']=='metal_'+axis) for axis in 'xyz']
                    task['active_indices']=metal_ids+task['active_indices']
                    require(len(task['active_indices'])==7 and len(set(task['active_indices']))==7 and np.array_equal(q[metal_ids],np.zeros(3)),'Original four donor modes or fixed metal differs')
                    require(np.array_equal(geom.full_q(task,q[task['active_indices']]),q),'Saved q outside active donor subspace')
                    pair={z:common.xyz(check(src['matrix'][z][base]['xyz'])) for z in METALS}
                    require(pair['Ca'][1:]==pair['La'][1:] and pair['Ca'][0][1:]==pair['La'][0][1:],'Saved paired geometry differs')
                    symbols=[x[0] for x in pair['Ca']];before=np.asarray([x[1:] for x in pair['Ca']]);physical_before,rendered,*_=kin.evaluate(q)
                    require(np.max(np.abs(rendered-before))<=1e-8 and np.array_equal(before[ids],origin[ids]),'Saved donor/PQQ geometry does not replay')
                    inherited=geom.final_geometry(kin,task,q[task['active_indices']],symbols)
                    require(inherited['physical_feasible'] and inherited['mapping_checks']['pass'] and inherited['all_other_physical_atoms_fixed'],'Original donor geometry is not admitted')
                    unionrows={z:s['source_union'][z][base] for z in METALS};unions={}
                    for z,ur in unionrows.items():
                        require(ur['status']=='prepared' and ur['state_signature']==src['state_signature'] and ur['source_preparation']==src['source_preparation'] and ur['archived_native_result']==src['matrix'][z][base]['native_frozen_ALPB_result'],'Wrong saved union source/state/cell')
                        check(ur['variants']['union']['mapping']);unions[z]=common.xyz(check(ur['variants']['union']['xyz']))
                        local_ids=ur['variants']['union']['original_local_indices']
                        require(unions[z][:192]==pair[z] and len(unions[z])==9474 and local_ids[:192]==list(range(192)) and all(i is None for i in local_ids[192:]),'Union first192/source identity differs')
                    require(unions['Ca'][1:]==unions['La'][1:] and unions['Ca'][0][1:]==unions['La'][0][1:],'Paired source union differs')
                    env=np.asarray([x[1:] for x in unions['Ca'][192:]]);env_symbols=[x[0] for x in unions['Ca'][192:]]
                    admission['source_unions'][base]={z:ur['variants']['union'] for z,ur in unionrows.items()}
                    for sign in SIGNS:
                        name=f'{base}__PQQ_{sign}1deg';angle=(-1 if sign=='minus' else 1)*np.pi/180
                        new=before.copy();v=before[ids]-center
                        new[ids]=center+v*np.cos(angle)+np.cross(e3,v)*np.sin(angle)+(v@e3)[:,None]*e3*(1-np.cos(angle))
                        delta=np.linalg.norm(new[ids,None]-new[None,ids],axis=2)-np.linalg.norm(before[ids,None]-before[None,ids],axis=2)
                        maxerr=float(np.max(np.abs(delta)));require(maxerr<=1e-10,'PQQ internal distances changed')
                        hg={h:{'parent':p,'residual_A':float(np.linalg.norm(new[named[h]]-new[named[p]])-np.linalg.norm(before[named[h]]-before[named[p]]))} for h,p in ctx['PQQ_hydrogen_parents'].items()}
                        require(all(abs(x['residual_A'])<=1e-10 for x in hg.values()),'Existing PQQ H-parent bond changed')
                        fixed=[i for i in range(192) if i not in ids];require(np.array_equal(new[fixed],before[fixed]),'Non-PQQ atom moved')
                        pnew=physical_before.copy()
                        for local,physical_i in physical.items():pnew[physical_i]=new[local]
                        berr=max(abs(float(np.linalg.norm(pnew[i]-pnew[j])-np.linalg.norm(physical_before[i]-physical_before[j]))) for i,j in kin.data['bonds'])
                        require(berr<=1e-9,'Physical source bond changed')
                        extent=float(np.max(np.linalg.norm(pnew[kin.heavy]-kin.positions[kin.heavy],axis=1)))
                        require(extent<=0.8+m['original_worker_settings']['final_physical_tolerance_A'],'Original physical bound exceeded')
                        local_contacts=overlaps(np,before,new,symbols,range(192));external_contacts=overlaps(np,before,new,symbols,ids,env,env_symbols)
                        gate={'pass':local_contacts['pass'] and external_contacts['pass'],'original_base':base,'angle_degrees':-1 if sign=='minus' else 1,
                          'inherited_original_donor_admission':inherited,'PQQ_internal_max_distance_residual_A':maxerr,
                          'PQQ_H_parent_checks':hg,'physical_source_bond_max_residual_A':berr,'maximum_total_physical_heavy_displacement_A':extent,
                          'maximum_PQQ_heavy_displacement_A':float(np.max(np.linalg.norm(new[heavy]-before[heavy],axis=1))),
                          'non_PQQ_fixed':True,'caps_and_anchors_fixed':True,'local_contact_check':local_contacts,
                          'full_source_contact_check':external_contacts,'full_source_energy_added':False,'PQQ_donor_kinematic_replay_claimed':False}
                        admission['variants'][name]=gate
                        if not gate['pass']:
                            for z in METALS:row['matrix'][z][name].update(reason='new severe contact',geometry_checks=gate)
                            continue
                        for z in METALS:
                            cell=src['matrix'][z][base];xp=directory/f'{z}_{name}.xyz';lines=check(cell['xyz']).read_text().splitlines(keepends=True);saved=list(lines)
                            for i in ids:lines[i+2]=symbols[i]+' '+' '.join(format(float(x),'.17g') for x in new[i])+'\n'
                            with xp.open('x') as f:f.write(''.join(lines))
                            coords=common.xyz(xp);require(all(lines[i+2]==saved[i+2] for i in fixed),'Excluded atom line changed')
                            require(np.array_equal(np.asarray([x[1:] for x in coords]),new),'Pose serialization changed')
                            row['matrix'][z][name]={'status':'prepared','xyz':pin(xp),'charge':cell['charge'],'multiplicity':1,'atom_count':192,
                              'valence_electrons':cell['GFN2_valence_electrons'],'base_candidate':base,'sign':sign,'original_xyz':cell['xyz'],
                              'original_MACE':cell['MACE'],'original_MACE_eV':cell['MACE_eV'],'original_native_ALPB':cell['native_frozen_ALPB_result'],
                              'atomic_check':atomic.check_atoms(coords,cell['charge']),'geometry_checks':gate}
                        ca=common.xyz(check(row['matrix']['Ca'][name]['xyz']));la=common.xyz(check(row['matrix']['La'][name]['xyz']))
                        require(ca[1:]==la[1:] and ca[0][1:]==la[0][1:],'New paired coordinates differ')
                complete=all(v['status']=='prepared' for cells in row['matrix'].values() for v in cells.values())
                row.update(status='prepared' if complete else 'unavailable',reason=None if complete else 'one or more required variants inadmissible')
                admission.update(status=row['status'],**{'pass':complete})
            except Exception as e:
                row.update(status='unavailable',reason=f'{type(e).__name__}: {e}');admission.update(reason=row['reason'],traceback=traceback.format_exc())
            write(directory/'GEOMETRY_ADMISSION.json',admission);row['geometry_admission']=pin(directory/'GEOMETRY_ADMISSION.json')
            for cells in row['matrix'].values():
                for v in cells.values():v['geometry_admission']=row['geometry_admission']
        receipt['status']='complete'
    except BaseException as e:receipt.update(error=str(e),traceback=traceback.format_exc());raise
    finally:
        receipt.update(elapsed_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(out/'EXECUTION.json',receipt)
        result={'protocol':PROTOCOL,'manifest':mp,'terminal':True,'status':'prepared' if len(rows)==2 and all(r['status']=='prepared' for r in rows) else 'unavailable',
          'inventory':m['inventory'],'rows':rows,'declared_sources':2,'declared_new_cells':24,'reused_original_cells':12,
          'prepared_cells':sum(v['status']=='prepared' for r in rows for cells in r['matrix'].values() for v in cells.values()),'execution':pin(out/'EXECUTION.json'),
          'new_energy_calls':0,'new_optimizations':0,'production_changed':False}
        write(out/'PREPARATION.json',result)
    print(json.dumps({'preparation':pin(out/'PREPARATION.json'),'prepared_cells':result['prepared_cells']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('freeze');a.add_argument('--output',required=True)
    a=s.add_parser('execute');a.add_argument('--manifest',required=True);a.add_argument('--manifest-sha',required=True)
    a=p.parse_args()
    if a.action=='freeze':freeze(a.output)
    else:execute(a.manifest,a.manifest_sha)
