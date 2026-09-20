"""Fixed whole-carboxylate OMOL/solvent profiles with preselected DFT validation."""
from __future__ import annotations
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import re
import resource
import shutil
import time
import numpy as np
from affordable_common import HA_TO_KCAL, InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL, check_atoms, write_xyz
from mace_site_kinematics import Kinematics
from compact_solvation import input_text, completed, diagnostics

PROTOCOL='physical_carboxylate_solvent_profiles_v1'
CASES=('1H4I','4MAE','PQQSEQ_83440678cbbd658047c9','PQQSEQ_07ab500e3df76b30d71c')
ANGLES=(-.4,-.2,0.,.2,.4)
ROLES={'anchor_glutamate':('GLU','chi3'),'extra_acidic_ligand_homolog':('ASP','chi2')}


def role_mode(parent,role):
    r=parent['fixed_core']['requested_roles'][role]
    if isinstance(r,str):
        ch,rn,n,ic=re.fullmatch(r'([^:]+):([A-Z]+)(-?\d+)([A-Za-z]?)',r).groups()
        r={'chain':ch,'resname':rn,'resnum':int(n),'icode':ic}
    rn,chi=ROLES[role]
    if r['resname']!=rn or r.get('icode',''):raise InvalidArtifact('unsupported role chemistry/insertion mode ambiguity')
    return f"{r['chain']}/{r['resnum']}/{chi}"


def geometry_check(kin,q,symbols):
    checks=kin.check(q);coords=kin.evaluate(q)[1];origin=kin.core
    pairs=[]
    for i in range(len(coords)):
        for j in range(i):
            limit=.55 if 'H' in (symbols[i],symbols[j]) else 1.
            old=float(np.linalg.norm(origin[i]-origin[j]));new=float(np.linalg.norm(coords[i]-coords[j]))
            if new<limit and old>=limit:pairs.append({'indices':[j,i],'before_A':old,'after_A':new,'threshold_A':limit})
    return {**checks,'maximum_heavy_displacement_A':kin.displacement(q),'new_severe_overlaps':pairs,
            'pass':checks['pass'] and not pairs}


def prepare(root,agreement,parent_plan,output):
    from second_shell_context import parent_state
    from environment_context_chemistry import complete_expansion
    from coordination_preparation_context import geometry
    from hydration_square import endpoint
    root=Path(root).resolve();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    compact_path=root/'workspaces/compact_solvation_20260920/full_v1/manifest.json';compact=read_json(compact_path)
    static_path=root/'workspaces/second_shell_20260919/prepared_v2/manifest.json';static=read_json(static_path)
    config=read_json(verify(static['configuration']));inv=read_json(verify(compact['inventory']))
    analytic=read_json(root/'workspaces/accommodation_response_20260920/prepared_v2/design.json')
    cases=[];points=[];low=[];dft=[];reused_dft={}
    for case in CASES:
        cd=out/'sources'/case;cd.mkdir(parents=True)
        if case in CASES[:2]:
            item=next(s for s in static['states'] if s['case']==case);source=item['source'];parent=read_json(verify(source['parent']))
            prep=read_json(verify(item['preparation']));context={};orig={}
            for metal in ('Ca','La'):
                task=next(t for t in compact['all_tasks'] if (t['case_id'],t['representation'],t['metal'],t['medium'])==(case,'context',metal,'vacuum'))
                context[metal]=xyz(verify(task['xyz']));orig[metal]=task['source_endpoint']
            prep_pin=item['preparation'];state=parent_state(source,config['topology'])
        else:
            pp=root/'workspaces/plm_xoxf_all_fixed_core_20260916/retry_1200785/prepared_batch/prepared'/(case+'_AF3_sample0_carve_manifest.json')
            parent=read_json(pp);eps={}
            for metal in ('Ca','La'):
                eps[metal]={'xyz':parent['outputs'][metal+'_xyz'],'input':parent['outputs'][metal+'_input'],
                            'charge':parent['charge_ledger'][metal+'_total'],'multiplicity':1}
            source={'case':case,'parent':record(pp),'endpoints':eps}
            state=parent_state(source,config['topology'],require_endpoint_receipts=False)
            context,prep=complete_expansion(state);pp=cd/'context_preparation.json';write_new(pp,prep);prep_pin=record(pp);orig={}
        roles=['anchor_glutamate']+(['extra_acidic_ligand_homolog'] if case!='1H4I' else [])
        modes={r:role_mode(parent,r) for r in roles};maps={};origins={}
        for metal in ('Ca','La'):
            rows=context[metal];charge=source['endpoints'][metal]['charge']+prep['added_formal_charge'];check_atoms(rows,charge)
            gp=geometry(state,prep,rows,state['original'][metal]);mp=cd/(metal+'_mapping.json');write_new(mp,gp);maps[metal]=record(mp)
            kin=Kinematics(gp['context']);names=[m['id'] for m in kin.modes]
            if any(v not in names for v in modes.values()):raise InvalidArtifact('missing physical role mode')
            xp=cd/(metal+'_origin.xyz');write_xyz(xp,rows);origins[metal]={'xyz':record(xp),'charge':charge,'multiplicity':1}
            if case in CASES[:2]:
                old=analytic['sources'][case+'__'+metal]
                if xyz(verify(old['xyz']))!=rows or charge!=old['charge']:raise InvalidArtifact('control origin changed')
            candidates=[('origin',None,0.)]+[(r+'_'+str(a).replace('-','m').replace('.','p'),r,a) for r in roles for a in ANGLES if a]
            for name,role,angle in candidates:
                q=np.zeros(len(kin.modes))
                if role:q[names.index(modes[role])]=angle
                coords=kin.evaluate(q)[1];tid=case+'__'+name+'__'+metal;pd=out/'points'/tid;pd.mkdir(parents=True)
                xp=pd/'core.xyz';write_xyz(xp,[(a[0],*p) for a,p in zip(rows,coords)])
                checks=geometry_check(kin,q,[a[0] for a in rows]);pt={'task_id':tid,'case_id':case,'metal':metal,'point':name,
                    'role':role,'mode_id':modes.get(role),'angle_radian':angle,'q':q.tolist(),'xyz':record(xp),'charge':charge,
                    'multiplicity':1,'mapping':maps[metal],'geometry_checks':checks,'atoms':len(rows)}
                points.append(pt)
                for medium in ('vacuum','alpb'):
                    ld=out/'low'/'tasks'/(tid+'__'+medium);ld.mkdir(parents=True);lx=ld/'core.xyz';shutil.copyfile(xp,lx)
                    ip=ld/'endpoint.inp';ip.write_text(input_text(charge,1,medium,'native'))
                    low.append({'task_id':tid+'__'+medium,'point_id':tid,'case':case,'case_id':case,'metal':metal,
                        'medium':medium,'charge':charge,'multiplicity':1,'xyz':record(lx),'input':record(ip),'output_path':str(ld/'endpoint.out'),
                        'geometry_supported':checks['pass']})
                native=(role=='extra_acidic_ligand_homolog' and abs(angle)==.2) or (case in CASES[2:] and name=='origin')
                if native:
                    dd=out/'dft'/'tasks'/tid;dd.mkdir(parents=True);dx=dd/'core.xyz';shutil.copyfile(xp,dx)
                    ip=dd/'endpoint.inp';ip.write_text('! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3\n'+f'* xyzfile {charge} 1 core.xyz\n')
                    dft.append({'task_id':tid,'case':case,'case_id':case,'metal':metal,'point_id':tid,'charge':charge,'multiplicity':1,
                        'xyz':record(dx),'input':record(ip),'output_path':str(dd/'endpoint.out'),'geometry_supported':checks['pass']})
            if case=='4MAE':
                oldtask=next(t for t in static['tasks'] if t['case']==case and t['metal']==metal)
                receipt=endpoint(record(oldtask['output_path']),record(oldtask['output_path']+'.execution.json'),oldtask['xyz'],oldtask['input'])
                if receipt is None or xyz(verify(oldtask['xyz']))!=rows or oldtask['charge']!=charge:raise InvalidArtifact('exact original 4MAE DFT missing')
                reused_dft[case+'__origin__'+metal]=receipt
        cases.append({'case_id':case,'source':source,'preparation':prep_pin,'modes':modes,'maps':maps,'origins':origins,
            'archived_control_endpoints':orig,'label_scope':'consumed_functional_PQQ_class' if case in CASES[:2] else 'unknown_PLM_prediction'})
    if (len(points),len(low),len(dft))!=(64,128,16):raise InvalidArtifact('fixed manifest denominator changed')
    design={'protocol_id':PROTOCOL,'agreement':record(agreement),'parent_plan':record(parent_plan),'angles_radian':list(ANGLES),
        'cases':cases,'points':points,'reused_DFT_centers':reused_dft,'compact_source':record(compact_path),
        'software':inv['software'],'model':inv['model'],'orca':compact['orca'],'topology':config['topology'],
        'expected_new_calls':{'MACE':64,'GFN2':128,'DFT':16},'affinity_correction':None,'baseline_changed':False}
    write_new(out/'design.json',design);impl=out/'implementation';impl.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        dp=impl/p.name;shutil.copyfile(p,dp);pins[p.name]=record(dp)
    for stage,tasks,ranks,concurrent in [('low',low,8,8),('dft',dft,16,4)]:
        m={'protocol_id':PROTOCOL,'stage':stage,'design':record(out/'design.json'),'agreement':record(agreement),'implementation':pins,
           'orca':compact['orca'],'tasks':[t for t in tasks if t['geometry_supported']],'all_tasks':tasks,
           'execution_policy':{'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},
           'execution_resources':{'mpi_ranks':ranks,'concurrent_tasks':concurrent},'baseline_changed':False}
        write_new(out/stage/'manifest.json',m)
    write_new(out/'manifest.json',{'protocol_id':PROTOCOL,'design':record(out/'design.json'),'implementation':pins,
        'low_manifest':record(out/'low/manifest.json'),'dft_manifest':record(out/'dft/manifest.json')})
    return validate(out/'manifest.json')


def validate(manifest):
    from affordable_workflow import dry_run
    m=read_json(manifest);d=read_json(verify(m['design']));verify(d['agreement']);verify(d['parent_plan'])
    if d['protocol_id']!=PROTOCOL or tuple(d['angles_radian'])!=ANGLES or [c['case_id'] for c in d['cases']]!=list(CASES):raise InvalidArtifact('frozen experiment changed')
    for p in m['implementation'].values():verify(p)
    for pt in d['points']:
        kin=Kinematics(read_json(verify(pt['mapping']))['context']);rows=xyz(verify(pt['xyz']))
        if not np.allclose(kin.evaluate(pt['q'])[1],[a[1:] for a in rows],atol=1e-9,rtol=0):raise InvalidArtifact('physical coordinates changed')
        check_atoms(rows,pt['charge'])
    return {'status':'validated','manifest':record(manifest),'points':len(d['points']),
        'unsupported_geometry':sum(not p['geometry_checks']['pass'] for p in d['points']),
        'low':dry_run(verify(m['low_manifest'])),'dft':dry_run(verify(m['dft_manifest']))}


def mace(manifest,output):
    import torch
    from ase import Atoms
    from mace.calculators import mace_omol
    from mace_omol import input_batch
    if not os.environ.get('SLURM_JOB_ID') or not torch.cuda.is_available():raise InvalidArtifact('allocated GPU required')
    m=read_json(manifest);d=read_json(verify(m['design']));out=Path(output);out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK']));torch.set_default_dtype(torch.float64);torch.cuda.reset_peak_memory_stats()
    start=time.monotonic();calc=mace_omol(model=str(verify(d['model']['checkpoint'])),device='cuda',default_dtype='float64');model=calc.models[0]
    if (type(model).__name__!='ScaleShiftMACE' or dict(model.embedding_specs)!=d['model']['embedding_specs']
        or list(model.heads)!=['omol'] or calc.energy_units_to_eV!=1.):raise InvalidArtifact('checkpoint state model differs')
    for p in model.parameters():p.requires_grad_(False)
    load=time.monotonic()-start;rows=[]
    for pt in d['points']:
        tick=time.monotonic();row={'task_id':pt['task_id'],'xyz':pt['xyz'],'status':'unsupported_geometry','energy_eV':None}
        if pt['geometry_checks']['pass']:
            atoms=xyz(verify(pt['xyz']));a=Atoms([x[0] for x in atoms],positions=[x[1:] for x in atoms],pbc=False);a.info.update(charge=pt['charge'],spin=1)
            batch=calc._atoms_to_batch(a);row['state_check']=input_batch(calc,a,pt['charge'],1,batch,0)
            with torch.no_grad():
                result=model(batch.to_dict(),training=False,compute_force=False,compute_virials=False,compute_stress=False,compute_displacement=False,
                    compute_hessian=False,compute_edge_forces=False,compute_atomic_stresses=False)
            if result['forces'] is not None:raise InvalidArtifact('unexpected force computation')
            e=float(result['energy'].detach().cpu().item())
            if not np.isfinite(e):raise InvalidArtifact('nonfinite native MACE energy')
            row.update(status='complete',energy_eV=e,evaluation_seconds=time.monotonic()-tick)
        rows.append(row);write_new(out/(pt['task_id']+'.json'),row);print(json.dumps({'point':pt['task_id'],'status':row['status']}),flush=True)
    result={'manifest':record(manifest),'design':m['design'],'model':d['model'],'rows':rows,'new_MACE_calls':sum(r['status']=='complete' for r in rows),
        'model_load_seconds':load,'wall_seconds':time.monotonic()-start,'peak_host_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'peak_cuda_allocated_bytes':torch.cuda.max_memory_allocated(),'job_id':os.environ['SLURM_JOB_ID'],
        'allocated_CPUs':os.environ['SLURM_CPUS_PER_TASK'],'device':torch.cuda.get_device_name(),'versions':{k:importlib.metadata.version(k) for k in ('torch','mace-torch','ase')}}
    write_new(out/'result.json',result);return {'complete':result['new_MACE_calls'],'denominator':len(rows)}


def execute(manifest):
    from affordable_workflow import execute as existing
    return existing(manifest)


def collect_native(manifest,output):
    m=read_json(manifest);rows=[]
    for t in m['all_tasks']:
        row={'task_id':t['task_id'],'point_id':t['point_id'],'status':'unavailable','energy_hartree':None}
        if not t['geometry_supported']:row.update(status='unsupported_geometry')
        else:
            pin=completed(manifest,t['task_id'])
            if pin:
                try:
                    audit=diagnostics(pin,t) if m['stage']=='low' else {}
                    if audit and audit['charge_sanity_status']!='pass':raise InvalidArtifact('charge sanity failed')
                    row.update(pin);row.update(audit);row['status']='complete'
                except (InvalidArtifact,KeyError,OSError,ValueError) as exc:row.update(reason=str(exc))
        rows.append(row)
    result={'manifest':record(manifest),'stage':m['stage'],'rows':rows,'complete':sum(r['status']=='complete' for r in rows),'denominator':len(rows)}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='operation',required=True)
    q=s.add_parser('prepare')
    for k in ('root','agreement','parent_plan','output'):q.add_argument('--'+k.replace('_','-'),required=True)
    for op in ('validate','mace','execute','collect_native'):
        q=s.add_parser(op);q.add_argument('--manifest',required=True)
        if op in ('mace','collect_native'):q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('operation');r=globals()[op](**a);print(json.dumps({k:v for k,v in r.items() if k!='rows'},indent=2))


if __name__=='__main__':main()
