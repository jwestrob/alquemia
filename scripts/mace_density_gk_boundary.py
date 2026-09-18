"""Versioned AMOEBA boundary/source preparation; no energy or response solve."""
from __future__ import annotations
import argparse
import copy
import json
import math
from pathlib import Path
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new
from mace_density_multipoles import parse_moments
from mace_tinker_framework_solver import NONPOLAR,native_call,validate as validate_frameworks

PROTOCOL='density_direct_common_Ca2018_GK_source_boundary_v1'
CASES=('GGR_extended','GGR_connected','ALPHA_1F6S','ALPHA_6IP9')
METALS=dict(Ca=dict(type=4998,Z=20,mass=40.078),La=dict(type=4999,Z=57,mass=138.90547))
TOL=dict(environment_charge_e=1e-9,residue_formal_e=1e-8,QM_charge_e=1e-5,native_moments=1e-10)
CHARGES_SHA='2689bcdfa9b3f2d0358c5e21fa8b5b38533b36502c7e5c3f9826ee7487041d10'
TRIAL_PROTOCOL='responsive_trial_common_Ca2018_GK_source_boundary_v1'
CHARGE_SOURCES={PROTOCOL:CHARGES_SHA,
    TRIAL_PROTOCOL:'fe1ab39ef2ef2ef45f43b49c57a5bce23561fa25ee6c060dfbca43b27df72145'}


def residue(pid):return '/'.join(pid.split('/')[:3])


def boundary(state,mapping,parent_params):
    support=set(state['projection_support_ids']);ids=mapping['system_atom_ids']
    if state['physical_atoms']!=mapping['physical_atoms'] or 'metal' not in support:
        raise InvalidArtifact('physical source framework differs')
    old=read_json(verify(mapping['source_mapping']));bonds=old['protein_bonds']+old['retained_water_bonds']
    neighbors={i:set() for i in ids}
    for a,b in bonds:neighbors[a].add(b);neighbors[b].add(a)
    q={pid:p['charge_e'] for pid,p in zip(ids,parent_params['atoms'])};env=dict(q)
    groups={}
    for pid in ids:groups.setdefault(residue(pid),[]).append(pid)
    source_formal=state['ligand_ledger'];legacy={l['residue']:l for l in state['boundary_ledger']}
    ledger=[];formals={}
    for resid,atoms in groups.items():
        total=math.fsum(q[i] for i in atoms);formal=round(total)
        if abs(total-formal)>TOL['residue_formal_e']:raise InvalidArtifact('noninteger AMOEBA residue charge: '+resid)
        formals[resid]=formal
        selected=sorted(set(atoms)&support)
        if not selected:continue
        fragment=source_formal.get(resid,0)
        if int(fragment)!=fragment:raise InvalidArtifact('noninteger source formal charge')
        if resid in legacy and (formal!=legacy[resid]['forcefield_formal_charge_e'] or fragment!=legacy[resid]['selected_fragment_formal_charge_e']):
            raise InvalidArtifact('chemical formal-charge ledger differs')
        outside=sorted(set(atoms)-support);target=formal-fragment
        recipients=sorted({b for a in selected for b in neighbors[a] if b in outside})
        recipient_bonds=sorted([sorted([a,b]) for a in selected for b in neighbors[a] if b in recipients])
        if resid in legacy and (recipients!=legacy[resid]['recipient_ids'] or selected!=legacy[resid]['projection_support_ids']):
            raise InvalidArtifact('source/recipient graph differs from pinned chemical boundary')
        original=math.fsum(q[i] for i in outside);delta=target-original
        if not outside and target!=0:raise InvalidArtifact('fully selected residue has nonzero exterior target')
        if not recipients and abs(delta)>TOL['environment_charge_e']:
            raise InvalidArtifact('charge residual without bonded exterior recipients')
        for i in selected:env[i]=0.
        if recipients:
            for i in recipients:env[i]+=delta/len(recipients)
        actual=math.fsum(env[i] for i in atoms)
        if abs(actual-target)>TOL['environment_charge_e']:raise InvalidArtifact('residue exterior charge closure failed')
        ledger.append(dict(residue=resid,original_formal_charge_e=formal,selected_fragment_formal_charge_e=fragment,
            source_support_ids=selected,exterior_ids=outside,exterior_target_e=target,original_exterior_charge_e=original,
            delta_e=delta,recipient_ids=recipients,recipient_bonds=recipient_bonds,
            original_recipient_charges_e=[q[i] for i in recipients],new_recipient_charges_e=[env[i] for i in recipients],
            actual_exterior_charge_e=actual))
    if not set(source_formal)<=groups.keys():raise InvalidArtifact('source formal residue absent from physical system')
    total=math.fsum(env.values());expected=math.fsum(formals.values())-math.fsum(source_formal.values())
    if expected!=state['environment_charge'] or abs(total-expected)>TOL['environment_charge_e']:
        raise InvalidArtifact('whole exterior formal-charge closure failed')
    return dict(protocol=PROTOCOL,physical_ids=ids+['metal'],environment_charge_by_id={**env,'metal':0.},
        original_charge_by_id={**q,'metal':0.},ledger=ledger,source_support_ids=sorted(support),
        environment_formal_charge_e=int(expected),actual_environment_charge_e=total,
        complete_residue_formal_charges_e=formals,full_physical_bonds=bonds,ligand_ledger_source=state['ligand_ledger_source'])


def source_key(task,m):
    return cache_key(dict(task={k:v for k,v in task.items() if k!='cache_key'},protocol=m['protocol'],
        sources=m['sources'],implementation=m['implementation'],tolerances=TOL,software=m['software'],plan=m['plan']))


def additional_key():
    lines=[]
    for metal,v in METALS.items():
        typ=v['type']
        lines.extend([f'atom {typ} 99 {metal}Q "Frozen QM {metal}, common Ca2018 cavity" {v["Z"]} {v["mass"]} 0',
            f'solute {typ} 3.6670 3.6670 3.6497 0.1350',f'polarize {typ} 0.0 0.0',
            f'multipole {typ} 0 0 0.0','0.0 0.0 0.0','0.0','0.0 0.0','0.0 0.0 0.0'])
    return lines


def verify_native(text,task,meta,base_params):
    result=parse_moments(text,None);params=result['parameters'];ids=meta['physical_ids'];n=len(ids)
    if params['inventory']!=[n,n,n-1,len(task['frozen_indices']),1]:raise InvalidArtifact('native source inventory changed')
    if params['global']!=[1.,78.3,2.455,.3,.09,1e-7,100.,1e12] or params['GK_flags']!=['GRYCUK','T','T'] or params['nonpolar']!=list(NONPOLAR):
        raise InvalidArtifact('native solver/cavity settings changed')
    original={};bonds={};indices={};identities={}
    for line in text.splitlines():
        f=line.split()
        if not f:continue
        if f[0]=='ORIGINAL_POLE':original[int(f[1])]=list(map(float,f[2:]))
        elif f[0]=='BONDS':bonds[int(f[1])]=list(map(int,f[2:]))
        elif f[0]=='POLE_INDEX':indices[int(f[1])]=list(map(int,f[2:]))
        elif f[0]=='IDENTITY':identities[int(f[1])]=[int(v) for v in f[2:5]]+[float(f[5])]
    expected_ids=set(range(1,n+1))
    if any(set(v)!=expected_ids for v in (original,bonds,indices,identities)):
        raise InvalidArtifact('missing native bond/index/original source rows')
    by_id={pid:i+1 for i,pid in enumerate(ids)};actual_bonds=set()
    frozen=set(task['frozen_indices']);expected_coords=meta['positions_A'];q=task['charge_e'];metal=METALS[task['metal']]
    for i,(atom,moment) in enumerate(zip(params['atoms'],result['atoms']),1):
        if atom['xyz_A']!=expected_coords[i-1] or atom['charge_e']!=q[i-1] or atom['response_allowed']!=(i not in frozen):
            raise InvalidArtifact('native source coordinates/charge/mask mismatch')
        if indices[i]!=[i,i,13]:raise InvalidArtifact('source removed from native multipole map')
        if bonds[i][0]!=len(bonds[i])-1:raise InvalidArtifact('native bond count mismatch')
        actual_bonds.update(tuple(sorted((i,j))) for j in bonds[i][1:])
        if i<n:
            base=copy.deepcopy(base_params['atoms'][i-1]);base['charge_e']=q[i-1];base['response_allowed']=i not in frozen
            if atom!=base:raise InvalidArtifact('nonmetal native cavity/response parameter changed')
        else:
            if identities[i]!=[metal['type'],99,metal['Z'],metal['mass']]:raise InvalidArtifact('source element/mass/type mismatch')
            for k,v in dict(atomic_number=metal['Z'],radius_A=1.82485,descreen_radius_A=1.795,scale=.72,neck=1.,
                           polarizability_A3=0.,damping_A_half=0.).items():
                if atom[k]!=v:raise InvalidArtifact('source common cavity parameter mismatch: '+k)
            if original[i]!=[0.]*13:raise InvalidArtifact('invented classical source moments')
        if i in frozen:
            if any(v!=0 for v in moment['global_'][1:]):raise InvalidArtifact('FF source moments retained')
        elif moment['raw_local'][1:]!=original[i][1:]:raise InvalidArtifact('exterior multipoles changed')
    expected_bonds={tuple(sorted((by_id[a],by_id[b]))) for a,b in meta['full_physical_bonds']}
    if actual_bonds!=expected_bonds:raise InvalidArtifact('physical covalent graph changed')
    result.update(status='native_density_boundary_prepared',native_identities=identities,
        native_original_poles=original,native_multipole_map=indices,covalent_graph_exact=True)
    return result


def prepare(frameworks,charges,states,software,plan,output,protocol=PROTOCOL):
    start=time.monotonic();cpu=time.process_time();fm=validate_frameworks(frameworks);cr=read_json(charges);sw=read_json(software)
    if protocol not in CHARGE_SOURCES or record(charges)['sha256']!=CHARGE_SOURCES[protocol] or cr['status']!='complete' or not cr['projection_gate_pass']:
        raise InvalidArtifact('declared projected charge source required')
    if sw['returncode'] or sw['library']!=read_json(verify(fm['software']))['library']:
        raise InvalidArtifact('unavailable/incompatible native source adapter')
    cm=read_json(verify(cr['manifest']));root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    impl=root/'implementation';impl.mkdir()
    # Reuse the existing density diagnostic's actual dependency snapshot, then
    # preserve this adapter and its current shared parsers without altering it.
    deps=Path(__file__).parent
    required=('affordable_common.py','affordable_solver.py','density_embedding.py','mace_hybrid.py',
        'mace_qm_field.py','mace_omol_charges.py','mace_tinker_framework_solver.py',
        'mace_tinker_capability.py','mace_amoeba_capability.py','mace_native_field_accounting.py',
        'mace_density_multipoles.py','mace_density_gk_boundary.py')
    # Charge utility snapshots contain transitive scientific parser dependencies.
    for name,pin in cm['implementation'].items():shutil.copyfile(verify(pin),impl/name)
    for name in required:shutil.copyfile(deps/name,impl/name)
    m=dict(protocol=protocol,sources=dict(frameworks=record(frameworks),charges=record(charges),charge_manifest=cr['manifest']),
        software=record(software),plan=record(plan),tolerances=TOL,implementation={p.name:record(p) for p in sorted(impl.glob('*.py'))},
        tasks=[],cases={},new_energy_calls=0,new_response_solves=0,new_DFT_calls=0,new_MACE_calls=0,numerical_score=None)
    for case in CASES:
        gid='GGR_1GLG' if case.startswith('GGR_') else case
        old=next(t for t in fm['tasks'] if t['case_id']==gid and t['variant']=='all_standard')
        mapping=read_json(verify(old['mapping']));base=read_json(verify(old['parameters']));sp=Path(states)/case/'state.json';state=read_json(sp)
        meta=boundary(state,mapping,base);physical={a['id']:a for a in state['physical_atoms']};ids=meta['physical_ids']
        meta.update(state=record(sp),base_parameters=old['parameters'],base_mapping=old['mapping'],
            positions_A=[physical[pid]['xyz_A'] for pid in ids],assembly=state['assembly'],microstate=state['microstate'],
            explicit_waters=state['explicit_waters'],evidence=state['evidence'],evidence_use=state['evidence_use'])
        directory=root/'cases'/case;directory.mkdir(parents=True);write_new(directory/'boundary.json',meta)
        m['cases'][case]=record(directory/'boundary.json');template=verify(old['xyz']).read_text().splitlines()
        for metal in ('Ca','La'):
            ctask=next(t for t in cm['tasks'] if t['task_id']==case+'_'+metal);row=cr['rows'][ctask['task_id']]
            proj=read_json(verify(ctask['projection']))
            if set(proj['physical_ids'])!=set(meta['source_support_ids']):raise InvalidArtifact('QM support differs')
            for pid,pos in zip(proj['physical_ids'],proj['coordinates_A']):
                if pos!=physical[pid]['xyz_A']:raise InvalidArtifact('projected source geometry changed')
            actual_Q=math.fsum(row['projected_charge_e']);formal=-1 if metal=='Ca' else 0
            if abs(actual_Q-formal)>TOL['QM_charge_e']:raise InvalidArtifact('QM projected charge sum unsupported')
            srcq=dict(zip(proj['physical_ids'],row['projected_charge_e']))
            for mode in ('source','zero_source'):
                tid=case+'_'+metal+'_'+mode;d=root/'tasks'/tid;d.mkdir(parents=True)
                q=[meta['environment_charge_by_id'][pid]+(srcq.get(pid,0.) if mode=='source' else 0.) for pid in ids]
                n=len(ids);typ=METALS[metal]['type'];pos=physical['metal']['xyz_A']
                lines=[f'{n} {tid}; real physical framework with declared common source cavity',*template[1:]]
                lines.append(f'{n} {metal} '+' '.join(format(v,'.17g') for v in pos)+f' {typ}')
                (d/'framework.xyz').write_text('\n'.join(lines)+'\n')
                opts=['parameters '+str(verify(fm['native_parameters'])),'multipoleterm ONLY','polarizeterm','solvateterm',
                    'solvate GK','dielectric 1.0','gk-radius SOLUTE','gkc 2.455','polarization MUTUAL','polar-iter 100','polar-eps 1e-7']
                (d/'framework.key').write_text('\n'.join(opts+additional_key())+'\n')
                frozen=[i+1 for i,pid in enumerate(ids) if pid in meta['source_support_ids']]
                (d/'framework.freeze').write_text(str(len(frozen))+'\n'+'\n'.join(map(str,frozen))+'\n')
                (d/'framework.charges').write_text(str(n)+'\n'+''.join(f'{i} {v:.17g}\n' for i,v in enumerate(q,1)))
                task=dict(task_id=tid,case_id=case,metal=metal,mode=mode,boundary=m['cases'][case],
                    source_projection=ctask['projection'],source_charge_execution=row['execution_receipt'],
                    QM_formal_charge=formal,QM_actual_charge=actual_Q,QM_charge_residual=actual_Q-formal,
                    source_QM_xyz=ctask['xyz'],source_quantum_receipt=ctask['source_receipt'],
                    charge_e=q,frozen_indices=frozen,total_actual_charge=math.fsum(q),
                    xyz=record(d/'framework.xyz'),key=record(d/'framework.key'),mask=record(d/'framework.freeze'),overrides=record(d/'framework.charges'))
                task['cache_key']=source_key(task,m)
                receipt=native_call(sw['executable'],task,str(verify(task['overrides'])),d/'native.log',1)
                result=dict(status='failed',task_id=tid,cache_key=task['cache_key'],receipt=receipt)
                try:
                    if receipt['returncode']:raise InvalidArtifact('native boundary preparation failed')
                    result.update(verify_native(verify(receipt['log']).read_text(),task,meta,base))
                except Exception as e:result['failure_reason']=f'{type(e).__name__}: {e}'
                write_new(d/'result.json',result);m['tasks'].append(task)
                if result['status']=='failed':
                    write_new(root/'preparation_failure.json',dict(task_id=tid,reason=result['failure_reason'],result=record(d/'result.json'),
                        native_initializations_attempted=len(m['tasks']),wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu))
                    raise InvalidArtifact(result['failure_reason'])
    m['preparation_receipt']=dict(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu)
    write_new(root/'manifest.json',m);return collect(root/'manifest.json',root/'result.json')


def validate(path):
    m=read_json(path)
    if m['protocol'] not in CHARGE_SOURCES or m['tolerances']!=TOL or len(m['tasks'])!=16:raise InvalidArtifact('undeclared native boundary configuration')
    if m['sources']['charges']['sha256']!=CHARGE_SOURCES[m['protocol']]:raise InvalidArtifact('boundary density protocol/source mismatch')
    for pin in [m['software'],m['plan'],*m['sources'].values(),*m['implementation'].values(),*m['cases'].values()]:verify(pin)
    sw=read_json(verify(m['software']))
    for k in ('executable','library','source','log','parent_software'):verify(sw[k])
    for t in m['tasks']:
        for k in ('boundary','source_projection','source_charge_execution','source_QM_xyz','source_quantum_receipt','xyz','key','mask','overrides'):verify(t[k])
        if t['cache_key']!=source_key(t,m):raise InvalidArtifact('changed source state cache key')
    return m


def comparable(r):
    p=copy.deepcopy(r['parameters']);p['atoms'][-1]['atomic_number']=0
    return dict(parameters=p,moments=[a['global_'] for a in r['atoms']],axes=[(a['axis'],a['axis_indices']) for a in r['atoms']])


def collect(path,output):
    m=validate(path);root=Path(path).resolve().parent;rows={};checks=[];results={}
    for t in m['tasks']:
        p=root/'tasks'/t['task_id']/'result.json';r=read_json(p);verify(r['receipt']['log'])
        if r['cache_key']!=t['cache_key'] or r['status']!='native_density_boundary_prepared':raise InvalidArtifact('missing/unqualified native preparation')
        rows[t['task_id']]=dict(result=record(p),status=r['status'],total_actual_charge=t['total_actual_charge'],QM_charge_residual=t['QM_charge_residual'])
        results[t['task_id']]=r
    for case in CASES:
        ca=results[case+'_Ca_zero_source'];la=results[case+'_La_zero_source']
        checks.append(dict(name=case+'_environment_reference_identity',pass_=comparable(ca)==comparable(la)))
        # All source-charge states retain the same physical cavity and response.
        ref=copy.deepcopy(ca['parameters']);ref['atoms'][-1]['atomic_number']=0
        for a in ref['atoms']:a['charge_e']=0.
        for metal in ('Ca','La'):
            p=copy.deepcopy(results[case+'_'+metal+'_source']['parameters']);p['atoms'][-1]['atomic_number']=0
            for a in p['atoms']:a['charge_e']=0.
            checks.append(dict(name=case+'_'+metal+'_cavity_response_identity',pass_=p==ref))
    result=dict(protocol=m['protocol'],manifest=record(path),cases=m['cases'],tasks=rows,checks=checks,complete=True,
        gates_pass=all(c['pass_'] for c in checks),actual_native_initializations=16,new_energy_calls=0,new_response_solves=0,
        native_wall_seconds=sum(r['receipt']['wall_seconds'] for r in results.values()),
        native_CPU_seconds=sum(r['receipt']['child_CPU_seconds'] for r in results.values()),
        numerical_score=None,full_model_qualified=False,baseline_changed=False)
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('prepare')
    for name in ('frameworks','charges','states','software','plan','output'):s.add_argument('--'+name,required=True)
    s.add_argument('--density-source',choices=('vacuum','responsive-trial'),default='vacuum')
    for name in ('dry-run','collect'):
        s=sub.add_parser(name);s.add_argument('--manifest',required=True)
        if name=='collect':s.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='prepare':
        r=prepare(a.frameworks,a.charges,a.states,a.software,a.plan,a.output,
            PROTOCOL if a.density_source=='vacuum' else TRIAL_PROTOCOL)
        print(json.dumps({'complete':r['complete'],'gates_pass':r['gates_pass']}))
    elif a.command=='dry-run':print(json.dumps({'tasks':len(validate(a.manifest)['tasks']),'new_native_calls':0}))
    else:
        r=collect(a.manifest,a.output);print(json.dumps({'complete':r['complete'],'gates_pass':r['gates_pass']}))


if __name__=='__main__':main()
