"""Four-cell, state-explicit EF3 reference scout on one consumed LanM source.

Reuses the existing native executor and allocation renderer. No whole-protein
vacuum subtraction, optimizer, classifier bands, or production default changes.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
import numpy as np
from affordable_common import (InvalidArtifact, HA_TO_KCAL, cache_key, read_json,
                               record, verify, write_new, xyz)
from affordable_workflow import dry_run, execute
from metal_environment_reference import ORCA, read_pcgrad
from metal_environment_electronic_state import METHOD, describe, input_text, parse
from metal_environment_force_checks import check_pins

PROTOCOL='nikasha_lanm_ef3_embedded_native_explicit_state_v1'
METALS=('La','Dy')
CONFIGURATIONS=('A','B')


def check_config(config):
    check_pins(config)
    for required in ('source_id','source_state','core_mapping','boundary_mapping','perturbation'):
        if required not in config:raise InvalidArtifact('missing preparation field: '+required)
    if set(config['configurations'])!=set(CONFIGURATIONS):
        raise InvalidArtifact('exactly source and declared perturbation required')
    original=None
    original_charges=None
    for label in CONFIGURATIONS:
        c=config['configurations'][label]
        if set(c['endpoints'])!=set(METALS):raise InvalidArtifact('La/Dy pair required')
        ca=[]
        for metal in METALS:
            e=c['endpoints'][metal]
            state=describe(verify(e['xyz']),metal,e['charge'],e['multiplicity'])
            for key,wanted in [('all_electron_count',state['all_electron_count']),('oxidation_state',state['oxidation_state_hypothesis']),('physical_spin_2S',state['alpha_minus_beta'])]:
                if key in e and e[key] != wanted:
                    raise InvalidArtifact('prepared electronic-state metadata differs: '+key)
            atoms=xyz(verify(e['xyz']));ca.append(atoms)
        target=[i for i,a in enumerate(ca[0]) if a[0]=='La']
        if len(target)!=1 or len(ca[0])!=len(ca[1]):raise InvalidArtifact('target inventory differs')
        index=target[0]
        if ca[1][index][0]!='Dy' or any(a[0]!=b[0] for i,(a,b) in enumerate(zip(*ca)) if i!=index):
            raise InvalidArtifact('paired nonmetal identity differs')
        if not np.array_equal(np.array([a[1:] for a in ca[0]]),np.array([a[1:] for a in ca[1]])):
            raise InvalidArtifact('paired coordinates differ')
        if c['endpoints']['La']['charge']!=c['endpoints']['Dy']['charge']:
            raise InvalidArtifact('isovalent pair has different charge')
        identity=([a[0] for a in ca[0]],c['endpoints']['La']['charge'])
        if original is None:original=identity
        elif identity!=original:raise InvalidArtifact('perturbation changed inventory or charge')
        pc=verify(c['pointcharges']);lines=pc.read_text().splitlines()
        values=np.array([[float(v) for v in line.split()] for line in lines[1:]])
        if values.shape!=(int(lines[0]),4) or not np.isfinite(values).all():
            raise InvalidArtifact('invalid prepared point charges')
        if original_charges is None:
            original_charges=values[:,0].copy()
        elif not np.array_equal(values[:,0], original_charges):
            raise InvalidArtifact('perturbation changed permanent-charge inventory or ordering')
    return config


def prepare(inputs,plan,output,*,mpi_ranks,workers):
    if mpi_ranks<1 or not 1<=workers<=4:raise InvalidArtifact('invalid four-cell resource layout')
    config=check_config(read_json(inputs));out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=out/'implementation';impl.mkdir();pins={}
    for source in sorted(Path(__file__).parent.glob('*.py')):
        target=impl/source.name;shutil.copyfile(source,target);pins[source.name]=record(target)
    shutil.copyfile(impl/'render_orca_runtime_input.py',impl/'_base_render_orca_runtime_input.py')
    shutil.copyfile(impl/'metal_environment_runtime.py',impl/'render_orca_runtime_input.py')
    for name in ('_base_render_orca_runtime_input.py','render_orca_runtime_input.py'):pins[name]=record(impl/name)
    tasks=[]
    for label in CONFIGURATIONS:
        c=config['configurations'][label]
        for metal in METALS:
            e=c['endpoints'][metal];td=out/f'{metal}_{label}';td.mkdir()
            shutil.copyfile(verify(e['xyz']),td/'core.xyz')
            shutil.copyfile(verify(c['pointcharges']),td/'environment.pc')
            (td/'endpoint.inp').write_text(input_text(e['charge'],e['multiplicity']))
            t=dict(task_id=td.name,metal=metal,configuration=label,charge=e['charge'],multiplicity=e['multiplicity'],
                   xyz=record(td/'core.xyz'),input=record(td/'endpoint.inp'),pointcharges=record(td/'environment.pc'),
                   output_path=str(td/'endpoint.out'),engrad_path=str(td/'endpoint.engrad'),task_type='analytic_gradient',
                   electronic_state=describe(td/'core.xyz',metal,e['charge'],e['multiplicity']))
            tasks.append(t)
    m=dict(protocol_id=PROTOCOL,method=METHOD,source_id=config['source_id'],inputs=record(inputs),agreement=record(plan),
           tasks=tasks,orca=record(ORCA),implementation=pins,
           execution_policy={'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},
           execution_resources={'mpi_ranks':mpi_ranks,'concurrent_tasks':workers},
           energy_scope='finite_embedded_electronic_component_only',classification=None,affinity=None,
           compute_budget=None,wall_time_limit=None,full_hybrid_qualified=False)
    for t in tasks:t['cache_key']=cache_key(dict(task=t,protocol=PROTOCOL,method=METHOD,inputs=m['inputs'],orca=m['orca'],implementation=pins))
    write_new(out/'manifest.json',m)
    return validate(out/'manifest.json')


def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['method']!=METHOD:raise InvalidArtifact('reference protocol differs')
    c=check_config(read_json(verify(m['inputs'])))
    if m['source_id'] != c['source_id']:
        raise InvalidArtifact('manifest source identity differs')
    layout=m['execution_resources']
    if type(layout['mpi_ranks']) is not int or layout['mpi_ranks'] < 1 or type(layout['concurrent_tasks']) is not int or not 1 <= layout['concurrent_tasks'] <= 4:
        raise InvalidArtifact('invalid explicit four-cell resource layout')
    for role,name in [('task_runner','run_orca_task_manifest.py'),('runtime_renderer','render_orca_runtime_input.py')]:
        if m['execution_policy'][role] != m['implementation'][name]:
            raise InvalidArtifact('execution policy differs from pinned implementation')
    for pin in m['implementation'].values():verify(pin)
    if len(m['tasks'])!=4 or {(t['metal'],t['configuration']) for t in m['tasks']}!={(a,b) for a in METALS for b in CONFIGURATIONS}:
        raise InvalidArtifact('four-cell inventory differs')
    for t in m['tasks']:
        conf=c['configurations'][t['configuration']];e=conf['endpoints'][t['metal']]
        if any(t[k]!=e[k] for k in ('charge','multiplicity')) or t['xyz']['sha256']!=e['xyz']['sha256']:
            raise InvalidArtifact('prepared state/geometry differs')
        if t['pointcharges']['sha256']!=conf['pointcharges']['sha256']:raise InvalidArtifact('field differs')
        verify(t['pointcharges'])
        if verify(t['input']).read_text()!=input_text(t['charge'],t['multiplicity']):raise InvalidArtifact('Hamiltonian differs')
        if describe(verify(t['xyz']),t['metal'],t['charge'],t['multiplicity'])!=t['electronic_state']:
            raise InvalidArtifact('state metadata differs')
        bare={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key(dict(task=bare,protocol=PROTOCOL,method=METHOD,inputs=m['inputs'],orca=m['orca'],implementation=m['implementation'])):
            raise InvalidArtifact('cache identity differs')
    return dry_run(manifest)



def torsion_tangent(config, label):
    """Physical derivative dXYZ/dtheta, Å/radian; caps follow their anchors."""
    mapping=read_json(verify(config['core_mapping'][label]))
    motion=config['perturbation'];axis=np.asarray(motion['axis_unit'],dtype=float)
    pivot=np.asarray(motion['pivot_xyz_A'],dtype=float)
    moving=set(motion['rotate_source_atoms'])
    if abs(np.linalg.norm(axis)-1)>1e-10 or not np.isfinite(pivot).all():
        raise InvalidArtifact('invalid physical torsion axis')
    found={a.get('source_id') for a in mapping if a.get('source_id') in moving}
    if found!=moving:raise InvalidArtifact('torsion atom mapping incomplete')
    tangent=np.zeros((len(mapping),3))
    for a in mapping:
        i=a['qm_index']
        if a['kind']=='cap':
            for side in ('retained','omitted'):
                if a[side+'_source_id'] in moving:
                    derivative=np.cross(axis,np.asarray(a[side+'_xyz_A'])-pivot)
                    tangent[i]+=np.asarray(a['jacobian_'+side])@derivative
        elif a.get('source_id') in moving:
            tangent[i]=np.cross(axis,np.asarray(a['xyz_A'])-pivot)
    return tangent


def collect(manifest):
    from ggr_sensitivity import executed
    validate(manifest);m,rows=executed(manifest)
    config=read_json(verify(m['inputs']))
    for t in m['tasks']:
        row=rows[t['task_id']]
        if row['status']!='complete':continue
        try:
            row.update(parse(t,verify(row['output']),t['engrad_path']))
            tangent=torsion_tangent(config,t['configuration'])
            gradient=np.asarray(row['gradient_kcal_mol_per_A'])
            row['torsion_gradient_kcal_mol_per_radian']=float(np.sum(gradient*tangent))
            row['torsion_force_kcal_mol_per_radian']=-row['torsion_gradient_kcal_mol_per_radian']
            p=Path(t['output_path']).with_name('endpoint.runtime.pcgrad')
            n=int(verify(t['pointcharges']).read_text().splitlines()[0]);g=read_pcgrad(p,n)
            row.update(pointcharge_gradient=record(p),pointcharge_gradient_units='Hartree/bohr',
                       pointcharge_gradient_max_abs=float(abs(g).max()))
        except (ValueError,OSError) as exc:
            row.update(status='invalid',reason=str(exc),energy_hartree=None)
    work={}
    for metal in METALS:
        a,b=[rows[f'{metal}_{label}'] for label in CONFIGURATIONS]
        work[metal]=(b['energy_hartree']-a['energy_hartree'])*HA_TO_KCAL if a['status']==b['status']=='complete' else None
    delta=work['Dy']-work['La'] if all(x is not None for x in work.values()) else None
    return dict(protocol_id=PROTOCOL,manifest=record(manifest),source_id=m['source_id'],rows=rows,
                status='complete' if all(x['status']=='complete' for x in rows.values()) else 'incomplete',
                per_metal_perturbation_work_kcal_mol=work,Dy_minus_La_perturbation_work_kcal_mol=delta,
                meaning='local fixed-state response, not affinity or protein-level preference',
                electronic_state_qualified=False,full_hybrid_qualified=False,classification=None,affinity=None)


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare')
    for k in ('inputs','plan','output'):a.add_argument('--'+k,required=True)
    a.add_argument('--mpi-ranks',type=int,required=True);a.add_argument('--workers',type=int,required=True)
    for op in ('dry-run','execute','collect'):
        a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
    a=p.parse_args()
    if a.op=='prepare':r=prepare(a.inputs,a.plan,a.output,mpi_ranks=a.mpi_ranks,workers=a.workers)
    elif a.op=='dry-run':r=validate(a.manifest)
    elif a.op=='collect':r=collect(a.manifest)
    else:
        validate(a.manifest);layout=read_json(a.manifest)['execution_resources']
        if int(os.environ.get('SLURM_NTASKS','0'))<layout['mpi_ranks']*layout['concurrent_tasks']:
            raise InvalidArtifact('MPI layout exceeds actual slots')
        # Renderer memory division must use the manifest's actual concurrency,
        # never an inherited/default count from another campaign.
        os.environ['METAL_ENV_WORKERS']=str(layout['concurrent_tasks'])
        r=execute(a.manifest)
    if a.op!='prepare' and a.output:write_new(a.output,r)
    print(json.dumps(r,indent=2))

if __name__=='__main__':main()
