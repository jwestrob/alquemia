"""Configured full/core MACE short components using the qualified worker."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from mace_hybrid import accepted_attempt,check_atoms,EV_TO_KCAL
from mace_global_benchmark import numerical_parent_gate

SCHEMA='alquemia.mace_density_short.v1'
PROTOCOL='source_graph_POLAR_medium_exact_short_context_v1'
CHECKPOINT='fab8b8713c832f31a2a853aaa22fd638be8a369cbf5095e6b3e982a18d10e93a'


def key(t,m):
    return cache_key(dict(task={k:v for k,v in t.items() if k!='cache_key'},
        **{k:m[k] for k in ('protocol_id','preparation','model','software','implementation','reference_manifest')}))


def prepare(preparation,quantum_manifest,reference_manifest,output):
    from mace_density_inputs import audit,validate_quantum
    audit(preparation);validate_quantum(quantum_manifest)
    p=read_json(preparation);qm=read_json(quantum_manifest);old=read_json(reference_manifest)
    ref=read_json(verify(old['short_reference']))
    if ref['status']!='complete' or not ref['numerical_checks_pass'] or old['model']['checkpoint']['sha256']!=CHECKPOINT:
        raise InvalidArtifact('qualified exact short readout reference required')
    if qm['preparation']!=record(preparation):raise InvalidArtifact('quantum/short physical sources differ')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    paths={k:verify(v) for parent in (qm,old) for k,v in parent['implementation'].items()}
    for name in ('mace_density_short.py','mace_hybrid.py','mace_density_inputs.py','mace_responsive_charges.py'):
        paths[name]=Path(__file__).with_name(name)
    for name,path in paths.items():shutil.copyfile(path,impl/name)
    pins={p.name:record(p) for p in impl.glob('*.py')};model=copy.deepcopy(old['model'])
    model['preparation_policy']='source_graph_exact_paired_normalized_coordinates_v1'
    m=dict(schema_version=SCHEMA,protocol_id=PROTOCOL,preparation=record(preparation),quantum_manifest=record(quantum_manifest),
        reference_manifest=record(reference_manifest),short_reference=old['short_reference'],numerical_reference=old['numerical_reference'],
        agreement=qm['agreement'],software=old['software'],model=model,implementation=pins,tasks=[],baseline_changed=False,
        reference=None,calibrated_class=None,combined_gradient=None)
    for name,pin in p['cases'].items():
        c=read_json(verify(pin));physical=read_json(verify(c['normalized_global_preparation']))
        for kind,ends in [('core',c['endpoints']),('full',physical['endpoints'])]:
            for metal,e in ends.items():
                t=dict(task_id=f'{name}_{metal}_{kind}',case_id=name,metal=metal,kind=kind,source_mapping=pin,
                    xyz=e['xyz'],charge=e['charge'],spin_multiplicity=1,energy_component='interaction_energy')
                if kind=='full' and physical.get('background_metals'):
                    t['background_calcium_indices']=[i for i,a in enumerate(physical['physical_atoms']) if a['kind']=='background_metal']
                t['cache_key']=key(t,m);m['tasks'].append(t)
    m['new_MACE_calls']=len(m['tasks']);write_new(root/'manifest.json',m);return validate(root/'manifest.json')


def validate(path):
    m=read_json(path);p=read_json(verify(m['preparation']));old=read_json(verify(m['reference_manifest']))
    if m['schema_version']!=SCHEMA or m['protocol_id']!=PROTOCOL:raise InvalidArtifact('short protocol differs')
    expected=copy.deepcopy(old['model']);expected['preparation_policy']='source_graph_exact_paired_normalized_coordinates_v1'
    if m['model']!=expected or m['software']!=old['software'] or m['model']['checkpoint']['sha256']!=CHECKPOINT:
        raise InvalidArtifact('short Hamiltonian/checkpoint/software differs')
    if numerical_parent_gate(m)['status']!='pass':raise InvalidArtifact('qualified numerical parent required')
    for pin in [m['agreement'],m['quantum_manifest'],m['short_reference'],m['software'],m['model']['checkpoint'],*m['implementation'].values()]:verify(pin)
    sw=read_json(verify(m['software']))
    for pin in [sw['python'],sw['requirements'],*read_json(verify(sw['backend_source_inventory']))['files']]:verify(pin)
    if len(m['tasks'])!=4*len(p['cases']) or {t['task_id'] for t in m['tasks']}!={f'{n}_{e}_{k}' for n in p['cases'] for e in ('Ca','La') for k in ('core','full')}:
        raise InvalidArtifact('short task inventory differs')
    for t in m['tasks']:
        c=read_json(verify(p['cases'][t['case_id']]));physical=read_json(verify(c['normalized_global_preparation']))
        e=(c if t['kind']=='core' else physical)['endpoints'][t['metal']]
        if t['xyz']!=e['xyz'] or t['charge']!=e['charge'] or t['spin_multiplicity']!=1 or t['energy_component']!='interaction_energy':
            raise InvalidArtifact('short physical state/component differs')
        background=[i for i,a in enumerate(physical['physical_atoms']) if a['kind']=='background_metal'] if t['kind']=='full' else []
        if t.get('background_calcium_indices',[])!=background:raise InvalidArtifact('background calcium mapping differs')
        check_atoms(xyz(verify(t['xyz'])),t['charge'],background_calcium_indices=background)
        if t['cache_key']!=key(t,m):raise InvalidArtifact('short scientific cache differs')
    return dict(status='pass',tasks=len(m['tasks']),manifest=record(path))


def collect(path):
    validate(path);m=read_json(path);mp=Path(path).resolve();rows={};attempts=[]
    for t in m['tasks']:
        accepted=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp)
            attempts.append(dict(task_id=t['task_id'],path=str(a),accepted=r is not None,
                receipt=record(a/'receipt.json') if (a/'receipt.json').exists() else None))
            if r is not None:accepted.append(r)
        rows[t['task_id']]=accepted[-1] if accepted else dict(status='unavailable',energy_eV=None)
    complete=all(r['status']=='computed' for r in rows.values());cases={}
    if complete:
        for name in read_json(verify(m['preparation']))['cases']:
            ends={}
            for metal in ('Ca','La'):
                full,core=(rows[f'{name}_{metal}_{kind}']['energy_eV'] for kind in ('full','core'))
                ends[metal]=dict(full_eV=full,core_eV=core,short_context_kcal=(full-core)*EV_TO_KCAL)
            cases[name]=dict(endpoints=ends,short_context_R_kcal=ends['Ca']['short_context_kcal']-ends['La']['short_context_kcal'])
    return dict(protocol_id=PROTOCOL,manifest=record(path),status='complete' if complete else 'incomplete',
        rows=rows,attempts=attempts,cases=cases,baseline_changed=False,full_hybrid_score=None,calibrated_class=None,
        force_scope='negative_gradient_of_learned_short_component_only',combined_gradient=None)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('prepare')
    for n in ('preparation','quantum-manifest','reference-manifest','output'):s.add_argument('--'+n,required=True)
    for op in ('dry-run','collect'):
        s=sub.add_parser(op);s.add_argument('--manifest',required=True);s.add_argument('--output')
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.preparation,a.quantum_manifest,a.reference_manifest,a.output)
    elif a.command=='dry-run':r=validate(a.manifest)
    else:r=collect(a.manifest)
    if a.command!='prepare' and a.output:write_new(a.output,r)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','attempts','cases')}))
