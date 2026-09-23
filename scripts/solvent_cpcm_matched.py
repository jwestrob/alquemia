"""Matched ordinary-SCF vacuum controls for the actual CPCM solver switch."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import re
import shutil

from affordable_common import InvalidArtifact, HA_TO_KCAL, read_json, record, verify, write_new, xyz
from affordable_workflow import dry_run, execute as run_existing
from compact_solvation import completed, diagnostics
from compact_solvation_qualification import restart_diagnostics
from solvent_cpcm_pilot import PILOT, SETTINGS as CAVITY, validate as validate_primary
from mace_hybrid import EV_TO_KCAL

PROTOCOL = 'native_OMOL_plus_ordinary_GFN2_CPCM_water_transfer_v1'
SETTINGS = {**CAVITY, 'UseXTBMixer': False, 'TolE_hartree': 1e-6,
            'vacuum_restart': 'explicit_MORead_archived_native_vacuum_GBW',
            'CPCM_source': 'actual_primary_ordinary_SCF_no_fresh_CPCM_calls'}


def recipe(charge, multiplicity):
    return ('! Native-GFN2-xTB\n%maxcore 2000\n'
            '%method\n WriteXTBParam true\n ReadXTBParam false\nend\n'
            '%scf\n SmearTemp 300\n UseXTBMixer false\n MaxIter 500\n'
            ' Guess MORead\n MOInp "seed_source.gbw"\nend\n'
            f'* xyzfile {charge} {multiplicity} core.xyz\n')


def prepare(primary_manifest, agreement, output):
    validate_primary(primary_manifest); parent = read_json(primary_manifest)
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    tasks = []
    for c in parent['cases']:
        for metal, e in c['endpoints'].items():
            tid = f'{c["case_id"]}__{metal}__vacuum_ordinary'
            td = out/'tasks'/tid; td.mkdir(parents=True)
            xp = td/'core.xyz'; shutil.copyfile(verify(e['xyz']), xp)
            ip = td/'endpoint.inp'; ip.write_text(recipe(e['charge'], e['multiplicity']))
            source = verify(e['GFN2']['vacuum']['endpoint']['output']).parent/'endpoint.runtime.gbw'
            seed = td/'seed_source.gbw'; shutil.copyfile(source, seed)
            tasks.append({'task_id': tid, 'case_id': c['case_id'], 'case': c['case_id'],
                          'metal': metal, 'medium': 'vacuum', 'representation': 'context',
                          'charge': e['charge'], 'multiplicity': e['multiplicity'],
                          'xyz': record(xp), 'input': record(ip), 'output_path': str(td/'endpoint.out'),
                          'restart': {'source': record(source), 'copy': record(seed)},
                          'CPCM_source_task_id': next(t['task_id'] for t in parent['tasks'] if (t['case_id'],t['metal'])==(c['case_id'],metal))})
    impl = out/'implementation'; impl.mkdir(); pins = {}
    for p in Path(__file__).parent.glob('*.py'):
        target=impl/p.name; shutil.copyfile(p,target); pins[p.name]=record(target)
    m={'protocol_id': PROTOCOL, 'settings': SETTINGS, 'agreement': record(agreement),
       'primary_manifest': record(primary_manifest), **{k:parent[k] for k in ('model','software','orca')},
       'tasks': tasks, 'implementation': pins,
       'execution_resources': {'mpi_ranks':8, 'concurrent_tasks':8},
       'execution_policy': {'task_runner':pins['run_orca_task_manifest.py'], 'runtime_renderer':pins['render_orca_runtime_input.py']},
       'new_GFN2_vacuum_calls':8, 'new_CPCM_calls':0, 'new_DFT_MACE_calls':0, 'reference':None,
       'failed_or_unfinished_primary_CPCM': 'unavailable; never substitute ALPB or vacuum', 'production_changed':False}
    mp=out/'manifest.json'; write_new(mp,m); v=validate(mp); write_new(out/'PREFLIGHT.json',v); return v


def validate(manifest):
    m=read_json(manifest); parent=read_json(verify(m['primary_manifest']))
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS: raise InvalidArtifact('matched numerical model changed')
    if m['new_GFN2_vacuum_calls']!=8 or m['new_CPCM_calls']!=0: raise InvalidArtifact('fresh-call scope changed')
    expected={(c,z) for c in PILOT for z in ('Ca','La')}
    if len(m['tasks'])!=8 or {(t['case_id'],t['metal']) for t in m['tasks']}!=expected: raise InvalidArtifact('matched task scope changed')
    for key in ('model','software','orca'):
        if m[key]!=parent[key]: raise InvalidArtifact('original electronic identity differs')
    verify(m['agreement']); verify(m['orca'])
    for pin in m['implementation'].values(): verify(pin)
    for t in m['tasks']:
        e=next(c for c in parent['cases'] if c['case_id']==t['case_id'])['endpoints'][t['metal']]
        if (t['charge'],t['multiplicity'])!=(e['charge'],e['multiplicity']): raise InvalidArtifact('matched state differs')
        if verify(t['xyz']).read_bytes()!=verify(e['xyz']).read_bytes(): raise InvalidArtifact('matched coordinates differ')
        if verify(t['input']).read_text()!=recipe(t['charge'],t['multiplicity']): raise InvalidArtifact('matched recipe differs')
        source=verify(e['GFN2']['vacuum']['endpoint']['output']).parent/'endpoint.runtime.gbw'
        if record(source)!=t['restart']['source'] or verify(t['restart']['copy']).read_bytes()!=source.read_bytes():
            raise InvalidArtifact('actual native orbital restart differs')
        old=next(p for p in parent['tasks'] if p['task_id']==t['CPCM_source_task_id'])
        if (old['case_id'],old['metal'])!=(t['case_id'],t['metal']): raise InvalidArtifact('CPCM reuse identity differs')
        if not Path(t['output_path']).is_relative_to(Path(manifest).resolve().parent/'tasks'):
            raise InvalidArtifact('matched task escapes workspace')
    dry_run(manifest)
    return {'status':'validated','manifest':record(manifest),'fresh_vacuum':8,'fresh_CPCM':0,'CPCM_source_slots':8}


def audit_ordinary(endpoint, task, source, restart=False):
    result=diagnostics(endpoint,task); text=verify(endpoint['output']).read_text()
    if result['charge_sanity_status']!='pass': raise InvalidArtifact('ordinary atomic-charge check failed')
    if 'INFO: Using special xTB SCF mixer' in text: raise InvalidArtifact('matched endpoint retained special mixer')
    tol=re.findall(r'Energy Change\s+TolE\s+\.{4}\s+([-+0-9.eE]+)',text)
    if len(tol)!=1 or float(tol[0])!=SETTINGS['TolE_hartree']: raise InvalidArtifact('ordinary convergence tolerance differs')
    if read_json(verify(result['parameter_export']))!=read_json(verify(source['GFN2']['vacuum']['audit']['parameter_export'])):
        raise InvalidArtifact('ordinary parameter set differs')
    result['effective_solver']='ordinary_ORCA_SCF'; result['TolE_hartree']=float(tol[0])
    result['warnings']=[s for s in text.splitlines() if 'WARNING' in s.upper() or 'IGNOR' in s.upper()]
    if restart:
        r=restart_diagnostics(text,'ordinary_explicit_gbw'); result['restart']=r
        if r['restart_status']!='orbital_restart_confirmed': raise InvalidArtifact('actual explicit MORead unconfirmed')
    if task['medium']=='cpcm': result['CPCM']=boundary_audit(text,task)
    elif 'CPCM SOLVATION MODEL' in text: raise InvalidArtifact('vacuum contains CPCM boundary')
    return result


def boundary_audit(text,task):
    if 'CPCM SOLVATION MODEL' not in text or 'GAUSSIAN VDW' not in text: raise InvalidArtifact('actual CPCM cavity not confirmed')
    block=text.split('CPCM SOLVATION MODEL',1)[1].split('Overall time for CPCM initialization',1)[0]
    scalar=lambda name: float(re.findall(r'^\s*'+re.escape(name)+r'\s+\.\.\.\s+([-+0-9.eE]+)',block,re.M)[0])
    eps,refrac=scalar('Epsilon'),scalar('Refrac')
    if abs(eps-CAVITY['epsilon'])>1e-4 or abs(refrac-CAVITY['refrac'])>1e-4: raise InvalidArtifact('water dielectric defaults differ')
    density=[float(v) for v in re.findall(r'Threshold for (?:non-H|H) atoms\s+\.\.\.\s+([-+0-9.eE]+)',block)]
    if density!=[5.,5.]: raise InvalidArtifact('surface discretization differs')
    radii={z:float(r) for z,r in re.findall(r'Radius for (\w+)\s+used is\s+[-+0-9.eE]+ Bohr \(=\s+([-+0-9.eE]+) Ang\.\)',block)}
    for z in {a[0] for a in xyz(verify(task['xyz']))}:
        if z not in radii: raise InvalidArtifact('element cavity radius unavailable')
        if z in CAVITY['expected_radii_A'] and abs(radii[z]-CAVITY['expected_radii_A'][z])>1e-4:
            raise InvalidArtifact('element cavity default differs')
    hits=re.findall(r'^\s*CPCM Dielectric\s*:\s*([-+0-9.eE]+) Eh',text,re.M)
    if not hits or not math.isfinite(float(hits[-1])) or abs(float(hits[-1]))<1e-8: raise InvalidArtifact('actual reaction-field contribution unavailable')
    return {'epsilon':eps,'refrac':refrac,'radii_A':radii,'surface_density_A_minus2':density,
            'printed_dielectric_energy_hartree':float(hits[-1]),'setup_text':block,
            'interpretation':'printed electrostatic component; complete transfer is the difference of matched total energies'}


def collect(manifest,output):
    validate(manifest); m=read_json(manifest); primary=verify(m['primary_manifest']); parent=read_json(primary)
    rows=[]; endpoints=[]
    for c in parent['cases']:
        components={}
        for metal,source in c['endpoints'].items():
            t=next(t for t in m['tasks'] if (t['case_id'],t['metal'])==(c['case_id'],metal))
            old=next(p for p in parent['tasks'] if p['task_id']==t['CPCM_source_task_id'])
            pair={}
            for medium,mp,task,restart in [('vacuum',manifest,t,True),('cpcm',primary,old,False)]:
                e=completed(mp,task['task_id']); op=Path(task['output_path'])
                r={'case_id':c['case_id'],'metal':metal,'medium':medium,'task_id':task['task_id'],
                   'status':'unavailable','energy_hartree':None,'raw_converged_endpoint':e,
                   'artifacts':[record(p) for p in (op,Path(str(op)+'.execution.json')) if p.exists()]}
                if e:
                    try: r.update(status='complete',**e,audit=audit_ordinary(e,task,source,restart))
                    except (InvalidArtifact,IndexError,KeyError,ValueError) as exc: r['reason']=str(exc)
                else: r['reason']='failed, nonconverged or unrun'
                pair[medium]=r; endpoints.append(r)
            good=all(p['status']=='complete' for p in pair.values())
            vac=pair['vacuum']['energy_hartree']; cpcm=pair['cpcm']['energy_hartree']
            transfer=(cpcm-vac)*HA_TO_KCAL if good else None
            components[metal]={'ordinary_vacuum_hartree':vac,'ordinary_CPCM_hartree':cpcm,
                'native_vacuum_hartree':source['GFN2']['vacuum']['endpoint']['energy_hartree'],
                'ordinary_minus_native_vacuum_kcal_mol':(vac-source['GFN2']['vacuum']['endpoint']['energy_hartree'])*HA_TO_KCAL if vac is not None else None,
                'MACE_vacuum_eV':source['MACE_energy_eV'],'CPCM_transfer_kcal_mol':transfer,
                'composite_energy_kcal_mol':source['MACE_energy_eV']*EV_TO_KCAL+transfer if good else None}
        good=all(x['composite_energy_kcal_mol'] is not None for x in components.values())
        R=components['Ca']['composite_energy_kcal_mol']-components['La']['composite_energy_kcal_mol'] if good else None
        old=c['old_result']['R0']['composite_R_model_kcal_mol']
        rows.append({'case_id':c['case_id'],'status':'available' if good else 'unavailable',
            'expected_class':c['old_result']['expected_class'],'label_scope':c['old_result']['label_scope'],
            'components':components,'CPCM_R':R,'released_ALPB_R':old,'CPCM_minus_ALPB_R':R-old if good else None,
            'new_reference':None,'new_decision':None,'affinity_claim':False})
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'endpoints':endpoints,'rows':rows,
        'available_pairs':sum(r['status']=='available' for r in rows),'case_denominator':4,
        'complete_endpoints':sum(e['status']=='complete' for e in endpoints),'endpoint_denominator':16,
        'fresh_vacuum_calls':8,'fresh_CPCM_calls':0,'reused_CPCM_complete':sum(e['status']=='complete' and e['medium']=='cpcm' for e in endpoints),
        'reference':None,'production_changed':False}
    write_new(output,result); return {k:v for k,v in result.items() if k not in ('rows','endpoints')}


def execute(manifest):
    validate(manifest); return run_existing(manifest)


def main():
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='operation',required=True)
    a=sub.add_parser('prepare'); a.add_argument('--primary-manifest',required=True); a.add_argument('--agreement',required=True); a.add_argument('--output',required=True)
    for name in ('validate','execute','collect'):
        a=sub.add_parser(name); a.add_argument('--manifest',required=True)
        if name=='collect': a.add_argument('--output',required=True)
    args=vars(p.parse_args()); op=args.pop('operation'); print(json.dumps(globals()[op](**args)))


if __name__=='__main__': main()
