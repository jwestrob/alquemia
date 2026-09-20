"""Fixed-geometry native GFN2 vacuum/ALPB transfer for archived compact MACE."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import re
import shutil

from affordable_common import InvalidArtifact, HA_TO_KCAL, cache_key, energy, read_json, record, verify, write_new, xyz
from run_orca_task_manifest import load_manifest_tasks, _completed_attempt_is_valid

PROTOCOL='native_OMOL_plus_native_GFN2_ALPB_transfer_v1'
METHOD='ORCA6.1.1_native_GFN2_ALPBwater_300K_v1'


def input_text(charge, multiplicity, medium, solver):
    if medium not in ('vacuum','alpb') or solver not in ('native','ordinary_tight'):
        raise InvalidArtifact('unsupported medium or numerical solver')
    header='! Native-GFN2-xTB NoAutostart'+(' ALPB(Water)' if medium=='alpb' else '')
    if solver=='ordinary_tight':header+=' TightSCF'
    return (header+'\n%maxcore 2000\n%method\n WriteXTBParam true\n ReadXTBParam false\nend\n'
        +'%scf\n SmearTemp 300\n UseXTBMixer '+('true' if solver=='native' else 'false')+'\nend\n'
        +f'* xyzfile {charge} {multiplicity} core.xyz\n')


def completed(manifest, task_id):
    m,ts=load_manifest_tasks(Path(manifest));t=next(t for t in ts if t['task_id']==task_id)
    op=t['output'];rp=Path(str(op)+'.execution.json')
    if not op.exists() or not rp.exists():return None
    if not _completed_attempt_is_valid(rp,op,manifest_sha256=record(manifest)['sha256'],task=t,
            runner_identity=m['execution_policy']['task_runner'],runtime_renderer_identity=m['execution_policy']['runtime_renderer']):return None
    return {'energy_hartree':energy(op),'output':record(op),'receipt':record(rp),
            'manifest':record(manifest),'task_id':task_id}


def reusable(task, manifests):
    for source in manifests:
        m=read_json(source)
        if m.get('protocol_id')!=PROTOCOL:raise InvalidArtifact('reuse protocol differs')
        for old in m['all_tasks']:
            if old['scientific_key']!=task['scientific_key']:continue
            result=m['reused'].get(old['task_id']) or completed(source,old['task_id'])
            if result and completed(verify(result['manifest']),result['task_id'])==result:return result
    return None


def prepare(inventory,agreement,output,case=None,representation=None,solver='native',reuse_manifest=None):
    inv=read_json(inventory);selected=set(case or [c['case_id'] for c in inv['cases']])
    if selected-{c['case_id'] for c in inv['cases']}:raise InvalidArtifact('unknown selected case')
    reps=set(representation or ['core','context'])
    if reps-{'core','context'}:raise InvalidArtifact('unsupported representation')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=out/'implementation';impl.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        dest=impl/p.name;shutil.copyfile(p,dest);pins[p.name]=record(dest)
    release=read_json(Path(__file__).resolve().parents[1]/'params/baseline_water_v1.json')['artifacts']
    tasks=[];reused={}
    for c in inv['cases']:
        if c['case_id'] not in selected:continue
        for rep in sorted(reps):
            representation_data=c['representations'][rep]
            for metal,e in representation_data['endpoints'].items():
                if metal not in ('Ca','La') or e['multiplicity']!=1:raise InvalidArtifact('unsupported metal/spin')
                original=verify(e['xyz']);coords=xyz(original)
                if coords[0][0]!=metal:raise InvalidArtifact('metal/source order differs')
                for medium in ('vacuum','alpb'):
                    tid=c['case_id']+'__'+rep+'__'+metal+'__'+medium+'__'+solver
                    if not re.fullmatch(r'[a-zA-Z0-9_.-]+',tid):raise InvalidArtifact('invalid task id')
                    directory=out/'tasks'/tid;directory.mkdir(parents=True)
                    xp=directory/'core.xyz';shutil.copyfile(original,xp)
                    ip=directory/'endpoint.inp';ip.write_text(input_text(e['charge'],e['multiplicity'],medium,solver))
                    task={'task_id':tid,'case_id':c['case_id'],'case':c['case_id'],'representation':rep,'metal':metal,
                          'medium':medium,'solver':solver,'charge':e['charge'],'multiplicity':e['multiplicity'],
                          'source_endpoint':e,'source_preparation':representation_data['preparation'],
                          'xyz':record(xp),'input':record(ip),'output_path':str(directory/'endpoint.out')}
                    task['scientific_key']=cache_key({'xyz_sha256':e['xyz']['sha256'],'charge':e['charge'],
                        'multiplicity':e['multiplicity'],'medium':medium,'solver':solver,'method':METHOD,
                        'input_body':ip.read_text(),'orca':release['orca']})
                    tasks.append(task);prior=reusable(task,reuse_manifest or [])
                    if prior:reused[tid]=prior
    m={'protocol_id':PROTOCOL,'method_id':METHOD,'inventory':record(inventory),'agreement':record(agreement),
        'orca':release['orca'],'implementation':pins,'execution_resources':{'mpi_ranks':8,'concurrent_tasks':8},
        'execution_policy':{'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},
        'selection':{'cases':sorted(selected),'representations':sorted(reps),'solver':solver},
        'all_tasks':tasks,'reused':reused,'tasks':[t for t in tasks if t['task_id'] not in reused],
        'reference':None,'new_high_level_evaluations':0,'baseline_changed':False}
    write_new(out/'manifest.json',m);return validate(out/'manifest.json')


def validate(manifest):
    m=read_json(manifest);inv=read_json(verify(m['inventory']));verify(m['agreement']);verify(m['orca'])
    if m['protocol_id']!=PROTOCOL or m['method_id']!=METHOD:raise InvalidArtifact('method differs')
    for pin in m['implementation'].values():verify(pin)
    expected={(c,r,z,s) for c in m['selection']['cases'] for r in m['selection']['representations'] for z in ('Ca','La') for s in ('vacuum','alpb')}
    actual={(t['case_id'],t['representation'],t['metal'],t['medium']) for t in m['all_tasks']}
    if actual!=expected or len(actual)!=len(m['all_tasks']):raise InvalidArtifact('endpoint coverage differs')
    for t in m['all_tasks']:
        c=next(c for c in inv['cases'] if c['case_id']==t['case_id']);source=c['representations'][t['representation']]
        e=source['endpoints'][t['metal']]
        if t['source_endpoint']!=e or t['source_preparation']!=source['preparation']:raise InvalidArtifact('source changes')
        verify(t['source_preparation']);verify(e['native_MACE_receipt']);verify(e['source_manifest'])
        if t['charge']!=e['charge'] or t['multiplicity']!=e['multiplicity'] or t['solver']!=m['selection']['solver']:
            raise InvalidArtifact('electronic state changes')
        if verify(t['xyz']).read_bytes()!=verify(e['xyz']).read_bytes():raise InvalidArtifact('coordinates changed')
        body=input_text(e['charge'],e['multiplicity'],t['medium'],t['solver'])
        if verify(t['input']).read_text()!=body:raise InvalidArtifact('input recipe differs')
        key=cache_key({'xyz_sha256':e['xyz']['sha256'],'charge':e['charge'],'multiplicity':e['multiplicity'],
            'medium':t['medium'],'solver':t['solver'],'method':METHOD,'input_body':body,'orca':m['orca']})
        if t['scientific_key']!=key:raise InvalidArtifact('scientific cache differs')
        if t['task_id'] in m['reused']:
            pin=m['reused'][t['task_id']]
            if completed(verify(pin['manifest']),pin['task_id'])!=pin:raise InvalidArtifact('reused result changed')
            old=read_json(verify(pin['manifest']));oldtask=next(x for x in old['all_tasks'] if x['task_id']==pin['task_id'])
            if oldtask['scientific_key']!=key:raise InvalidArtifact('reused result describes different calculation')
    if m['tasks']!=[t for t in m['all_tasks'] if t['task_id'] not in m['reused']]:raise InvalidArtifact('reuse/task partition changed')
    if m['tasks']:
        from affordable_workflow import dry_run
        dry_run(manifest)
    return {'status':'validated','manifest':record(manifest),'tasks':len(m['tasks']),'reused':len(m['reused'])}


def diagnostics(pin,task):
    """Parse real native state/charges/parameters; expose unprinted solvent terms."""
    op=verify(pin['output']);text=op.read_text();receipt=read_json(verify(pin['receipt']))
    if receipt['orca_version']!='6.1.1':raise InvalidArtifact('unsupported ORCA version')
    if not re.search(r'XTB-Hamiltonian\s+Method\s+\.{4} GFN2-xTB',text):
        raise InvalidArtifact('native GFN2 Hamiltonian not confirmed')
    parameter=op.parent/'endpoint.runtime.xtb.json';params=read_json(parameter)
    if params['meta']['name']!='GFN2-xTB':raise InvalidArtifact('native parameter set differs')
    prop=op.parent/'endpoint.runtime.property.txt';pt=prop.read_text()
    def integer(name):
        hits=re.findall(r'&'+name+r' \[&Type "Integer"\] (-?\d+)',pt.split('$Calculation_Info',1)[1].split('$End',1)[0])
        if len(hits)!=1:raise InvalidArtifact('missing/ambiguous native '+name)
        return int(hits[0])
    atoms=xyz(verify(task['xyz']))
    count=sum(sum(params['element'][a[0]]['refocc']) for a in atoms)-task['charge']
    if (integer('Charge'),integer('Mult'),integer('NumOfAtoms'),integer('NumOfElectrons'))!=(task['charge'],task['multiplicity'],len(atoms),count):
        raise InvalidArtifact('native valence/state accounting differs')
    block=pt.split('&AtomicCharges',1)[1].split('$End',1)[0]
    charges=[(int(i),float(q)) for i,q in re.findall(r'^\s*(\d+)\s+([-+0-9.eE]+)\s*$',block,re.M)]
    if [i for i,q in charges]!=list(range(len(atoms))):raise InvalidArtifact('native charge ordering differs')
    values=[q for i,q in charges]
    if not all(math.isfinite(q) for q in values) or abs(sum(values)-task['charge'])>5e-4:
        raise InvalidArtifact('native charge closure failed')
    max_charge=max(abs(q) for q in values)
    components={}
    for name in ('Total Energy','Nuclear Repulsion','Electronic Energy','One Electron Energy','Two Electron Energy'):
        matches=re.findall(r'^'+re.escape(name)+r'\s*:\s*([-+0-9.eE]+) Eh',text,re.M)
        components[name]=float(matches[-1]) if matches else None
    if components['Total Energy'] is None or abs(components['Total Energy']-pin['energy_hartree'])>1e-10:
        raise InvalidArtifact('native total/property energy differs')
    runtime=verify(receipt['artifacts']['runtime_input']).read_text()
    body=verify(task['input']).read_text()
    if ('ALPB(Water)' in runtime)!=(task['medium']=='alpb') or ('ALPB(Water)' in body)!=(task['medium']=='alpb'):
        raise InvalidArtifact('runtime solvent selection differs')
    if 'SmearTemp 300' not in runtime:raise InvalidArtifact('runtime electronic temperature differs')
    return {'strict_component_audit':'native_state_and_parameters_verified',
        'parameter_export':record(parameter),'property_file':record(prop),
        'native_valence_electrons':count,'metal_reference_occupation':params['element'][task['metal']]['refocc'],
        'mulliken_charges_e':values,'charge_sum_e':sum(values),'max_absolute_charge_e':max_charge,
        'charge_sanity_status':'pass' if max_charge<=4 else 'outside_declared_range',
        'printed_SCF_components_hartree':components,'ALPB_component_hartree':None,
        'solvent_parameter_export':None,'solvent_detail_status':'native installed defaults; no separate ALPB term/radii printed',
        'electronic_temperature_K':300,'electronic_entropy_component_hartree':None,
        'PMIX_diagnostic_count':text.count('PMIX ERROR:'),'normal_termination_confirmed':True}


def collect(manifest,output):
    validate(manifest);m=read_json(manifest);endpoints={};rows=[]
    for t in m['all_tasks']:
        pin=m['reused'].get(t['task_id']) or completed(manifest,t['task_id'])
        if pin:
            try:r={'status':'complete',**pin,**diagnostics(pin,t)}
            except InvalidArtifact as exc:r={'status':'unsupported','reason':str(exc),'energy_hartree':None}
        else:
            op=Path(t['output_path']);rp=Path(str(op)+'.execution.json')
            r={'status':'unavailable' if not op.exists() else 'failed_or_incomplete','energy_hartree':None,
               'available_artifacts':[record(p) for p in (op,rp) if p.exists()]}
        endpoints[(t['case_id'],t['representation'],t['metal'],t['medium'])]=r
    for case in m['selection']['cases']:
        for rep in m['selection']['representations']:
            for metal in ('Ca','La'):
                eps={s:endpoints[(case,rep,metal,s)] for s in ('vacuum','alpb')}
                good=all(e['status']=='complete' for e in eps.values())
                rows.append({'case_id':case,'representation':rep,'metal':metal,'solver':m['selection']['solver'],
                    'status':'complete' if good else 'unavailable','endpoints':eps,
                    'vacuum_hartree':eps['vacuum']['energy_hartree'],'alpb_hartree':eps['alpb']['energy_hartree'],
                    'delta_solv_hartree':eps['alpb']['energy_hartree']-eps['vacuum']['energy_hartree'] if good else None})
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'inventory':m['inventory'],
        'status':'complete' if all(r['status']=='complete' for r in rows) else 'incomplete','rows':rows,
        'component_audit':'successful endpoints audited; failed endpoints unavailable; separate ALPB components unprinted','S_kcal_mol':None,'baseline_changed':False}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='operation',required=True)
    q=s.add_parser('prepare')
    for k in ('inventory','agreement','output'):q.add_argument('--'+k,required=True)
    q.add_argument('--case',action='append');q.add_argument('--representation',action='append');q.add_argument('--reuse-manifest',action='append')
    q.add_argument('--solver',choices=('native','ordinary_tight'),default='native')
    for name in ('validate','collect','execute'):
        q=s.add_parser(name);q.add_argument('--manifest',required=True)
        if name=='collect':q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('operation')
    if op=='execute':
        validate(a['manifest'])
        from affordable_workflow import execute
        result=execute(a['manifest'])
    else:result=globals()[op](**a)
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))

if __name__=='__main__':main()
