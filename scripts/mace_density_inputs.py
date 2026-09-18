"""Source-backed inputs for the frozen responsive-density hybrid; no fitting."""
from __future__ import annotations
import argparse
import copy
import json
import math
from pathlib import Path
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact,cache_key,paired,read_json,record,verify,write_new,xyz
from affordable_response import source_key
from affordable_peptide import SourceGraph
from mace_hybrid import check_atoms,write_xyz
from mace_mechanics import original_full,physical_id
from mace_omol_matched_h import normalized_physical
from mace_omol_charges import projection
from mace_omol_solvent import forcefield_background
from mace_omol_vacuum import METHOD,embedded_input,parse_endpoint

PROTOCOL='source_graph_normalized_responsive_density_inputs_v1'
QUANTUM_PROTOCOL='source_graph_r2scan3c_fixed_ff19SB_trial_density_v1'
NULLS=dict(reference=None,calibrated_class=None,combined_gradient=None,
           relaxation_correction=None,baseline_changed=False)


def normalize(physical_path,core_path):
    """Reuse the archived physical H rule; preserve heavy atoms and real caps."""
    p=read_json(physical_path);c=read_json(core_path)
    if p['status']!='prepared' or c['status']!='prepared':
        raise InvalidArtifact('prepared physical and source-core records required')
    if p['cofactor_charge_e'] or p.get('background_metals') or p['explicit_waters']:
        raise InvalidArtifact('this input version supports dry single-metal standard proteins only')
    if any(a['kind'] not in ('protein_source','selected_metal') for a in p['physical_atoms']):
        raise InvalidArtifact('unsupported nonstandard physical component')
    if c['explicit_water_inventory'] or c['protocol_id']!='generic_peptide_alpha_caps_native_r2scan3c_dev_v1':
        raise InvalidArtifact('declared dry peptide alpha-cap source required')
    if c['source_structure']!=p['source']:raise InvalidArtifact('core and protein source differ')
    graph=SourceGraph(verify(c['source_structure']),verify(c['topology_definition']))
    coords={source_key(graph.meta[k]):np.array(tuple(a.pos)) for k,a in graph.atoms.items()}
    original=original_full(p,graph,coords,np.array(c['selected_site']['xyz_A']),None)
    checks=normalized_physical(p,original)
    byid={a['id']:a for a in p['physical_atoms']};mapping=[];endpoints={}
    for metal in ('Ca','La'):
        old=xyz(verify(c['outputs'][metal]['xyz']));new=list(old);seen={0};moves=[]
        if old[0]!=(metal,*c['selected_site']['xyz_A']):raise InvalidArtifact('selected metal geometry differs')
        for item in c['atom_graph']['source_to_qm']:
            i=item['qm_index']
            if i in seen or i>=len(old):raise InvalidArtifact('duplicate or invalid core mapping')
            seen.add(i)
            if item['kind']=='source':
                pid=physical_id(item['source']);a=byid[pid]
                if old[i][0]!=a['element'] or not np.allclose(old[i][1:],coords[source_key(item['source'])],atol=1e-9,rtol=0):
                    raise InvalidArtifact('source core atom differs')
                if a['element']=='H':
                    new[i]=('H',*a['xyz_A'])
                    if new[i]!=old[i]:moves.append(dict(qm_index=i,physical_id=pid,before_A=list(old[i][1:]),after_A=a['xyz_A']))
                elif not np.allclose(old[i][1:],a['xyz_A'],atol=1e-9,rtol=0):raise InvalidArtifact('heavy coordinate changed')
                mapped=dict(qm_index=i,kind='source',physical_id=pid,source=item['source'])
            elif item['kind']=='sigma_link_H':
                x,y=(np.array(byid[physical_id(item[k])]['xyz_A']) for k in ('retained','omitted'))
                expected=x+item['length_A']*(y-x)/np.linalg.norm(y-x)
                if old[i][0]!='H' or not np.allclose(old[i][1:],expected,atol=1e-9,rtol=0):
                    raise InvalidArtifact('source cap geometry differs')
                mapped=copy.deepcopy(item)
            else:raise InvalidArtifact('unsupported source mapping')
            if metal=='Ca':mapping.append(mapped)
        charge=c['outputs'][metal]['charge']
        if seen!=set(range(len(old))) or charge!=c['ligand_formal_charge']+(2 if metal=='Ca' else 3):
            raise InvalidArtifact('core coverage/formal-charge ledger differs')
        endpoints[metal]=dict(atoms=new,charge=charge,spin_multiplicity=1,state=check_atoms(new,charge),
                              original_xyz=c['outputs'][metal]['xyz'],H_moves=moves)
    if endpoints['Ca']['atoms'][1:]!=endpoints['La']['atoms'][1:]:raise InvalidArtifact('paired ligand geometry differs')
    return dict(global_id=p['case_id'],normalized_global_preparation=record(physical_path),
        source_core=record(core_path),mapping=mapping,endpoints=endpoints,original_physical_atoms=original,
        geometry_checks=checks,evidence=p['evidence'],explicit_waters=p['explicit_waters'],cap_coordinates_unchanged=True)


def environment(case,proj,background):
    """The existing local residue-closure rule, before any quantum charge fit."""
    prep=read_json(verify(case['normalized_global_preparation']));core=read_json(verify(case['source_core']))
    physical=prep['physical_atoms'];ids=[a['id'] for a in physical];support=set(proj['physical_ids'])
    source_ids={a['physical_id'] for a in case['mapping'] if a['kind']=='source'}
    if not support<=set(ids) or set(ids)-set(background['charges_e'])!={'metal'}:
        raise InvalidArtifact('unaccounted physical/source atom')
    ledger={};seen=set()
    for item in core['charge_ledger']:
        s=item['source'];rid=f"{s['chain']}/{s['resnum']}/{s['insertion_code']}";q=item['formal_charge']
        if q:
            required={'ASP':{'CG','OD1','OD2'},'GLU':{'CD','OE1','OE2'}}.get(s['canonical_resname'])
            if q!=-1 or item['kind']!='sidechain' or required is None or rid in seen:
                raise InvalidArtifact('unsupported charged source functional group')
            if not {rid+'/'+n for n in required}<=source_ids:raise InvalidArtifact('incomplete charged source group')
            seen.add(rid)
        ledger[rid]=ledger.get(rid,0)+q
    if sum(ledger.values())!=core['ligand_formal_charge']:raise InvalidArtifact('core formal ledger fails')
    ff=background['charges_e'];env={i:0. if i in support else q for i,q in ff.items()};rows=[]
    for rid in sorted({i.rsplit('/',1)[0] for i in support if i in ff}):
        local=sorted(i for i in ff if i.rsplit('/',1)[0]==rid);outside=[i for i in local if i not in support]
        formal=round(math.fsum(ff[i] for i in local))
        if abs(math.fsum(ff[i] for i in local)-formal)>1e-6:raise InvalidArtifact('noninteger forcefield residue charge')
        target=formal-ledger.get(rid,0);before=math.fsum(env[i] for i in outside);delta=target-before
        recipients=set();bonds=[]
        for a,b in background['bonds']:
            if a in support and b in outside:recipients.add(b);bonds.append([a,b])
            elif b in support and a in outside:recipients.add(a);bonds.append([b,a])
        recipients=sorted(recipients)
        if not recipients and abs(delta)>1e-6:raise InvalidArtifact('no bonded exterior charge recipient')
        increment=delta/len(recipients) if recipients else 0.
        for i in recipients:env[i]+=increment
        if abs(math.fsum(env[i] for i in local)-target)>1e-6:raise InvalidArtifact('local charge closure fails')
        rows.append(dict(residue=rid,forcefield_formal_charge_e=formal,selected_fragment_formal_charge_e=ledger.get(rid,0),
            projection_support_ids=sorted(set(local)&support),exterior_target_e=target,original_exterior_charge_e=before,
            delta_e=delta,recipient_ids=recipients,recipient_bonds=bonds,increment_per_recipient_e=increment,
            original_recipient_charges_e=[ff[i] for i in recipients],new_recipient_charges_e=[env[i] for i in recipients]))
    values=[env.get(i,0.) for i in ids];expected=prep['protein_charge_e']-sum(ledger.values())
    if abs(math.fsum(values)-expected)>1e-6:raise InvalidArtifact('whole exterior charge closure fails')
    return dict(physical_atoms=physical,projection_support_ids=sorted(support),environment_charges_e=values,
        environment_charge=expected,boundary_ledger=rows,ligand_ledger=ledger,ligand_ledger_source=case['source_core'],
        forcefield_background=background,normalized_preparation=case['normalized_global_preparation'],
        source_case=case['source_core'],assembly=prep['assembly'],microstate=prep['microstate'],
        explicit_waters=prep['explicit_waters'],evidence=prep['evidence'],evidence_use=prep['evidence_use'],
        QM_charges=None,environment_correction=None)


def pointcharges(state):
    rows=[(q,*a['xyz_A']) for a,q in zip(state['physical_atoms'],state['environment_charges_e']) if q]
    return str(len(rows))+'\n'+''.join(' '.join(format(float(x),'.17g') for x in r)+'\n' for r in rows)


def prepare(config,output):
    started=time.monotonic();cpu=time.process_time();cfg=read_json(config)
    if cfg['protocol']!=PROTOCOL or not cfg['cases']:raise InvalidArtifact('explicit supported input configuration required')
    old=read_json(verify(cfg['quantum_implementation_parent']));root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    impl=root/'implementation';impl.mkdir();paths={k:verify(v) for k,v in old['implementation'].items()}
    for name in ('mace_density_inputs.py','mace_mechanics.py','mace_omol_matched_h.py','mace_omol_charges.py','mace_omol_solvent.py','mace_responsive_charges.py','mace_omol_vacuum.py'):
        paths[name]=Path(__file__).with_name(name)
    pins={}
    for name,src in paths.items():shutil.copyfile(src,impl/name);pins[name]=record(impl/name)
    m=dict(protocol_id=PROTOCOL,config=record(config),agreement=cfg['plan'],implementation=pins,cases={},**NULLS)
    for row in cfg['cases']:
        name=row['case_id']
        if not name.replace('_','').isalnum() or name in m['cases']:raise InvalidArtifact('invalid or duplicate case ID')
        pp=verify(row['physical']);cp=verify(row['core']);case=normalize(pp,cp)
        if case['global_id']!=name:raise InvalidArtifact('configuration/source identity differs')
        d=root/name;d.mkdir();write_new(d/'original_physical.json',case.pop('original_physical_atoms'))
        case['physical_atoms']=record(d/'original_physical.json')
        for metal,e in case['endpoints'].items():
            xp=d/(metal+'.xyz');write_xyz(xp,e.pop('atoms'));e['xyz']=record(xp)
        proj=projection(case,xyz(verify(case['endpoints']['Ca']['xyz'])))
        if projection(case,xyz(verify(case['endpoints']['La']['xyz'])))!=proj:raise InvalidArtifact('paired projection differs')
        write_new(d/'projection.json',proj);case['projection']=record(d/'projection.json')
        state=environment(case,proj,forcefield_background(read_json(pp)))
        write_new(d/'state.json',state);case['state']=record(d/'state.json')
        (d/'environment.pc').write_text(pointcharges(state));case['pointcharges']=record(d/'environment.pc')
        nearest=min(np.linalg.norm(np.array(a['xyz_A'])-np.array(q[1:])) for a,v in zip(state['physical_atoms'],state['environment_charges_e']) if v for q in xyz(verify(case['endpoints']['Ca']['xyz'])))
        if nearest<1e-6:raise InvalidArtifact('environment charge overlaps a QM nucleus/cap')
        case['nearest_QM_charge_separation_A']=float(nearest)
        write_new(d/'mapping.json',case);m['cases'][name]=record(d/'mapping.json')
    m['preparation_receipt']=dict(wall_seconds=time.monotonic()-started,CPU_seconds=time.process_time()-cpu)
    write_new(root/'preparation.json',m);return audit(root/'preparation.json')


def audit(path):
    m=read_json(path);cfg=read_json(verify(m['config']))
    if m['protocol_id']!=PROTOCOL or cfg['protocol']!=PROTOCOL:raise InvalidArtifact('input protocol differs')
    for pin in [m['agreement'],*m['implementation'].values()]:verify(pin)
    if set(m['cases'])!={r['case_id'] for r in cfg['cases']}:raise InvalidArtifact('case inventory differs')
    results={}
    for row in cfg['cases']:
        case=read_json(verify(m['cases'][row['case_id']]));fresh=normalize(verify(row['physical']),verify(row['core']))
        if read_json(verify(case['physical_atoms']))!=fresh['original_physical_atoms']:raise InvalidArtifact('original physical mapping changed')
        for key in ('mapping','evidence','normalized_global_preparation','source_core'):
            if case[key]!=fresh[key]:raise InvalidArtifact('source mapping/state differs: '+key)
        # normalized_physical already enforces the frozen 1e-12 A replay bound.
        # Its measured roundoff can differ between BLAS hosts; it is not an input.
        for value in (case['geometry_checks'],fresh['geometry_checks']):
            if value['maximum_replay_error_A']>1e-12:raise InvalidArtifact('H replay error exceeds declared bound')
        if ({k:v for k,v in case['geometry_checks'].items() if k!='maximum_replay_error_A'}!=
            {k:v for k,v in fresh['geometry_checks'].items() if k!='maximum_replay_error_A'}):
            raise InvalidArtifact('H replay inventory/state differs')
        for metal in ('Ca','La'):
            e=case['endpoints'][metal];f=fresh['endpoints'][metal]
            if xyz(verify(e['xyz']))!=f['atoms'] or {k:v for k,v in e.items() if k!='xyz'}!={k:v for k,v in f.items() if k!='atoms'}:
                raise InvalidArtifact('prepared core differs from source replay')
        proj=projection(case,xyz(verify(case['endpoints']['Ca']['xyz'])))
        from mace_responsive_charges import same_projection
        if not same_projection(read_json(verify(case['projection'])),proj):raise InvalidArtifact('source projection differs')
        state=environment(case,proj,forcefield_background(read_json(verify(row['physical']))))
        if state!=read_json(verify(case['state'])) or pointcharges(state)!=verify(case['pointcharges']).read_text():
            raise InvalidArtifact('permanent generating field differs')
        paired(verify(case['endpoints']['La']['xyz']),verify(case['endpoints']['Ca']['xyz']),case['endpoints']['La']['charge'],case['endpoints']['Ca']['charge'])
        results[row['case_id']]=dict(physical_atoms=len(state['physical_atoms']),core_atoms=len(fresh['endpoints']['Ca']['atoms']),
            environment_charge=state['environment_charge'],projection_support=len(proj['physical_ids']),status='source_replay_pass')
    return dict(status='pass',preparation=record(path),cases=results,new_scientific_calls=0)


def quantum_key(t,m):
    return cache_key(dict(task={k:v for k,v in t.items() if k!='cache_key'},
        **{k:m[k] for k in ('protocol_id','method','preparation','orca','implementation')}))


def prepare_quantum(preparation,output):
    audit(preparation);p=read_json(preparation);cfg=read_json(verify(p['config']));old=read_json(verify(cfg['quantum_implementation_parent']))
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir();pins={}
    for name,pin in p['implementation'].items():
        src=Path(__file__).with_name(name) if name in ('mace_density_inputs.py','mace_omol_vacuum.py') else verify(pin)
        shutil.copyfile(src,impl/name);pins[name]=record(impl/name)
    name='mace_responsive_charges.py'
    if name not in pins:
        shutil.copyfile(Path(__file__).with_name(name),impl/name);pins[name]=record(impl/name)
    m=dict(protocol_id=QUANTUM_PROTOCOL,preparation=record(preparation),agreement=p['agreement'],implementation=pins,
        method=METHOD,orca=old['orca'],tasks=[],compute_budget=None,wall_time_limit=None,
        execution_policy={k:pins[n] for k,n in [('task_runner','run_orca_task_manifest.py'),('runtime_renderer','render_orca_runtime_input.py')]},
        energy_scope=old['energy_scope'],variational_response_check=None,**NULLS)
    for name,pin in p['cases'].items():
        c=read_json(verify(pin))
        for metal,e in c['endpoints'].items():
            tid=name+'_'+metal;d=root/tid;d.mkdir()
            shutil.copyfile(verify(e['xyz']),d/'core.xyz');shutil.copyfile(verify(c['pointcharges']),d/'environment.pc')
            (d/'endpoint.inp').write_text(embedded_input(e['charge']))
            t=dict(task_id=tid,case_id=name,metal=metal,charge=e['charge'],multiplicity=1,task_type='analytic_gradient',
                input=record(d/'endpoint.inp'),xyz=record(d/'core.xyz'),pointcharges=record(d/'environment.pc'),
                output_path=str(d/'endpoint.out'),engrad_path=str(d/'endpoint.engrad'),source_mapping=pin)
            t['cache_key']=quantum_key(t,m);m['tasks'].append(t)
    m['new_DFT_calls']=len(m['tasks']);write_new(root/'manifest.json',m);return validate_quantum(root/'manifest.json')


def validate_quantum(path):
    from affordable_workflow import dry_run
    m=read_json(path);audit(verify(m['preparation']));p=read_json(verify(m['preparation']))
    parent=read_json(verify(read_json(verify(p['config']))['quantum_implementation_parent']))
    if m['protocol_id']!=QUANTUM_PROTOCOL or m['method']!=METHOD or m['orca']!=parent['orca']:
        raise InvalidArtifact('quantum method/source differs')
    if len(m['tasks'])!=2*len(p['cases']) or {t['task_id'] for t in m['tasks']}!={c+'_'+metal for c in p['cases'] for metal in ('Ca','La')}:
        raise InvalidArtifact('quantum task inventory differs')
    for pin in m['implementation'].values():verify(pin)
    for t in m['tasks']:
        c=read_json(verify(p['cases'][t['case_id']]));e=c['endpoints'][t['metal']]
        if (t['cache_key']!=quantum_key(t,m) or t['charge']!=e['charge'] or t['multiplicity']!=1 or
            t['xyz']['sha256']!=e['xyz']['sha256'] or t['pointcharges']['sha256']!=c['pointcharges']['sha256'] or
            verify(t['input']).read_text()!=embedded_input(e['charge'])):raise InvalidArtifact('quantum prepared state/field differs')
        verify(t['xyz']);verify(t['pointcharges'])
    return dry_run(path)


def collect_quantum(path):
    from ggr_sensitivity import executed
    validate_quantum(path);m,rows=executed(path)
    for t in m['tasks']:
        r=rows[t['task_id']]
        if r['status']=='complete':
            try:r.update(parse_endpoint(t,verify(r['output']),t['engrad_path'],permanent_field=True))
            except (ValueError,OSError) as exc:r.update(status='invalid',reason=str(exc),energy_hartree=None)
    return dict(protocol_id=QUANTUM_PROTOCOL,manifest=record(path),rows=rows,
        status='complete' if all(r['status']=='complete' for r in rows.values()) else 'incomplete',
        variational_response_check=None,variational_status='no_matching_vacuum_endpoint_requested',
        collection_implementation=record(__file__),parser_implementation=record(Path(__file__).with_name('mace_omol_vacuum.py')),**NULLS)


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare');p.add_argument('--config',required=True);p.add_argument('--output',required=True)
    for command in ('audit','prepare-quantum'):
        p=sub.add_parser(command);p.add_argument('--preparation',required=True)
        if command=='prepare-quantum':p.add_argument('--output',required=True)
    for command in ('dry-run','execute','collect'):
        p=sub.add_parser(command);p.add_argument('--manifest',required=True);p.add_argument('--output')
    a=parser.parse_args()
    if a.command=='prepare':r=prepare(a.config,a.output)
    elif a.command=='audit':r=audit(a.preparation)
    elif a.command=='prepare-quantum':r=prepare_quantum(a.preparation,a.output)
    elif a.command=='dry-run':r=validate_quantum(a.manifest)
    elif a.command=='collect':r=collect_quantum(a.manifest)
    else:
        from affordable_workflow import execute
        validate_quantum(a.manifest);r=execute(a.manifest)
    if a.command in ('dry-run','execute','collect') and a.output:write_new(a.output,r)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','tasks')},indent=2))


if __name__=='__main__':main()
