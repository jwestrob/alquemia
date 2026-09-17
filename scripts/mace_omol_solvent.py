"""Full-boundary frozen QM/ff19SB OBC2 challenger; no calibrated affinity claim."""
from __future__ import annotations
import argparse
import copy
import json
import math
import os
from pathlib import Path
import resource
import shutil
import time
import traceback
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_file_checks import cached_file_checks
from mace_gb import MODEL as OLD_GB_MODEL, RADII, SCALES, EV_TO_KJ, build_system
from mace_hybrid import accepted_attempt, rotation, write_xyz

SCHEMA='alquemia.mace_QMFF_gb.v1'
PROTOCOL='normalized_QM_projected_ff19SB_full_OBC2_vacuum_hybrid_v1'
COMPONENT='OBC2_reaction_energy_QMFF_charges'
CHARGE_MANIFEST_SHA='31caf15249cbfab1ff381f1185c7243e77c69a003c735f52ca80e5ad63bc3fc6'
CHARGE_RESULT_SHA='2689bcdfa9b3f2d0358c5e21fa8b5b38533b36502c7e5c3f9826ee7487041d10'
CASES=('GGR_extended','GGR_connected','ALPHA_1F6S','ALPHA_6IP9')
MODEL={**OLD_GB_MODEL,'charge_representation':'native_CHELPG_barycentric_QM_plus_local_ff19SB_exterior',
       'boundary_rule':'residue_formal_ledger_direct_bond_neighbor_equal_redistribution_v1',
       'energy_definition':'full_physical_boundary_GB_reaction_only; no_bare_Coulomb; no_core_CPCM',
       'vacuum_parent':'masked_omol_matched_normalized_H_vacuum_hybrid_v1'}
TOL={'charge_e':5e-5,'boundary_charge_e':1e-6,'numerical_kcal_mol':.01,'identity_kcal_mol':1e-6,
     'algebra_kcal_scale':1e-7,'partition_kcal_scale':2.,'ordering_kcal_scale':.02}
NULLS={'aqueous_affinity_score':None,'calibrated_class':None,'reference':None,'combined_gradient':None,
       'response_status':'response_model_not_validated','relaxation_correction':None,'baseline_changed':False}


def rid(pid):return pid.rsplit('/',1)[0]


def forcefield_background(prep):
    """Reproduce the already frozen protein topology, without energy evaluation."""
    import openmm as mm
    from openmm import app,unit
    from mace_global_prepare import STANDARD,atom_id
    source=verify(prep['source']);ff=verify(prep['forcefield']);verify(prep['water_forcefield'])
    if prep['preparation_details']['terminal_additions']:raise InvalidArtifact('unimplemented added terminal topology')
    pdb=app.PDBFile(str(source));model=app.Modeller(pdb.topology,pdb.positions)
    model.delete([a for a in model.topology.atoms() if a.residue.chain.id!='A' or a.residue.name not in STANDARD])
    atoms=list(model.topology.atoms());coords=np.array(model.positions.value_in_unit(unit.angstrom))
    original=prep['preparation_details']['original_protein_atoms'];physical={a['id']:a for a in prep['physical_atoms']}
    ids=[atom_id(a) for a in atoms]
    if ids!=[a['id'] for a in original] or len(set(ids))!=len(ids):raise InvalidArtifact('ff protein source inventory differs')
    for a,old in zip(atoms,original):
        if a.element.symbol!=old['element'] or not np.allclose(coords[a.index],old['xyz_A'],atol=1e-12,rtol=0):raise InvalidArtifact('ff source geometry differs')
        if a.element.symbol!='H' and old['xyz_A']!=physical[atom_id(a)]['xyz_A']:raise InvalidArtifact('normalized heavy geometry differs')
    bonds=sorted(tuple(sorted((atom_id(a),atom_id(b)))) for a,b in model.topology.bonds())
    expected=sorted(tuple(sorted((b['atom_a_id'],b['atom_b_id']))) for b in prep['preparation_details']['bonds'])
    if bonds!=expected:raise InvalidArtifact('protein connectivity/disulfides differ')
    system=app.ForceField(str(ff)).createSystem(model.topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
    nb=next(f for f in system.getForces() if isinstance(f,mm.NonbondedForce))
    charges={atom_id(a):float(nb.getParticleParameters(a.index)[0].value_in_unit(unit.elementary_charge)) for a in atoms}
    if abs(math.fsum(charges.values())-prep['protein_charge_e'])>TOL['boundary_charge_e']:raise InvalidArtifact('forcefield protein charge differs')
    if set(charges)!={a['id'] for a in prep['physical_atoms'] if a['kind']=='protein_source'}:
        raise InvalidArtifact('unaccounted protein/cofactor atoms')
    return {'charges_e':charges,'bonds':[list(b) for b in bonds],'forcefield':prep['forcefield'],'source':prep['source']}


def ligand_ledger(case,physical,source_ids):
    mechanics=read_json(verify(case['source_mapping']));repair=read_json(verify(mechanics['source_preparation']))
    ledger={};seen=set()
    for item in repair['charge_ledger']:
        s=item['source'];key=f"{s['chain']}/{s['resnum']}/{s['insertion_code']}";q=item['formal_charge']
        if q:
            # Current declared data contain only complete anionic Asp/Glu groups.
            required={'ASP':{'CG','OD1','OD2'},'GLU':{'CD','OE1','OE2'}}.get(s['canonical_resname'])
            if q!=-1 or item['kind']!='sidechain' or required is None or key in seen:
                raise InvalidArtifact('unsupported or duplicate charged functional group')
            if not {key+'/'+a for a in required}<=source_ids:raise InvalidArtifact('incomplete QM charged functional group')
            seen.add(key)
        ledger[key]=ledger.get(key,0)+q
    if sum(ledger.values())!=repair['ligand_formal_charge'] or repair['ligand_formal_charge']!=-3:
        raise InvalidArtifact('selected ligand formal-charge ledger differs')
    return ledger,mechanics['source_preparation']


def assemble(case,charge_rows,background):
    prep=read_json(verify(case['normalized_global_preparation']));physical=prep['physical_atoms'];ids=[a['id'] for a in physical];lookup={k:i for i,k in enumerate(ids)}
    if len(lookup)!=len(ids) or prep['cofactor_charge_e']!=0:raise InvalidArtifact('unsupported full physical inventory/cofactor')
    first=charge_rows['Ca'];proj=read_json(verify(first['projection']));support=set(proj['physical_ids'])
    source_ids={a['physical_id'] for a in case['mapping'] if a['kind']=='source'}
    ledger,ledger_pin=ligand_ledger(case,physical,source_ids)
    if not support<=lookup.keys():raise InvalidArtifact('projected QM atom absent from full cavity')
    if not {a['id'] for a in physical if a['kind']!='protein_source'}<=support:raise InvalidArtifact('nonprotein atom lacks QM representation')
    if any(a['element'] not in (*RADII,'M') or a['radius_A']!=RADII['Ca' if a['element']=='M' else a['element']] for a in physical):
        raise InvalidArtifact('physical radius inventory differs')
    ff=background['charges_e'];env={k:(0. if k in support else q) for k,q in ff.items()};ledgers=[]
    affected=sorted({rid(k) for k in support if k in ff})
    for key in affected:
        local=sorted(k for k in ff if rid(k)==key);outside=[k for k in local if k not in support]
        formal=round(math.fsum(ff[k] for k in local))
        if abs(math.fsum(ff[k] for k in local)-formal)>TOL['boundary_charge_e']:raise InvalidArtifact('noninteger forcefield residue charge')
        target=formal-ledger.get(key,0);before=math.fsum(env[k] for k in outside);delta=target-before
        recipients=set();bonds=[]
        for a,b in background['bonds']:
            if a in support and b in outside:recipients.add(b);bonds.append([a,b])
            elif b in support and a in outside:recipients.add(a);bonds.append([b,a])
        recipients=sorted(recipients)
        if not recipients and abs(delta)>TOL['boundary_charge_e']:raise InvalidArtifact('local boundary cannot close: '+key)
        increment=delta/len(recipients) if recipients else 0.
        for k in recipients:env[k]+=increment
        if abs(math.fsum(env[k] for k in local)-target)>TOL['boundary_charge_e']:raise InvalidArtifact('residue exterior charge closure failed')
        ledgers.append({'residue':key,'forcefield_formal_charge_e':formal,'selected_fragment_formal_charge_e':ledger.get(key,0),
            'projection_support_ids':sorted(set(local)&support),'exterior_target_e':target,'original_exterior_charge_e':before,
            'delta_e':delta,'recipient_ids':recipients,'recipient_bonds':bonds,'increment_per_recipient_e':increment,
            'original_recipient_charges_e':[ff[k] for k in recipients],'new_recipient_charges_e':[env[k] for k in recipients]})
    qenv=np.array([env.get(k,0.) for k in ids]);expected_env=prep['protein_charge_e']+3
    if any(qenv[lookup[k]]!=0 for k in support) or abs(float(qenv.sum())-expected_env)>TOL['boundary_charge_e']:
        raise InvalidArtifact('environment closure or QM overlap failed')
    ends={}
    for metal in ('Ca','La'):
        row=charge_rows[metal];p=read_json(verify(row['projection']));qm=np.zeros(len(ids))
        if p!=proj:raise InvalidArtifact('paired charge projection differs')
        fit=np.array(row['charge_e']);projected=np.array(p['weights'])@fit
        if not np.allclose(projected,row['projected_charge_e'],atol=1e-12,rtol=0):raise InvalidArtifact('actual projected charges differ')
        # Preserve the qualified serialized distribution exactly across BLAS hosts.
        projected=np.array(row['projected_charge_e'])
        for key,q in zip(p['physical_ids'],projected):qm[lookup[key]]+=q
        core_charge=-1 if metal=='Ca' else 0;total=qm+qenv;formal=prep['endpoints'][metal]['charge']
        if abs(float(qm.sum())-core_charge)>TOL['charge_e'] or abs(float(total.sum())-formal)>TOL['charge_e']:
            raise InvalidArtifact('native printed charge closure failed')
        expected=[(metal if a['element']=='M' else a['element'],*a['xyz_A']) for a in physical]
        if xyz(verify(prep['endpoints'][metal]['xyz']))!=expected:raise InvalidArtifact('full exact paired coordinates differ')
        ends[metal]={'QM_charges_e':qm.tolist(),'full_charges_e':total.tolist(),'charge':formal,
                     'QM_charge':core_charge,'xyz':prep['endpoints'][metal]['xyz'],'actual_charge_receipt':row['execution_receipt']}
    return {'physical_atoms':physical,'projection_support_ids':sorted(support),'environment_charges_e':qenv.tolist(),
       'environment_charge':expected_env,'boundary_ledger':ledgers,'ligand_ledger':ledger,'ligand_ledger_source':ledger_pin,
       'forcefield_background':background,'endpoints':ends,'source_case':case['source_mapping'],
       'normalized_preparation':case['normalized_global_preparation'],'assembly':prep['assembly'],'microstate':prep['microstate'],
       'explicit_waters':prep['explicit_waters'],'evidence':prep['evidence'],'evidence_use':'consumed_method_development'}


@cached_file_checks
def source_data(charges,hybrid):
    from mace_omol_charges import validate as validate_charges,accepted as accepted_charge
    if record(charges)['sha256']!=CHARGE_RESULT_SHA:raise InvalidArtifact('declared actual charge report differs')
    r=read_json(charges);mp=verify(r['manifest'])
    if r['manifest']['sha256']!=CHARGE_MANIFEST_SHA or not(r['fit_quality_gate_pass'] and r['projection_gate_pass']):raise InvalidArtifact('charge qualification unavailable')
    validate_charges(mp);cm=read_json(mp);p=read_json(verify(cm['preparation']));h=read_json(hybrid)
    if h['preparation']!=cm['preparation'] or h['status']!='complete' or not h['numerical_gate_pass'] or h['protocol_id']!=MODEL['vacuum_parent']:
        raise InvalidArtifact('matched vacuum parent differs')
    for t in cm['tasks']:
        row=r['rows'][t['task_id']];receipt=verify(row['execution_receipt'])
        if accepted_charge(receipt.parent,t,mp) is None:raise InvalidArtifact('actual charge utility receipt unavailable')
    systems={};backgrounds={}
    for name in CASES:
        case=read_json(verify(p['cases'][name]));gp=verify(case['normalized_global_preparation']);key=record(gp)['sha256']
        if key not in backgrounds:backgrounds[key]=forcefield_background(read_json(gp))
        systems[name]=assemble(case,{metal:r['rows'][name+'_'+metal] for metal in ('Ca','La')},backgrounds[key])
    if systems['GGR_extended']['physical_atoms']!=systems['GGR_connected']['physical_atoms']:
        raise InvalidArtifact('GGR physical boundary changes with partition')
    return systems,cm,h


def specs():
    out=[]
    for name in CASES:
        for component in ('full','QM'):
            for metal in ('Ca','La'):out.append((name,metal,component,'primary','CUDA'))
        out.append((name,'Ca','environment','primary','CUDA'))
        for variant in ('identity','rotate','translate'):
            for metal in ('Ca','La'):out.append((name,metal,'full',variant,'CUDA'))
    for variant,platform in [('repeat','CUDA'),('reference','Reference')]:
        for metal in ('Ca','La'):out.append(('GGR_connected',metal,'full',variant,platform))
    assert len(out)==48
    return out


def task_payload(spec,system,state_pin,xyz_pin,charge_pin):
    name,metal,component,variant,platform=spec
    charge=system['environment_charge'] if component=='environment' else system['endpoints'][metal]['charge' if component=='full' else 'QM_charge']
    return {'task_id':'_'.join(spec[:4]),'case_id':name,'metal':metal,'charge_component':component,'variant':variant,'platform':platform,
        'kind':'full','energy_component':COMPONENT,'xyz':xyz_pin,'charges':charge_pin,'charge':charge,'physical_state':state_pin,
        'solvent_dielectric':1. if variant=='identity' else MODEL['solvent_dielectric'],'solver':'native',
        'state_semantics':'component charge is a GB source, not a separate electronic/protonation state',
        'force_definition':'fixed-charge reaction-field gradient only; not combined model gradient'}


def geometry(system,metal,variant):
    atoms=xyz(verify(system['endpoints'][metal]['xyz']));coords=np.array([a[1:] for a in atoms]);center=coords[next(i for i,a in enumerate(system['physical_atoms']) if a['id']=='metal')]
    if variant=='rotate':coords=(coords-center)@rotation().T+center
    elif variant=='translate':coords+=np.array([10.,-7.,3.])
    return [(a[0],*x) for a,x in zip(atoms,coords)]


def prepare(charges,hybrid,software,agreement,output):
    systems,parent,h=source_data(charges,hybrid);sm=read_json(software)
    if sm['openmm_version']!='8.5.1':raise InvalidArtifact('unsupported solver version')
    for k in ('python','requirements','backend_source_inventory'):verify(sm[k])
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir();pins={}
    sources={name:verify(pin) for name,pin in parent['implementation'].items()}
    for name in ('mace_hybrid.py','mace_omol_solvent.py','mace_gb.py','mace_global_prepare.py','mace_file_checks.py'):sources[name]=Path(__file__).with_name(name)
    for name,src in sources.items():shutil.copyfile(src,impl/name);pins[name]=record(impl/name)
    stored={};qfiles={}
    for name,system in systems.items():
        d=out/'states'/name;d.mkdir(parents=True);write_new(d/'state.json',system);stored[name]=record(d/'state.json')
        for metal in ('Ca','La'):
            for component in ('full','QM'):
                f=d/(metal+'_'+component+'.npy');np.save(f,np.array(system['endpoints'][metal][component+'_charges_e']));qfiles[(name,metal,component)]=record(f)
        f=d/'environment.npy';np.save(f,np.array(system['environment_charges_e']));qfiles[(name,'Ca','environment')]=record(f)
    tasks=[];gd=out/'coordinates';gd.mkdir()
    for spec in specs():
        name,metal,component,variant,_=spec;f=gd/('_'.join(spec[:4])+'.xyz');write_xyz(f,geometry(systems[name],metal,variant))
        task=task_payload(spec,systems[name],stored[name],record(f),qfiles[(name,metal,component)])
        task['cache_key']=cache_key({'task':task,'model':MODEL,'software':record(software),'implementation':pins});tasks.append(task)
    m={'schema_version':SCHEMA,'protocol_id':PROTOCOL,'model':MODEL,'tolerances':TOL,'software':record(software),'agreement':record(agreement),
       'charges':record(charges),'vacuum_hybrid':record(hybrid),'states':stored,'tasks':tasks,'implementation':pins,
       'compute_budget':None,'wall_time_limit':None,'new_GB_calls':48,'new_DFT_calls':0,'new_MACE_calls':0,**NULLS}
    write_new(out/'manifest.json',m);return validate(out/'manifest.json')


@cached_file_checks
def validate(manifest):
    m=read_json(manifest)
    if m['schema_version']!=SCHEMA or m['protocol_id']!=PROTOCOL or m['model']!=MODEL or m['tolerances']!=TOL or len(m['tasks'])!=48:
        raise InvalidArtifact('solvent scientific settings/inventory differ')
    for pin in [m['agreement'],m['software'],*m['implementation'].values()]:verify(pin)
    sm=read_json(verify(m['software']))
    for pin in [sm['python'],sm['requirements'],sm['backend_source_inventory'],*read_json(verify(sm['backend_source_inventory']))['files']]:verify(pin)
    systems,_,_=source_data(verify(m['charges']),verify(m['vacuum_hybrid']))
    for name,pin in m['states'].items():
        if read_json(verify(pin))!=systems[name]:raise InvalidArtifact('full charge state/boundary ledger changed')
    for spec,t in zip(specs(),m['tasks']):
        name,metal,component,variant,_=spec;s=systems[name];atoms=xyz(verify(t['xyz']));q=np.load(verify(t['charges']))
        wanted_q=s['environment_charges_e'] if component=='environment' else s['endpoints'][metal][component+'_charges_e']
        wanted=geometry(s,metal,variant)
        if [a[0] for a in atoms]!=[a[0] for a in wanted] or not np.allclose([a[1:] for a in atoms],[a[1:] for a in wanted],atol=1e-12,rtol=0) or not np.array_equal(q,wanted_q):
            raise InvalidArtifact('solver geometry/charges differ')
        payload=task_payload(spec,s,m['states'][name],t['xyz'],t['charges'])
        if {k:v for k,v in t.items() if k!='cache_key'}!=payload or t['cache_key']!=cache_key({'task':payload,'model':MODEL,'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('solver state/cache changed')
    return {'status':'pass','tasks':48,'manifest':record(manifest),'new_DFT_calls':0,'new_MACE_calls':0}


def worker(manifest,task_id,output,memory_mode):
    import openmm as mm
    from openmm import unit
    m=read_json(manifest);t=next(t for t in m['tasks'] if t['task_id']==task_id);out=Path(output);start=time.monotonic()
    r={'status':'failed','task_id':task_id,'cache_key':t['cache_key'],'energy_component':COMPONENT,'manifest':record(manifest),
       'memory_mode':memory_mode,'peak_cuda_allocated_bytes':None,'peak_cuda_reserved_bytes':None,**NULLS}
    try:
        if not os.environ.get('SLURM_JOB_ID') or mm.__version__!='8.5.1':raise InvalidArtifact('approved allocated solver required')
        atoms=xyz(verify(t['xyz']));q=np.load(verify(t['charges']))
        if q.shape!=(len(atoms),) or not np.isfinite(q).all() or abs(float(q.sum())-t['charge'])>TOL['charge_e']:raise InvalidArtifact('solver charge inventory/closure invalid')
        system=build_system(atoms,q,'native',t['solvent_dielectric']);sp=out/'system.xml';sp.write_text(mm.XmlSerializer.serialize(system))
        integrator=mm.VerletIntegrator(.001);platform=mm.Platform.getPlatformByName(t['platform'])
        props={'Precision':'double','DeterministicForces':'true'} if t['platform']=='CUDA' else {}
        tick=time.monotonic();context=mm.Context(system,integrator,platform,props);r['context_seconds']=time.monotonic()-tick
        r['platform_properties']={k:platform.getPropertyValue(context,k) for k in platform.getPropertyNames()}
        context.setPositions(np.array([a[1:] for a in atoms])*.1*unit.nanometer);tick=time.monotonic()
        state=context.getState(getEnergy=True,getForces=True,groups={0});kj=float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        force=np.array(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer))/(EV_TO_KJ*10.)
        r['evaluation_seconds']=time.monotonic()-tick
        if not math.isfinite(kj) or force.shape!=(len(atoms),3) or not np.isfinite(force).all():raise InvalidArtifact('nonfinite solver result')
        fp=out/'fixed_charge_forces_eV_A.npy';np.save(fp,force)
        r.update(status='computed',energy_eV=kj/EV_TO_KJ,GB_reaction_kcal_mol=kj/4.184,GB_reaction_kJ_mol=kj,
            forces=record(fp),charges=t['charges'],charge_check=True,charge_sum_e=float(q.sum()),charge_error_e=float(q.sum())-t['charge'],
            serialized_system=record(sp),extracted_force_groups=[0],direct_Coulomb_included=False,
            force_definition=t['force_definition'],platform=t['platform'])
        del context,integrator
    except Exception as exc:r.update(reason=str(exc),exception_type=type(exc).__name__);traceback.print_exc()
    finally:
        r.update(wall_seconds=time.monotonic()-start,peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            slurm_job_id=os.environ.get('SLURM_JOB_ID'),allocated_cpus=os.environ.get('SLURM_CPUS_PER_TASK'),
            allocated_host_mem_MiB=os.environ.get('SLURM_MEM_PER_NODE'),versions={'openmm':mm.__version__,'numpy':np.__version__})
        write_new(out/'result.json',r)
    print(json.dumps({k:r.get(k) for k in ('task_id','status','reason','evaluation_seconds')}),flush=True);return r


@cached_file_checks
def collect(manifest):
    validate(manifest);m=read_json(manifest);mp=Path(manifest).resolve();rows={};attempts=[]
    for t in m['tasks']:
        good=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp);attempts.append({'task_id':t['task_id'],'directory':str(a),'accepted':r is not None,
                'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
            if r:good.append(r)
        rows[t['task_id']]=good[-1] if good else {'status':'unavailable','energy_eV':None}
    complete=all(r['status']=='computed' for r in rows.values());checks=[];cases={};contrasts=[];partition=None
    def energy(name,metal,component='full',variant='primary'):return rows['_'.join((name,metal,component,variant))]['GB_reaction_kcal_mol']
    def check(name,error,tol):checks.append({'name':name,'error_kcal_scale':error,'tolerance':tol,'pass':abs(error)<=tol})
    if complete:
        h=read_json(verify(m['vacuum_hybrid']))
        for name in CASES:
            for metal in ('Ca','La'):
                base=energy(name,metal)
                check(name+'_'+metal+'_identity',energy(name,metal,variant='identity'),TOL['identity_kcal_mol'])
                for variant in ('rotate','translate')+ (('repeat','reference') if name=='GGR_connected' else ()):
                    check(name+'_'+metal+'_'+variant,energy(name,metal,variant=variant)-base,TOL['numerical_kcal_mol'])
            baseR=energy(name,'Ca')-energy(name,'La')
            for variant in ('rotate','translate')+(('repeat','reference') if name=='GGR_connected' else ()):
                check(name+'_'+variant+'_contrast',energy(name,'Ca',variant=variant)-energy(name,'La',variant=variant)-baseR,TOL['numerical_kcal_mol'])
            genv=energy(name,'Ca','environment');ends={}
            for metal in ('Ca','La'):
                total=energy(name,metal);qm=energy(name,metal,'QM');cross=total-qm-genv
                ends[metal]={'GB_total_kcal_mol':total,'GB_QM_self_and_internal_kcal_mol':qm,'GB_env_kcal_mol':genv,
                    'GB_cross_kcal_mol':cross,'vacuum_hybrid_kcal_scale':h['cases'][name]['endpoints'][metal]['hybrid_kcal_scale'],
                    'solvent_hybrid_kcal_scale':h['cases'][name]['endpoints'][metal]['hybrid_kcal_scale']+total}
            qdiff=ends['Ca']['GB_QM_self_and_internal_kcal_mol']-ends['La']['GB_QM_self_and_internal_kcal_mol'];xdiff=ends['Ca']['GB_cross_kcal_mol']-ends['La']['GB_cross_kcal_mol']
            corrected=h['cases'][name]['hybrid_R_kcal_scale']+baseR
            check(name+'_component_closure',baseR-qdiff-xdiff,TOL['algebra_kcal_scale'])
            check(name+'_direct_algebra',corrected-(ends['Ca']['solvent_hybrid_kcal_scale']-ends['La']['solvent_hybrid_kcal_scale']),TOL['algebra_kcal_scale'])
            cases[name]={'endpoints':ends,'vacuum_R_kcal_scale':h['cases'][name]['hybrid_R_kcal_scale'],
                'GB_R_kcal_mol':baseR,'GB_QM_R_kcal_mol':qdiff,'GB_cross_R_kcal_mol':xdiff,
                'solvent_hybrid_R_kcal_scale':corrected,'evidence':read_json(verify(m['states'][name]))['evidence']}
        partition=cases['GGR_connected']['solvent_hybrid_R_kcal_scale']-cases['GGR_extended']['solvent_hybrid_R_kcal_scale']
        for a in ('ALPHA_1F6S','ALPHA_6IP9'):
            for b in ('GGR_extended','GGR_connected'):
                delta=cases[a]['solvent_hybrid_R_kcal_scale']-cases[b]['solvent_hybrid_R_kcal_scale']
                contrasts.append({'alpha':a,'GGR':b,'difference_kcal_scale':delta,'pass':delta>TOL['ordering_kcal_scale']})
    numeric=complete and bool(checks) and all(c['pass'] for c in checks)
    return {'protocol_id':PROTOCOL,'manifest':record(mp),'status':'complete' if complete else 'incomplete','rows':rows,'attempts':attempts,
        'checks':checks,'cases':cases,'contrasts':contrasts,'numerical_gate_pass':numeric,'partition_shift_kcal_scale':partition,
        'partition_gate_pass':partition is not None and abs(partition)<=TOL['partition_kcal_scale'],
        'ordering_gate_pass':numeric and len(contrasts)==4 and all(c['pass'] for c in contrasts),
        'tolerances':TOL,'collection_implementation':record(__file__),**NULLS}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare')
    for k in ('charges','hybrid','software','agreement','output'):a.add_argument('--'+k,required=True)
    for op in ('dry-run','collect'):
        a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
    a=p.parse_args()
    if a.op=='prepare':r=prepare(a.charges,a.hybrid,a.software,a.agreement,a.output)
    elif a.op=='dry-run':r=validate(a.manifest)
    else:r=collect(a.manifest)
    if a.op!='prepare' and a.output:write_new(a.output,r)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','attempts')},indent=2))
