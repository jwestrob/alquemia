#!/usr/bin/env python3
"""Prepare the two selected XoxF sites; audit H before any endpoint admission."""
import copy
import json
from pathlib import Path
import shutil
import time
import numpy as np
import openmm as mm
from openmm import app,unit
import gemmi
import select_inputs as inputs
from geometry_checks import check,unchanged_non_H

A=inputs.A;OUT=inputs.OUT;HERE=Path(__file__).resolve().parent
helper=inputs.load('xoxf_H_diagnostic_helpers',A/'diagnostics/plm_adh9_af3_20260916/hydrogen_repair/repair_hydrogens.py')
record=inputs.record;read=inputs.read;write=inputs.write

def capture_original(w,target,pins,original,cp):
    cid=target['case_id'];audit=OUT/'hydrogen_audit'/cid;audit.mkdir(parents=True,exist_ok=False)
    original_minimizer=mm.LocalEnergyMinimizer.minimize;calls=[]
    def observe(context,tolerance=10,maxIterations=0,reporter=None):
        system=context.getSystem();n=len(calls);mask=np.asarray([system.getParticleMass(i).value_in_unit(unit.dalton)>0 for i in range(system.getNumParticles())])
        xml=audit/f'system_{n}.xml';xml.write_text(mm.XmlSerializer.serialize(system))
        before,bm=helper.state(context,mask);np.save(audit/f'initial_positions_{n}_nm.npy',before)
        trace=helper.Reporter();t=time.monotonic();original_minimizer(context,tolerance,maxIterations,trace)
        after,am=helper.state(context,mask);np.save(audit/f'final_positions_{n}_nm.npy',after)
        if not np.array_equal(before[~mask],after[~mask]):raise ValueError('Original minimizer moved heavy atoms')
        calls.append({'system':record(xml),'tolerance':tolerance,'max_iterations':maxIterations,'initial':bm,'final':am,
          'iterations_reported':len(trace.entries),'elapsed_seconds':time.monotonic()-t,'trace':trace.entries,'heavy_displacement_A':0.0})
    mm.LocalEnergyMinimizer.minimize=observe
    try:carve_record=w.prepare_one(target,OUT/'original_preparation',pins,original,cp)
    finally:
        mm.LocalEnergyMinimizer.minimize=original_minimizer
        write(audit/'minimization_observation.json',{'observational_only':True,'original_protonator_unchanged':True,
          'script':record(__file__),'calls':calls})
    if len(calls)!=1:raise ValueError('Expected exactly one standard-protonation minimizer call')
    return carve_record,audit

def identity_coordinates(path):
    s=gemmi.read_structure(str(path));result={}
    for c in s[0]:
        for r in c:
            for a in r:
                key=(c.name,r.name,r.seqid.num,r.seqid.icode,a.name,a.element.name)
                if key in result:raise ValueError('Duplicate atom identity')
                result[key]=(a.pos.x,a.pos.y,a.pos.z)
    return s,result

def recover(w,target,original_carve,audit):
    """Same captured objective, failed PDB restart, at most1000 numerical steps."""
    old=read(original_carve['path']);cid=target['case_id'];folder=OUT/'hydrogen_recovery'/cid;folder.mkdir(parents=True,exist_ok=False)
    source=Path(old['source_structure']['path']);pdb=app.PDBFile(str(source));sys=mm.XmlSerializer.deserialize((audit/'system_0.xml').read_text())
    h=np.asarray([a.element==app.element.hydrogen for a in pdb.topology.atoms()]);ctx=mm.Context(sys,mm.VerletIntegrator(0),mm.Platform.getPlatformByName('CPU'),{'Threads':'1'})
    if sys.getNumParticles()!=len(h):raise ValueError('Captured System atom count mismatch')
    if any((sys.getParticleMass(i).value_in_unit(unit.dalton)>0)!=v for i,v in enumerate(h)):raise ValueError('Captured H identities differ')
    ctx.setPositions(pdb.positions);before,bm=helper.state(ctx,h);reporter=helper.Reporter();start=time.monotonic()
    mm.LocalEnergyMinimizer.minimize(ctx,1.0,1000,reporter);after,am=helper.state(ctx,h)
    if not np.array_equal(before[~h],after[~h]):raise ValueError('Recovery moved heavy atoms')
    if am['hydrogen_force_rms_kJ_mol_nm']>1.0:raise ValueError('H recovery did not converge')
    np.save(folder/'positions_nm.npy',after);recovered=folder/'protonated_recovered.pdb'
    with recovered.open('x') as f:app.PDBFile.writeFile(pdb.topology,after*unit.nanometer,f,keepIds=True)
    serialized=app.PDBFile(str(recovered));ctx.setPositions(serialized.positions);rounded,rm=helper.state(ctx,h)
    result={'status':'PASS','source':record(source),'original_minimization':record(audit/'minimization_observation.json'),
       'system':record(audit/'system_0.xml'),'objective':'Exact captured original standard-only forcefield=None objective',
       'deviation':'Restart failed serialized PDB H; iteration ceiling50→1000; same objective/tolerance/CPU/thread and all heavy positions',
       'new_random_draws':0,'original_seed':20260914,'initial':bm,'final':am,'iterations_reported':len(reporter.entries),
       'elapsed_seconds':time.monotonic()-start,'heavy_displacement_A':0.0,'platform':'CPU','threads':1,
       'post_PDB_serialization':rm,'max_H_rounding_displacement_A':float(np.max(np.linalg.norm(rounded[h]-after[h],axis=1))*10),
       'convergence_evaluated_before_serialization':True,'output':record(recovered),'script':record(__file__),'trace':reporter.entries}
    write(folder/'minimization.json',result)
    structure,newcoords=identity_coordinates(recovered);_,oldcoords=identity_coordinates(source)
    if newcoords.keys()!=oldcoords.keys() or any(v!=newcoords[k] for k,v in oldcoords.items() if k[-1]!='H'):raise ValueError('Full protein heavy atoms or H states changed')
    fixed=w.fixed;keys,residues,partner=w.role_state(structure[0],target);site,pkey,pqq,prepared,contacts,*_=w.site_state(structure[0],target)
    fragments=[copy.deepcopy(old['qm_fragments'][0])];atoms=list(prepared.atoms)
    if fixed.pqq_atom_records(pqq,prepared)!=fragments[0]['atom_records']:raise ValueError('PQQ atom state changed')
    for role in fixed.ROLE_ORDER:
        if role not in old['fixed_core']['included_roles']:continue
        if role=='catalytic_asp_cationic_partner':aa,ff=fixed.cationic_sidechain_fragment(residues[role],keys[role],site.atom.pos)
        else:aa,ff=fixed.canonical_sidechain_fragment(residues[role],keys[role],role,site.atom.pos)
        atoms.extend(aa);fragments.append(ff)
    changes=unchanged_non_H(old['qm_fragments'],fragments);geometry=check(fragments)
    if geometry['status']!='PASS':raise ValueError('Recovered core fails geometry:'+str(geometry['errors']))
    new=copy.deepcopy(old);new['qm_fragments']=fragments;new['charge_ledger']['fragments']=fragments
    new['source_structure']=record(recovered)
    protonation=folder/'protonation_manifest.json';write(protonation,{'protocol_id':'standard_only_PDBFixer_with_same_objective_H_recovery',
       'original':old['protonation_manifest'],'source':old['normalized_source_structure'],'output':record(recovered),
       'minimization':record(folder/'minimization.json'),'new_protonation_states':False,'changed_numerical_preparation':True})
    heavy=folder/'heavy_coordinate_check.json';write(heavy,{'passes':True,'source':old['source_cif'],'protonated':record(recovered),
       'original_check':old['heavy_coordinate_check'],'full_protein_atom_identity_preserved':True,'heavy_displacement_A':0.0})
    new.update(protonation_manifest=record(protonation),heavy_coordinate_check={**record(heavy),'passes':True},
       hydrogen_recovery={'approval':record(inputs.APPROVAL),'original_failed_carve':original_carve,
       'minimization':record(folder/'minimization.json'),'changed_source_hydrogens':changes,
       'calibration_compatibility':'Same fixed core/state/heavy geometry/electronic method; changed numerical H convergence explicitly recorded; no calibration rerun/refit.'})
    outputs={};tasks=[]
    for metal in ('La','Ca'):
        xyz=folder/f'{cid}_{metal}_qm.xyz';inp=folder/f'sp_{cid}_{metal}.inp';charge=old['charge_ledger'][metal+'_total']
        aa=[(metal,site.atom.pos.x,site.atom.pos.y,site.atom.pos.z),*atoms];w.base._validate_singlet(metal,aa,charge)
        fixed.write_xyz(xyz,stem=cid,label=metal,atoms=aa,charge=charge);fixed.write_orca_input(inp,xyz_name=xyz.name,charge=charge,stem=cid,label=metal)
        outputs[metal+'_xyz']=record(xyz);outputs[metal+'_input']=record(inp)
        tasks.append({'task_id':metal,'xyz':record(xyz),'input':record(inp),'output_path':str(inp.with_suffix('.out'))})
    new.update(outputs=outputs,tasks=tasks)
    return new,geometry,record(folder/'minimization.json')

def main():
    w=inputs.wrapper();pins,original,cp=w.verify_pins(OUT/'candidate_implementation_pins.json');source=read(OUT/'candidate_manifest.json')
    if len(source['targets'])!=2 or set(t['target_id'] for t in source['targets'])!=set(inputs.EXPECTED):raise ValueError('Two-target scope changed')
    (OUT/'original_preparation').mkdir(exist_ok=False);final=OUT/'prepared';final.mkdir(exist_ok=False)
    cases=[];case_checks={}
    for target in source['targets']:
        cid=target['case_id'];case={k:target[k] for k in ('case_id','target_id','rank','source_cif')}
        try:
            raw,audit=capture_original(w,target,pins,original,cp);carve=read(raw['path']);geometry=check(carve['qm_fragments']);original_geometry=copy.deepcopy(geometry)
            write(audit/'original_core_geometry.json',original_geometry);recovery=None
            if geometry['status']!='PASS':carve,geometry,recovery=recover(w,target,raw,audit)
            else:
                # Original frozen50-step preparation is retained when its actual core passes.
                # Its observed residual force is recorded without pretending it met tolerance.
                for rec in carve['outputs'].values():rec['path']=str(Path(raw['path']).parent/rec['path'])
            details={'status':'PASS','target_id':target['target_id'],'gene_id':target['gene_id'],'selected_sample':target['rank'],
                'source_cif':target['source_cif'],'geometry':geometry,'original_geometry':original_geometry,
                'original_minimizer_observation':record(audit/'minimization_observation.json'),
                'hydrogen_recovery':recovery,'preparation_mode':'recovered_same_objective' if recovery else 'original_frozen_protonation',
                'calibration_implementation_changed':bool(recovery),'all_heavy_atoms_preserved':True,'pqq_microstate':'pqq_ox_3minus_v1'}
            checkpath=final/f'{cid}_geometry_validation.json';write(checkpath,details)
            carve.update(geometry_validation=record(checkpath),hydrogen_geometry_audit=details,
                metadata_packager=record(__file__),original_preparation_carve=raw)
            mp=final/f'{cid}_carve_manifest.json';write(mp,carve)
            case.update(status='ready_for_orca',reason=None,carve_manifest=record(mp));case_checks[cid]=details
        except Exception as exc:
            reason=f'{type(exc).__name__}:{exc}';case.update(status='unsupported',reason=reason,carve_manifest=None)
            case_checks[cid]={'status':'UNSUPPORTED','target_id':target['target_id'],'reason':reason}
        cases.append(case);print(cid,case['status'],case.get('reason'),flush=True)
    validation=final/'geometry_validation.json';write(validation,{'status':'PASS' if all(c['status']=='ready_for_orca' for c in cases) else 'PARTIAL',
      'target_count':2,'ready_count':sum(c['status']=='ready_for_orca' for c in cases),'cases':case_checks,'approval':record(inputs.APPROVAL)})
    implementation=final/'implementation_pins.json';write(implementation,{'protocol_id':w.PROTOCOL,'approval':record(inputs.APPROVAL),
       'candidate_pins':record(OUT/'candidate_implementation_pins.json'),'candidate_manifest':record(OUT/'candidate_manifest.json'),
       'selection_review':record(OUT/'selection_review.json'),'geometry_validation':record(validation),
       'scripts':{p.name:record(p) for p in HERE.glob('*.py')},'H_observer_helpers':record(helper.__file__),
       'original_calibration_implementation_pins':record(cp)})
    write(final/'prepared_pairs.json',{'schema_version':'plm.adh9.prepared_pairs.v1','protocol_id':w.PROTOCOL,
       'candidate_manifest':record(OUT/'candidate_manifest.json'),'implementation_pins':record(implementation),
       'geometry_validation':record(validation),'cases':cases,'prepared_pair_count':sum(c['status']=='ready_for_orca' for c in cases),
       'unsupported_count':sum(c['status']!='ready_for_orca' for c in cases),'orca_executed':False})
    print(final/'prepared_pairs.json',flush=True)

if __name__=='__main__':main()
