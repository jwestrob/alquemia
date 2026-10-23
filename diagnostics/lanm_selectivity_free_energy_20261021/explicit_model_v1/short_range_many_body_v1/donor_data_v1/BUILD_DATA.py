#!/usr/bin/env python3
"""Materialize admitted saved donor data. No native, quantum or fitting interface."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shlex
import time

HERE = Path(__file__).resolve().parent
PINS = {}
CHECKS = []

def check(ok, name):
    if not ok:
        raise ValueError(name)
    CHECKS.append(name)

def pin(path):
    p = Path(path).resolve()
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def verify(p):
    if p['path'] not in PINS:
        check(pin(p['path']) == p, 'pin ' + p['path'])
        PINS[p['path']] = p
    else:
        check(PINS[p['path']] == p, 'consistent pin ' + p['path'])
    return Path(p['path'])

def read(p):
    return json.loads(verify(p).read_text())

def save(p, x):
    p.parent.mkdir(parents=True, exist_ok=True)
    check(not p.exists(), 'fresh output ' + str(p))
    p.write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')
    return pin(p)

def load(p):
    spec = importlib.util.spec_from_file_location('saved_parameter_resolver', verify(p))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

def key_records(text):
    return [shlex.split(line, comments=True) for line in text.splitlines()
            if shlex.split(line, comments=True)]

def nonpair_body(text):
    # Include all multipole continuation lines, bonded terms and atom records.
    excluded = {'parameters', 'vdwpr', 'polpair', 'digits', 'openmp-threads',
                'polar-eps', 'response-matrix-export', 'polar-response-readback',
                'polar-print', 'save-uinduce'}
    return [r for r in key_records(text) if r[0].lower() not in excluded]

def typed_geometry(inputs, kind, resolver):
    xyz = verify(inputs['xyz']).read_text().splitlines()
    rows = [s.split() for s in xyz[1:] if s.strip()]
    n = int(xyz[0].split()[0])
    check(len(rows) == n and [int(r[0]) for r in rows] == list(range(1,n+1)), 'typed atom order')
    prm = verify(inputs['parameters']).read_text()
    key = verify(inputs['key']).read_text()
    check(not any(r[0].lower().startswith(('polar-local','polar-donor','polar-radial','exchpol')) for r in key_records(key)), 'source has no added response')
    parameters = resolver.parse_parameters(prm,key)
    types = [int(r[5]) for r in rows]
    z = [parameters['atom'][t]['atomic_number'] for t in types]
    bonds = [[int(v)-1 for v in r[6:]] for r in rows]
    check(types[0] in (500,501) and not bonds[0], 'isolated metal first')
    for i, bs in enumerate(bonds):
        check(all(0 <= j < n and i in bonds[j] for j in bs), 'symmetric covalent topology')
    unseen = set(range(1,n)); ligands = []
    while unseen:
        seed=min(unseen); todo=[seed]; component=set()
        while todo:
            i=todo.pop()
            if i in component: continue
            component.add(i); todo.extend(bonds[i])
        check(0 not in component, 'intact ligand excludes metal')
        unseen -= component
        inds=sorted(component); elems=sorted(z[i] for i in inds)
        water=elems==[1,1,8]
        ligand_kind='water' if water else kind
        check(water or ligand_kind in ('acetate','amide'), 'known intact ligand kind')
        oxy=[i for i in inds if z[i]==8]
        frames={}
        for o in oxy:
            frames[str(o)]=sorted(j for j in bonds[o] if j in component)
            check([z[j] for j in frames[str(o)]] == ([1,1] if water else [6]), 'covalent oxygen frame')
        check(water or len(oxy)==(2 if kind=='acetate' else 1), 'all donor oxygen atoms retained')
        ligands.append(dict(id='ligand_'+str(len(ligands)),kind=ligand_kind,
                            atom_indices=inds,donor_atoms=oxy,frame_atoms_by_donor=frames))
    check(sorted(i for l in ligands for i in l['atom_indices'])==list(range(1,n)), 'exact fragment partition')
    return dict(xyz_A=[[float(v)for v in r[2:5]]for r in rows],atomic_numbers=z,
                atom_ids=list(range(1,n+1)),metal_atom=0,ligands=ligands,
                tinker_types=types,bonds_zero_based=bonds,index_convention='zero-based arrays; atom_ids native one-based'),parameters,key

def compatibility(parameters,key,geometry,shared,shared_key,prm_sha):
    types=sorted(set(geometry['tinker_types'])); mt=geometry['tinker_types'][0]
    body_same=nonpair_body(key)==nonpair_body(shared_key)
    # Full common non-pair parsed settings plus literal remaining key body.
    common=all(parameters[k]==shared[k] for k in ('settings','atom','vdw','polarize','metal_multipoles'))
    differences=[]
    for category, domain in [('polpair',set(types)),('vdwpr',{parameters['atom'][t]['class']for t in types})]:
        for pair in sorted(set(parameters[category]) | set(shared[category])):
            if set(pair) <= domain and parameters[category].get(pair)!=shared[category].get(pair):
                differences.append(dict(category=category,pair=list(pair),saved=parameters[category].get(pair),shared=shared[category].get(pair)))
    return dict(compatible=common and body_same and not differences and prm_sha,
                same_parameter_file=prm_sha,common_nonpair_values_equal=common,
                complete_nonpair_key_body_equal=body_same,applicable_pair_differences=differences,
                present_types=types,metal_type=mt,old_h_H_zero=True,
                rule='Same full key body and common parameters; only pair records whose two types/classes occur can affect this one-metal fixture.')

def operator_fixture(rec,geometry,units):
    check(rec['status']=='complete' and rec['numerical_checks_passed'] and rec['operator']['passed'], 'admitted absent operator')
    check(rec['response']['active'] is False and rec['job']['profile']=='absent', 'zero added h/H')
    a,b=rec['operator']['exports']
    check(a['A']==b['A'] and a['sites']==b['sites'], 'two native operator identities')
    check([v[:2]for v in a['vectors']]==[v[:2]for v in b['vectors']], 'two raw dual field identities')
    for p in rec['exports']:verify(p)
    for s in a['sites']:
        i=s['atom']-1
        check(s['xyz_A']==geometry['xyz_A'][i] and s['type']==geometry['tinker_types'][i], 'native site geometry/type map')
    n=len(a['sites']);dim=3*n
    check(len(a['A'])==dim and all(len(r)==dim for r in a['A']), 'operator shape')
    check(len(a['vectors'])==dim and all(len(v)==4 for v in a['vectors']), 'dual vector shape')
    return dict(site_atom_indices=[s['atom']-1 for s in a['sites']],alpha_A3=[s['alpha_A3']for s in a['sites']],
                A0_Aminus3=a['A'],Ed_e_per_A2=[v[0]for v in a['vectors']],Ep_e_per_A2=[v[1]for v in a['vectors']],
                native_u_eA=[v[2]for v in a['vectors']],native_v_eA=[v[3]for v in a['vectors']],
                K_kcal_A_per_mol_e2=units['electric'],permanent_mu_eA=rec['response']['permanent_lab_eA'][0],
                nonpolar_energy_kcal=rec['value']['energy']-rec['value']['components']['Polarization']['energy_kcal_mol'],
                nonpolar_energy_definition='Saved total minus printed native polarization; printed-energy precision retained.',
                charge_e=rec['moments']['charge_e'],origin_A=[0.,0.,0.],raw_exports=rec['exports'],
                coordinate_derivatives=None,full_force_provider_available=False,
                scope='Admitted finite vacuum operator; both native channels retained. No total force from frozen arrays.')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True);args=ap.parse_args()
    t0=time.monotonic();c0=time.process_time();out=Path(args.out).resolve()
    check(not out.exists(),'fresh materialization directory');out.mkdir(parents=True)
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
        check(os.environ.get(name)=='1','one thread '+name)
    sources=read(pin(HERE/'INPUTS.json'))
    # Never import a molecular runner or the optimizer.
    resolver=load(sources['resolver']);wb=read(sources['water_baseline']);wf=read(sources['water_fixture'])
    review=read(sources['donor_review']);training=read(sources['donor_training'])
    dm=read(sources['donor_derivative_manifest']);qm=read(sources['qualification_manifest'])
    qr=read(sources['qualification']);qa=read(sources['qualification_review'])
    sens=read(sources['sensitivity']);ds=read(sources['design_sources'])
    inv=read(ds['donor_inventory']);rp=read(ds['donor_authorities']['reference_panel'])
    field=read(ds['La_quantum_field_collection']);read(ds['La_quantum_field_actual_review'])
    mech=read(sources['mechanical_comparison']);cc=read(sources['canonical_rotation'])
    check(review['status']=='passed' and review['frozen_water_unchanged'],'selected donor admission')
    check(qa['qualification_passed'] and qr['numerical_checks_passed'],'operator admission')
    check(field['numerical_checks_passed'] and field['complete']==5,'qualified La field mode')
    vals=dm['selected_parameters'];check(len(vals)==8,'eight donor parameters')
    prm=verify(wf['inputs']['parameters']).read_text();water_key=verify(wf['inputs']['key']).read_text()
    water=wb['fixed_six_parameters'];extra=[]
    for metal,mt,k in [('La',500,0),('Dy',501,2)]:
        extra += [f'vdwpr 34 {mt} {water[k]:.17g} {water[4+(mt==501)]:.17g}',f'polpair 36 {mt} {water[k+1]:.17g}']
    for metal,mt,start in [('La',500,0),('Dy',501,4)]:
        for family,cl,typ,offset,eps in [('acetate',403,403,0,.228 if mt==500 else .3795688850865599),('amide',55,105,2,.3019 if mt==500 else .3795688850865599)]:
            extra += [f'vdwpr {cl} {mt} {vals[start+offset]:.17g} {eps:.17g}',f'polpair {typ} {mt} {vals[start+offset+1]:.17g}']
    shared_key=water_key.rstrip()+'\n\n# Explicit common-water plus donor17621 baseline; no added response.\n'+'\n'.join(extra)+'\n'
    check(not any(r[0].lower().startswith(('polar-local','polar-donor','polar-radial','exchpol')) for r in key_records(shared_key)), 'no historical response correction')
    shared=resolver.parse_parameters(prm,shared_key)
    (out/'SHARED_BASELINE.key').write_text(shared_key)
    baseline=dict(id='common_water17676_plus_disjoint_donor17621_zero_h_H',water_authority=sources['water_baseline'],
                  donor_authority=sources['donor_review'],donor_parameters=vals,donor_parameter_order=['La_acetate_R','La_acetate_damping','La_amide_R','La_amide_damping','Dy_acetate_R','Dy_acetate_damping','Dy_amide_R','Dy_amide_damping'],
                  parameter_file=wf['inputs']['parameters'],key=pin(out/'SHARED_BASELINE.key'),explicit_appended_records=extra,
                  old_h_H_zero=True,not_fitted=True,scope='Explicit disjoint parameter union. Fixture compatibility is evaluated by actual chemical applicability, not full key hash equality.')
    baseline_pin=save(out/'BASELINE.json',baseline)
    fixtures={};compat={};missing=[]
    for gid,src in review['selected_sources'].items():
        if not ('acetate' in gid or 'acetamide' in gid):continue
        cell=read(src['authority']);inp=src['effective_input'];kind='acetate' if 'acetate' in gid else 'amide'
        geom,p,key=typed_geometry(inp,kind,resolver);tr=training['rows'][gid]
        check(geom['xyz_A']==tr['geometry']['positions_A'],'training exact geometry '+gid)
        check(cell['status']=='complete','selected cell complete')
        co=compatibility(p,key,geom,shared,shared_key,inp['parameters']['sha256']==wf['inputs']['parameters']['sha256']);compat[gid]=co
        q=read(tr['reference_result'])
        target=dict(method='ORCA ECP PBE0-D4, original saved protocol; consumed development',
                    energy_hartree=tr['reference']['energy_hartree'],energy_kcal_mol=tr['reference']['energy'],
                    forces_kcal_mol_A=tr['reference']['forces'],source=tr['reference_result'],
                    printed_scf_criteria=q.get('printed_scf_criteria'),native_admission=q['checks'],
                    derivative_scope='Only separately attached sampled derivative qualifications; full force accuracy not established.',
                    physical_state_qualified=False,method_accuracy_qualified=False)
        check(len(target['forces_kcal_mol_A'])==len(geom['xyz_A']),'reference force dimension')
        f=dict(schema='short_range_many_body_fixture_v1',id=gid,baseline_id=baseline['id'],baseline=baseline_pin,
               geometry=geom,inputs=inp,baseline_identity=co,compatible=co['compatible'],
               baseline_observables=dict(energy_kcal_mol=cell['energy'],forces_kcal_mol_A=cell['forces'],
                   components=cell['native'].get('components'),source=src['authority'],native_outputs=src['native_outputs']),
               historical_baseline_label='donor17621 with older water6 in hydrated states',
               reference_targets=target,native_operator=None,consumed_development=True,
               full_force_from_frozen_operator_supported=False,
               missing=['baseline operator/field coordinate derivatives','native operator export for this geometry'],
               fragment_cycle_residual=None,
               absolute_energy_warning='QM and classical absolute energies have different zeros; do not subtract them as a fragment interaction target.',
               original_fragment_metadata=dict(bodies=tr['bodies'],geometry=tr['geometry']))
        if gid.startswith('matched_dry_acetate_mono_central_'):
            metal=gid.rsplit('_',1)[1];row=qr['rows']['absent__'+metal+'_central__analyze'];g=qr['rows']['absent__'+metal+'_central__testgrad']
            check(inp==qm['rows'][metal+'_central']['input'],'operator exact input '+metal)
            check(row['value']['energy']==cell['energy'] and [[-v for v in r]for r in g['value']['gradient']]==cell['forces'],'saved operator energy/force join')
            f['native_operator']=operator_fixture(row,geom,qm['units']);f['native_qualification']=dict(result=sources['qualification'],review=sources['qualification_review'],row='absent__'+metal+'_central__analyze')
            f['baseline_observables']['conjugate_mu_lab_au']=row['mu_lab_au']
            f['baseline_observables']['saved_directional_values']=sens['central'][metal]['baseline']['values']
            f['reference_targets']['moment']=dict(total_dipole_lab_au=inv['rows'][metal+'_central']['total_dipole_lab_au'],
                precise_native_origin_Bohr=inv['rows'][metal+'_central']['precise_origin_Bohr'],charge_e=2,
                common_origin_A=[0.,0.,0.],source=ds['donor_inventory'],raw=inv['rows'][metal+'_central'])
            f['reference_targets']['mechanical_modes']=dict(analytic=rp['central'][metal]['analytic'],
                qualified_two_step=[r for r in mech['rows']if r['comparison']==metal],
                qualification=ds['donor_authorities']['reference_actual_review'],source=sources['mechanical_comparison'])
            f['reference_targets']['susceptibility']=dict(value=field['response']['per_metal']['La'] if metal=='La' else None,
                direction_lab=qm['field_direction'],units='a0^3',qualified=metal=='La',source=ds['La_quantum_field_collection'] if metal=='La' else None,
                limitation='One La field projection only; failed Dy and paired susceptibility remain unavailable.')
            f['reference_targets']['canonical_rotation']=dict(value=cc['La_only'] if metal=='La' else None,
                source=sources['canonical_rotation'],qualified=metal=='La',limitation='La four endpoints only; Dy canonical convergence failed.')
            f['missing'].remove('native operator export for this geometry')
        fixtures[gid]=save(out/'fixtures'/ (gid+'.json'),f)
    check(len(fixtures)==24,'exact 24 donor geometry rows')
    check(sum(c['compatible'] for c in compat.values())==6,'six dry compatible; eighteen wet incompatible')
    watergeom,wp,wkey=typed_geometry(wf['inputs'],'water',resolver)
    watercompat=compatibility(wp,wkey,watergeom,shared,shared_key,True)
    check(watercompat['compatible'],'explicit union preserves water example')
    # Materialize all ten saved lab-origin dipoles and their literal geometry rows.
    moments=[]
    for mid,row in inv['rows'].items():
        gi,_,_=typed_geometry(row['model_input'],'acetate',resolver)
        moments.append(dict(id=mid,geometry=gi,reference=row,
                            method='PBE0 electronic moment; precise charge/origin translation; no D4 moment'))
    moments_pin=save(out/'MOMENT_GEOMETRIES.json',dict(source=ds['donor_inventory'],rows=moments,
                       saved_geometry_derivatives=inv['saved_quantum_dipole_geometry_derivatives']))
    # Complete signed wet/dry exchange cycles are copied, never reinterpreted as independent observations.
    matched=read(sources['matched_comparison'])
    cycles=save(out/'HISTORICAL_CYCLES.json',dict(source=sources['matched_comparison'],
          reference_cycles=[r for r in matched['matched_water_addition_contrasts']if r['method']=='ecp_pbe0'],
          historical_model_cycles=[r for r in matched['matched_water_addition_contrasts']if r['method']=='AMOEBA_transfer'],
          model_compatible_with_shared_baseline=False,independent_observations=False,
          limitation='Historical transfer model differs from both donor17621 and shared baseline; retained descriptively only.'))
    cp=save(out/'COMPATIBILITY.json',dict(baseline=baseline_pin,rows=compat,water_example=watercompat,
          compatible_dry=6,incompatible_wet=18,shared_baseline_wet_model_values_available=False,
          reason='Wet rows contain the older water6 radius/epsilon/damping. No baseline residual is generated for them.'))
    save(out/'DATASET.json',dict(schema='short_range_many_body_donor_dataset_v1',baseline=baseline_pin,
        fixtures=fixtures,moment_geometries=moments_pin,historical_cycles=cycles,compatibility=cp,
        counts=dict(donor_states=24,dry_compatible=6,wet_incompatible=18,native_operator_fixtures=2,
                    qualified_La_susceptibility_modes=1,qualified_Dy_susceptibility_modes=0),
        force_provider=None,coefficients_fitted=False,physical_accuracy_qualified=False,
        missing=['Shared-baseline wet E/F/operators','Dry amide and bidentate operator exports','Full operator/field coordinate derivatives','Qualified Dy susceptibility and canonical torque','Complete donor fragment-energy cycles for residual U2/U3']))
    save(out/'RESULT.json',dict(status='materialized_saved_donor_dataset',implementation=pin(__file__),inputs=pin(HERE/'INPUTS.json'),
        dataset=pin(out/'DATASET.json'),checks=len(CHECKS),verified_pins=list(PINS.values()),
        cost=dict(wall_seconds=time.monotonic()-t0,CPU_seconds=time.process_time()-c0,threads=1),
        new_native_calls=0,new_quantum_calls=0,new_fits=0,new_parameter_vectors=0))
    print(json.dumps(dict(result=pin(out/'RESULT.json'),dataset=pin(out/'DATASET.json'),checks=len(CHECKS))))

if __name__=='__main__':main()
