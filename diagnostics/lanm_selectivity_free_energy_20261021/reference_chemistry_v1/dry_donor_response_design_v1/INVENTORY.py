#!/usr/bin/env python3
"""Saved-only authority/dipole inventory for the proposed dry-donor diagnostic."""
import hashlib, importlib.util, json, math, re, time
from pathlib import Path
import numpy as np

ROOT=Path.cwd();REL=Path('lanm_selectivity_free_energy_20261021');D=ROOT/'diagnostics'/REL;W=ROOT/'workspaces'/REL
HERE=Path(__file__).resolve().parent;OUT=W/'reference_chemistry_v1/dry_donor_response_design_v1'
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def read(p):return json.loads(Path(p).read_text())
def verify(a):
    assert pin(a['path'])==a;return Path(a['path'])
def half(token):
    mantissa,_,exponent=token.lower().partition('e');return .5*10**(int(exponent or 0)-(len(mantissa.partition('.')[2]) if '.' in mantissa else 0))

def main():
    tick,cpu=time.monotonic(),time.process_time();OUT.mkdir(exist_ok=False)
    panel_path=W/'reference_chemistry_v1/dry_acetate_force_qualification_v1/preparation_v1/PANEL.json';panel=read(panel_path)
    qm_path=W/'reference_chemistry_v1/dry_acetate_force_qualification_v1/remaining_submission_v1/collection_v1/RESULT.json';qm=read(qm_path)
    mb=W/'explicit_model_v1/paired_donor_fit_v1/dry_acetate_component_derivatives_v1'
    mp=mb/'preparation_v1/MANIFEST.json';model=read(mp);mr=mb/'run_v1/execution/RESULT.json';model_result=read(mr)
    assert model['reference_panel']==pin(panel_path) and model_result['complete_cells']==16
    helper_path=D/'reference_chemistry_v1/saved_electronic_redistribution_v1/independent_actual_review_v1/REVIEW_v3.py'
    spec=importlib.util.spec_from_file_location('qualified_property_format',helper_path);helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    ids={m+'_central':(panel['central'][m]['saved_result'],model['central'][m]['source']['effective_input'],model['central'][m]['source']['authority']) for m in ['La','Dy']}
    for c in panel['cells']:
        if c['direction']=='rotation':
            cid=c['cell_id'];ids[cid]=(qm['slots'][cid]['result'],model['rows'][cid]['input'],pin(Path(model_result['rows'][cid]['receipts'][0]['output']['path']).parent/'RESULT.json'))
    rows={}
    for cid,(result_pin,model_input,model_authority) in ids.items():
        q=read(verify(result_pin));raw=verify(q['output']).read_text();property_path=Path(q['output']['path']).with_name('endpoint.property.txt');prop=property_path.read_text()
        geo=helper.section(prop,'Geometry');assert '&Units "Bohr"' in geo
        geom=re.findall(r'^\s*([A-Za-z]+)\s+('+helper.N+r')\s+('+helper.N+r')\s+('+helper.N+r')\s*$',geo,re.M)
        assert len(geom)==8;coords=np.array([[float(x) for x in g[1:]] for g in geom]);z=np.array([11,6,6,8,8,1,1,1]);assert z.sum()==42
        native=raw.split('CARTESIAN COORDINATES (A.U.)');assert len(native)==2
        atoms=[];mass_tokens=[]
        pattern=r'\s*(\d+)\s+([A-Za-z]+)\s+('+helper.N+r')\*?\s+\d+\s+('+helper.N+r')\s+('+helper.N+r')\s+('+helper.N+r')\s+('+helper.N+r')\s*'
        for line in native[1].splitlines():
            m=re.fullmatch(pattern,line)
            if m:atoms.append((int(m[1]),m[2],float(m[3]),float(m[4])));mass_tokens.append(m[4])
            elif atoms:break
        assert len(atoms)==8 and [a[1] for a in atoms]==[g[0] for g in geom] and [a[2] for a in atoms]==z.tolist()
        dip=helper.section(prop,'SCF_Dipole_Moment');tot=np.array(helper.array(dip,'dipoleTotal',3));nuc=np.array(helper.array(dip,'dipoleNucContrib',3));elec=np.array(helper.array(dip,'dipoleElecContrib',3))
        assert np.max(np.abs(tot-nuc-elec))<1e-11
        origin=(z@coords-nuc)/42.;printed_origin=np.array([float(v) for v in re.findall(r'Position of electric origin\s*\.\.\.\s+('+helper.N+r')\s+('+helper.N+r')\s+('+helper.N+r')',raw)[0]])
        assert np.max(np.abs(origin-printed_origin))<=.5e-6+1e-11
        masses=np.array([a[3] for a in atoms]);dm=np.array([half(t) for t in mass_tokens]);com=masses@coords/masses.sum()
        com_bound=(dm@np.abs(coords-com))/(masses.sum()-dm.sum())+1e-11
        assert np.all(np.abs(origin-com)<=com_bound)
        printed=np.array([float(v) for v in re.findall(r'Total Dipole Moment\s*:\s+('+helper.N+r')\s+('+helper.N+r')\s+('+helper.N+r')',raw)[0]])
        assert np.max(np.abs(tot-printed))<=.5e-9+1e-12
        lab=tot+2*origin;assert np.max(np.abs(lab-(printed+2*printed_origin)))<=1.0005e-6+1e-11
        for a in model_input.values():verify(a)
        rows[cid]=dict(quantum_result=result_pin,quantum_output=q['output'],property_file=pin(property_path),model_input=model_input,model_zero_energy_authority=model_authority,
            effective_nuclear_charge=42,electrons=40,charge=2,precise_origin_Bohr=origin.tolist(),rounded_mass_COM_Bohr=com.tolist(),COM_roundoff_bound_Bohr=com_bound.tolist(),
            total_dipole_native_origin_au=tot.tolist(),total_dipole_lab_au=lab.tolist(),printed_origin_total_bound_au=1.0005e-6,
            quantum_static_polarizability_available=False,model_conjugate_dipole_available=False)
        if cid.endswith('_central'):
            rows[cid]['population_metal_charge_e']={scheme:helper.array(helper.section(prop,'SCF_'+scheme+'_Population_Analysis'),'AtomicCharges',8)[0] for scheme in ['Mulliken','Loewdin']}
    assert len(rows)==10
    g=panel['central']['La']['original_cell']['geometry'];x=np.array(g['positions_A']);u=x[3]-x[0];u/=np.linalg.norm(u);center=x[1:].mean(axis=0);axis=np.cross(u,x[4]-x[3]);axis/=np.linalg.norm(axis)
    tangent=np.cross(axis,x[3]-center);tangent/=np.linalg.norm(tangent)
    qderiv={}
    for metal in ['La','Dy']:
        qderiv[metal]={}
        for step,suffix in [(.01,'010'),(.005,'005')]:
            plus=np.array(rows[f'pbe0__{metal}__rotation__plus_{suffix}__complex']['total_dipole_lab_au']);minus=np.array(rows[f'pbe0__{metal}__rotation__minus_{suffix}__complex']['total_dipole_lab_au'])
            value=(plus-minus)/(2*step);qderiv[metal][str(step)]=dict(vector_ea0_per_A=value.tolist(),tangent_projection_ea0_per_A=float(value@tangent))
    qderiv['paired']={s:dict(vector_ea0_per_A=(np.array(qderiv['Dy'][s]['vector_ea0_per_A'])-qderiv['La'][s]['vector_ea0_per_A']).tolist(),
        tangent_projection_ea0_per_A=qderiv['Dy'][s]['tangent_projection_ea0_per_A']-qderiv['La'][s]['tangent_projection_ea0_per_A']) for s in ['0.01','0.005']}
    authorities=dict(reference_panel=pin(panel_path),reference_derivative_result=pin(qm_path),reference_actual_review=pin(W/'reference_chemistry_v1/dry_acetate_force_qualification_v1/independent_actual_review_v1/REVIEW.json'),
        model_manifest=pin(mp),model_derivative_result=pin(mr),selected_donor_review=model['selected_review'],selected_donor_result=model['selected_fit_result'],original_analyze=model['analyze'],
        property_origin_helper=pin(helper_path),origin_evidence=pin(D/'reference_chemistry_v1/electrical_response_design_v1/ORIGIN_ECP_EVIDENCE.json'),
        previous_field_spec=pin(D/'reference_chemistry_v1/conditional_susceptibility_v1/SPEC.md'),previous_field_model_ready=pin(D/'explicit_model_v1/conditional_susceptibility_v1/model/READY.json'))
    q=dict(status='saved_response_evidence_and_design_inventory',implementation=pin(__file__),authorities=authorities,rows=rows,
        geometry_defined_field_direction=tangent.tolist(),definition='normalize(n cross (O3 minus ligand centroid)); same fixed lab vector for all ten geometries and both metals',
        saved_quantum_dipole_geometry_derivatives=qderiv,
        proposed_first_stage=dict(quantum_calls=0,original_model_energy_calls=40,geometries=10,fields_au=[-.001,-.0005,.0005,.001],zero_energies_reused=10,model_parameters_changed=False,geometry_changes=False),
        conditional_QM_susceptibility_extension=dict(quantum_calls=10,compatible_central_SP_zeros=2,central_fields=8,model_calls=0,explanation='Only if separately chosen after first-stage evidence; center model field quartet already present in the forty-call design.'),
        new_native_calls=0,new_quantum_calls=0,new_model_calls=0,held_outcomes_opened=False,wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1)
    f=OUT/'SOURCE_INVENTORY.json'
    with f.open('x') as stream:json.dump(q,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(inventory=pin(f),central={m:rows[m+'_central']['total_dipole_lab_au'] for m in ['La','Dy']},tangent=tangent.tolist(),derivatives=qderiv,cost=dict(wall=q['wall_seconds'],CPU=q['CPU_seconds']))))

if __name__=='__main__':main()
