#!/usr/bin/env python3
"""Materialize admitted saved data. No native calls, solvers or model evaluation."""
import argparse, hashlib, json, os, platform, time, re
from pathlib import Path

def pin(p):
    p=Path(p).resolve(); return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    q=Path(p['path']); assert pin(q)==p, ('changed source',q); return q
def read(p): return json.loads(verify(p).read_text())
def save(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    s=json.dumps(x,indent=2,allow_nan=False)+'\n'
    if p.exists(): assert p.read_text()==s, ('immutable output differs',p)
    else: p.write_text(s)
    return pin(p)
def geometry(inp):
    rows=[s.split() for s in verify(inp['xyz']).read_text().splitlines()[1:] if s.strip()]
    ty=[int(r[5]) for r in rows]; bonds=[[int(v)-1 for v in r[6:]] for r in rows]
    zs={500:57,501:66,36:8,37:1}
    assert all(t in zs for t in ty)
    lig=[]
    for i,t in enumerate(ty):
        if t!=36:continue
        hs=sorted(bonds[i]);assert len(hs)==2 and all(ty[h]==37 and i in bonds[h] for h in hs)
        lig.append({'id':f'water_{len(lig)+1:02d}','kind':'water','atom_indices':[i]+hs,'donor_atoms':[i],'frame_atoms_by_donor':{str(i):hs}})
    metals=[i for i,t in enumerate(ty) if t in (500,501)];assert len(metals)<=1
    return {'xyz_A':[[float(v) for v in r[2:5]] for r in rows],'atomic_numbers':[zs[t] for t in ty],'atom_ids':[int(r[0]) for r in rows],'metal_atom':metals[0] if metals else None,'ligands':lig,'types':ty,'bonds_zero_based':bonds,'coordinate_origin_A':[0.,0.,0.]}
def key_identity(inp):
    verify(inp['parameters']);s=verify(inp['key']).read_text()
    tokens=[x.split('#')[0].split() for x in s.splitlines()];tokens=[t for t in tokens if t]
    forbidden={'polar-local-response','polar-radial-response','polar-pair-response','exchpol','exchange-polar','external-field'}
    assert not any(t[0].lower() in forbidden for t in tokens), 'not zero-response baseline'
    records=[t for t in tokens if t[0].lower() in ('vdwpr','polpair','polarize')]
    return {'parameters':inp['parameters'],'effective_parameter_records_in_order':records,'old_h_H_zero':True,'added_response_keywords_absent':True}
def raw_export(saved):
    lines=[s.split() for s in verify(saved['source']).read_text().splitlines() if s.strip()]
    assert lines[0]==['TINKER_RESPONSE_MATRIX_V1'] and lines[-1]==['END_TINKER_RESPONSE_MATRIX_V1']
    n=len(saved['sites']);dim=3*n
    assert list(map(int,lines[1][1:]))==[saved['sequence'],n,n,dim]
    assert len(lines)==8+n+2*dim
    def nums(xs): return [float(v.replace('D','E')) for v in xs]
    for i,s in enumerate(saved['sites']):
        t=lines[7+i];assert list(map(int,t[1:4]))==[i+1,s['atom'],s['type']]
        assert nums(t[5:])==[s['alpha_A3']]+s['xyz_A']
    for i in range(dim):
        assert nums(lines[7+n+2*i][2:])==saved['vectors'][i]
        assert nums(lines[8+n+2*i][2:])==saved['A'][i]
    return saved

def reference_model_id(id):
    if not id.startswith('pbe0__'):return id
    t=id.split('__');parent=f'{t[1]}_central_{t[2]}'
    if t[3]=='ion_ghost_shell':return f'ion_{t[2]}'
    if t[3]=='shell_ghost_ion':return parent+'__shell'
    if t[3].startswith('water_'):
        return parent+'__'+('dimer_'+t[3][6:] if t[4]=='metal_water' else t[3])
    assert t[3] in ('near','far') and t[4] in ('metal_two_water','two_water'),id
    return parent+'__'+t[3]+'__'+t[4]

def references_and_cycles(cfg,out,m,targets,fixtures):
    art=targets['artifacts'];dev=read(art['development_data']);qm=read(art['quantum_full_precision_review'])
    common=read(qm['collection']);panel=read(common['panel']);cp=read(panel['cp_panel'])
    tr=read(art['trimer_comparison']);trmeta=read(art['trimer_native_admission']);trgeom=read(art['trimer_geometry'])
    fields=read(art['field_reference_collection']);model_field_source=read(m['configuration']['inventory'])['artifacts']['field_model']
    # Resolve the admitted zero-field baseline and finite-field source only; do not re-evaluate it.
    fm=read(model_field_source)
    static={}
    for id,s in m['static'].items():
        g=geometry(s['input']);identity=key_identity(s['input'])
        assert g['metal_atom'] is None
        static[id]={'id':id,'geometry':g,'inputs':s['input'],'baseline_identity':identity,'baseline_id':'water_fixed_six_selected17676_zero_h_H','energy_kcal_mol':s['energy_kcal_mol'],'mu_Debye':s['mu_Debye'],'mu_print_bound':s['mu_print_bound'],'permanent_lab_Debye':s['permanent_lab_Debye'],'native_operator':None,'operator_missing_reason':'No raw operator exported for this metal-free cycle term; saved baseline observables retained unchanged. Metal-centered correction is identically zero when metal is absent.','consumed_development':True,'source':s['source'],'source_authority':cfg['screen_manifest']}
    allgeom={id:read(p)['geometry'] for id,p in fixtures.items()};allgeom.update({id:s['geometry'] for id,s in static.items()})
    cells={c['cell_id']:c for c in panel['cells']+cp['cells']+trgeom['reference_cells']}
    for id,g in panel['geometries'].items():
        cells[id]={'cell_id':id,'positions_A':g['positions_A'],'symbols':g['symbols'],'ghost_mask':[False]*len(g['symbols']),'real_parent_atom_indices':list(range(len(g['symbols']))),'charge':g['charge'],'physical_multiplicity':g['physical_multiplicity'],'method':'pbe0','d4':True,'role':'full_complex','source_geometry':g['xyz']}
    refs={};joins={};symZ={'La':57,'Dy':66,'O':8,'H':1}
    raw={**common['new_fragments'],**common['reused_terms'],**trmeta['reference_rows']}
    moments={**qm['rows'],**tr['quantum_moment_rows']}
    for id,r in raw.items():
        cell=cells[id];mid=reference_model_id(id);assert mid in allgeom,(id,mid)
        indices=cell.get('real_parent_atom_indices',[i for i,g in enumerate(cell['ghost_mask']) if not g])
        xyz=[cell['positions_A'][i] for i in indices];zs=[symZ[cell['symbols'][i]] for i in indices]
        assert xyz==allgeom[mid]['xyz_A'] and zs==allgeom[mid]['atomic_numbers'],('reference/model geometry differs',id,mid)
        mom=moments[id]
        # Individual labels copied at original precision. No cycle reconstruction here.
        refs[id]={'id':id,'model_fixture_id':mid,'energy_electronic_hartree':r['electronic_energy_hartree'],'energy_D4_hartree':r.get('D4_energy_hartree'),'energy_final_hartree':r.get('final_energy_hartree',r.get('energy_hartree')),'dipole_total_origin0_au':mom['total_dipole_au'],'dipole_origin_A':[0.,0.,0.],'raw_moment_record':mom,'reference_cell':cell,'native_metadata':r,'energy_recipe':'common parent orbital/auxiliary basis electronic PBE0; D4 separate, excluded from E3','consumed_development':True,'physical_accuracy_qualified':False,'provenance':{'energy':qm['collection'] if id in common['new_fragments'] or id in common['reused_terms'] else art['trimer_native_admission'],'moments':art['quantum_full_precision_review'] if id in qm['rows'] else art['trimer_comparison']}}
        joins[id]={'fixture_id':mid,'real_geometry_exact':True,'reference_ghost_basis_retained':True,'operator_available':mid in fixtures}
    force={}
    for id,r in dev['rows'].items():
        baseline=dev['baseline_model'][id]
        assert baseline['provenance']['requested_vector']==targets['fixed_six_parameters']
        assert r['geometry']['positions_A']==allgeom[id]['xyz_A']
        native=read(fixtures[id])['baseline_observables'];assert baseline['energy']==native['energy_kcal_mol']
        force[id]={'id':id,'model_fixture':fixtures[id],'reference':r['reference'],'reference_result':r['reference_result'],'reference_geometry':r['geometry'],'baseline_model':baseline,'bodies':r['bodies'],'force_units':'kcal/mol/Angstrom','energy_recipe':'original total PBE0+D4 development energy/gradient, not electronic-only fragment labels','reference_accuracy_qualified':False,'consumed_development':True,'source':art['development_data'],'full_kernel_force_supported':False}
    field_records={}
    aliases=m['field_aliases'];stateids=targets['groups']['conditional_susceptibility']['state_ids'];inverse={v:k for k,v in stateids.items()}
    for id,r in fields['rows'].items():
        baseid=id.split('__Fz_')[0];alias=inverse[baseid];mid=aliases[alias]['fixture'];cell=r['cell'];idx=cell['real_parent_atom_indices']
        assert [cell['positions_A'][i] for i in idx]==allgeom[mid]['xyz_A']
        suffix=id.split('__Fz_')[1];modelid=alias+'__'+suffix
        field_records[id]={'id':id,'zero_operator_fixture_id':mid,'zero_reference_id':baseid,'baseline_zero_operator':fixtures[mid],'operator_at_nonzero_field_exported':False,'reference_cell':cell,'parsed':r['parsed'],'properties':r['properties'],'reference_input':r['input'],'reference_output':r['output'],'reference_result':r['result'],'model_finite_field_record':fm['rows'][modelid],'model_source':model_field_source,'consumed_development':True,'field_protocol':art['field_specification'],'not_a_seven_term_susceptibility':True}
    cycles={'pair':targets['groups']['pair_and_collective_dipoles']['pair_cycles'],'parent':targets['groups']['pair_and_collective_dipoles']['parents'],'trimer':tr['cycles'],'conditional_susceptibility':targets['groups']['conditional_susceptibility'],'common_basis_decomposition':common['decomposition'],'original_energy_force_checks':targets['groups']['original_energy_force'],'reference_limitations':targets['policy']['quantum_limitations'],'original_targets_retained':True,'new_tighter_reference_labels_automatically_substituted':False}
    for cycle in cycles['pair']+list(cycles['trimer'].values()):
        assert all(id in refs for id in cycle['reference_terms'])
        assert all(id in allgeom for id in cycle['model_terms'])
    return {'static':save(out/'STATIC.json',static),'references':save(out/'REFERENCES.json',refs),'reference_model_joins':save(out/'REFERENCE_JOINS.json',joins),'forces':save(out/'FORCES.json',force),'fields':save(out/'FIELDS.json',field_records),'cycles':save(out/'CYCLES.json',cycles),'counts':{'operators':len(fixtures),'static_terms':len(static),'individual_common_basis_reference_states':len(refs),'individual_force_reference_states':len(force),'nonzero_field_reference_states':len(field_records),'pair_cycles':len(cycles['pair']),'trimer_cycles':len(cycles['trimer'])}}
def materialize(cfg,out,one=None):
    tick=time.monotonic();cpu=time.process_time()
    m=read(cfg['screen_manifest']);em=read(cfg['export_manifest']);result=read(cfg['export_result']);review=read(cfg['export_review']);targets=read(cfg['targets'])
    assert result['status']=='baseline_matrix_inventory_complete' and result['error'] is None
    assert result['complete_fixtures']==54 and result['admitted_snapshots']==108
    baseline={'id':'water_fixed_six_selected17676_zero_h_H','candidate':targets['fixed_model_candidate'],'fixed_six_parameters':targets['fixed_six_parameters'],'parameter_order':['La_R','La_POLPAIR','Dy_R','Dy_POLPAIR','La_epsilon','Dy_epsilon'],'export_result':cfg['export_result'],'export_review':cfg['export_review'],'native_analyze':result['config']['analyze'],'native_source':result['config']['export_source_manifest'],'old_h_H_zero':True,'units':m['units'],'not_radial14_terminal':True}
    save(out/'BASELINE.json',baseline)
    fixtures={}; ids=[one] if one else list(m['rows'])
    for id in ids:
        row=m['rows'][id];rec=result['rows'][id]
        if rec.get('inherited'):
            q=rec['admission']['qualification'];native=rec['admission']['on_row'];producer=rec['producer_result']
        else:q=rec['numerical'];native=rec;producer=cfg['export_result']
        assert q['passed'];exports=[raw_export(s) for s in q['exports']]
        assert len(exports)==2 and exports[0]['A']==exports[1]['A']
        geom=geometry(row['input']);identity=key_identity(row['input'])
        assert geom['xyz_A']==row['positions_A']
        for s in exports[0]['sites']: assert geom['xyz_A'][s['atom']-1]==s['xyz_A']
        p=exports[0];v=p['vectors'];mom=native['moments'];K=m['units']['electric'];D=m['units']['debye_per_eA']
        op={'site_atom_indices':[s['atom']-1 for s in p['sites']],'alpha_A3':[s['alpha_A3'] for s in p['sites']],'A0_Aminus3':p['A'],'Ed_e_per_A2':[x[0] for x in v],'Ep_e_per_A2':[x[1] for x in v],'u_native_eA':[x[2] for x in v],'v_native_eA':[x[3] for x in v],'K_kcal_A_per_mol_e2':K,'permanent_mu_eA':[x/D for x in row['permanent_lab_Debye']],'permanent_mu_origin_A':[0.,0.,0.],'permanent_energy_kcal':None,'nonpolar_energy_kcal':mom['energy_kcal_mol']-q['polarization_energy']['printed'],'nonpolar_energy_meaning':'all native energy other than polarization, not just permanent electrostatics','snapshot':p['source'],'coordinate_derivatives':None,'field_coupling':'uniform field added to both native RHS; exact conversion in dataset units'}
        fixture={'schema':'short_range_many_body_cached_fixture_v1','id':id,'baseline_id':baseline['id'],'geometry':geom,'native_operator':op,'inputs':row['input'],'baseline_identity':identity,'baseline_observables':mom,'native_qualification':{'producer':producer,'review':cfg['export_review'],'polarization_energy':q['polarization_energy'],'two_solve_exact':q['two_solve_exact'],'snapshots':[{'source':x['source'],'sequence':x['sequence'],'vectors':x['vectors'],'numerical_bounds':x['numerical_bounds']} for x in exports]},'compatible':True,'consumed_development':True,'full_force_from_frozen_operator_supported':False,'missing':['dA0/dR','dEd/dR','dEp/dR','permanent multipole frame derivatives'],'reference_ids':[]}
        fixtures[id]=save(out/'fixtures'/f'{id}.json',fixture)
    dataset=None
    if not one:
        save(out/'FIXTURES.json',{'baseline':pin(out/'BASELINE.json'),'fixtures':fixtures})
        data=references_and_cycles(cfg,out,m,targets,fixtures)
        data.update(schema='short_range_many_body_water_dataset_v1',baseline=pin(out/'BASELINE.json'),fixtures=pin(out/'FIXTURES.json'),implementation=pin(__file__),inputs=cfg)
        dataset=save(out/'DATASET.json',data)
    return {'status':'materialized_saved_water_operators','dataset':dataset,'fixtures':fixtures,'baseline':pin(out/'BASELINE.json'),'source':pin(__file__),'configuration':cfg,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'host':platform.node()},'new_native_calls':0,'new_quantum_calls':0,'operator_solves':0,'fits':0}
def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--one');a=p.parse_args()
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): assert os.environ.get(k,'1')=='1'
    a.out.mkdir(parents=True,exist_ok=True);r=materialize(json.loads(a.inputs.read_text()),a.out,a.one);save(a.out/'RESULT.json',r);print(json.dumps({'status':r['status'],'fixtures':len(r['fixtures']),'cost':r['cost']}))
if __name__=='__main__':main()
