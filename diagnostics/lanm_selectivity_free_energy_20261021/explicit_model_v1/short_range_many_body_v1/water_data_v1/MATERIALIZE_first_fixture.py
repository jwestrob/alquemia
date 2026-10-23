#!/usr/bin/env python3
"""Materialize admitted saved data. No native calls, solvers or model evaluation."""
import argparse, hashlib, json, os, platform, time
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
    if not one:save(out/'FIXTURES.json',{'baseline':pin(out/'BASELINE.json'),'fixtures':fixtures})
    return {'status':'materialized_saved_water_operators','fixtures':fixtures,'baseline':pin(out/'BASELINE.json'),'source':pin(__file__),'configuration':cfg,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'host':platform.node()},'new_native_calls':0,'new_quantum_calls':0,'operator_solves':0,'fits':0}
def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--one');a=p.parse_args()
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): assert os.environ.get(k,'1')=='1'
    a.out.mkdir(parents=True,exist_ok=True);r=materialize(json.loads(a.inputs.read_text()),a.out,a.one);save(a.out/'RESULT.json',r);print(json.dumps({'status':r['status'],'fixtures':len(r['fixtures']),'cost':r['cost']}))
if __name__=='__main__':main()
