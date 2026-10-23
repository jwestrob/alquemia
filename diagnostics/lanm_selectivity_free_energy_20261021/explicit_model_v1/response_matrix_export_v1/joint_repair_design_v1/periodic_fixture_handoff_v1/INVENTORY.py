#!/usr/bin/env python3
"""Metadata only; no engine import, subprocess, geometry writing or tail sum."""
from collections import Counter
from decimal import Decimal, localcontext, ROUND_CEILING
import hashlib
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p,'Changed pin '+p['path']
    return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def sqrt_upper(x):
    with localcontext() as ctx:
        ctx.prec=60
        y=x.sqrt().quantize(Decimal('.000000000001'),rounding=ROUND_CEILING)
        while y*y<x:y+=Decimal('.000000000001')
        assert y*y>=x
        return str(y)
def xyz(p):
    lines=verify(p).read_text().splitlines();n=int(lines[0].split()[0]);box=[Decimal(s) for s in lines[1].split()];atoms=[]
    assert len(lines)==n+2
    for i,line in enumerate(lines[2:],1):
        t=line.split();assert int(t[0])==i;atoms.append({'id':i,'name':t[1],'xyz':[Decimal(s) for s in t[2:5]],'type':int(t[5]),'bonds':list(map(int,t[6:]))})
    return box,atoms

def main():
    tick,cpu=time.monotonic(),time.process_time();inputs=json.loads((HERE/'INPUTS.json').read_text())
    for p in inputs.values():verify(p)
    mapping=read(inputs['mapping_manifest']);paired=read(inputs['paired_preparation']);admission=read(inputs['BIO_actual_review']);native=read(inputs['BIO_actual_result']);custom=read(inputs['custom_qualification_manifest'])
    assert mapping['status']=='aqueous_water_mapping_prepared_no_calls' and mapping['periodic_numerical_policy_bound'] is False
    assert admission['status']=='actual_mapping_review_passed' and admission['actual_result']==inputs['BIO_actual_result']
    outputs=[];coordinates=[]
    periodic_keywords={'ewald','ewald-alpha','pewald-alpha','pme-grid','pme-order','ppme-order','ewald-cutoff','vdw-cutoff','vdw-taper','ewald-boundary','lrc','vdw-correction','cutoff','neighbor-list','polar-iter'}
    for row in mapping['aqueous_endpoints']:
        metal=row['metal'];box,atoms=xyz(row['xyz']);original=paired['aqueous_boxes'][metal]
        assert verify(row['xyz']).read_bytes()==verify(original).read_bytes() and row['source']==original
        assert box==[Decimal('49.324')]*3+[Decimal(90)]*3
        for p in row['files'].values():verify(p)
        assert row['parameters']['sha256']==mapping['BIO_parameters_preserved_byte_for_byte']['sha256']
        counts=Counter(a['type'] for a in atoms);mt=900 if metal=='La' else 901
        assert counts==Counter({mt:1,402:3980,403:7960,405:7,414:10})
        assert atoms[0]['type']==mt and not atoms[0]['bonds']
        water=[];covered=set()
        for a in atoms:
            if a['type']==402:
                hs=a['bonds'];assert len(hs)==2 and len(set(hs))==2
                assert all(atoms[h-1]['type']==403 and atoms[h-1]['bonds']==[a['id']] for h in hs)
                ids=[a['id']]+hs;assert not covered.intersection(ids);covered.update(ids);water.append(ids)
            elif a['type'] in (mt,405,414):assert not a['bonds']
        assert len(covered)==11940 and len(water)==3980
        charges={mt:Decimal(3),402:Decimal('-.51966'),403:Decimal('.25983'),405:Decimal(1),414:Decimal(-1)}
        charge=sum(charges[t]*n for t,n in counts.items());assert charge==0
        mins=[min(a['xyz'][k] for a in atoms) for k in range(3)];maxs=[max(a['xyz'][k] for a in atoms) for k in range(3)]
        origin=atoms[0]['xyz'];oxygen=[a for a in atoms if a['type']==402];delta=[[a['xyz'][k]-origin[k] for k in range(3)] for a in oxygen];r2=[sum(v*v for v in d) for d in delta];d2=max(r2);maxaxis=[max(abs(d[k]) for d in delta) for k in range(3)];rect2=sum(v*v for v in maxaxis);L=box[0]
        key=verify(row['key']).read_text();tokens=[s.split('#',1)[0].split() for s in key.splitlines() if s.split('#',1)[0].split()];bound=[t for t in tokens if t[0].lower() in periodic_keywords]
        assert not bound,'Unexpected actually bound periodic policy'
        output={'id':row['id'],'metal':metal,'geometry':row['xyz'],'original_geometry':original,'geometry_byte_identical_to_original':True,'key':row['key'],'parameters':row['parameters'],'atoms':len(atoms),'type_counts':dict(sorted(counts.items())),'water_molecules':3980,'salt':{'Na':7,'Cl':10},'net_monopole_charge_e':str(charge),'topology':{'intact_reciprocal_water_graphs':3980,'isolated_metal_and_salt':18,'water_triplet_map_sha256':hashlib.sha256(json.dumps(water,separators=(',',':')).encode()).hexdigest(),'zero_metal_water_bonds':True},'cell_A_degrees':[str(v) for v in box],'metal_coordinate_A':[str(v) for v in origin],'literal_coordinate_min_A':[str(v) for v in mins],'literal_coordinate_max_A':[str(v) for v in maxs],'lattice_bounds':{'H':'49.324 I Angstrom; exact declared orthogonal cube','singular_lower_A':'49.324','singular_upper_A':'49.324','volume_A3':str(L**3),'fixed_literal_M_to_O_D_squared_A2':str(d2),'fixed_literal_M_to_O_D_upper_A':sqrt_upper(d2),'rectangular_M_to_O_D_upper_A':sqrt_upper(rect2),'canonical_cell_pair_D_squared_A2':str(3*L*L),'canonical_cell_pair_D_upper_A':sqrt_upper(3*L*L),'minimum_image_radius_upper_A':sqrt_upper(3*L*L/4),'minimum_image_insphere_radius_A':str(L/2),'future_configuration_bound':False},'actually_bound_periodic_keywords':bound,'key_polar_eps':'1.0e-8','periodic_status':'prepared-only; no periodic numerical admission','custom_baseline_key':custom['rows']['aquo9_central_'+metal]['input']['key']}
        outputs.append(output);coordinates.append(atoms)
    assert all(a['xyz']==b['xyz'] and a['bonds']==b['bonds'] and (a['id']==1 or (a['name']==b['name'] and a['type']==b['type'])) for a,b in zip(*coordinates))
    result={'status':'periodic_fixture_metadata_handoff_only','implementation':pin(__file__),'declaration':pin(HERE/'DECLARATION.md'),'inputs':pin(HERE/'INPUTS.json'),'fixtures':outputs,'La_Dy_coordinates_topology_nonmetal_types_exactly_identical':True,'prior_BIO_actual_status':native['status'],'prior_BIO_review_status':admission['status'],'prior_BIO_scope':'two nonperiodic 25-atom aquo8 controls; four native programs','old_mapped_candidate':'17605 six-water candidate, not current custom local-response candidate','proposed_unbound_periodic_policy':{'authority':inputs['paired_report'],'fixed_cube_A':'49.324','surface':'conducting/TINFOIL; omit EWALD-BOUNDARY','Ewald_real_cutoff_A':9,'permanent_Ewald_alpha_A_minus1':.4,'polar_Ewald_alpha_A_minus1':.4,'PME_grid':[64,64,64],'permanent_PME_order':5,'polar_PME_order':5,'polar_eps':1e-7,'polar_maxiter':100,'vdw_cutoff_A':9,'vdw_taper_A':8.1,'long_range_vdw_correction':False,'actually_bound':False},'new_native_calls':0,'new_geometry':False,'tail_evaluations':0,'parameter_replacement':False,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    dest=HERE/'INVENTORY.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({'inventory':pin(dest),'cost':result['cost'],'fixture_cell_bounds':[{'id':r['id'],'bounds':r['lattice_bounds'],'min':r['literal_coordinate_min_A'],'max':r['literal_coordinate_max_A']} for r in outputs]}))

if __name__=='__main__':main()
