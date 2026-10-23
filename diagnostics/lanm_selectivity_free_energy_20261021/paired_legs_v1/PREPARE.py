#!/usr/bin/env python3
"""Pinned bound audit and neutral aqueous geometry, no native programs."""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import shlex
import time


def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    assert pin(p['path'])==p,p['path'];return Path(p['path'])


def load(p):return json.loads(verify(p).read_text())


def save(p,x):
    with Path(p).open('x')as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')


def xyz(p):
    lines=Path(p).read_text().splitlines();n=int(lines[0].split()[0]);box=None
    if len(lines)==n+2:box=list(map(float,lines.pop(1).split()))
    assert len(lines)==n+1
    rows=[]
    for i,line in enumerate(lines[1:]):
        t=line.split();assert int(t[0])==i+1
        rows.append({'label':t[1],'xyz':list(map(float,t[2:5])),'type':int(t[5]),'bonds':[int(x)-1 for x in t[6:]]})
    return rows,box


def write_xyz(path,rows,L):
    text=f'{len(rows)} AQUEOUS GEOMETRY ONLY; unqualified reserved metal types\n{L:.8f} {L:.8f} {L:.8f} 90.0 90.0 90.0\n'
    for i,r in enumerate(rows):
        text+=f'{i+1:6d} {r["label"]:3s} '+''.join(f'{v:18.10f}'for v in r['xyz'])+f' {r["type"]:5d}'+''.join(f' {v+1:6d}'for v in r['bonds'])+'\n'
    path.write_text(text)


def d2(a,b,L):return sum((x-y-L*round((x-y)/L))**2 for x,y in zip(a,b))


def topology(rows,defs,metal_index):
    assert len({tuple(r['xyz'])for r in rows})==len(rows),'Duplicate coordinates'
    for i,row in enumerate(rows):
        assert all(0<=b<len(rows)and i in rows[b]['bonds']for b in row['bonds']),'Asymmetric bond'
        if i==metal_index:assert row['type']in (900,901)and not row['bonds']
        else:assert row['type']in defs
        if row['type']==402:
            assert len(row['bonds'])==2 and all(rows[b]['type']==403 and rows[b]['bonds']==[i]for b in row['bonds'])
        if row['type']==403:assert len(row['bonds'])==1 and rows[row['bonds'][0]]['type']==402
        if row['type']in (405,414):assert not row['bonds']


def contacts(rows,L,metal_index):
    waterO=[r['xyz']for r in rows if r['type']==402];waterH=[r['xyz']for r in rows if r['type']==403]
    ions=[r for r in rows if r['type']in (405,414)];metal=rows[metal_index]['xyz']
    result={'metal_water_O_min_A':math.sqrt(min(d2(metal,v,L)for v in waterO)),
            'metal_water_H_min_A':math.sqrt(min(d2(metal,v,L)for v in waterH)),
            'ion_ion_min_A':math.sqrt(min(d2(a['xyz'],b['xyz'],L)for i,a in enumerate(ions)for b in ions[i+1:])),
            'ion_metal_min_A':math.sqrt(min(d2(metal,r['xyz'],L)for r in ions))}
    assert result['metal_water_O_min_A']>=2.4-1e-8 and result['metal_water_H_min_A']>=1.5-1e-8
    assert result['ion_ion_min_A']>=5-1e-8
    for typ,name in [(405,'Na'),(414,'Cl')]:
        chosen=[r['xyz']for r in ions if r['type']==typ]
        result[name+'_water_O_min_A']=math.sqrt(min(d2(a,b,L)for a in chosen for b in waterO))
        result[name+'_water_H_min_A']=math.sqrt(min(d2(a,b,L)for a in chosen for b in waterH))
        assert result[name+'_water_O_min_A']>=2 and result[name+'_water_H_min_A']>=1.2
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--out',required=True)
    args=ap.parse_args();start=time.monotonic();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False)
    result={'status':'started','implementation':pin(__file__),'config':pin(args.config),'local_threads':1,
            'native_calls':0,'energy_calls':0,'force_calls':0,'dynamics_steps':0,'physical_accuracy_qualified':False}
    try:
        cfg=load(result['config']);assert cfg['implementation']==pin(__file__);verify(cfg['specification'])
        parent=load(cfg['bound_result']);verify(parent['implementation']);verify(parent['source_cif'])
        review=load(cfg['bound_charge_review']);verify(review['native_parameter_output']);assert review['sum_printed_monopoles_e']=='-4.00000'
        la,box=xyz(verify(parent['boxes']['La']));dy,box2=xyz(verify(parent['boxes']['Dy']))
        L=parent['box_length_A'];assert box==box2==[L,L,L,90.,90.,90.]and L==49.324
        assert len(la)==len(dy)==11955 and parent['protein_atoms']==257
        assert la[:257]==dy[:257]and la[258:]==dy[258:]and la[257]['xyz']==dy[257]['xyz']
        assert la[257]['type']==900 and dy[257]['type']==901
        prm=verify(parent['protein_parameters']);defs={}
        for line in prm.read_text().splitlines():
            if line.split()[:1]==['atom']:
                a=shlex.split(line);defs[int(a[1])]={'class':int(a[2]),'atomic_number':int(a[5]),'mass':float(a[6]),'label':a[3]}
        assert 900 not in defs and 901 not in defs
        assert {t:defs[t]['class']for t in [402,403,405,414]}=={402:73,403:74,405:76,414:85}
        topology(la,defs,257);bound_contacts=contacts(la,L,257)
        counts=Counter(r['type']for r in la[258:]);assert counts=={402:3894,403:7788,405:8,414:7}
        # Rebuild the complete water fill; do not reuse the peptide-carved cavity.
        source,_=xyz(verify(parent['water_source']));assert len(source)==1500
        L0=L/2;metal=la[257]['xyz'];candidates=[];tiled={}
        for ix in (0,1):
          for iy in (0,1):
           for iz in (0,1):
            shift=[(v-.5)*L0 for v in (ix,iy,iz)]
            for j in range(0,len(source),3):
                water=source[j:j+3];assert [r['type']for r in water]==[1,2,2]
                o=water[0]['xyz'];pts=[[o[k]+shift[k]for k in range(3)]]
                for h in water[1:]:pts.append([pts[0][k]+h['xyz'][k]-o[k]-L0*round((h['xyz'][k]-o[k])/L0)for k in range(3)])
                rec={'xyz':pts,'source_water_index_zero_based':j//3,'tile':[ix,iy,iz],
                     'minimum_metal_distance_squared_A2':d2(pts[0],metal,L)}
                tiled[tuple(round(x,8)for x in pts[0])]=rec
                if d2(pts[0],metal,L)<2.4**2 or any(d2(h,metal,L)<1.5**2 for h in pts[1:]):continue
                candidates.append(rec)
        assert len(tiled)==4000
        selected=[]
        for idx in sorted(range(len(candidates)),key=lambda i:(-candidates[i]['minimum_metal_distance_squared_A2'],i)):
            if all(d2(candidates[idx]['xyz'][0],candidates[k]['xyz'][0],L)>=25 for k in selected):selected.append(idx)
            if len(selected)==17:break
        assert len(selected)==17
        rows=[{'label':'La','type':900,'xyz':metal,'bonds':[]}];mapping=[{'index_1based':1,'role':'metal','bound_index_1based':258,'La_type':900,'Dy_type':901,'physical_charge':3}]
        for i,w in enumerate(candidates):
            if i in selected:continue
            n=len(rows);rows.extend([{'label':'O','type':402,'xyz':w['xyz'][0],'bonds':[n+1,n+2]},
                                    {'label':'H','type':403,'xyz':w['xyz'][1],'bonds':[n]},
                                    {'label':'H','type':403,'xyz':w['xyz'][2],'bonds':[n]}])
            for k in range(3):mapping.append({'index_1based':n+k+1,'role':'water','water_atom':('O','H1','H2')[k],
                                               'source_water_index_zero_based':w['source_water_index_zero_based'],'tile':w['tile'],
                                               'type':rows[n+k]['type'],'class':defs[rows[n+k]['type']]['class']})
        for i,idx in enumerate(selected):
            w=candidates[idx];t=405 if i<7 else 414;n=len(rows)
            rows.append({'label':'Na' if i<7 else 'Cl','type':t,'xyz':w['xyz'][0],'bonds':[]})
            mapping.append({'index_1based':n+1,'role':'Na' if i<7 else 'Cl','type':t,'class':defs[t]['class'],
                            'source_water_index_zero_based':w['source_water_index_zero_based'],'tile':w['tile']})
        topology(rows,defs,0);aqueous_contacts=contacts(rows,L,0)
        aqueous_count=Counter(r['type']for r in rows);assert aqueous_count[405]==7 and aqueous_count[414]==10
        assert 3+aqueous_count[405]-aqueous_count[414]==0
        outputs={}
        for metal_name,metal_type in [('La',900),('Dy',901)]:
            copy=[dict(r)for r in rows];copy[0]=dict(copy[0],label=metal_name,type=metal_type)
            path=out/(metal_name+'_aqueous.xyz');write_xyz(path,copy,L);actual,actualbox=xyz(path)
            topology(actual,defs,0);assert actualbox==box and len(actual)==len(rows);outputs[metal_name]=pin(path)
        actual_la,_=xyz(outputs['La']['path']);actual_dy,_=xyz(outputs['Dy']['path'])
        assert actual_la[1:]==actual_dy[1:]and actual_la[0]['xyz']==actual_dy[0]['xyz']
        bound_mapping=[]
        by_index={r['tinker_index']:r for r in parent['source_joins']};assert len(by_index)==257
        for i,row in enumerate(la):
            if i<257:bound_mapping.append({'index_1based':i+1,'role':'peptide','type':row['type'],'class':defs[row['type']]['class'],'source_join':by_index[i+1]})
            elif i==257:bound_mapping.append({'index_1based':258,'role':'metal','La_type':900,'Dy_type':901,'physical_charge':3})
            else:
                role='water'if row['type']in (402,403)else row['label'];oxygen_index=i if row['type']==402else(row['bonds'][0]if row['type']==403else i)
                origin=tiled[tuple(round(x,8)for x in la[oxygen_index]['xyz'])]
                bound_mapping.append({'index_1based':i+1,'role':role,'type':row['type'],'class':defs[row['type']]['class'],
                                      'source_water_index_zero_based':origin['source_water_index_zero_based'],'tile':origin['tile']})
        save(out/'ATOM_MAPPING.json',{'indexing':'one_based_native; source_water_index explicitly zero_based',
              'within_each_leg_La_to_Dy':'identity for every atom; only metal label/type changes',
              'between_legs':'metal bound258 to aqueous1; no atomwise alchemical map is required between separate thermodynamic legs',
              'bound':bound_mapping,'aqueous':mapping})
        c0=1/(6.02214076e23*L**3*1e-27)
        result.update(status='geometry_prepared_hamiltonian_unassigned',parent_bound_result=cfg['bound_result'],
                      source_water=parent['water_source'],protein_parameters=parent['protein_parameters'],bound_boxes=parent['boxes'],
                      aqueous_boxes=outputs,atom_mapping=pin(out/'ATOM_MAPPING.json'),box_length_A=L,
                      bound={'atoms':len(la),'peptide_atoms':257,'waters':3894,'Na':8,'Cl':7,'formal_net_charge':0,'contacts':bound_contacts},
                      aqueous={'atoms':len(rows),'peptide_atoms':0,'waters':aqueous_count[402],'Na':7,'Cl':10,'formal_net_charge':0,'contacts':aqueous_contacts},
                      one_particle_concentration_M=c0,background_salt_pairs=7,background_salt_pair_concentration_M=7*c0,
                      bound_mobile_monovalent_ionic_strength_M=7.5*c0,aqueous_mobile_monovalent_ionic_strength_M=8.5*c0,
                      fixed_charge_path=True,metal_parameters_selected=False,executable_periodic_key_prepared=False,
                      water_geometry_origin='TIP3P packing assigned flexible AMOEBA water; not equilibrated',
                      root_convention=cfg['root_convention'],type_definitions={str(t):defs[t]for t in [402,403,405,414]})
    except Exception as exc:result.update(status='failed',error_type=type(exc).__name__,error=str(exc))
    result['elapsed_seconds']=time.monotonic()-start;save(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json'),'error':result.get('error')}))
    return 0 if result['status']=='geometry_prepared_hamiltonian_unassigned'else 1


if __name__=='__main__':raise SystemExit(main())
