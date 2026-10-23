#!/usr/bin/env python3
"""Geometric mixed-donor fixtures; zero electronic calls."""
import argparse
from collections import Counter
import importlib.util
import json
import math
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'REFERENCE.py'
SOURCE_SHA='eec0e60d982bda4e1a8685c6120f25c6ab74cc57e114cd3d8f370b4baccba129'
spec=importlib.util.spec_from_file_location('frozen_chemistry',SOURCE)
ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
ref.check(ref.pin(SOURCE)['sha256']==SOURCE_SHA,'Frozen helper differs')
PANEL=ref.ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/prepared_v2/PANEL.json'
PANEL_SHA='d5929526cf70ac957c847a5936ed229e9edf1322764b0cc287a59292e01204f9'


def project_perp(v,u): return ref.sub(v,ref.mul(u,ref.dot(v,u)))


def basis_for(s,family):
    x=s['positions']; oxy=[i for i,a in enumerate(s['symbols']) if a=='O']
    if family=='acetate_bi':
        center=[sum(x[i][k] for i in oxy)/2 for k in range(3)]
        e1=ref.unit(ref.sub(x[oxy[1]],x[oxy[0]]))
        e0=ref.unit(project_perp(ref.sub(center,x[0]),e1))
        half=ref.norm(ref.sub(x[oxy[1]],x[oxy[0]]))/2
        radius=math.sqrt(2.5**2-half**2)
    else:
        donor=min(oxy,key=lambda i:ref.norm(ref.sub(x[i],x[0])))
        center=x[donor]; e0=ref.unit(ref.sub(center,x[0]))
        candidates=[project_perp(ref.sub(x[i],center),e0) for i in range(1,len(x)) if i!=donor]
        e1=ref.unit(max(candidates,key=ref.norm));radius=2.5;oxy=[donor]
    e2=ref.cross(e0,e1)
    return center,[e0,e1,e2],radius,oxy


def transformed(s,center,basis,u,tangent,radius,angle,shift):
    f0=u; tt=ref.cross(u,tangent)
    f1=ref.add(ref.mul(tangent,math.cos(angle)),ref.mul(tt,math.sin(angle)))
    f2=ref.cross(f0,f1);target=ref.mul(u,radius+shift)
    out=[]
    for x in s['positions'][1:]:
        y=ref.sub(x,center); coeff=[ref.dot(y,e) for e in basis]
        out.append(ref.add(target,[sum(coeff[j]*f[k] for j,f in enumerate((f0,f1,f2))) for k in range(3)]))
    det=ref.dot(f0,ref.cross(f1,f2))*ref.dot(basis[0],ref.cross(basis[1],basis[2]))
    return out,det


def contacts(ligand,ls,waters,ws):
    pairs=[(ref.norm(ref.sub(a,b)),i,j) for i,a in enumerate(ligand) for j,b in enumerate(waters)]
    heavy=[p for p in pairs if ls[p[1]]!='H' and ws[p[2]]!='H']
    withh=[p for p in pairs if ls[p[1]]=='H' or ws[p[2]]=='H']
    return {'all':min(pairs),'heavy_heavy':min(heavy),'H_containing':min(withh)}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    tick=time.monotonic(); p=json.loads(PANEL.read_text());ref.check(ref.pin(PANEL)['sha256']==PANEL_SHA,'Panel changed')
    aq={m:next(g for g in p['geometries'] if g['geometry_id']=='aquo8_central_'+m) for m in ('La','Dy')}
    for g in aq.values():ref.verify(g['xyz'])
    ref.check(aq['La']['positions_A']==aq['Dy']['positions_A'],'Aquo states differ')
    base=aq['La']['positions_A']; old=ref.published_structures();args.out.mkdir(parents=True,exist_ok=False)
    (args.out/'xyz').mkdir();geometries=[];reactions=[];construction=[];failures=[]
    for family,remove in [('acetate_mono',[0]),('acetate_bi',[0,1]),('acetamide',[0])]:
        s=old[family];center,basis,radius,donors=basis_for(s,family)
        removed=[base[1+3*j] for j in remove];u=ref.unit([sum(x[k] for x in removed) for k in range(3)])
        tangent=ref.unit(project_perp(ref.sub(removed[-1],removed[0]),u)) if len(remove)==2 else ref.unit(ref.cross(u,[0.,0.,1.]))
        kept=[j for j in range(8) if j not in remove]
        waters=[base[1+3*j+k] for j in kept for k in range(3)];ws=['O','H','H']*len(kept);ls=s['symbols'][1:]
        candidates=[]
        for turn in range(8):
            angle=turn*math.pi/4;cs=[]
            for shift in (-.15,0.,.15):
                xx,det=transformed(s,center,basis,u,tangent,radius,angle,shift)
                cs.append(contacts(xx,ls,waters,ws)['all'][0])
            candidates.append({'turn':turn,'angle_radians':angle,'minimum_contact_over_three_positions_A':min(cs)})
        chosen=max(candidates,key=lambda r:r['minimum_contact_over_three_positions_A'])
        construction.append({'family':family,'removed_water_indices':remove,'retained_water_indices':kept,'candidates':candidates,'chosen':chosen,
            'selection':'largest minimum ligand-water separation across all three translations; no energies; identical for both metals'})
        for label,shift in [('minus',-.15),('central',0.),('plus',.15)]:
            ligand,det=transformed(s,center,basis,u,tangent,radius,chosen['angle_radians'],shift)
            cc=contacts(ligand,ls,waters,ws)
            if cc['heavy_heavy'][0]<1.8 or cc['H_containing'][0]<1.2:failures.append({'family':family,'distortion':label,'contacts':cc})
            original=s['positions'][1:]
            error=max(abs(ref.norm(ref.sub(ligand[i],ligand[j]))-ref.norm(ref.sub(original[i],original[j]))) for i in range(len(ligand)) for j in range(i+1,len(ligand)))
            ref.check(error<1e-10 and abs(det-1)<1e-12,'Rigid ligand transformation failed')
            for metal in ('La','Dy'):
                symbols=[metal]+ls+ws;positions=[[0.,0.,0.]]+ligand+waters
                charge=2 if family.startswith('acetate') else 3;mult=1 if metal=='La' else 6
                ref.check((sum(ref.Z[a] for a in symbols)-charge-mult+1)%2==0,'Electron parity')
                gid=f'hydrated_{family}_{label}_{metal}';xyz=ref.xyz_write(args.out/'xyz'/(gid+'.xyz'),symbols,positions,gid)
                ss,xx=ref.xyz_read(xyz['path']);nlig=len(ligand);bonds=[[]]+[[a for a in neighbors if a!=0] for neighbors in s['bonds'][1:]]
                for j in range(len(kept)):
                    start=1+nlig+3*j;bonds += [[start+1,start+2],[start],[start]]
                row={'geometry_id':gid,'case_id':'hydrated_'+family,'donor_family':family,'distortion':label,
                    'metal':metal,'metal_index':0,'charge':charge,'physical_multiplicity':mult,'atoms':len(ss),
                    'symbols':ss,'positions_A':xx,'xyz':xyz,'water_count':len(kept),'displaced_water_count':len(remove),
                    'donor_indices':donors,'ligand_indices':list(range(1,nlig+1)),'radial_direction':u,
                    'radial_shift_A':shift,'contacts_A':cc,'rigid_internal_distance_max_error_A':error,
                    'source_bonds_zero_based':bonds,'source_tinker_types':s['tinker_types']+[v for j in kept for v in (36,37,37)],
                    'source':'Rigid exact published ligand plus exact retained aquo8 waters; constructed nonequilibrium fixture'}
                ref.check(len(bonds)==len(ss)==len(row['source_tinker_types']),'Topology length')
                for i,bb in enumerate(bonds):
                    for j in bb:ref.check(i in bonds[j],'Nonsymmetric covalent graph')
                ref.check(xx[-3*len(kept):]==waters,'Retained water coordinates changed')
                left=Counter(aq[metal]['symbols'])+Counter(ls);right=Counter(ss)+Counter({'O':len(remove),'H':2*len(remove)})
                ref.check(left==right and 3+(-1 if family.startswith('acetate') else 0)==charge,'Reaction imbalance')
                geometries.append(row)
        reactions.append({'family':family,'water_displaced':len(remove),'left':['M(W)8',family+'_free'],
            'right':['hydrated_'+family,str(len(remove))+'*water_free'],'charge_balanced':True,'atom_counts_balanced':True,
            'contrast':'E(Dy mixed)-E(Dy aquo8)-E(La mixed)+E(La aquo8); free fragments cancel'})
    ref.check(len(geometries)==18,'Wrong mixed denominator')
    for g in geometries:
        peer=next(q for q in geometries if q['case_id']==g['case_id'] and q['distortion']==g['distortion'] and q['metal']!=g['metal'])
        ref.check(peer['positions_A']==g['positions_A'],'Paired geometry differs')
    status='prepared' if not failures else 'failed_contact_gate'
    free_water_positions=base[1:4]
    free_water_xyz=ref.xyz_write(args.out/'xyz/water_free_aquo_geometry.xyz',['O','H','H'],free_water_positions,'Exact aquo8 water internal geometry; no optimization')
    free_water={'geometry_id':'water_free_aquo_geometry','case_id':'water','charge':0,'physical_multiplicity':1,'atoms':3,
        'symbols':['O','H','H'],'positions_A':free_water_positions,'xyz':free_water_xyz,
        'source_bonds_zero_based':[[1,2],[0],[0]],'source_tinker_types':[36,37,37],
        'source':'Exact first water of reused aquo8; all aquo waters are rigid copies'}
    result={'status':status,'protocol':'lady_hydrated_mixed_preparation_v1','specification':ref.pin(HERE/'SPEC.md'),'implementation':ref.pin(__file__),
        'source_panel':ref.pin(PANEL),'helper':ref.pin(SOURCE),'geometries':geometries,'reused_controls':list(aq.values()),
        'free_ligands':[f for f in p['free_fragments'] if f['case_id']!='water'],'free_water':free_water,
        'construction':construction,'balanced_reactions':reactions,
        'contact_failures':failures,'molecular_calls':0,'local_numeric_threads':1,'elapsed_seconds':time.monotonic()-tick,
        'interpretation':'Geometric stress fixtures only; scoring and physical qualification not yet released'}
    ref.write(args.out/'PANEL.json',result)
    print(json.dumps({'status':status,'new_geometries':len(geometries),'reused_controls':2,'failures':failures}))
    return 0 if not failures else 1


if __name__=='__main__':raise SystemExit(main())
