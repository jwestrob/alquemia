#!/usr/bin/env python3
"""Prepare four fixed W9 radial geometries; metadata only, no model imports."""
import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
METALS=('La','Dy');SIGNS=(('minus',-.10),('plus',.10));TOLERANCE_A=2e-12

def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def require(ok,message):
    if not ok:raise ValueError(message)
def verify(p):
    require(pin(p['path'])==p,'Changed source: '+p['path']);return Path(p['path'])
def save(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def distance(a,b):return math.sqrt(sum((x-y)**2 for x,y in zip(a,b)))

def default_config():
    return {'implementation':pin(__file__),'root_declaration':pin(HERE/'ROOT_DECLARATION.md'),
        'original_typed_manifest':pin(W/'explicit_model_v1/prepared/MANIFEST.json'),
        'admitted17670_typed_xyz':{metal:pin(W/f'aqueous_qualification_design_v1/coordination_screen_v1/preparation_v2/inputs/aquo9_central_{metal}/system.xyz') for metal in METALS},
        'radial_steps_A':dict(SIGNS),'coordinate_decimal_places':12,'arithmetic_tolerance_A':TOLERANCE_A,
        'new_cells':4,'atoms_per_cell':28,'waters_per_cell':9,'new_native_calls':0,'candidate_binding':None}

def parse_typed(path):
    lines=Path(path).read_text().splitlines();require(len(lines)==29 and int(lines[0].split()[0])==28,'Typed source must have28 atoms')
    symbols=[];positions=[];types=[];bonds=[]
    for i,line in enumerate(lines[1:]):
        t=line.split();require(int(t[0])==i+1,'Typed index order differs')
        symbols.append(t[1]);positions.append([float(v) for v in t[2:5]]);types.append(int(t[5]));bonds.append([int(v)-1 for v in t[6:]])
    require(all(math.isfinite(v) for x in positions for v in x),'Nonfinite source coordinate')
    return {'symbols':symbols,'positions_A':positions,'types':types,'bonds_zero_based':bonds}

def partition(atoms):
    graph=atoms['bonds_zero_based'];require(len(graph)==28,'Graph count differs')
    for i,neighbors in enumerate(graph):
        require(len(set(neighbors))==len(neighbors) and all(0<=j<28 and j!=i and i in graph[j] for j in neighbors),'Invalid or asymmetric graph')
    unseen=set(range(28));groups=[]
    while unseen:
        todo=[min(unseen)];group=[]
        while todo:
            i=todo.pop()
            if i not in unseen:continue
            unseen.remove(i);group.append(i);todo.extend(graph[i])
        groups.append(sorted(group))
    require(groups[0]==[0] and len(groups)==10,'Expected one unbonded metal and nine waters')
    waters=[]
    for group in groups[1:]:
        require(len(group)==3 and sorted(atoms['symbols'][i] for i in group)==['H','H','O'],'Nonwater fragment')
        oxygen=next(i for i in group if atoms['symbols'][i]=='O');hydrogens=[i for i in group if i!=oxygen]
        require(atoms['types'][oxygen]==36 and all(atoms['types'][i]==37 for i in hydrogens) and set(graph[oxygen])==set(hydrogens) and all(graph[h]==[oxygen] for h in hydrogens),'Water graph/type differs')
        waters.append({'indices':group,'oxygen_index':oxygen,'hydrogen_indices':hydrogens})
    return waters

def contacts(x,waters):
    membership={i:j for j,w in enumerate(waters) for i in w['indices']};oi=[w['oxygen_index'] for w in waters];hi=[h for w in waters for h in w['hydrogen_indices']]
    def minimum(pairs):
        values=[(distance(x[i],x[j]),i,j) for i,j in pairs];a=min(values)
        return {'distance_A':a[0],'atom_indices_zero_based':[a[1],a[2]]}
    return {'metal_oxygen_minimum':minimum([(0,i) for i in oi]),'metal_oxygen_maximum_A':max(distance(x[0],x[i]) for i in oi),
        'metal_hydrogen_minimum':minimum([(0,i) for i in hi]),
        'interwater_oxygen_minimum':minimum([(i,j) for a,i in enumerate(oi) for j in oi[a+1:]]),
        'interwater_any_atom_minimum':minimum([(i,j) for i in range(1,28) for j in range(i+1,28) if membership[i]!=membership[j]])}

def write_xyz(directory,gid,metal,positions,atoms):
    directory.mkdir(parents=True,exist_ok=False);symbols=atoms['symbols'];ordinary=directory/'geometry.xyz';typed=directory/'system.xyz'
    def xyz(x):return ' '.join(f'{v:.12f}' for v in x)
    ordinary.write_text('28\n'+gid+'\n'+''.join(symbol+' '+xyz(x)+'\n' for symbol,x in zip(symbols,positions)))
    typed.write_text('28  '+gid+'\n'+''.join(f'{i+1:5d} {symbol:<3s} '+xyz(x)+f' {atoms["types"][i]:5d}'+''.join(f' {j+1:5d}' for j in atoms['bonds_zero_based'][i])+'\n' for i,(symbol,x) in enumerate(zip(symbols,positions))))
    parsed=parse_typed(typed);plain=ordinary.read_text().splitlines()[2:]
    require(parsed['positions_A']==positions and parsed['symbols']==symbols and parsed['types']==atoms['types'] and parsed['bonds_zero_based']==atoms['bonds_zero_based'],'Written typed geometry differs')
    require([[float(v) for v in line.split()[1:]] for line in plain]==positions and [line.split()[0] for line in plain]==symbols,'Written ordinary geometry differs')
    return {'ordinary_xyz':pin(ordinary),'typed_xyz':pin(typed)}

def prepare(config_path,out):
    start,cpu=time.monotonic(),time.process_time();cfg=json.loads(Path(config_path).read_text());require(cfg==default_config(),'Exact declared geometry config required');require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread required')
    source=json.loads(verify(cfg['original_typed_manifest']).read_text());source={r['id']:r for r in source['rows']};original={};checks=[]
    def check(ok,name):
        checks.append({'name':name,'passed':bool(ok)});require(ok,name)
    for metal in METALS:
        gid='aquo9_central_'+metal;r=source[gid];atoms=parse_typed(verify(r['xyz']));g=r['source_geometry'];plain=verify(g['xyz']).read_text().splitlines();verify(cfg['admitted17670_typed_xyz'][metal])
        check(r['xyz']['sha256']==cfg['admitted17670_typed_xyz'][metal]['sha256'],metal+' original XYZ equals admitted17670 input bytes')
        check(atoms['positions_A']==g['positions_A'] and atoms['symbols']==g['symbols'] and atoms['types']==r['types'] and atoms['bonds_zero_based']==r['bonds_zero_based'],metal+' original geometry/type/graph metadata exactly replay')
        check(len(plain)==30 and int(plain[0])==28 and [[float(v) for v in line.split()[1:]] for line in plain[2:]]==g['positions_A'] and [line.split()[0] for line in plain[2:]]==g['symbols'],metal+' source ordinary XYZ matches typed XYZ and metadata')
        check(g['charge']==3 and g['physical_multiplicity']==(1 if metal=='La' else 6) and g['metal_index']==0 and g['atoms']==28 and g['water_count']==9 and atoms['types'][0]==(500 if metal=='La' else 501),metal+' exact physical state')
        waters=partition(atoms);original[metal]={'row':r,'geometry':g,'atoms':atoms,'waters':waters}
    check(original['La']['atoms']['positions_A']==original['Dy']['atoms']['positions_A'],'Central La/Dy coordinates identical')
    check(original['La']['waters']==original['Dy']['waters'] and original['La']['atoms']['bonds_zero_based']==original['Dy']['atoms']['bonds_zero_based'],'Central La/Dy graphs identical')
    x=original['La']['atoms']['positions_A'];waters=original['La']['waters'];directions=[]
    for w in waters:
        vector=[a-b for a,b in zip(x[w['oxygen_index']],x[0])];radius=math.sqrt(sum(v*v for v in vector));require(radius>0,'Zero metal-O distance');directions.append([v/radius for v in vector])
    source_contacts=contacts(x,waters);rows={}
    for sign,step in SIGNS:
        ideal=copy.deepcopy(x)
        for water,u in zip(waters,directions):
            for i in water['indices']:ideal[i]=[a+step*v for a,v in zip(x[i],u)]
        positions=[[float(f'{v:.12f}') for v in p] for p in ideal]
        check(positions[0]==x[0],sign+' metal fixed exactly')
        maximum_internal=0.;maximum_translation=0.;maximum_radius=0.
        for water,u in zip(waters,directions):
            for i in water['indices']:maximum_translation=max(maximum_translation,max(abs(positions[i][k]-x[i][k]-step*u[k]) for k in range(3)))
            for a,i in enumerate(water['indices']):
                for j in water['indices'][a+1:]:maximum_internal=max(maximum_internal,abs(distance(positions[i],positions[j])-distance(x[i],x[j])))
            oi=water['oxygen_index'];maximum_radius=max(maximum_radius,abs(distance(positions[oi],positions[0])-distance(x[oi],x[0])-step))
        check(maximum_internal<=TOLERANCE_A,sign+' all27 intrawater pair distances preserved')
        check(maximum_translation<=TOLERANCE_A,sign+' all intact waters receive exact declared translation')
        check(maximum_radius<=TOLERANCE_A,sign+' all nine O radii change by declared step')
        target_contacts=contacts(positions,waters)
        check(all(v['distance_A']>0 for v in target_contacts.values() if isinstance(v,dict)),sign+' no coincident distinct fragments')
        for metal in METALS:
            gid=f'aquo9_{sign}_{metal}';old=original[metal];files=write_xyz(out/'cells'/gid,gid,metal,positions,old['atoms']);g=copy.deepcopy(old['geometry'])
            g.update(geometry_id=gid,distortion=sign,positions_A=copy.deepcopy(positions),xyz=files['ordinary_xyz'],kind='fresh_coordination_radial_validation',radial_step_A=step,water_radial_directions=copy.deepcopy(directions),water_groups=copy.deepcopy(waters),source_geometry=copy.deepcopy(old['geometry']),source_central_id='aquo9_central_'+metal)
            inp={k:copy.deepcopy(old['row'][k]) for k in ('metal','case_id','charge','atoms','types','bonds_zero_based')};inp.update(id=gid,kind=g['kind'],xyz=files['typed_xyz'],source_geometry=g,derivative_test=False)
            row={'geometry':g,'input':inp,'bodies':[{'indices':[0],'category':'ion'}]+[{'indices':w['indices'],'category':'water'} for w in waters],
                'source':{'central_id':g['source_central_id'],'central_typed_xyz':old['row']['xyz'],'central_ordinary_xyz':old['geometry']['xyz'],'admitted17670_typed_xyz':cfg['admitted17670_typed_xyz'][metal]},
                'geometry_checks':{'maximum_internal_pair_distance_change_A':maximum_internal,'maximum_translation_component_error_A':maximum_translation,'maximum_radial_step_error_A':maximum_radius,'coordinate_quantization_A':1e-12,'graph_unchanged':True,'metal_fixed':True,'water_orientations_unchanged_to_coordinate_precision':True},
                'contacts':{'central':source_contacts,'displaced':target_contacts,'policy':'Geometry report only; no predicted energy or contact-based fixture selection.'}}
            rows[gid]=row
        check(rows[f'aquo9_{sign}_La']['geometry']['positions_A']==rows[f'aquo9_{sign}_Dy']['geometry']['positions_A'],sign+' prepared La/Dy coordinates identical')
    for record in [cfg['implementation'],cfg['root_declaration'],cfg['original_typed_manifest'],*cfg['admitted17670_typed_xyz'].values()]:verify(record)
    manifest={'status':'geometry_metadata_prepared','implementation':pin(__file__),'configuration':pin(config_path),'root_declaration':cfg['root_declaration'],'source_manifest':cfg['original_typed_manifest'],'rows':rows,'cell_order':list(rows),'checks':checks,'checks_passed':len(checks),'count':4,'atoms_per_cell':28,'water_count_per_cell':9,'candidate_binding':None,'reference_binding':None,'fit_trajectory_or_outcomes_read':False,'validation_energies_or_forces_read':False,'new_native_calls':0,'eligible_for_molecular_execution':False,'required_gate':'Terminal17676 all joint-development controls pass plus independent selected-ten actual review, then root source/config release.','wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}
    save(out/'MANIFEST.json',manifest);return manifest

def main():
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False)
    try:r=prepare(a.config,out)
    except Exception as exc:
        save(out/'FAILURE.json',{'status':'geometry_preparation_failed','error':repr(exc),'implementation':pin(__file__),'native_calls':0});raise
    print(json.dumps({'status':r['status'],'manifest':pin(out/'MANIFEST.json'),'checks':r['checks_passed'],'wall_seconds':r['wall_seconds']}))

if __name__=='__main__':main()
