"""Pinned saved-coordinate q9 audit. File/geometry inspection only; no energies."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','OMP_THREAD_LIMIT']:os.environ[k]='1'
import csv, hashlib, json, math, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from coordination_policy import donor_type, COORDINATION_CUTOFF_A, DIAGNOSTIC_CUTOFF_A
from affordable_peptide import SourceGraph
PREP=ROOT/'workspaces/orca_parameter_audit_20261012/context_factorial_preparation_v1/PREPARATION.json'
PINS={}

def pin(p):
    p=Path(p).resolve();v={'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()};PINS[str(p)]=v;return v
def load(p):
    if isinstance(p,dict):assert pin(p['path'])['sha256']==p['sha256'];p=p['path']
    else:pin(p)
    return json.loads(Path(p).read_text())
def xyz(p):
    if isinstance(p,dict):assert pin(p['path'])['sha256']==p['sha256'];p=p['path']
    else:pin(p)
    rows=Path(p).read_text().splitlines();a=[l.split() for l in rows[2:] if l.strip()];assert len(a)==int(rows[0])
    return [x[0] for x in a],np.array([[float(y) for y in x[1:4]] for x in a])
def key(k):return json.dumps(k,separators=(',',':'))
def residue(a):return (a['chain'],a['resnum'],a['icode'],a['resname'])
def label(a):return f"{a['chain']}:{a['resname']}{a['resnum']}{a['icode']}:{a['atom']}"
def distance_bin(d):return '<=3.5' if d<=3.5 else '3.5-4.0' if d<=4 else '>4.0'
def angle(a,h,c):
    u=a-h;v=c-h;return float(np.degrees(np.arccos(np.clip(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)),-1,1))))
def dihedral(a,b,c,d):
    b0=a-b;b1=c-b;b2=d-c;b1=b1/np.linalg.norm(b1)
    v=b0-np.dot(b0,b1)*b1;w=b2-np.dot(b2,b1)*b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1,v),w),np.dot(v,w))))
def write(name,value):
    p=OUT/name
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def main():
    prep=load(PREP);case=next(c for c in prep['cases'] if c['protein']=='q9z4j7-pqq-la_model')
    parent=case['source'];assert pin(parent['path'])['sha256']==parent['sha256']
    atoms=[];positions=[]
    for line in Path(parent['path']).read_text().splitlines():
        if line[:6].strip() not in ['ATOM','HETATM']:continue
        a={'chain':line[21].strip(),'resnum':int(line[22:26]),'icode':line[26].strip(),'resname':line[17:20].strip(),
           'atom':line[12:16].strip(),'altloc':line[16].strip(),'element':line[76:78].strip(),'serial':int(line[6:11])}
        a['id']=['source',a['chain'],a['resnum'],a['icode'],a['atom'],a['altloc'],a['element']]
        atoms.append(a);positions.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
    pos=np.array(positions);lookup={key(a['id']):i for i,a in enumerate(atoms)};assert len(lookup)==len(atoms)
    graph_config=load(ROOT/'workspaces/accommodation_goal_20260920/folds_v1/preparation_reconciled_v1.json')['config']['topology']
    assert pin(graph_config['path'])['sha256']==graph_config['sha256']
    graph=SourceGraph(parent['path'],graph_config['path'])
    graph_to_parent={k:lookup[key(['source',m['chain'],m['resnum'],m.get('insertion_code',''),m['atom'],m.get('altloc',''),m['element']])] for k,m in graph.meta.items()}
    adjacency={i:set() for i in range(len(atoms))}
    for a,b in list(graph.edges)+list(graph.hparents.items()):
        i,j=graph_to_parent[a],graph_to_parent[b];adjacency[i].add(j);adjacency[j].add(i)
    def bond_separation(i,j,limit=3):
        seen={i};front={i}
        for step in range(1,limit+1):
            front=set().union(*(adjacency[k] for k in front))-seen
            if j in front:return step
            seen|=front
        return None
    metal=next(i for i,a in enumerate(atoms) if a['element'].lower()=='la')
    included={};qmaps={};origin_qm={};metadata={}
    for context,c in case['contexts'].items():
        included[context]=set();qmaps[context]={key(a['id']):a['index'] for a in c['source_id_to_qm']}
        sy,coords=xyz(c['origin']);origin_qm[context]=(sy,coords)
        m=load(c['mapping_pins']['adaptive_La']);p=load(c['proposal_receipts']['adaptive_La'])
        metadata[context]={'origin':c['origin'],'atom_count':c['atom_count'],'cap_count':c['cap_count'],
            'selection':c['selection'],'final_geometry':p['final_geometry'],'boundary_flag':p['boundary_flag'],
            'residual':p['residual'],'optimizer':p['optimizer'],'proposal_receipt':c['proposal_receipts']['adaptive_La']}
        for k,qi in qmaps[context].items():
            identity=json.loads(k)
            if identity[0]=='source':
                pi=lookup[k];included[context].add(pi);assert np.max(np.abs(coords[qi]-pos[pi]))<=1e-7
            elif identity[0]=='opaque_cofactor':
                matches=[i for i,a in enumerate(atoms) if a['resname']=='PQQ' and a['atom']==identity[2] and a['element']==identity[3]]
                if matches:
                    assert len(matches)==1;pi=matches[0];assert np.max(np.abs(coords[qi]-pos[pi]))<=1e-7;included[context].add(pi)
            elif identity[0]=='metal':included[context].add(metal)
    states={'origin':pos.copy()};source_rows={}
    for geometry in ['original_adaptive_La','expanded_adaptive_La']:
        rows=[r for r in prep['rows'] if r['protein']==case['protein'] and r['geometry_id']==geometry]
        mapped=[]
        for row in rows:
            context=row['target_context'];sy,coords=xyz(row['endpoints']['La']['xyz']);p=pos.copy()
            for ident,qi in qmaps[context].items():
                if json.loads(ident)[0]=='source':p[lookup[ident]]=coords[qi]
            mapped.append(p)
        assert np.array_equal(mapped[0],mapped[1]);states[geometry]=mapped[0];source_rows[geometry]=rows
    moved=sorted(set(i for s in states.values() for i in np.where(np.linalg.norm(s-pos,axis=1)>1e-8)[0]))
    moved_heavy=[i for i in moved if atoms[i]['element']!='H']
    def record_atom(i):
        return {'id':atoms[i]['id'],'label':label(atoms[i]),'element':atoms[i]['element'],
                'included_original':i in included['original'],'included_expanded':i in included['expanded']}
    movements=[{**record_atom(i),'displacement_A':{g:float(np.linalg.norm(s[i]-pos[i])) for g,s in states.items()}} for i in moved]
    coordination=[]
    for i,a in enumerate(atoms):
        typ=donor_type(a['resname'],a['atom'],a['element'])
        if not typ:continue
        ds={g:float(np.linalg.norm(s[i]-s[metal])) for g,s in states.items()}
        if min(ds.values())<=DIAGNOSTIC_CUTOFF_A:
            coordination.append({**record_atom(i),'donor_type':typ,'distances_A':ds,'direct':{g:d<=COORDINATION_CUTOFF_A for g,d in ds.items()}})
    cn={g:sum(r['direct'][g] for r in coordination) for g in states}
    denticity={g:{} for g in states}
    for r in coordination:
        group=':'.join(r['label'].split(':')[:-1])
        for g in states:
            if r['direct'][g]:denticity[g][group]=denticity[g].get(group,0)+1
    torsion_specs={('A',213,'chi1'):['N','CA','CB','CG'],('A',213,'chi2'):['CA','CB','CG','CD'],('A',213,'chi3'):['CB','CG','CD','OE1'],
                   ('A',350,'chi1'):['N','CA','CB','CG'],('A',350,'chi2'):['CA','CB','CG','OD1'],('A',300,'chi1'):['N','CA','CB','CG'],('A',300,'chi2'):['CA','CB','CG','OD1']}
    torsions=[]
    for (chain,resnum,name),names in torsion_specs.items():
        ids=[next(i for i,a in enumerate(atoms) if a['chain']==chain and a['resnum']==resnum and a['atom']==n) for n in names]
        angles={g:dihedral(*(s[i] for i in ids)) for g,s in states.items()}
        torsions.append({'id':f'{chain}/{resnum}/{name}','atoms':[label(atoms[i]) for i in ids],'angles_degrees':angles,
            'changes_degrees':{g:(a-angles['origin']+180)%360-180 for g,a in angles.items()},
            'declared_rotation_degrees':{g:float(np.degrees(source_rows[g][0]['physical_q'][source_rows[g][0]['mode_ids'].index(f'{chain}/{resnum}/{name}')])) for g in source_rows}})
    heavy=[i for i,a in enumerate(atoms) if a['element'] not in ['H','La','LA','Ca']]
    contacts=[];seen=set()
    for i in moved_heavy:
        for j in heavy:
            if i==j or residue(atoms[i])==residue(atoms[j]) or tuple(sorted([i,j])) in seen:continue
            seen.add(tuple(sorted([i,j])));ds={g:float(np.linalg.norm(s[i]-s[j])) for g,s in states.items()}
            if min(ds.values())>4.:continue
            contacts.append({'moving':record_atom(i),'partner':record_atom(j),'distances_A':ds,'bins':{g:distance_bin(d) for g,d in ds.items()},
                'delta_A':{g:d-ds['origin'] for g,d in ds.items()},'either_partner_omitted_from_both':any(k not in included['original']|included['expanded'] for k in [i,j])})
    omitted=[]
    for i,a in enumerate(atoms):
        if i in included['original']|included['expanded']:continue
        nearest={g:min((float(np.linalg.norm(s[i]-s[j])),j) for j in moved) for g,s in states.items()}
        if min(x[0] for x in nearest.values())<=4:
            omitted.append({**record_atom(i),'resname':a['resname'],'backbone_atom':a['atom'] in ['N','H','CA','HA','C','O','OXT'],
                'typed_metal_donor':donor_type(a['resname'],a['atom'],a['element']),
                'nearest_moved':{g:{'distance_A':d,'atom':label(atoms[j])} for g,(d,j) in nearest.items()}})
    # Same declared H geometry bins, now using the actual existing source graph.
    dh=[];ambiguous=[]
    for gh,gd in graph.hparents.items():
        h,d=graph_to_parent[gh],graph_to_parent[gd]
        if atoms[d]['element'] in ['N','O'] and np.linalg.norm(pos[h]-pos[d])<=1.3:dh.append((d,h))
    protonated={d for d,h in dh}
    acceptors=[i for i,a in enumerate(atoms) if a['element']=='O' or (a['element']=='N' and i not in protonated and ((a['resname']=='HIS' and a['atom'] in ['ND1','NE2']) or(a['resname']=='PQQ' and a['atom']=='N2')))]
    hb=[]
    for d,h in dh:
        for a in acceptors:
            if a==d or bond_separation(a,d) in [1,2] or not(set([a,d,h])&set(moved)):continue
            measures={g:{'D_A_A':float(np.linalg.norm(s[d]-s[a])),'H_A_A':float(np.linalg.norm(s[h]-s[a])),
                         'D_H_A_degrees':angle(s[d],s[h],s[a])} for g,s in states.items()}
            if min(v['D_A_A'] for v in measures.values())>4:continue
            for v in measures.values():v['geometry_flag']=v['D_A_A']<=3.5 and v['H_A_A']<=2.5 and v['D_H_A_degrees']>=120
            hb.append({'donor':record_atom(d),'hydrogen':record_atom(h),'acceptor':record_atom(a),'geometry':measures})
    omitted_pairs=[]
    for i in moved:
        for j in range(len(atoms)):
            if j in included['original']|included['expanded']:continue
            ds={g:float(np.linalg.norm(s[i]-s[j])) for g,s in states.items()}
            if min(ds.values())>4:continue
            separation=bond_separation(i,j)
            omitted_pairs.append({'moving':record_atom(i),'omitted':record_atom(j),'distances_A':ds,
                'bins':{g:distance_bin(d) for g,d in ds.items()},'source_bond_separation_up_to3':separation,
                'excluded_bond_or_angle_pair':separation in [1,2],
                'one_four_pair':separation==3,'both_heavy':atoms[i]['element']!='H' and atoms[j]['element']!='H'})
    graft_bonds=[]
    for i in moved:
        for j in adjacency[i]:
            graft_bonds.append({'a':record_atom(i),'b':record_atom(j),'lengths_A':{g:float(np.linalg.norm(s[i]-s[j])) for g,s in states.items()}})
    max_bond_error=max(abs(d-r['lengths_A']['origin']) for r in graft_bonds for d in r['lengths_A'].values())
    assert max_bond_error<=1e-7
    summaries={g:{'coordination_number':cn[g],'denticity':denticity[g],
        'maximum_real_heavy_displacement_A':max(float(np.linalg.norm(s[i]-pos[i])) for i in moved_heavy),
        'contact_changes_at_3_5_A':{'new':[r for r in contacts if r['distances_A']['origin']>3.5 and r['distances_A'][g]<=3.5],
                                     'lost':[r for r in contacts if r['distances_A']['origin']<=3.5 and r['distances_A'][g]>3.5]},
        'hydrogen_bond_geometry_changes':{'new':[r for r in hb if not r['geometry']['origin']['geometry_flag'] and r['geometry'][g]['geometry_flag']],
                                           'lost':[r for r in hb if r['geometry']['origin']['geometry_flag'] and not r['geometry'][g]['geometry_flag']]}} for g,s in states.items()}
    result={'question':'Missing environment resistance versus local model response; geometric diagnosis only','source_parent':parent,'preparation':pin(PREP),
        'definitions':pin(OUT/'PLAN.md'),'typed_coordination_policy':pin(ROOT/'scripts/coordination_policy.py'),'implementation':pin(__file__),
        'parent_atom_count':len(atoms),'same_saved_motion_in_both_contexts_verified':True,'context_metadata':metadata,
        'movements':movements,'coordination':coordination,'torsions':torsions,'contacts':contacts,'omitted_source_atoms_near_movement':omitted,
        'hydrogen_bond_geometry':hb,'ambiguous_parent_hydrogens':ambiguous,'omitted_pairs':omitted_pairs,'source_graft_bond_checks':{'all_preserved':True,'maximum_bond_error_A':max_bond_error,'bonds':graft_bonds},'source_graph_topology':graph_config,'source_graph_implementation':pin(ROOT/'scripts/affordable_peptide.py'),'revision':'v2 includes same-residue backbone H-bond candidates and real source bond exclusions, with all prior distance/angle bins unchanged.','summaries':summaries,'pins':list(PINS.values()),
        'numerical_threads':1,'new_molecular_calls':0,'new_optimizations':0,'interpretation_limit':'Distances and contact flags identify hypotheses, not interaction energies or unique causes.'}
    write('RESULT_V2.json',result)
    print(json.dumps({'CN':cn,'moved_atoms':len(moved),'moved_heavy':len(moved_heavy),'contacts':len(contacts),'omitted_atoms_within4A':len(omitted),'hydrogen_bond_candidates':len(hb)}))

if __name__=='__main__':main()
