"""Hydrogen and fixed-core integrity checks on the actual prepared atoms."""
import collections
import math

def check(fragments):
    atoms=[(f['id'],a) for f in fragments for a in f['atom_records']]
    lookup={(rid,a['name']):a for rid,a in atoms};errors=[];parent_distances=[];hh=[];all_distances=[]
    if len(lookup)!=len(atoms):errors.append('duplicate core atom identity')
    for rid,a in atoms:
        if not all(math.isfinite(x) for x in a['xyz_A']):errors.append(f'nonfinite:{rid}:{a["name"]}')
        if a['origin']=='source_protonation_hydrogen':
            p=lookup.get((rid,a['parent_atom']))
            if p is None:errors.append(f'unknown H parent:{rid}:{a["name"]}');continue
            d=math.dist(a['xyz_A'],p['xyz_A']);parent_distances.append(d)
            if not .8<=d<=1.35:errors.append(f'H-parent outside0.8–1.35A:{rid}:{a["name"]}:{d}')
    for i,(rid,a) in enumerate(atoms):
        for sid,b in atoms[i+1:]:
            d=math.dist(a['xyz_A'],b['xyz_A']);all_distances.append(d)
            if a['element']==b['element']=='H':hh.append(d)
            if d<.5:errors.append(f'overlap:{rid}:{a["name"]}/{sid}:{b["name"]}:{d}')
    return {'status':'FAIL' if errors else 'PASS','errors':errors,'nonmetal_atom_count':len(atoms),
      'total_core_atom_count':len(atoms)+1,'formula_excluding_metal':dict(collections.Counter(a['element'] for _,a in atoms)),
      'minimum_atom_distance_A':min(all_distances),'minimum_HH_distance_A':min(hh),
      'source_H_parent_distance_range_A':[min(parent_distances),max(parent_distances)],
      'minimum_allowed_atom_distance_A':.5,'source_H_parent_allowed_range_A':[.8,1.35]}

def unchanged_non_H(old,new):
    aa=[(f['id'],a) for f in old for a in f['atom_records']];bb=[(f['id'],a) for f in new for a in f['atom_records']]
    ident=lambda x:(x[0],x[1]['name'],x[1]['element'],x[1]['origin'],x[1]['parent_atom'])
    if [ident(x) for x in aa]!=[ident(x) for x in bb]:raise ValueError('Core atom identities/order changed')
    changes=0
    for (rid,a),(_,b) in zip(aa,bb):
        if a['xyz_A']!=b['xyz_A']:
            if a['origin']!='source_protonation_hydrogen':raise ValueError('Frozen heavy/PQQ-H/cap coordinate changed')
            changes+=1
    if [f['formal_charge'] for f in old]!=[f['formal_charge'] for f in new]:raise ValueError('Fragment protonation changed')
    return changes
