"""Fail-closed checks for the approved, H-only ADH9 preparation repair."""
import collections
import math

def validate(original_fragments, repaired_fragments, minimization):
    old=[(f['id'],a) for f in original_fragments for a in f['atom_records']]
    new=[(f['id'],a) for f in repaired_fragments for a in f['atom_records']]
    if len(old)!=77 or len(new)!=77:raise ValueError('Expected77 nonmetal core atoms')
    identity=lambda x:(x[0],x[1]['name'],x[1]['element'],x[1]['origin'],x[1]['parent_atom'])
    if [identity(x) for x in old]!=[identity(x) for x in new]:raise ValueError('Core atom identity/order/state changed')
    lookup={(rid,a['name']):a for rid,a in new};changes=0;parent_distances=[]
    for (rid,a),(_,b) in zip(old,new):
        if not all(math.isfinite(v) for v in b['xyz_A']):raise ValueError('Nonfinite coordinates')
        allowed=a['origin']=='source_protonation_hydrogen'
        if not allowed and a['xyz_A']!=b['xyz_A']:raise ValueError('Frozen heavy/PQQ/cap atom moved')
        if allowed:
            changes+=a['xyz_A']!=b['xyz_A']
            d=math.dist(b['xyz_A'],lookup[rid,b['parent_atom']]['xyz_A'])
            if not .8<=d<=1.35:raise ValueError(f'Invalid hydrogen-parent distance:{rid}:{b["name"]}={d}')
            parent_distances.append(d)
    for f,g in zip(original_fragments,repaired_fragments):
        if (f['id'],f['formal_charge'],f['atom_count'],f['role'])!=(g['id'],g['formal_charge'],g['atom_count'],g['role']):raise ValueError('Fragment identity/charge changed')
    distances=[];hh=[]
    for i,(rid,a) in enumerate(new):
        for sid,b in new[i+1:]:
            d=math.dist(a['xyz_A'],b['xyz_A']);distances.append(d)
            if a['element']==b['element']=='H':hh.append(d)
            if d<.5:raise ValueError(f'Overlapping atoms:{rid}:{a["name"]}/{sid}:{b["name"]}={d}')
    formula=dict(collections.Counter(a['element'] for _,a in new))
    if formula!={'C':27,'H':31,'N':4,'O':15}:raise ValueError('Core formula changed')
    if minimization['status']!='PASS' or minimization['final']['hydrogen_force_rms_kJ_mol_nm']>1.0:raise ValueError('Hydrogen minimizer has not converged')
    if minimization['heavy_displacement_A']!=0:raise ValueError('Minimizer moved heavy atoms')
    return {'status':'PASS','nonmetal_core_atom_count':77,'total_core_atom_count':78,
       'formula_excluding_metal':formula,'changed_source_hydrogens':changes,'core_atom_identity_and_order_preserved':True,
       'all_heavy_PQQ_hydrogen_and_cap_coordinates_unchanged':True,
       'minimum_core_atom_distance_A':min(distances),'minimum_core_HH_distance_A':min(hh),
       'source_H_parent_distance_range_A':[min(parent_distances),max(parent_distances)],
       'minimum_allowed_atom_distance_A':.5,'source_H_parent_allowed_range_A':[.8,1.35],
       'hydrogen_minimization':{k:v for k,v in minimization.items() if k!='trace'}}
