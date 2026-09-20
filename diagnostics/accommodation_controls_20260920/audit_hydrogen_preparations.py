#!/usr/bin/env python3
"""Coordinate/topology-only audit of saved failures plus four fixed controls."""
import argparse,json
from collections import Counter
from pathlib import Path
import numpy as np
from openmm.app import PDBFile
from audit_failed_fold_geometry import atoms,overlaps,closest,pin,distance

p=argparse.ArgumentParser()
p.add_argument('--prepared',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
failures={p.stem:json.loads(p.read_text())['reason'] for p in sorted((a.prepared/'failures').glob('*.json'))}
controls=['mmol_1770-pqq-la_model__conditioned_La__seed-1_sample-0',
          'mmol_1770-pqq-la_model__conditioned_La__seed-1_sample-3',
          'a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-0',
          'a0a3f2yly8-pqq-la_model__conditioned_La__seed-1_sample-0']
result={'scope':'all_saved_failures_plus_four_prespecified_successes_coordinate_only',
        'failure_categories':dict(Counter(failures.values())),
        'control_selection':'same-protein La sample0 and canonical sample3, first alphabetical panel member Ca0 and La0',
        'rows':[]}
for case in list(failures)+controls:
    path=a.prepared/'cases'/case/'protonated.pdb';rows=atoms(path);pairs=overlaps(rows)
    heavy=[r for r in rows if r['element'] not in ('H','D')]
    metal=next(r for r in rows if r['element'] in ('La','Ca'))
    pdb=PDBFile(str(path));aa=list(pdb.topology.atoms());neighbors=[[] for r in aa]
    for x,y in pdb.topology.bonds():neighbors[x.index].append(y.index);neighbors[y.index].append(x.index)
    assert len(rows)==len(aa)
    hidx=[r.index for r in aa if r.element.symbol=='H'];bdist=[];bad=[]
    for i in hidx:
        ns=[j for j in neighbors[i] if aa[j].element.symbol!='H']
        if len(ns)!=1:bad.append({'atom':rows[i],'bonded_heavy':len(ns)})
        else:bdist.append(distance(rows[i],rows[ns[0]]))
    item={'case_id':case,'status':'failure' if case in failures else 'prepared_success',
          'failure':failures.get(case),'protonated':pin(path),'atom_count':len(rows),'H_count':len(hidx),
          'pairs_below_0p45A':len(pairs),'pair_element_counts':dict(Counter('/'.join(sorted([x['a']['element'],x['b']['element']])) for x in pairs)),
          'shortest_pairs':pairs[:5],'H_without_exactly_one_heavy_bond':bad,
          'H_parent_bond_distance_percentiles_A':np.percentile(bdist,[0,5,50,95,100]).tolist(),
          'H_parent_distance_within_0p002A_of_1A':sum(abs(d-1)<.002 for d in bdist),
          'minimum_PQQ_distance_A':closest(metal,[r for r in heavy if r['resname']=='PQQ'])[0]['distance_A'],
          'nearest_protein_ON':closest(metal,[r for r in heavy if r['chain']=='A' and r['element'] in ('O','N')],3)}
    result['rows'].append(item)
with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(a.output)
