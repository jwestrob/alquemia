"""Export existing six-source geometry JSON; no molecular calculations."""
import csv,json,sys
from pathlib import Path
p=Path(sys.argv[1]); out=Path(sys.argv[2]); a=json.loads(p.read_text()); rows=[];works=[]
for r in a['rows']:
    for name,g in r['points'].items():
        for role,v in g['roles'].items():
            for atom in v['atoms']:
                rows.append(dict(case_id=r['case_id'],source_sample=r['source_evidence']['source_provenance']['source_sample'],source_conditioning_metal=r['source']['source_conditioning_metal'],known_biological_label='',geometry=name,selected_by_Ca=name==r['selected_candidates']['Ca'],selected_by_La=name==r['selected_candidates']['La'],role=role,chain=v['source_selector']['chain'],resnum=v['source_selector']['resnum'],atom=atom['atom'],source_atom_id=json.dumps(atom['source_id']),metal_O_distance_A=atom['distance_A'],atom_displacement_A=atom['displacement_A'],typed_role_contact_count_3p1A=v['typed_contacts_within_3p1_A'],maximum_heavy_displacement_A=g['maximum_heavy_displacement_A'],boundary=g['boundary_flag']))
    w=dict(case_id=r['case_id'],known_biological_label='',origin_R=r['origin_R'],selected_R=r['pool_R'],delta_R=r['pool_R']-r['origin_R'])
    for z in ('Ca','La'):
        w[z+'_candidate']=r['selected_candidates'][z]
        w.update({z+'_'+k:v for k,v in r['selected_work'][z].items()})
    w.update({'differential_'+k:v for k,v in r['differential_work'].items()}); works.append(w)
assert len(rows)==126 and len(works)==6
for name,records in [('DISTANCES.csv',rows),('COMPONENT_WORKS.csv',works)]:
    with (out/name).open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    assert len(list(csv.DictReader((out/name).open())))==len(records)
print('Exported all 126 donor-distance rows and all six component-work rows.')
