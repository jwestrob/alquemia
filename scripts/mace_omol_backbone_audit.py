"""Independent peptide-connectivity geometry audit of real whole-chain inputs."""
from __future__ import annotations
import argparse
import math
from pathlib import Path
from affordable_common import read_json,record,verify,write_new

LIMITS={'minimum_C_N_A':1.0,'maximum_C_N_A':1.8}


def inspect(data):
    rows=[]
    for row in data['rows']:
        entry={'case_id':row['case_id'],'preparation':row.get('preparation')}
        if row['status']!='prepared':
            entry.update(status='unsupported',reason=row.get('reason','source preparation unavailable'),bonds=[])
        else:
            p=read_json(verify(row['preparation']))
            atoms={a['id']:a for a in p['physical_atoms']};bonds=[]
            for b in p['preparation_details']['bonds']:
                a,c=b['atom_a_id'],b['atom_b_id']
                if {a.split('/')[-1],c.split('/')[-1]}!={'C','N'} or a.rsplit('/',1)[0]==c.rsplit('/',1)[0]:continue
                distance=math.dist(atoms[a]['xyz_A'],atoms[c]['xyz_A'])
                bonds.append({'atom_a_id':a,'atom_b_id':c,'distance_A':distance,
                              'forcefield_equilibrium_A':b['forcefield_equilibrium_A'],
                              'pass':LIMITS['minimum_C_N_A']<=distance<=LIMITS['maximum_C_N_A']})
            bad=[b for b in bonds if not b['pass']]
            entry.update(status='pass' if not bad else 'unsupported',bonds=bonds,invalid_bonds=bad,
                         reason=None if not bad else 'inferred peptide connectivity crosses unresolved/malformed source geometry')
        rows.append(entry)
    return {'limits':LIMITS,'rows':rows,'pass':all(r['status']=='pass' for r in rows),
            'unsupported_cases':[r['case_id'] for r in rows if r['status']!='pass'],
            'scope':'peptide_connectivity_geometry; not exhaustive chemical validation',
            'new_energy_evaluations':0}


def audit(preparation,agreement,output):
    result=inspect(read_json(preparation))
    result.update(preparation_manifest=record(preparation),agreement=record(agreement),implementation=record(__file__))
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    write_new(out/'result.json',result)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('preparation','agreement','output'):parser.add_argument('--'+key,required=True)
    result=audit(**vars(parser.parse_args()))
    print(result['unsupported_cases'])
