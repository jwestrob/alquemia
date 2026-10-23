#!/usr/bin/env python3
"""Independent Decimal arithmetic from raw endpoints; no molecular execution."""
import argparse
from decimal import Decimal as D
import hashlib
import json
from pathlib import Path
import re


def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--result',required=True);ap.add_argument('--out',required=True)
    a=ap.parse_args();r=json.loads(Path(a.result).read_text());cfg=json.loads(Path(r['config']['path']).read_text())
    panel=json.loads(Path(cfg['prepared_inventory']['path']).read_text());geoms={g['geometry_id']:g for g in panel['geometries']}
    wet=json.loads(Path(json.loads(Path(cfg['wet_native_panel']['path']).read_text())['source_inventory']['path']).read_text())
    geoms.update({g['geometry_id']:g for g in wet['geometries']+wet['reused_controls']})
    checks=0;residuals=[];raw={m:{} for m in r['endpoints']}
    for p in r['reviewed_pins']:
        assert pin(p['path'])==p;checks+=1
    for method,ends in r['endpoints'].items():
        for gid,e in ends.items():
            src=json.loads(Path(e['result']['path']).read_text(),parse_float=D)
            if method=='ecp_pbe0':
                lines=[v.strip() for v in Path(src['engrad']['path']).read_text().splitlines() if v.strip() and not v.startswith('#')]
                raw[method][gid]=D(lines[1])*D('627.509474')
            elif method=='MACE':
                row=next(x for x in src['rows'] if x['id']==gid)
                raw[method][gid]=row['energy_eV']*D('23.060547830619')
                state=row['input_state'];assert state['charge']==e['charge'] and state['spin_multiplicity']==e['physical_multiplicity']
                assert state['atoms']==len(e['force_kcal_mol_A']) and state['graph_count']==1 and state['selected_metal_index']==0
                checks+=5
            else:
                row=next(x for x in src['rows'] if x['id']==gid)
                text=Path(row['native']['analyze_output']['path']).read_text()
                raw[method][gid]=D(re.search(r'Total Potential Energy\s*:\s*([-+0-9.Ee]+)',text)[1])
            residuals.append(abs(float(raw[method][gid])-e['energy_kcal_mol']));checks+=1
    for row in r['matched_water_addition_contrasts']:
        p=row['donor']+'_central_';v=raw[row['method']]
        independent=(v['hydrated_'+p+'Dy']-v['matched_dry_'+p+'Dy'])-(v['hydrated_'+p+'La']-v['matched_dry_'+p+'La'])
        residuals.append(abs(float(independent)-row['value_kcal_mol']));checks+=1
    radial_residual=[]
    for row in r['ligand_radial_response']:
        for context in ('wet','dry'):
            for metal in ('La','Dy'):
                gid=('hydrated_' if context=='wet' else 'matched_dry_')+row['donor']+'_central_'+metal
                g=geoms[gid];f=r['endpoints'][row['method']][gid]['force_kcal_mol_A']
                force_sum=[sum(f[i][axis] for i in g['ligand_indices']) for axis in range(3)]
                derivative=-sum(x*u for x,u in zip(force_sum,g['radial_direction']))
                radial_residual.append(abs(derivative-row[context+'_'+metal+'_kcal_mol_A']));checks+=1
    assert max(residuals)<2e-8 and max(radial_residual)<1e-10
    out={'status':'passed','checks':checks,'pins':len(r['reviewed_pins']),
         'max_energy_or_contrast_residual_kcal_mol':max(residuals),
         'max_ligand_radial_residual_kcal_mol_A':max(radial_residual),
         'source':pin(__file__),'result':pin(a.result),'new_molecular_calls':0,'local_threads':1}
    with Path(a.out).open('x') as h:json.dump(out,h,indent=2);h.write('\n')
    print(json.dumps(out))


if __name__=='__main__':main()
