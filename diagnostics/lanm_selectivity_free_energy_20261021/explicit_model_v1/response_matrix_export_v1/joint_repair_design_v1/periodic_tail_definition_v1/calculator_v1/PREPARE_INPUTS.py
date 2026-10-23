"""Copy exactly the two pinned fixture metadata rows; no tail/geometry evaluation."""
from decimal import Decimal
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
J=HERE.parents[1]
INVENTORY=J/'periodic_fixture_handoff_v1/INVENTORY.json'
EXPECTED='ae407a66e4ec2041f4f28cb7f5f391acfc358f6b15b372457226d9183ad469a7'

def pin(path):
    return {'path':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def main():
    ip=pin(INVENTORY)
    assert ip['sha256']==EXPECTED
    inventory=json.loads(INVENTORY.read_text())
    assert inventory['new_native_calls']==0 and inventory['tail_evaluations']==0
    assert [r['id'] for r in inventory['fixtures']]==['La_aqueous_BIO','Dy_aqueous_BIO']
    rows=[]
    outdir=HERE/'inputs_v1'
    outdir.mkdir(exist_ok=True)
    for row in inventory['fixtures']:
        cell=row['cell_A_degrees'];bounds=row['lattice_bounds']
        assert len(cell)==6 and all(isinstance(s,str) for s in cell)
        assert all(Decimal(angle)==90 for angle in cell[3:])
        assert min(map(Decimal,cell[:3]))==Decimal(bounds['singular_lower_A'])
        assert max(map(Decimal,cell[:3]))==Decimal(bounds['singular_upper_A'])
        cfg={'schema':'orthorhombic_literal_tail_input_v1','box_id':row['id'],'precision':80,
             'cell_sides_A':cell[:3],'water_count':row['water_molecules'],
             'representatives':{'convention':'literal_primary_coordinates',
                 'displacement_upper_A':bounds['fixed_literal_M_to_O_D_upper_A'],
                 'squared_displacement_upper_A2':bounds['fixed_literal_M_to_O_D_squared_A2']},
             'authorities':{'inventory':ip,'geometry':row['geometry']},
             'scope':'geometry-only; older key/parameter/periodic boundary choices are not admitted',
             'real_box_evaluation_requires_root_release':True}
        path=outdir/(row['id']+'.json')
        with path.open('x') as f:json.dump(cfg,f,indent=2,sort_keys=True);f.write('\n')
        rows.append(pin(path))
    result={'status':'prepared_metadata_only','implementation':pin(Path(__file__)),
            'inventory':ip,'configurations':rows,'tail_evaluations':0,'native_calls':0,
            'geometry_files_read':0,'fit_outputs_read':0}
    with (outdir/'MANIFEST.json').open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(result))
if __name__=='__main__':main()
