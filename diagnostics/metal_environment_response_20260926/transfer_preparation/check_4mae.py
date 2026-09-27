"""Real prepared fixture and existing six-cell dry-run, no molecular execution."""
import argparse,json,tempfile,sys
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--repository',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();sys.path.insert(0,str(a.repository/'scripts'))
from affordable_common import read_json,verify,write_new
from metal_environment_reference import prepare,validate
from metal_environment_force_checks import check_pins,xyz_data,pc_data
x=read_json(a.inputs);check_pins(x);checks=[]
_,ca,qa=xyz_data(verify(x['endpoints']['Ca']['xyz']));_,la,ql=xyz_data(verify(x['endpoints']['La']['xyz']));assert np.array_equal(qa,ql) and ca[1:]==la[1:];assert len(qa)==80
nums={'H':1,'C':6,'N':7,'O':8,'Ca':20,'La':57};electrons={z:sum(nums[e] for e in symbols)-x['endpoints'][z]['charge'] for z,symbols in [('Ca',ca),('La',la)]};assert all(n%2==0 for n in electrons.values());checks.append({'check':'paired_80_atom_geometry_state_and_parity','status':'pass','physical_electrons':electrons,'explicit_La_electrons':electrons['La']-46})
ea,eb=[read_json(verify(x['environments'][k]['atoms'])) for k in ('A','B')];diff=[i for i,(u,v) in enumerate(zip(ea,eb)) if u!=v];assert diff==[x['perturbation']['pointcharge_index']];assert [v['charge_e'] for v in ea]==[v['charge_e'] for v in eb];assert abs(sum(v['charge_e'] for v in ea)+1)<1e-10;checks.append({'check':'only_declared_hydroxyl_H_moves_same_charges','status':'pass','changed_index':diff[0]})
for label,atoms in [('A',ea),('B',eb)]:
 _,q,y=pc_data(verify(x['environments'][label]['pointcharges']));assert np.array_equal(q,np.array([v['charge_e'] for v in atoms]));assert np.array_equal(y,np.array([v['xyz_A'] for v in atoms]));assert len(atoms)==8774
checks.append({'check':'field_serialization_exact_source_mapping','status':'pass'})
with tempfile.TemporaryDirectory() as d:
 r=prepare(a.inputs,a.repository/'diagnostics/metal_environment_response_20260926/TRANSFER_4MAE_PLAN.md',Path(d)/'six_cells',workers=6,mpi_ranks=57);validate(verify(r['manifest']));checks.append({'check':'existing_six_cell_runner_prepare_and_dry_run','status':'pass','molecular_calls':0})
write_new(a.output,{'status':'pass','checks':checks,'molecular_evaluations':0,'numerical_qualification':'not_established'});print(json.dumps(checks,indent=2))
