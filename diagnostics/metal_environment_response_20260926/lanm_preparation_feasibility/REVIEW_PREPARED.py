"""Independent read-only checks of prepared Hans EF3 atoms, field and motion."""
import argparse,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new
p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
c=read_json(args.inputs);ex=read_json(verify(c['export']));native=read_json(verify(ex['atoms']));byid={a['id']:a for a in native}
expect=read_json(Path(__file__).with_name('BOUNDARY_EXPECTATIONS.json'))
env=read_json(verify(c['environment_atoms']));ed={a['id']:a for a in env};ba=read_json(verify(c['boundary_mapping']))
expectedq={a['id']:a['charge_e']for a in native}
for row in expect['local_boundary_expectations']:
 for ident,charge in row['expected_shifted_recipients'].items():expectedq[ident]=charge
for ident in ba['all_selected_QM_source_ids']+ba['removed_MM1_CA_ids']:expectedq.pop(ident)
for a in ba['spectators']:expectedq[a['id']]=a['charge_e']
assert set(expectedq)==set(ed) and len(ed)==1696
charge_error=max(abs(expectedq[k]-ed[k]['charge_e'])for k in expectedq)
assert charge_error<1e-13
assert abs(sum(a['charge_e']for a in env)-7)<1e-12
assert {(a['site'],a['element'],a['charge_e'])for a in ba['spectators']}=={('EF1','La',3.),('EF2','La',3.),('EF4','Na',1.)}
for a in env:
 if a['id']in byid:assert a['xyz_A']==byid[a['id']]['xyz_A']
for a in ba['spectators']:assert ed[a['id']]['xyz_A']==a['xyz_A']
for config in ('A','B'):
 pc=verify(c['configurations'][config]['pointcharges'])
 assert pc==verify(c['configurations']['A']['pointcharges'])
 values=np.loadtxt(pc,skiprows=1);assert values.shape==(1696,4)
 np.testing.assert_allclose(values,np.array([[a['charge_e'],*a['xyz_A']]for a in env]),rtol=0,atol=2e-13)
ma,mb=[read_json(verify(c['configurations'][k]['core_mapping']))for k in ('A','B')]
assert len(ma)==len(mb)==195
changed=[a['id']for a,b in zip(ma,mb)if a['xyz_A']!=b['xyz_A']]
assert set(changed)==set(c['perturbation']['changed_source_ids'])=={'A/85//CG','A/85//OD1','A/85//OD2','A/85//HB2','A/85//HB3'}
source_proposal=read_json(Path(__file__).with_name('RESULT_v2.json'))['sources'][0]['proposal']
assert {a['id']for a in ma if a['kind']=='source'}=={a['id']for a in source_proposal['source_atoms']}
axis=np.array(c['perturbation']['axis_unit']);pivot=np.array(c['perturbation']['pivot_xyz_A']);ang=np.deg2rad(2.)
max_mode=0.;max_jac=0.
for a,b in zip(ma,mb):
 assert a['id']==b['id']
 if a['kind']=='source':
  np.testing.assert_array_equal(a['xyz_A'],byid[a['id']]['xyz_A'])
  v=np.array(a['xyz_A'])-pivot
  expected=pivot+v*np.cos(ang)+np.cross(axis,v)*np.sin(ang)+axis*(axis@v)*(1-np.cos(ang)) if a['id']in changed else a['xyz_A']
  max_mode=max(max_mode,float(np.max(abs(np.asarray(expected)-b['xyz_A']))))
 if a['kind']=='cap':
  assert a==b
  rr=np.asarray(a['retained_xyz_A']);ro=np.asarray(a['omitted_xyz_A']);d=ro-rr;norm=np.linalg.norm(d);u=d/norm
  np.testing.assert_array_equal(rr,byid[a['retained_source_id']]['xyz_A']);np.testing.assert_array_equal(ro,byid[a['omitted_source_id']]['xyz_A'])
  np.testing.assert_allclose(rr+a['length_A']*u,a['xyz_A'],atol=1e-13,rtol=0)
  j=a['length_A']/norm*(np.eye(3)-np.outer(u,u));max_jac=max(max_jac,float(np.max(abs(j-a['jacobian_omitted']))))
  np.testing.assert_allclose(j+np.array(a['jacobian_retained']),np.eye(3),atol=1e-14,rtol=0)
assert max_mode<1e-13 and max_jac<1e-14
xyz_residual=0.
for config,mapping in [('A',ma),('B',mb)]:
 allrows={}
 for metal in ('La','Dy'):
  task=c['configurations'][config]['endpoints'][metal];lines=verify(task['xyz']).read_text().splitlines();rows=[line.split()for line in lines[2:]]
  assert int(lines[0])==len(rows)==195 and rows[0][0]==metal
  coords=np.array([[float(z)for z in r[1:]]for r in rows]);xyz_residual=max(xyz_residual,float(np.max(abs(coords-np.array([a['xyz_A']for a in mapping])))))
  assert (task['all_electron_count']-(task['multiplicity']-1))%2==0
  allrows[metal]=coords
 np.testing.assert_array_equal(allrows['La'],allrows['Dy'])
assert xyz_residual<2e-10
result=dict(status='independent_prepared_boundary_geometry_review_pass',inputs=record(args.inputs),review_implementation=record(__file__),source_real_nonmetal_atoms=190,QM_atoms=195,field_rows=len(env),field_charge_e=sum(a['charge_e']for a in env),full_charge_e=sum(a['charge_e']for a in env)-1,maximum_field_charge_difference_e=charge_error,changed_source_ids=changed,maximum_mode_coordinate_difference_A=max_mode,maximum_cap_jacobian_difference=max_jac,maximum_XYZ_serialization_difference_A=xyz_residual,identical_field_all_endpoints=True,physical_state_qualification='electron parity only; actual Dy Hamiltonian/spin localization not evaluated',energy_force_calls=0)
write_new(args.output,result);print(result)
