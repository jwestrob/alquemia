"""Real prepared/archived fixtures only; no synthetic scientific results."""
import copy,json,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import metal_environment_partition as p
from affordable_common import read_json,record,verify,InvalidArtifact
BASE=ROOT/'workspaces/metal_environment_response_20260926';INPUTS=BASE/'preparation/expanded_region_v1/INPUTS.json';SOURCE=BASE/'reference_scout_v1/FINAL_COLLECTION_1219207.json'
class PartitionTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  if not INPUTS.exists() or not SOURCE.exists():raise unittest.SkipTest('Actual prepared/computed fixtures absent')
  cls.config=read_json(INPUTS)
 def test_exact_endpoint_geometries_and_fields(self):
  x=p.validate_inputs(self.config)
  self.assertFalse((x['Ca_A']==x['Ca_B']).all());self.assertTrue((x['La_B']==x['Ca_B']).all())
 def test_wrong_B_geometry_rejected(self):
  c=copy.deepcopy(self.config);c['endpoints']['Ca_B']['xyz']=c['endpoints']['Ca_A']['xyz']
  with self.assertRaises(InvalidArtifact):p.validate_inputs(c)
 def test_state_corruption_rejected(self):
  c=copy.deepcopy(self.config);c['endpoints']['La_A']['charge']=-1
  with self.assertRaisesRegex(InvalidArtifact,'state'):p.validate_inputs(c)
 def test_actual_response_algebra_and_missing_not_zero(self):
  s=read_json(SOURCE);per,d=p.responses(s['rows']);self.assertAlmostEqual(d,s['delta_env_el_kcal_mol'],places=10)
  r=copy.deepcopy(s['rows']);r['La_B']['status']='unavailable';per,d=p.responses(r);self.assertIsNone(d);self.assertIsNone(per['La']);self.assertIsNotNone(per['Ca'])
 def test_prepare_validate_and_missing_matrix_cell(self):
  plan=ROOT/'diagnostics/metal_environment_response_20260926/PARTITION_DIAGNOSTIC_PLAN.md'
  with tempfile.TemporaryDirectory() as d:
   result=p.prepare(INPUTS,SOURCE,plan,Path(d)/'prepared',workers=4,mpi_ranks=86);mp=verify(result['manifest']);p.validate(mp);m=read_json(mp)
   self.assertEqual(len(m['tasks']),4);self.assertNotEqual(m['tasks'][0]['xyz']['sha256'],m['tasks'][1]['xyz']['sha256']);m['tasks'].pop();mp.write_text(json.dumps(m))
   with self.assertRaisesRegex(InvalidArtifact,'four-cell'):p.validate(mp)
if __name__=='__main__':unittest.main()
