"""Real normalized core/field fixtures; no fabricated native outputs."""
from pathlib import Path
import sys,copy,tempfile,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from mace_responsive_field import validate,pointcharge_text,key
from mace_omol_vacuum import parse_endpoint,embedded_input
W=ROOT/'workspaces/mace_omol_20260917';M=W/'responsive_quantum_v2/manifest.json'

@unittest.skipUnless(M.exists(),'real embedded preparation unavailable')
class ResponsiveField(unittest.TestCase):
    def test_exact_fields_and_electronic_pairs(self):
        m=read_json(M);self.assertEqual(validate(M)['tasks'],8)
        for t in m['tasks']:
            pc=np.loadtxt(verify(t['pointcharges']),skiprows=1);s=read_json(verify(t['source_state']))
            weights=read_json(verify(t['source_task']['weights']));idx=weights['physical_indices']
            np.testing.assert_array_equal(pc[:,0],np.array(s['environment_charges_e'])[idx])
            np.testing.assert_array_equal(pc[:,1:],[s['physical_atoms'][i]['xyz_A'] for i in idx])
            self.assertEqual(verify(t['input']).read_text(),embedded_input(t['charge']))
            self.assertEqual(t['xyz']['sha256'],t['source_task']['xyz']['sha256'])

    def test_corrupted_real_field_and_input_rejected(self):
        m=read_json(M)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);t=m['tasks'][0];pc=verify(t['pointcharges']).read_text();lines=pc.splitlines();cols=lines[1].split();cols[0]=str(float(cols[0])+.1);lines[1]=' '.join(cols)
            cp=p/'environment.pc';cp.write_text('\n'.join(lines)+'\n');t['pointcharges']=record(cp)
            mp=p/'manifest.json';write_new(mp,m)
            with self.assertRaisesRegex(InvalidArtifact,'field/cache changed'):validate(mp)
        m=read_json(M);t=m['tasks'][0]
        altered=copy.deepcopy(t);altered['source_state']['sha256']='corrupted_real_fixture'
        self.assertNotEqual(key(t,m),key(altered,m))

    def test_actual_vacuum_output_cannot_pass_as_embedded(self):
        m=read_json(M);t=m['tasks'][0];q=read_json(W/'matched_H_quantum_v1/manifest.json')
        old=next(v for v in q['tasks'] if v['task_id']==t['task_id'])
        with self.assertRaisesRegex(InvalidArtifact,'external-charge inventory'):
            parse_endpoint(t,old['output_path'],old['engrad_path'],permanent_field=True)

    def test_real_conditional_warning_and_corrupted_numerical_header(self):
        t=read_json(M)['tasks'][0];p=Path(t['output_path'])
        if not p.exists():self.skipTest('actual embedded output unavailable')
        row=parse_endpoint(t,p,t['engrad_path'],permanent_field=True)
        self.assertTrue(row['native_components']['analytic_scf'])
        with tempfile.TemporaryDirectory() as tmp:
            bad=Path(tmp)/'corrupted_real.out';bad.write_text(p.read_text()+'\nORCA NUMERICAL GRADIENT CALCULATION\n')
            with self.assertRaisesRegex(InvalidArtifact,'numerical gradients unsupported'):
                parse_endpoint(t,bad,t['engrad_path'],permanent_field=True)

    @unittest.skipUnless((W/'responsive_quantum_result_v2.json').exists(),'actual embedded outputs unavailable')
    def test_actual_embedded_native_states(self):
        r=read_json(W/'responsive_quantum_result_v2.json');self.assertEqual(r['status'],'complete')
        self.assertEqual(len(r['rows']),8)
        for row in r['rows'].values():
            self.assertEqual(row['energy_scope'],'embedded_permanent_field_endpoint')
            self.assertTrue(row['core_gradient_includes_permanent_field']);self.assertIsNone(row['combined_gradient'])
            self.assertTrue(all(row['native_components'].values()))

if __name__=='__main__':unittest.main()
