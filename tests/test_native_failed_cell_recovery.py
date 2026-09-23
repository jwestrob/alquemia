import copy,json,re,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import native_failed_cell_recovery as n
from affordable_common import read_json,verify,InvalidArtifact,xyz
ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'workspaces/native_failed_cell_recovery_20260923/run_v1'
INV=ROOT/'diagnostics/native_failed_cell_recovery_20260923/SEEDS.json'

class ActualRecovery(unittest.TestCase):
    def test_actual_seeds_order_state(self):
        d=read_json(INV);self.assertEqual([r['seed_kind']for r in d['seeds']],list(n.SEEDS))
        order=[a[0]for a in xyz(verify(d['target']['task']['xyz']))]
        self.assertEqual(len(order),204)
        for row in d['seeds']:
            self.assertEqual([a[0]for a in xyz(verify(row['task']['xyz']))],order)
            self.assertEqual(row['bytes'],{'gbw':3855516,'xtbw':17368})
            self.assertEqual((row['SCF']['charge'],row['SCF']['multiplicity'],row['electron_count']),(-3,1,636))
            for pin in row['seed_source'].values():verify(pin)
        self.assertFalse(d['target']['failure_receipt']['normal_termination'])
        self.assertFalse(d['target']['failure_receipt']['scf_converged'])
    def test_actual_two_task_preflight_and_recipe(self):
        p=RUN/'manifest.json';m=read_json(p);self.assertEqual(n.validate(p)['tasks'],2)
        target=read_json(INV)['target']['task']
        for t in m['tasks']:
            self.assertEqual(verify(t['xyz']).read_bytes(),verify(target['xyz']).read_bytes())
            self.assertEqual(verify(t['input']).read_text(),verify(target['input']).read_text().replace(' NoAutostart',''))
    def test_wrong_seed_or_state_rejected(self):
        m=read_json(RUN/'manifest.json')
        for field,value in [('seed_source',m['tasks'][1]['seed_source']),('charge',-2)]:
            bad=copy.deepcopy(m);bad['tasks'][0][field]=value
            with tempfile.NamedTemporaryFile(mode='w',suffix='.json',dir=RUN) as f:
                json.dump(bad,f);f.flush()
                with self.assertRaises(InvalidArtifact):n.validate(Path(f.name))
    def test_real_result_or_incomplete_preserves_failure(self):
        p=RUN/'COLLECTION.json'
        if p.exists():
            d=read_json(p);self.assertEqual(d['primary_status'],'unavailable_unchanged')
            self.assertFalse(d['primary_result_changed'])
            for r in d['rows']:
                if r['status']=='confirmed_restart':
                    text=verify(r['actual']['output']).read_text()
                    self.assertIn('INITIAL GUESS: XTBRESTART',text)
                    block=text.split('CARTESIAN COORDINATES (ANGSTROEM)',1)[1].split('CARTESIAN COORDINATES (A.U.)',1)[0]
                    observed=[(e,*map(float,(x,y,z)))for e,x,y,z in re.findall(r'^\s*([A-Z][a-z]?)\s+([-+0-9.]+)\s+([-+0-9.]+)\s+([-+0-9.]+)\s*$',block,re.M)]
                    target=xyz(verify(read_json(INV)['target']['task']['xyz']))
                    self.assertEqual([a[0]for a in observed],[a[0]for a in target])
                    self.assertLessEqual(max(abs(a-b)for arow,brow in zip(observed,target)for a,b in zip(arow[1:],brow[1:])),5.1e-7)
            if d['agreement_pass']:
                self.assertEqual(d['recovery_energy_hartree'],d['rows'][0]['energy_hartree'])
                self.assertLessEqual(abs(d['seed_difference_kcal_mol']),.1)
            else:self.assertIsNone(d['recovery_energy_hartree'])
        else:
            with tempfile.TemporaryDirectory() as d:
                p=Path(d)/'INCOMPLETE.json';v=n.collect(RUN/'manifest.json',p)
                self.assertEqual(v['status'],'recovery_unavailable');self.assertIsNone(v['recovery_energy_hartree'])
                self.assertTrue(all(r['status']=='unavailable'for r in read_json(p)['rows']))

if __name__=='__main__':unittest.main()
