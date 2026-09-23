"""Pinned real pool/seed/state checks; no synthetic molecular results."""
import sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import strict_native_pool as s
from affordable_common import read_json,verify,HA_TO_KCAL
from mace_hybrid import EV_TO_KCAL
from native_pool_continuation import recipe as old_recipe
D=ROOT/'workspaces/strict_native_pool_20260923/run_v1'

class StrictPoolTests(unittest.TestCase):
    def test_population_reuse_and_paths(self):
        for branch,n,reuse in [('fresh',384,0),('cold_seed',374,10)]:
            mp=D/branch/'manifest.json';v=s.validate(mp);self.assertEqual((v['fresh_tasks'],v['reused'],v['missing']),(n,reuse,0))
            m=read_json(mp);self.assertEqual(len(m['tasks'])+len(m['reused']),384)
            for t in m['tasks']:
                self.assertTrue(Path(t['output_path']).is_relative_to(mp.parent))
                self.assertEqual(verify(t['xyz']).read_bytes(),verify(t['source']['xyz']).read_bytes())
    def test_only_tolerance_and_initialization_change(self):
        for branch in s.BRANCHES:
            m=read_json(D/branch/'manifest.json')
            for t in m['tasks']:
                body=verify(t['input']).read_text();base=body.replace(' TolE 1e-10\n','').replace(' NoAutostart','')
                self.assertEqual(base,old_recipe(t['charge'],1,t['medium']))
                self.assertEqual('NoAutostart' in body,branch=='fresh')
                if branch=='fresh':self.assertIsNone(t['seed_source']);self.assertNotIn('immutable_seed',t)
                else:
                    self.assertEqual(t['immutable_seed']['sha256'],t['source']['xtbw']['sha256'])
                    self.assertEqual(t['immutable_gbw']['sha256'],t['source']['gbw']['sha256'])
    def test_actual_final_initialization_and_denominators(self):
        if not all((D/b/'COLLECTION.json').exists() for b in s.BRANCHES):self.skipTest('strict32 molecular calls unrun')
        for branch in s.BRANCHES:
            c=read_json(D/branch/'COLLECTION.json');self.assertEqual(len(c['rows']),384)
            for r in c['rows']:
                if r['status']!='complete':self.assertIsNone(r['energy_hartree']);continue
                self.assertEqual(r['observed_TolE_hartree'],1e-10)
                text=verify(r['actual']['output']).read_text()
                self.assertIn('INITIAL GUESS: '+('SAD' if branch=='fresh' else 'XTBRESTART'),text)
                self.assertTrue(r['details']['native_mixer_observed'])


    def test_actual_pools_signs_independent_references_and_failed_gates(self):
        path=D/'COMPARISON.json'
        if not path.exists():self.skipTest('strict32 comparison unrun')
        from nikasha_pool import choose_rows
        from accommodation_folds_compare import decision
        c=read_json(path);refs=read_json(D/'REFERENCES.json')
        self.assertEqual((len(c['rows']),len(c['cells'])),(32,384))
        self.assertEqual(sum(r['role']=='calibration' for r in c['rows']),25)
        for row in c['rows']:
            for b in s.BRANCHES:
                v=row['branches'][b];self.assertEqual(choose_rows(v['matrix'],['origin','adaptive_Ca','adaptive_La']),v['pool'])
                for q in ['origin','adaptive_Ca','adaptive_La']:
                    cells={z:v['matrix'][z][q] for z in ['Ca','La']}
                    if all(t['status']=='complete' for t in cells.values()):
                        energies={z:t['components']['MACE_eV']*EV_TO_KCAL+(t['components']['GFN2_ALPB_hartree']-t['components']['GFN2_vacuum_hartree'])*HA_TO_KCAL for z,t in cells.items()}
                        self.assertTrue(all(isinstance(e,float) for e in energies.values()))
                        # Direct endpoint difference preserves E_Ca-E_La and the additive ALPB-vac correction.
                        direct=energies['Ca']-energies['La']
                        components=(cells['Ca']['components']['MACE_eV']-cells['La']['components']['MACE_eV'])*EV_TO_KCAL
                        components+=((cells['Ca']['components']['GFN2_ALPB_hartree']-cells['Ca']['components']['GFN2_vacuum_hartree'])-(cells['La']['components']['GFN2_ALPB_hartree']-cells['La']['components']['GFN2_vacuum_hartree']))*HA_TO_KCAL
                        self.assertAlmostEqual(direct,components,places=7)
                for variant in ['mathematical','operational']:
                    bands=refs['branches'][b]['variants'][variant]['bands']
                    expected=decision(v['qualified_R'][variant],bands) if bands else 'unavailable'
                    self.assertEqual(v['own_reference'][variant]['decision'],expected)
                    if row['qualified_status']!='available':self.assertIsNone(v['qualified_R'][variant])
        for cell in c['cells']:
            if cell['status']=='complete':
                delta=(cell['strict_hartree']['cold_seed']-cell['strict_hartree']['fresh'])*HA_TO_KCAL
                self.assertAlmostEqual(delta,cell['cold_seed_minus_fresh_kcal_mol'],places=10)
                self.assertEqual(abs(delta)<=.1,cell['initialization_agreement_pass'])
            else:self.assertIsNone(cell['cold_seed_minus_fresh_kcal_mol'])
        for b in s.BRANCHES:
            self.assertFalse(refs['branches'][b]['crystals_or_noncanonical_used_for_fit'])
            if any(r['qualified_status']!='available' for r in c['rows'] if r['role']=='calibration'):
                self.assertTrue(all(v['bands'] is None for v in refs['branches'][b]['variants'].values()))

if __name__=='__main__':unittest.main()
