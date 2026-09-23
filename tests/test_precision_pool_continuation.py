"""Uniform32 preparation and eventual real-output analysis checks."""
import json,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,InvalidArtifact,HA_TO_KCAL
import precision_pool_continue_run as run
from precision_pool_continuation import fingerprint,CANDIDATES,data
from nikasha_pool import choose_rows
from nikasha_pool_compare import extrema_reference
D=ROOT/'workspaces/precision_pool_continuation_20260923';M=D/'run_v1/stage1/manifest.json'

class PoolContinuationTests(unittest.TestCase):
    def test_actual_exact_population_and_reuses(self):
        inv,rows=run.source_rows(D/'INVENTORY_v2.json')
        self.assertEqual((inv['available_seed_pairs'],inv['two_pass_reuse_cells'],inv['maximum_new_scalar_calls']),(384,6,756))
        self.assertEqual(sum(r['role']=='calibration' for r in inv['cases']),25)
        for r in rows:
            self.assertLessEqual(r['pool_source_max_coordinate_component_difference_A'],1e-12)
            if r['continued_reuse']:
                reuse=r['continued_reuse']
                self.assertEqual(fingerprint(r),fingerprint(reuse['stage1']['source']))
                a=read_json(verify(reuse['stage1_collection']));b=read_json(verify(reuse['stage2_collection']))
                self.assertEqual(read_json(verify(b['manifest']))['previous'],reuse['stage1_collection'])
                self.assertEqual(a['stage'],1);self.assertEqual(b['stage'],2)

    def test_actual_finite_staged_manifest(self):
        v=run.validate(M)
        self.assertEqual((v['fresh_tasks'],v['reused'],v['missing'],v['cell_denominator']),(378,6,0,384))
        m=read_json(M)
        self.assertEqual(m['execution_resources'],{'mpi_ranks':1,'concurrent_tasks':32})
        self.assertEqual(m['settings']['reported_stage'],2)
        self.assertEqual(m['new_MACE_DFT_optimization_calls'],0)

    def test_cold_canonical_rule_replays_frozen_reference(self):
        inv=read_json(D/'INVENTORY_v2.json');p34=read_json(verify(inv['precision34']));ref=read_json(verify(p34['new_reference']))
        for variant in ('mathematical','operational'):
            rows=[]
            for group in inv['cases']:
                if group['role']!='calibration':continue
                case=next(c for c in data(group['collection'])['cases'] if c['case_id']==group['case_id'])
                pool=choose_rows(case['matrix'],list(CANDIDATES))
                rows.append({'case_id':group['case_id'],'expected_class':group['expected_class'],'R_model_kcal_mol':pool[variant]['composite_R_model_kcal_mol']})
            actual=extrema_reference(rows,variant,'archived_replay_only')
            self.assertEqual(actual['bands'],ref['variants'][variant]['bands'])
            self.assertEqual(actual['gap_model_kcal_mol'],ref['variants'][variant]['gap_model_kcal_mol'])

    def test_actual_terminal_384_cells_and32_qualification(self):
        path=D/'run_v1/COMPARISON.json'
        if not path.exists():self.skipTest('new uniform molecular qualification unrun')
        c=read_json(path);self.assertEqual((len(c['cells']),len(c['rows'])),(384,32))
        for cell in c['cells']:
            if cell['status']=='complete':
                delta=(cell['stage2_hartree']-cell['stage1_hartree'])*HA_TO_KCAL
                self.assertAlmostEqual(delta,cell['pass2_minus_pass1_kcal_mol'],places=10)
                self.assertEqual(cell['energy_settling_pass'],abs(delta)<=.1)
        for row in c['rows']:
            for stage in ('stage1','stage2'):
                self.assertEqual(choose_rows(row[stage+'_matrix'],list(CANDIDATES)),row[stage+'_pool'])
            stable=(row['cells_stable']==12 and all(x['pass'] for x in row['same_geometry_R_checks'].values()) and all(x['pass'] for x in row['pool_R_checks'].values()))
            self.assertEqual(row['qualified_status']=='available',stable)
        ref=read_json(verify(c['reference']));self.assertFalse(ref['crystals_or_noncanonical_used_for_fit'])
        for variant in ('mathematical','operational'):
            selected=[{'case_id':r['case_id'],'expected_class':r['expected_class'],'R_model_kcal_mol':r['qualified_R'][variant]} for r in c['rows'] if r['role']=='calibration']
            expected=extrema_reference(selected,variant,'precision_two_native_continuations_canonical25_v1')
            self.assertEqual(ref['variants'][variant]['bands'],expected['bands'])
            self.assertEqual(ref['variants'][variant]['status'],expected['status'])

if __name__=='__main__':unittest.main()
