"""Real whole-chain preparation, typed reuse and corruption regressions."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from mace_global_prepare import prepare_case,protein
from mace_omol_panel_prepare import validate
from mace_omol_panel import expected_tasks,reuse
from mace_omol_backbone_audit import inspect as inspect_backbone

W=ROOT/'workspaces/mace_omol_20260917'
P=W/'intact_panel_prepared_v1/preparation_manifest.json'
C=W/'intact_benchmark_v1/collection_job_1200816.json'


@unittest.skipUnless(P.exists(),'requires actual frozen whole-chain canonical preparation')
class PanelTests(unittest.TestCase):
    def test_archived_preparation_states_and_domain_flags(self):
        p=validate(P)
        self.assertEqual(p['status'],'complete');self.assertEqual(len(p['rows']),28)
        self.assertEqual(sum(r['outside_reported_training_charge_range'] for r in p['rows']),4)
        self.assertTrue(all(r['outside_reported_training_size_range'] for r in p['rows']))
        tasks=expected_tasks(p);self.assertEqual(len(tasks),112)
        self.assertEqual(sum(t['case_id'] not in ('1H4I','4MAE') for t in tasks),104)
        self.assertEqual(sum(r['evaluation_role']=='calibration' for r in p['rows']),25)
        for row in p['rows']:
            self.assertEqual(row['endpoint_charges']['La']-row['endpoint_charges']['Ca'],1)

    def test_existing_crystal_reuse_has_actual_matching_receipts(self):
        p=validate(P);refs,values,_=reuse(p,C)
        self.assertEqual(len(refs),8);self.assertEqual(set(refs),set(values))
        self.assertTrue(all(v['status']=='computed' and v['energy_only'] for v in values.values()))
        self.assertTrue(all(v['execution_adapter']['id']=='omol_nonlinear_exact_edge_batches_v1' for v in values.values()))

    def test_corrupted_real_saved_energy_cannot_be_reused(self):
        # Explicitly corrupted copy of an actual collection, not simulated output.
        p=validate(P);c=read_json(C)
        c['rows']['PQQ_1H4I_La_bound_primary']['energy_eV']+=1
        with tempfile.TemporaryDirectory() as temp:
            f=Path(temp)/'corrupted_collection.json';f.write_text(json.dumps(c))
            with self.assertRaises(InvalidArtifact):reuse(p,f)

    def test_actual_crystal_algebra_and_missing_calibration(self):
        from mace_omol_panel_report import summarize
        p=validate(P);_,values,_=reuse(p,C);result=summarize(p,values,inspect_backbone(p))
        self.assertIsNone(result['bands']);self.assertEqual(result['calibration_valid_count'],0)
        self.assertEqual(result['transfer_valid_count'],2)
        self.assertEqual(result['unavailable_score_count'],26)
        self.assertFalse(result['canonical_decision_gate_pass'])
        prior=read_json(W/'intact_report_v1/result.json')
        scores={r['case_id']:r for r in result['scores']}
        for name in ('1H4I','4MAE'):
            self.assertEqual(scores[name]['R_coord_kcal_mol'],prior['scores']['PQQ_'+name]['R_coord_kcal_mol'])
            self.assertIsNone(scores[name]['calibrated_class'])

    def test_real_1KB0_gaps_override_earlier_template_success(self):
        result=inspect_backbone(validate(P))
        self.assertEqual(result['unsupported_cases'],['1KB0'])
        row=next(r for r in result['rows'] if r['case_id']=='1KB0')
        self.assertEqual(len(row['invalid_bonds']),2)
        self.assertTrue(all(b['distance_A']>4 for b in row['invalid_bonds']))

    @unittest.skipUnless(importlib.util.find_spec('openmm'),'requires existing OpenMM preparation environment')
    def test_new_guard_rejects_real_gapped_input_before_template_matching(self):
        p=read_json(P);row=next(r for r in p['rows'] if r['case_id']=='1KB0')
        state=read_json(verify(row['preparation']))
        with self.assertRaisesRegex(InvalidArtifact,'unsupported peptide connection'):
            protein(state['source_audit_row'],allow_terminal_completion=True,check_peptide_connectivity=True)

    def test_one_recorded_1KB0_terminal_addition_and_preserved_heavy_coordinates(self):
        p=read_json(P);row=next(r for r in p['rows'] if r['case_id']=='1KB0')
        state=read_json(verify(row['preparation']));details=state['preparation_details']
        self.assertEqual(len(details['terminal_additions']),1)
        self.assertTrue(details['terminal_additions'][0]['id'].endswith('/OXT'))
        original={a['id']:a for a in details['original_protein_atoms']}
        for atom in state['physical_atoms']:
            if atom['kind']=='protein_source' and atom['element']!='H':
                self.assertEqual(atom['xyz_A'],original[atom['id']]['xyz_A'])
        self.assertEqual(row['endpoint_charges'],{'La':5,'Ca':4})

    @unittest.skipUnless(importlib.util.find_spec('openmm'),'requires existing OpenMM preparation environment')
    def test_default_helper_reproduces_real_old_GGR_and_PQQ_coordinates(self):
        old=read_json(ROOT/'workspaces/mace_global_benchmark_20260916/prepared_v1/preparation_manifest.json')
        with tempfile.TemporaryDirectory() as temp:
            for row in old['rows']:
                if row['case_id'] not in ('GGR_1GLG','PQQ_1H4I'):continue
                before=read_json(verify(row['preparation']))
                ref=prepare_case(before['source_audit_row'],Path(temp)/row['case_id'])
                after=read_json(verify(ref))
                self.assertEqual(after['physical_atoms'],before['physical_atoms'])
                self.assertEqual(after['preparation_details'],before['preparation_details'])
                for metal in ('La','Ca'):
                    # The archived writer's comment called every structure 1H4I.
                    # That comment was corrected earlier; all actual atom lines
                    # and the atom count must still reproduce byte for byte.
                    old_lines=verify(before['endpoints'][metal]['xyz']).read_text().splitlines()
                    new_lines=verify(after['endpoints'][metal]['xyz']).read_text().splitlines()
                    self.assertEqual(new_lines[:1]+new_lines[2:],old_lines[:1]+old_lines[2:])
                    self.assertEqual(after['endpoints'][metal]['state'],before['endpoints'][metal]['state'])


if __name__=='__main__':unittest.main()
