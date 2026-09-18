"""Real expansion preparations and qualified worker dispatch; no mock science."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,xyz
from mace_density_frameworks import validate as frameworks
from mace_density_short import validate as short,collect
from mace_tinker_framework_solver import check_parameters
from mace_hybrid import dry_run

BASE=ROOT/'workspaces/mace_omol_20260917'
FRAME=BASE/'trial_gk_expansion_frameworks_v1/manifest.json'
SHORT=BASE/'trial_gk_expansion_short_v1/manifest.json'


class ExpansionTests(unittest.TestCase):
    def test_actual_native_parameter_reads_preserve_full_physical_inventory(self):
        m=frameworks(FRAME)
        self.assertEqual({t['case_id'] for t in m['tasks']},{'GGR_2FW0','GGR_2FVY'})
        self.assertEqual(m['new_energy_calls'],0)
        for t in m['tasks']:
            mapping=read_json(verify(t['mapping']));params=read_json(verify(t['parameters']))
            self.assertEqual(len(params['atoms'])+1,len(mapping['physical_atoms']))
            self.assertEqual(len(set(mapping['system_atom_ids'])),len(params['atoms']))
            self.assertTrue(check_parameters(params,t,mapping))

    def test_corrupted_real_native_radius_fails(self):
        m=frameworks(FRAME);t=m['tasks'][0];params=read_json(verify(t['parameters']))
        params['atoms'][0]['radius_A']=0.
        with self.assertRaisesRegex(InvalidArtifact,'radius'):
            check_parameters(params,t,read_json(verify(t['mapping'])))

    def test_eight_exact_source_tasks_use_existing_qualified_dispatch(self):
        self.assertEqual(dry_run(SHORT)['tasks'],8);m=read_json(SHORT)
        for t in m['tasks']:
            self.assertEqual(t['energy_component'],'interaction_energy')
            case=read_json(verify(t['source_mapping']))
            p=case if t['kind']=='core' else read_json(verify(case['normalized_global_preparation']))
            self.assertEqual(xyz(verify(t['xyz'])),xyz(verify(p['endpoints'][t['metal']]['xyz'])))

    def test_changed_real_short_checkpoint_is_rejected(self):
        m=read_json(SHORT);m['model']['dtype']='float32' # explicitly corrupted real configuration
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_manifest.json';p.write_text(json.dumps(m))
            with self.assertRaisesRegex(InvalidArtifact,'Hamiltonian'):short(p)

    def test_missing_real_results_never_become_a_hybrid_score(self):
        m=read_json(SHORT)
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'copied_real_manifest_without_executions.json';p.write_text(json.dumps(m))
            r=collect(p)
            self.assertEqual(r['status'],'incomplete');self.assertIsNone(r['full_hybrid_score'])
            self.assertTrue(all(v['energy_eV'] is None for v in r['rows'].values()))

    def test_real_completed_quantum_and_short_outputs(self):
        path=BASE/'trial_gk_expansion_collection_v3/result.json'
        if not path.exists():self.skipTest('actual quantum collection not available')
        q=read_json(path);verify(q['collection_implementation']);verify(q['parser_implementation'])
        self.assertEqual(q['status'],'complete');self.assertEqual(len(q['rows']),4)
        for row in q['rows'].values():
            self.assertEqual(row['status'],'complete');self.assertTrue(all(row['native_components'].values()))
            verify(row['receipt']);verify(row['output']);verify(row['engrad'])
        self.assertIsNone(q['variational_response_check'])
        result=collect(SHORT)
        if result['status']!='complete':self.skipTest('actual MACE outputs not complete')
        self.assertEqual(len(result['rows']),8);self.assertIsNone(result['full_hybrid_score'])


if __name__=='__main__':unittest.main()
