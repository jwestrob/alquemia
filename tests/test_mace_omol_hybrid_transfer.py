"""Real source-backed structural transfer; no fabricated scientific outputs."""
from pathlib import Path
import copy,json,sys,tempfile,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,xyz,InvalidArtifact
from mace_omol_hybrid_transfer import normalized_core,prepared_states,validate,validate_quantum
BASE=ROOT/'workspaces/mace_omol_hybrid_transfer_20260918'
PREP=BASE/'prepared_v2/preparation.json'
OLD=ROOT/'workspaces/mace_omol_20260917/matched_H_prepared_v1'

@unittest.skipUnless((OLD/'preparation.json').exists(),'requires real archived normalized GGR cores')
class TransferTests(unittest.TestCase):
    def test_uniform_mapping_reproduces_both_archived_1GLG_cores(self):
        for representation in ('extended','connected'):
            c=read_json(OLD/('GGR_'+representation)/'mapping.json')
            original=read_json(verify(c['source_mapping']))
            repair=read_json(verify(original['source_preparation']))
            ends,mapping,checks=normalized_core(repair,read_json(verify(c['normalized_global_preparation'])))
            self.assertTrue(checks['heavy_coordinates_unchanged'])
            for metal,e in ends.items():
                archive=xyz(verify(c['endpoints'][metal]['xyz']))
                self.assertEqual([a[0] for a in e['atoms']],[a[0] for a in archive])
                np.testing.assert_allclose([a[1:] for a in e['atoms']],[a[1:] for a in archive],atol=5e-10,rtol=0)
                self.assertEqual(e['charge'],c['endpoints'][metal]['charge'])

    @unittest.skipUnless(PREP.exists(),'requires real new transfer preparations')
    def test_transfer_counts_pairs_and_source_graphs(self):
        p,cases=prepared_states(PREP);self.assertEqual(len(cases),4)
        for name,c in cases.items():
            repair=read_json(verify(c['source_preparation']));w=read_json(verify(c['whole_preparation']))
            replay,mapping,checks=normalized_core(repair,w)
            self.assertEqual(mapping,c['mapping']);self.assertEqual(checks,c['normalization_replay'])
            expected=111 if name.endswith('connected') else 58
            for metal,e in replay.items():
                rows=xyz(verify(c['endpoints'][metal]['xyz']));self.assertEqual(len(rows),expected)
                np.testing.assert_allclose([a[1:] for a in rows],[a[1:] for a in e['atoms']],atol=5e-10,rtol=0)
            self.assertEqual(replay['Ca']['atoms'][1:],replay['La']['atoms'][1:])

    @unittest.skipUnless((BASE/'initial_v1/initial.json').exists(),'requires actual finite manifests')
    def test_finite_tasks_and_corrupted_real_charge_rejected(self):
        initial=read_json(BASE/'initial_v1/initial.json');total=0
        for ref in initial['MACE_manifests']:total+=validate(verify(ref))['tasks']
        self.assertEqual(total,156);self.assertEqual(validate_quantum(verify(initial['quantum_manifest']))['status'],'dry_run_pass')
        m=read_json(verify(initial['MACE_manifests'][0]));m['tasks'][0]['charge']+=2
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'corrupted_actual_charge.json';p.write_text(json.dumps(m))
            with self.assertRaises(InvalidArtifact):validate(p)

if __name__=='__main__':unittest.main()
