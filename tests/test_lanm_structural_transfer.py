"""Physical source mapping and finite-call checks on both real Hans crystals."""
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import lanm_structural_transfer as transfer
from affordable_common import read_json,verify,xyz

MANIFEST=ROOT/'workspaces/lanm_series_followup_20260923/dy_transfer_v1/manifest.json'


class RealHansTransfer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.m=read_json(MANIFEST)

    def test_finite_source_scope(self):
        r=transfer.validate(MANIFEST)
        self.assertEqual((r['supported_sites'],r['MACE_tasks'],r['GFN_initial_tasks']),(3,6,12))
        self.assertEqual(r['mapping_matches'],[True,True,True])
        self.assertEqual(self.m['finite_new_calls'],{'MACE':6,'GFN_initial':12,'GFN_self_continuation':12})
        old=read_json(verify(self.m['original_collection']))
        self.assertEqual(len([e for e in old['endpoints'] if e['endpoint_id'].startswith('Mex_') and e['status']=='complete']),6)

    def test_whole_carboxylate_is_retained_despite_denticity_change(self):
        for c,resnum in zip(self.m['cases'],[42,66,91]):
            r=read_json(verify(c['repair_manifest']))
            source=[a['source'] for a in r['atom_graph']['source_to_qm'] if a['kind']=='source']
            self.assertTrue({'OE1','OE2'} <= {a['atom'] for a in source if a['chain']=='A' and a['resnum']==resnum})
            self.assertEqual(r['paired_invariants']['atom_count'],50)
            self.assertEqual(r['explicit_water_inventory'],[])
            self.assertEqual(r['outputs']['La']['charge'],-1)

    def test_source_graph_caps_and_elements_match_without_coordinate_edits(self):
        for c in self.m['cases']:
            new=read_json(verify(c['repair_manifest']));old=read_json(verify(c['original_repair']))
            self.assertEqual(transfer.physical_map(new),transfer.physical_map(old))
            a,b=xyz(verify(new['outputs']['La']['xyz'])),xyz(verify(old['outputs']['La']['xyz']))
            self.assertEqual([x[0] for x in a],[x[0] for x in b])
            self.assertNotEqual(a,b)  # Real different crystals, not a reused/altered q0 label.
            self.assertEqual(new['heavy_coordinate_max_displacement_A'],0.)
        self.assertEqual(transfer.original.heavy_atoms(verify(self.m['source'])),
                         transfer.original.heavy_atoms(Path(MANIFEST).parent/'hans_Dy_protonated.pdb'))

    def test_actual_result_reuses_Mex_and_closes_the_exchange_cycle(self):
        result=read_json(MANIFEST.parent/'final_collection.json')
        old=read_json(verify(self.m['original_collection']))
        self.assertEqual(result['reused_Mex'],{e['endpoint_id']:e for e in old['endpoints'] if e['endpoint_id'].startswith('Mex_')})
        self.assertEqual(result['fresh_complete'],6)
        for v in result['ordered_vector']:
            work=v['Dy_source_minus_La_source_work']
            self.assertAlmostEqual(v['change_in_D_kcal_mol'],work['Dy']['composite_kcal_mol']-work['La']['composite_kcal_mol'],places=8)
            self.assertAlmostEqual(v['D_Dy_crystal_kcal_mol'],v['D_native_MACE_kcal_mol']+v['D_solvent_transfer_kcal_mol'],places=8)

    def test_actual_physical_batches_and_native_restarts(self):
        result=read_json(MANIFEST.parent/'final_collection.json')
        for e in result['fresh_endpoints']:
            mult=6 if e['endpoint_id'].endswith('__Dy') else 1
            self.assertEqual(e['MACE_receipt']['input_state_check']['spin_multiplicity'],mult)
            for cell in e['cells'].values():
                self.assertEqual(cell['audit']['native_valence_electrons'],150)
                self.assertIn('INITIAL GUESS: XTBRESTART',verify(cell['continuation']['output']).read_text())


if __name__=='__main__':unittest.main()
