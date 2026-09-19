"""Real pinned preparation fixtures; no invented energies or scientific mocks."""
import sys
import unittest
import copy
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import second_shell_context as ss
from affordable_common import read_json, InvalidArtifact

ROOT=Path(__file__).resolve().parents[1]

class TestSecondShell(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config=read_json(ROOT/'diagnostics/second_shell_20260919/CONFIG.json')
        cls.data={}
        for item in cls.config['cases']:
            state=ss.parent_state(item,cls.config['topology'])
            cls.data[item['case']]=(state,*ss.expansion(state))

    def test_real_source_coordinates_and_waters(self):
        for case,(state,rows,audit) in self.data.items():
            for metal in ('Ca','La'):
                for old,new in audit['core_to_context'].items():
                    self.assertEqual(state['original'][metal][old],rows[metal][new],case)
            allowed={audit['core_to_context'][i] for k,i in state['old_sources'].items()
                     if k in state['water_keys'] and state['graph'].meta[k]['element']=='H'}
            self.assertLessEqual(set(audit['paired_differing_nonmetal_indices']),allowed)

    def test_complete_unchanged_pqq(self):
        for case in ('1H4I','4MAE'):
            state,rows,audit=self.data[case]
            self.assertEqual(len(state['opaque']),27)
            self.assertEqual(audit['added_formal_charge'],0)
            for entry in state['opaque']:
                i=entry['qm_index'];j=audit['core_to_context'][i]
                self.assertEqual(state['original']['La'][i],rows['La'][j])

    def test_connectivity_overlap_and_charge(self):
        for case,(state,rows,audit) in self.data.items():
            entries=audit['mapping']['source_to_qm']
            self.assertEqual({a['qm_index'] for a in entries},set(range(1,len(rows['La']))))
            ids=[ss.atom_key(a['source']) for a in entries if 'source' in a]
            self.assertEqual(len(ids),len(set(ids)))
            for cut in audit['mapping']['cut_bonds_and_caps']:
                self.assertFalse({cut['retained']['atom'],cut['omitted']['atom']}=={'C','N'})
            for metal in ('Ca','La'):
                item=next(i for i in self.config['cases'] if i['case']==case)
                ss.check_atoms(rows[metal],item['endpoints'][metal]['charge']+audit['added_formal_charge'])

    def test_source_graph_rejects_corrupted_real_neighbor(self):
        item=next(i for i in self.config['cases'] if i['case']=='1H4I')
        state=ss.parent_state(item,self.config['topology'])
        _,fragments,_,_=ss.discover(state)
        rkey=next(r for kind,r in fragments if kind=='sidechain' and state['graph'].residues[r].resname=='TRP')
        key=state['graph'].key(rkey,'NE1')
        del state['graph'].atoms[key]
        with self.assertRaises((InvalidArtifact,KeyError)):
            ss.sidechain(state['graph'],rkey)

    def test_actual_mace_receipt_and_corrupted_charge(self):
        from mace_omol import accepted_state
        path=ROOT/'workspaces/second_shell_20260919/prepared_v2/mace_collection_1202083.json'
        if not path.exists():self.skipTest('actual GPU integration outputs unavailable; no substitute')
        collection=read_json(path)
        manifest=read_json(collection['manifest']['path'])
        for row in collection['rows']:
            task=next(t for t in manifest['tasks'] if t['task_id']==row['task_id'])
            self.assertEqual(row['status'],'complete')
            self.assertTrue(accepted_state(row['accepted'],task))
            corrupted=copy.deepcopy(row['accepted']);corrupted['input_state_check']['charge']+=1
            self.assertFalse(accepted_state(corrupted,task))

    def test_actual_mace_sign_and_once_only_conversion(self):
        from mace_hybrid import EV_TO_KCAL
        path=ROOT/'workspaces/second_shell_20260919/prepared_v2/mace_collection_1202083.json'
        if not path.exists():self.skipTest('actual GPU integration outputs unavailable; no substitute')
        collection=read_json(path)
        for case in collection['cases']:
            energies={(r['variant'],r['metal']):r['accepted']['energy_eV'] for r in collection['rows'] if r['case']==case['case']}
            change=((energies['expanded','Ca']-energies['expanded','La'])-(energies['core','Ca']-energies['core','La']))*EV_TO_KCAL
            self.assertAlmostEqual(change,case['delta_R_kcal_mol'],places=7)

    def test_actual_DFT_context_component_accounting(self):
        path=ROOT/'diagnostics/second_shell_20260919/result_v1/result.json'
        if not path.exists():self.skipTest('actual completed context DFT unavailable; no substitute')
        report=read_json(path)
        for row in report['rows']:
            comp=row['component_delta_R_kcal_mol']
            # CPCM is already inside SCF and must not be added for a second time.
            self.assertAlmostEqual(comp['SCF']+comp['dispersion']+comp['gCP'],row['DFT_delta_R_kcal_mol'],places=5)

if __name__=='__main__':unittest.main()
