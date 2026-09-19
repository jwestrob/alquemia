"""Independent parvalbumin source geometry and archived native-energy fixtures."""
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import water_reference_transfer as wr
from affordable_common import read_json,verify,xyz
ROOT=Path(__file__).resolve().parents[1]
MP=ROOT/'workspaces/water_reference_validation_20260919/proposals_v1/manifest.json'


class WaterReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.m=read_json(MP)

    def test_ordered_vector_and_actual_baseline_receipts(self):
        self.assertEqual([p['site'] for p in self.m['parents']],['CD','EF'])
        for p in self.m['parents']:wr.checked_parent(p)
        self.assertEqual(wr.validate(MP)['tasks'],2)
        cd,ef=[read_json(verify(p['preparation'])) for p in self.m['parents']]
        self.assertEqual(cd['explicit_water_inventory'],[])
        self.assertEqual([w['resnum'] for w in ef['explicit_water_inventory']],[166])

    def test_actual_source_context_water_roles(self):
        context=read_json(verify(self.m['contexts']['PARV_4CPV_EF']))
        self.assertEqual([(w['source']['resnum'],w['role']) for w in context['water_groups']],[(166,'variable'),(257,'frozen_outer')])
        for t in self.m['tasks']:
            source=xyz(verify(t['starts']['source']));alternate=xyz(verify(t['starts']['radial_away']))
            for i in range(len(source)):
                if i not in t['mobile_indices']:self.assertEqual(source[i],alternate[i])
            self.assertEqual(len(t['mobile_indices']),2)

    def test_no_unrequested_endpoints_or_protonation_change(self):
        self.assertEqual({(t['case_id'],t['metal']) for t in self.m['tasks']},{('PARV_4CPV_EF','Ca'),('PARV_4CPV_EF','La')})
        ca=next(t for t in self.m['tasks'] if t['metal']=='Ca');la=next(t for t in self.m['tasks'] if t['metal']=='La')
        self.assertEqual(la['charge']-ca['charge'],1)
        self.assertEqual(xyz(verify(ca['xyz']))[1:],xyz(verify(la['xyz']))[1:])

    def test_original_core_hamiltonian_and_only_water_H_transfer(self):
        m=read_json(ROOT/'workspaces/water_reference_validation_20260919/dft_v1/manifest.json')
        self.assertEqual(len(m['tasks']),2)
        self.assertEqual([r['site'] for r in m['identities']],['CD'])
        for t in m['tasks']:
            p=read_json(verify(t['parent']));before=xyz(verify(p['outputs'][t['metal']]['xyz']));after=xyz(verify(t['xyz']))
            self.assertEqual(len(before),len(after));self.assertEqual(t['charge'],p['outputs'][t['metal']]['charge'])
            for i in range(len(before)):
                if i not in t['mobile_H_indices']:self.assertEqual(before[i],after[i])
                else:self.assertEqual(before[i][0],'H')
            old=verify(p['outputs'][t['metal']]['input']).read_text().split('* xyzfile')[0]
            new=verify(t['input']).read_text().split('* xyzfile')[0]
            self.assertEqual(old,new)

    def test_actual_completed_proposal_and_native_receipt_replay(self):
        root=MP.parent.parent
        self.assertEqual(wr.collect(MP),read_json(root/'proposals_v1/collection_job_1202475.json'))
        manifest=root/'dft_v1/manifest.json'
        expected=read_json(root/'dft_v1/collection_1202476.json')
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'collected.json'
            self.assertEqual(wr.collect_dft(manifest,output),{'status':'complete','sites':2})
            self.assertEqual(read_json(output),expected)

    def test_actual_ordered_contrast_and_explicit_missing_decisions(self):
        result=read_json(MP.parent.parent/'dft_v1/collection_1202476.json')
        cd,ef=result['rows']
        self.assertEqual(cd['endpoints'],cd['baseline_endpoints'])
        self.assertEqual(cd['delta_R_kcal_mol'],0.0)
        delta={m:(ef['endpoints'][m]['energy_hartree']-ef['baseline_endpoints'][m]['energy_hartree'])*wr.HA_TO_KCAL for m in ('Ca','La')}
        self.assertAlmostEqual(ef['delta_R_kcal_mol'],delta['Ca']-delta['La'],places=8)
        self.assertAlmostEqual(ef['delta_R_kcal_mol'],-3.3215921674855053,places=8)
        for row in result['rows']:
            self.assertIsNone(row['calibrated_decision'])
            self.assertIsNone(row['occupancy_probability'])
        # Preserve the actual pre-execution collection; missing computations were
        # unavailable, never a zero correction or a successful fallback.
        partial=read_json(MP.parent.parent/'dft_v1/collection_before_execution.json')
        self.assertEqual(partial['rows'][0]['delta_R_kcal_mol'],0.0)
        self.assertEqual(partial['rows'][1]['status'],'unavailable')
        self.assertIsNone(partial['rows'][1]['prepared_R_kcal_mol'])
        self.assertIsNone(partial['rows'][1]['delta_R_kcal_mol'])


if __name__=='__main__':unittest.main()
