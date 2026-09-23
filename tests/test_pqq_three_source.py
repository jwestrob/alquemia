"""Real source/graph fixtures; no invented scientific results or energy calls."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import pqq_three_source as p

BASE=ROOT/'workspaces/pqq_three_source_20260923'
FIX=BASE/'fixtures_v2'
PLM=('PQQSEQ_07ab500e3df76b30d71c','PQQSEQ_83440678cbbd658047c9')


class ThreeSource(unittest.TestCase):
    def test_real_requests_need_no_historical_member_or_label(self):
        paths=[FIX/(t+'_REQUEST.json') for t in ('triple_000','triple_096')]+[BASE/'plm_v1'/(n+'_REQUEST.json') for n in PLM]
        for path in paths:
            req=p.check_request(p.read_json(path))
            self.assertFalse(req['anchor_is_biological_label_or_calibration_member'])
            self.assertEqual(req['physical_state_anchor'],min(c['case_id'] for c in req['cases']))
            self.assertTrue(all('canonical_coordinate_match' not in c and 'expected_class' not in c for c in req['cases']))

    def test_actual_reference_union_state_and_geometry_replay(self):
        info=p.read_json(FIX/'FIXTURES.json');old=p.read_json(p.verify(info['original_preparation']))
        for t in ('triple_000','triple_096'):
            actual=p.read_json(FIX/t/'PREPARATION.json');check=p.dry_run(FIX/t/'PREPARATION.json')
            self.assertEqual(check['supported_sources'],3);self.assertEqual(check['new_molecular_calls'],0)
            triple=next(x for x in old['triples'] if x['triple_id']==t)
            for row in actual['cases']:
                previous=next(x for x in old['cases'] if (x['case_id'],x['selection_id'])==(row['case_id'],triple['selection_id']))
                self.assertEqual(p.read_json(p.verify(row['state_signature'])),p.read_json(p.verify(previous['state_signature'])))
                for z in ('Ca','La'):
                    a,b=[p.xyz(p.verify(r['representations']['context']['endpoints'][z]['xyz'])) for r in (row,previous)]
                    self.assertEqual([v[0] for v in a],[v[0] for v in b])
                    np.testing.assert_allclose([v[1:] for v in a],[v[1:] for v in b],atol=1e-12,rtol=0)

    def test_paired_maps_and_full_state_for_actual_plm_sources(self):
        for name in PLM:
            path=BASE/'plm_v1'/name/'PREPARATION.json';data=p.read_json(path);check=p.dry_run(path)
            self.assertEqual(check['status'],'prepared');self.assertEqual(len(data['tasks']),6)
            self.assertFalse(data['new_protonation']);self.assertFalse(data['future_execution']['launched'])
            req=p.read_json(p.verify(data['request']))
            self.assertTrue(all(x['expected_class'] is None and x['evidence_stratum']=='unresolved_PLM_prediction' for x in req['source_evidence'].values()))
            for row in data['cases']:
                kin=p.Kinematics(p.read_json(p.verify(row['physical_mapping']))['context']);x=kin.evaluate(np.zeros(len(kin.modes)))[1]
                for z in ('Ca','La'):
                    ep=row['representations']['context']['endpoints'][z]
                    np.testing.assert_allclose(x,[a[1:] for a in p.xyz(p.verify(ep['xyz']))],atol=1e-12,rtol=0)
                signature=p.read_json(p.verify(row['state_signature']))
                self.assertEqual(set(signature),{'source_atoms_and_cofactor','H_covalent_parents','caps','retained_bonds','formal_charges','endpoint_charges','multiplicities','core_microstates','water_inventory'})

    def test_real_failed_source_invalidates_full_declared_triple(self):
        path=FIX/'triple_037/PREPARATION.json';data=p.read_json(path);check=p.dry_run(path)
        self.assertEqual(check['status'],'unavailable');self.assertEqual(data['denominator'],3)
        self.assertEqual(len(data['cases']),3);self.assertEqual(data['tasks'],[]);self.assertIsNone(data['union'])
        self.assertEqual(data['missing_members'],['mmol_1770-pqq-la_model__conditioned_La__seed-1_sample-4'])
        self.assertTrue(all(x['status']=='group_source_unavailable' for x in data['cases']))

    def test_corrupted_real_request_rejects_anchor_member_and_reference_changes(self):
        original=p.read_json(FIX/'triple_000_REQUEST.json')
        for kind in ('anchor','conditioning','duplicate'):
            value=copy.deepcopy(original)
            if kind=='anchor':value['physical_state_anchor']=value['cases'][1]['case_id']
            elif kind=='conditioning':value['cases'][0]['source_conditioning_metal']='Ca'
            else:value['cases'][1]=copy.deepcopy(value['cases'][0])
            with self.assertRaises(p.InvalidArtifact):p.check_request(value)
        with tempfile.TemporaryDirectory() as d:
            ref=p.read_json(p.verify(original['reference']));ref['variants']['operational']['bands']['La_min']+=1
            path=Path(d)/'explicitly_corrupted_real_reference.json';p.write_new(path,ref)
            value=copy.deepcopy(original);value['reference']=p.record(path)
            with self.assertRaises(p.InvalidArtifact):p.check_request(value)

    def test_no_executable_handoff_on_corrupted_state_or_membership(self):
        original=p.read_json(FIX/'triple_000/PREPARATION.json')
        for field in ('denominator','physical_state_anchor'):
            value=copy.deepcopy(original);value[field]=2 if field=='denominator' else value['cases'][1]['case_id']
            with tempfile.TemporaryDirectory() as d:
                path=Path(d)/'explicitly_corrupted_real_preparation.json';p.write_new(path,value)
                with self.assertRaises(p.InvalidArtifact):p.dry_run(path)

if __name__=='__main__':unittest.main()
