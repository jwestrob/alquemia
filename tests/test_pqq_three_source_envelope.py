"""Pinned real source preparations only; no fabricated successful molecular output."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import pqq_three_source_envelope as p

BASE=ROOT/'workspaces/pqq_three_source_envelope_20260923/plm_v1'
NAMES=('PQQSEQ_07ab500e3df76b30d71c','PQQSEQ_83440678cbbd658047c9')


class ThreeSourceEnvelope(unittest.TestCase):
    def test_actual_calibration_and_unknown_source_identity(self):
        for name in NAMES:
            req=p.check_request(p.read_json(BASE/(name+'_REQUEST.json')))
            self.assertEqual(req['reference']['sha256'],p.REFERENCE_SHA)
            self.assertEqual(req['reference_status']['status'],'available')
            self.assertEqual(req['physical_state_anchor'],min(c['case_id'] for c in req['cases']))
            self.assertFalse(req['anchor_is_biological_label_or_calibration_member'])
            self.assertTrue(all(v['expected_class'] is None for v in req['source_evidence'].values()))
            self.assertTrue(all('canonical_coordinate_match' not in v for v in req['cases']))

    def test_real_preparation_whole_state_pairing_and_strict_scalar_recipe(self):
        for name in NAMES:
            path=BASE/name/'PREPARATION.json';data=p.read_json(path);check=p.dry_run(path)
            self.assertEqual(check['status'],'prepared');self.assertEqual(check['supported_sources'],3)
            self.assertEqual(len(data['origin_MACE_tasks']),6);self.assertEqual(len(data['origin_scalar_tasks']),12)
            self.assertEqual(check['prospective_calls'],p.finite_calls(True))
            self.assertFalse(check['execution_implemented']);self.assertFalse(check['new_protonation'])
            self.assertEqual(check['new_molecular_calls'],0)
            for row in data['cases']:
                kin=p.Kinematics(p.read_json(p.verify(row['maps']['Ca']))['context'])
                coord=kin.evaluate(np.zeros(len(kin.modes)))[1]
                for z in ('Ca','La'):
                    xyz=p.xyz(p.verify(row['representations']['context']['endpoints'][z]['xyz']))
                    np.testing.assert_allclose(coord,[a[1:] for a in xyz],rtol=0,atol=1e-12)
            for task in data['origin_scalar_tasks']:
                text=p.verify(task['input']).read_text()
                self.assertEqual(text,p.scalar.recipe(task['charge'],task['medium'],'fresh'))
                self.assertIn('TolE 1e-10',text);self.assertIn('NoAutostart',text)

    def test_exact_discovery_union_and_immutable_old_preparations(self):
        for name in NAMES:
            data=p.read_json(BASE/name/'PREPARATION.json');discovery=p.read_json(p.verify(data['discovery']))
            sources=p.read_json(p.verify(data['source_preparation']))['cases']
            actual=[p.envelope.discover_source(s,data['config']) for s in sources]
            self.assertEqual(discovery['cases'],json.loads(json.dumps(actual)))
            self.assertTrue(all(x['anchor_identity_replayed'] and x['original_radius_A']==3.5 and x['envelope_radius_A']==4.3 for x in actual))
            old=ROOT/'workspaces/pqq_three_source_20260923/plm_v1'/name/'PREPARATION.json'
            self.assertEqual(p.original.dry_run(old)['status'],'prepared')

    def test_missing_reference_stops_before_any_molecular_task(self):
        data=p.read_json(BASE/(NAMES[0]+'_REQUEST.json'))
        value=copy.deepcopy(data);value['reference']=None;value['reference_status']=p.reference_status(None,value['config'])
        self.assertEqual(p.check_request(value)['reference_status']['status'],'unsupported_reference')
        self.assertIsNone(value['reference_status']['bands'])
        self.assertTrue(all(n==0 for n in p.finite_calls(False).values()))
        prepared=p.read_json(BASE/NAMES[0]/'PREPARATION.json')
        with tempfile.TemporaryDirectory(dir=ROOT/'workspaces',prefix='envelope_missing_reference_') as td:
            request=Path(td)/'explicit_missing_reference_request.json';p.write_new(request,value)
            result=p.prepare(request,p.verify(prepared['source_archive']),Path(td)/'blocked_preparation')
            blocked=p.read_json(p.verify(result['preparation']))
            self.assertEqual(result['status'],'unsupported_reference')
            self.assertEqual(blocked['origin_MACE_tasks'],[]);self.assertEqual(blocked['origin_scalar_tasks'],[])
            self.assertIsNone(blocked['union']);self.assertEqual(blocked['new_molecular_calls'],0)
        old=p.read_json(ROOT/'workspaces/pqq_three_source_20260923/plm_v1'/(NAMES[0]+'_REQUEST.json'))
        with self.assertRaises(p.InvalidArtifact):p.reference_status(old['reference'],value['config'])

    def test_corrupted_real_request_state_anchor_and_membership_rejected(self):
        original=p.read_json(BASE/(NAMES[0]+'_REQUEST.json'))
        for mutation in ('anchor','condition','member','config'):
            value=copy.deepcopy(original)
            if mutation=='anchor':value['physical_state_anchor']=value['cases'][1]['case_id']
            elif mutation=='condition':value['cases'][0]['source_conditioning_metal']='Ca'
            elif mutation=='member':value['cases'][1]=copy.deepcopy(value['cases'][0])
            else:value['config']['model']['dtype']='float32'
            with self.assertRaises(p.InvalidArtifact):p.check_request(value)

    def test_corrupted_real_preparation_cannot_emit_incomplete_or_changed_tasks(self):
        original=p.read_json(BASE/NAMES[0]/'PREPARATION.json')
        for mutation in ('denominator','state','missing','recipe'):
            value=copy.deepcopy(original)
            if mutation=='denominator':value['denominator']=2
            elif mutation=='state':value['cases'][1]['state_key']='explicitly_corrupted'
            elif mutation=='missing':value['status']='unavailable';value['missing_members']=[value['cases'][0]['case_id']]
            else:value['origin_scalar_tasks'][0]['charge']+=1
            with tempfile.TemporaryDirectory() as td:
                path=Path(td)/'explicitly_corrupted_real_preparation.json';path.write_text(json.dumps(value))
                with self.assertRaises(p.InvalidArtifact):p.dry_run(path)


if __name__=='__main__':unittest.main()
