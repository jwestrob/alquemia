"""Real author/model fixtures only; no fabricated scientific outputs."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, record, write_new
from mace_khoury_benchmark import CASES, normalize
from mace_omol_multisite import author_mapping, build
from mace_omol_prepared import audit_preparation

WORK=ROOT/'workspaces/mace_omol_20260917'


class AuthorPreparation(unittest.TestCase):
    def test_real_author_atom_map_and_corrupted_coordinate_rejection(self):
        source=WORK/'d5sc02315g_reading_v1/SC-016-D5SC02315G-s001.pdb'
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp); normalized,sites=normalize(source,CASES['A0A7'],d)
            cfg={'author_normalization':record(d/'normalization.json')}
            mapping=author_mapping(cfg,record(normalized))
            self.assertEqual(len(mapping['source_to_normalized']),659)
            self.assertEqual([s['site_key'] for s in sites],list('BCDEFG'))
            # Explicit malformed copy of actual author fixture; no scientific energy.
            lines=normalized.read_text().splitlines(); idx=next(i for i,l in enumerate(lines) if l.startswith('HETATM'))
            l=lines[idx];lines[idx]=l[:30]+f'{float(l[30:38])+1:8.3f}'+l[38:]
            normalized.write_text('\n'.join(lines)+'\n')
            mapping['normalized']=record(normalized);write_new(d/'corrupted_receipt.json',mapping)
            with self.assertRaisesRegex(InvalidArtifact,'changed chemistry or source coordinates'):
                author_mapping({'author_normalization':record(d/'corrupted_receipt.json')},record(normalized))

    def test_legacy_parvalbumin_policy_replays_without_changes(self):
        p=read_json(WORK/'multisite_prepared_v2/PARV_4CPV_CD/preparation.json')
        actual=build(p['recipe']['path'],p['selected_site_key'])
        self.assertEqual(actual,{k:v for k,v in p.items() if k not in ('implementation','agreement','endpoints')})

    def test_real_pH6_pair_and_complete_background_inventory(self):
        root=WORK/'khoury_author_domains_recovery_v1'
        manifest=root/'manifest.json'
        if not manifest.exists(): self.skipTest('actual author preparation is not yet complete')
        p=read_json(manifest)
        self.assertEqual(p['prepared_endpoint_forwards'],44)
        self.assertFalse(p['source_failures'])
        for name,cfg in CASES.items():
            rows=[r for r in p['rows'] if r['domain']==name]
            self.assertEqual(len(rows),cfg['sites'])
            for row in rows:
                state=read_json(row['preparation']['path']);audit=audit_preparation(row['preparation']['path'])
                self.assertEqual(audit['status'],'pass')
                self.assertEqual(len(state['background_metals']),cfg['sites']-1)
                self.assertEqual(state['endpoints']['La']['charge']-state['endpoints']['Ca']['charge'],1)
                self.assertFalse(state['explicit_waters'])


if __name__=='__main__': unittest.main()
