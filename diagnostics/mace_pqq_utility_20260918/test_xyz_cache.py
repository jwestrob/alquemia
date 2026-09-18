"""Strict-parser cache tests on actual PQQ XYZ files and corrupted copies."""
from pathlib import Path
import json
import os
import shutil
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'workspaces/mace_pqq_utility_20260918/interface_source_v5/implementation'))
import affordable_common as common
from mace_file_checks import cached_file_checks
PANEL=ROOT/'workspaces/mace_omol_20260917/intact_panel_prepared_v2/preparation_manifest.json'


class ActualXYZCacheTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.preparations=[common.read_json(common.verify(r['preparation'])) for r in common.read_json(PANEL)['rows'] if r['status']=='prepared']
        cls.xyz=common.verify(cls.preparations[0]['endpoints']['La']['xyz'])

    @cached_file_checks
    def test_all_real_paired_coordinates_are_bitwise_equal(self):
        count=0
        for p in self.preparations:
            for e in p['endpoints'].values():
                path=common.verify(e['xyz'])
                self.assertEqual(common.xyz(path),common._xyz_uncached(path))
                self.assertEqual(common.xyz(path),common._xyz_uncached(path))
                count+=1
        self.assertEqual(count,54)

    def test_parser_runs_once_per_operation_and_returns_independent_lists(self):
        original=common._xyz_uncached;expected=original(self.xyz)
        with patch.object(common,'_xyz_uncached',wraps=original) as counted:
            @cached_file_checks
            def operation():
                first=common.xyz(self.xyz);first.pop()
                self.assertEqual(common.xyz(self.xyz),expected)
            operation();self.assertEqual(counted.call_count,1)
            operation();self.assertEqual(counted.call_count,2)
        self.assertEqual(common.xyz(self.xyz),expected)

    def test_actual_file_corruption_with_restored_mtime_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'corrupted_real_xyz.xyz';shutil.copyfile(self.xyz,path)
            time.sleep(2.1)
            @cached_file_checks
            def operation():
                common.xyz(path);before=path.stat();lines=path.read_text().splitlines()
                fields=lines[2].split();fields[1]='nan';lines[2]=' '.join(fields)
                path.write_text('\n'.join(lines)+'\n')
                os.utime(path,ns=(before.st_atime_ns,before.st_mtime_ns))
                with self.assertRaisesRegex(common.InvalidArtifact,'nonfinite coordinates'):
                    common.xyz(path)
            operation()

    def test_deleted_real_file_is_not_served_from_cache(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'real_xyz_copy.xyz';shutil.copyfile(self.xyz,path)
            @cached_file_checks
            def operation():
                common.xyz(path);path.unlink()
                with self.assertRaises(FileNotFoundError):common.xyz(path)
            operation()


if __name__=='__main__':unittest.main()
