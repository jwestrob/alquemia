"""File-cache integrity tests use explicitly corrupted copies of a real manifest."""
from pathlib import Path
import os
import shutil
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import affordable_common as common
from mace_file_checks import cached_file_checks
SOURCE=ROOT/'workspaces/mace_omol_20260917/intact_core_v1/manifest.json'


@unittest.skipUnless(SOURCE.exists(),'requires real prepared OMOL manifest')
class FileCheckTests(unittest.TestCase):
    def test_nested_digest_reuse_is_scoped_and_preserves_real_hash(self):
        ref=common.record(SOURCE);original=common.digest
        with patch.object(common,'digest',wraps=original) as counted:
            @cached_file_checks
            def inner():return common.record(SOURCE)
            @cached_file_checks
            def outer():
                self.assertEqual(common.record(SOURCE),ref)
                self.assertEqual(inner(),ref)
                common.verify(ref)
            outer();self.assertEqual(counted.call_count,1)
            outer();self.assertEqual(counted.call_count,2)
            self.assertIs(common.digest,counted)
        self.assertIs(common.digest,original)

    def test_same_size_corruption_with_restored_mtime_is_rejected(self):
        for delay in (0,2.1):
            with self.subTest(delay=delay),tempfile.TemporaryDirectory() as directory:
                p=Path(directory)/'corrupted_real_manifest.json';shutil.copyfile(SOURCE,p)
                if delay:time.sleep(delay)
                ref=common.record(p);before=p.stat()
                @cached_file_checks
                def check():
                    common.verify(ref);data=bytearray(p.read_bytes());data[-2]^=1;p.write_bytes(data)
                    os.utime(p,ns=(before.st_atime_ns,before.st_mtime_ns))
                    with self.assertRaises(common.InvalidArtifact):common.verify(ref)
                check()

    def test_replacement_and_deletion_do_not_reuse_old_digest(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'copy.json';shutil.copyfile(SOURCE,p);ref=common.record(p)
            @cached_file_checks
            def check():
                common.verify(ref);other=p.with_name('corrupted_copy.json')
                data=bytearray(p.read_bytes());data[-2]^=1;other.write_bytes(data);os.replace(other,p)
                with self.assertRaises(common.InvalidArtifact):common.verify(ref)
                p.unlink()
                with self.assertRaises(common.InvalidArtifact):common.verify(ref)
            check()

    def test_exception_restores_original_verifier(self):
        original=common.digest
        @cached_file_checks
        def fail():
            common.verify(common.record(SOURCE));raise RuntimeError('intentional verification abort')
        with self.assertRaises(RuntimeError):fail()
        self.assertIs(common.digest,original)


if __name__=='__main__':unittest.main()
