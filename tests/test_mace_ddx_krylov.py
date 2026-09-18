"""Real native-source preservation and frozen GMRES qualification inputs."""
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,cache_key,read_json,verify,write_new
from mace_ddx_krylov import preflight
W=ROOT/'workspaces/mace_omol_20260917'


class DDXKrylov(unittest.TestCase):
    def test_actual_built_adapter_preserves_all_native_numerical_modules(self):
        m=read_json(W/'ddx_kernel_software_v1/manifest.json')
        parent=read_json(verify(m['parent_build_manifest']))
        with tarfile.open(parent['source_archive']) as f:
            self.assertEqual(len(m['unchanged_native_modules']),19)
            for name,pin in m['unchanged_native_modules'].items():
                self.assertEqual(verify(pin).read_bytes(),f.extractfile('pyddx-0.9.0/src/'+name).read())
        receipt=read_json(W/'ddx_kernel_software_v1/attempt_0001/receipt.json')
        self.assertEqual(receipt['status'],'built_and_imported');self.assertEqual(receipt['scientific_calculations'],0)

    def test_real_pair_uses_completed_native_coefficients(self):
        m=preflight(W/'ddx_krylov_v1/manifest.json');r=read_json(verify(m['source_group']))
        self.assertEqual(m['states'],['Ca','La']);self.assertEqual(r['status'],'complete')
        self.assertEqual(r['actual_new_forward_solve_starts'],0)
        self.assertGreater(r['reciprocity_error_kcal'],.05) # actual coarse failure must remain visible
        self.assertEqual(m['solver']['rtol'],1e-10);self.assertEqual(m['tolerances']['energy_kcal'],1e-6)

    def test_corrupted_real_solver_settings_rejected_even_with_updated_digest(self):
        m=read_json(W/'ddx_krylov_v1/manifest.json');m['solver']['rtol']=1e-4
        m['cache_key']=cache_key({k:v for k,v in m.items() if k!='cache_key'})
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_solver_scope.json';write_new(p,m)
            with self.assertRaisesRegex(InvalidArtifact,'settings changed'):preflight(p)


if __name__=='__main__':unittest.main()
