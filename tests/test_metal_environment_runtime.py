"""Resource algebra and immutable-template rendering; no molecular executions."""
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from metal_environment_runtime import allocation_policy, render_runtime_input, exclusive_node_memory

ENV = {'SLURM_CPUS_ON_NODE':'64', 'SLURM_NTASKS':'64',
       'SLURM_MEM_PER_NODE':'262144'}


class AllocationTests(unittest.TestCase):
    def test_actual_exclusive_node_ignores_smaller_request(self):
        captured=json.loads((ROOT/'diagnostics/metal_environment_response_20260926/RESOURCE_POLICY_20260927.json').read_text())
        self.assertEqual(exclusive_node_memory(captured['job'],captured['node'],captured['partition']),8256990)
        env=dict(ENV,SLURM_JOB_ID='1219207',SLURM_CPUS_ON_NODE='344',SLURM_NTASKS='344')
        with patch('metal_environment_runtime.subprocess.check_output',side_effect=[captured['job'],captured['node'],captured['partition']]):
            result=allocation_policy(env,workers=6,ranks=57)
        self.assertEqual(result['allocated_memory_MiB'],8256990)
        self.assertEqual(result['workers']*result['ranks_per_worker'],342)
        self.assertGreater(result['maxcore_MB'],18000)
        shared=captured['partition'].replace('OverSubscribe=EXCLUSIVE','OverSubscribe=NO')
        self.assertIsNone(exclusive_node_memory(captured['job'],captured['node'],shared))

    def test_shared_zero_request_uses_observed_available_memory(self):
        # Resource-only scheduler fixtures; no scientific data or energies.
        env=dict(ENV,SLURM_JOB_ID='1219782',SLURM_CPUS_ON_NODE='24',SLURM_NTASKS='24',SLURM_MEM_PER_NODE='0')
        job='NumNodes=1 NumCPUs=24 NodeList=node-48-384g-1 Partition=standard-shared'
        node='CPUTot=48 RealMemory=385557'
        partition='OverSubscribe=YES'
        with patch('metal_environment_runtime.subprocess.check_output',side_effect=[job,node,partition]), patch('metal_environment_runtime.Path.read_text',return_value='MemAvailable: 83886080 kB\n'):
            result=allocation_policy(env,workers=2,ranks=12)
        self.assertEqual(result['allocated_memory_MiB'],81920)
        self.assertIn('not exclusive entitlement',result['memory_source'])
        self.assertEqual(result['workers']*result['ranks_per_worker'],24)

    def test_approved_256_gib_allocation(self):
        p = allocation_policy(ENV)
        self.assertEqual((p['workers'],p['ranks_per_worker'],p['maxcore_MB']), (4,16,3200))
        self.assertAlmostEqual(p['per_rank_decimal_MB_before_rounding'],3221.225472)
        p = allocation_policy(dict(ENV,SLURM_MEM_PER_NODE='256000'))
        self.assertEqual(p['maxcore_MB'],3100)
        self.assertAlmostEqual(p['per_rank_decimal_MB_before_rounding'],3145.728)

    def test_per_cpu_memory_and_mpi_slots(self):
        env = dict(ENV); del env['SLURM_MEM_PER_NODE']; env['SLURM_MEM_PER_CPU']='4096'
        self.assertEqual(allocation_policy(env)['maxcore_MB'],3200)
        p = allocation_policy(dict(ENV,SLURM_NTASKS='1'))
        self.assertEqual((p['workers'],p['ranks_per_worker']),(1,1))
        p = allocation_policy(dict(ENV,SLURM_NTASKS='32',SLURM_CPUS_ON_NODE='32'))
        self.assertEqual((p['workers'],p['ranks_per_worker']),(4,8))

    def test_missing_zero_and_multinode_evidence_fail(self):
        for change in ({'SLURM_MEM_PER_NODE':'0'}, {'SLURM_CPUS_ON_NODE':'0'},
                       {'SLURM_JOB_NUM_NODES':'2'}):
            with self.assertRaises(ValueError):
                allocation_policy(dict(ENV,**change))
        with self.assertRaises(ValueError):
            allocation_policy({'SLURM_CPUS_ON_NODE':'64','SLURM_NTASKS':'64'})

    def test_render_real_archived_scientific_input(self):
        # Use actual immutable archived ORCA input, no invented molecular fixture.
        source = ROOT/'diagnostics/pqq_balanced_embedding_20260914/prepared/mxaf/qm33/Ca.inp'
        if not source.exists():
            self.skipTest('Archived real ORCA template unavailable')
        archived = source.read_bytes()
        with tempfile.TemporaryDirectory() as directory, patch.dict(os.environ,ENV,clear=True):
            # Resource-policy adaptation of actual archived input; preserve all
            # method, charge, spin, embedding and geometry-reference lines.
            template = Path(directory)/'source.inp'
            original = b''.join(line for line in archived.splitlines(keepends=True)
                                if not line.lower().lstrip().startswith(b'%maxcore'))
            template.write_bytes(original)
            runtime = Path(directory)/'endpoint.runtime.inp'
            record = render_runtime_input(template,runtime,nprocs=16)
            self.assertEqual(source.read_bytes(),archived)
            self.assertEqual(template.read_bytes(),original)
            self.assertEqual(record['template']['sha256'],hashlib.sha256(original).hexdigest())
            self.assertEqual(record['runtime_input']['sha256'],hashlib.sha256(runtime.read_bytes()).hexdigest())
            self.assertEqual(record,json.loads(runtime.with_suffix('.json').read_text()))
            self.assertTrue(runtime.read_text().startswith('%maxcore 3200\n'))
            self.assertEqual(record['schema_version'],'alchemical_bvs.orca_runtime_input.v1')
            with patch.dict(os.environ,{'SLURM_NTASKS':'1'}):
                with self.assertRaisesRegex(ValueError,'layout exceeds'):
                    render_runtime_input(template,Path(directory)/'bad.inp',nprocs=16)


if __name__ == '__main__':
    unittest.main()
