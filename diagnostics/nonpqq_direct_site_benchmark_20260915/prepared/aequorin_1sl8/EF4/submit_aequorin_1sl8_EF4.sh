#!/bin/bash
#SBATCH -p memory
#SBATCH -J orca_aequorin_1sl8_EF4
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/nonpqq_direct_site_benchmark_20260915/prepared/aequorin_1sl8/EF4/slurm_%j.out
#SBATCH -e /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/nonpqq_direct_site_benchmark_20260915/prepared/aequorin_1sl8/EF4/slurm_%j.err

set -euo pipefail
ORCA_PATH=/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export ORCA_PATH
export PATH=$ORCA_PATH:$PATH
export LD_LIBRARY_PATH=$ORCA_PATH:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/nonpqq_direct_site_benchmark_20260915/prepared/aequorin_1sl8/EF4
python3 /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/scripts/run_orca_task_manifest.py aequorin_1sl8_EF4_carve_manifest.json     --expected-runner-sha256 b9278e74e5317e859e8cb05541e1ec57cb335b26dcdf65a83b89675407343b9a
