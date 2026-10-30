#!/bin/bash
#SBATCH --job-name=lady-almo-memory-aquo8_central_La
#SBATCH --account=banfield
#SBATCH --partition=standard
#SBATCH --qos=standard_expand
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=256G
#SBATCH --no-requeue
set -euo pipefail
exec /home/jwestrob/.pyenv/versions/3.12.0/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lady_electronic_decomposition_20261023/cp2k_memory_recovery_v1/RUN.py /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lady_electronic_decomposition_20261023/cp2k_memory_recovery_v1/CONFIG.json /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lady_electronic_decomposition_20261023/cp2k_quartet_memory_v1/tasks/aquo8_central_La aquo8_central_La
