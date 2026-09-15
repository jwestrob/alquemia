#!/bin/bash
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
while [ ! -s /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/global_representation_20260915/terminal_1199003.json ]; do sleep 30; done
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/affordable_global_collect.py --manifest /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/global_representation_20260915/retry_64rank_v1/global_manifest.json --output /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/global_representation_20260915/collection_1199003.json
