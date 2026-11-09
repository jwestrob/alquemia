# Exact failing input and recomputation recipe

`minimal_start0_A.extxyz.gz` is the exact 7,288-atom input used by the native
scout, compressed without changing its coordinates or atom order. Its hashes
and state are in MINIMAL_CASE.json. It contains no experimental affinity values.
This is a new solvated development system; no archived fixed-panel fixture is
replaced by it.

The default `reproduce.sh` performs saved-data analyses only. The commands below
recompute the recorded failure and are not submitted automatically. Run from
the repository root with the existing pinned environment. They change only
input/output paths in the original settings.

```bash
workspaces/so3lr_tauri_screen_20261106/software/env/bin/python - <<'PY'
from pathlib import Path
import gzip, hashlib, json, yaml
d = Path('diagnostics/so3lr_neutral_pair_20261107')
w = Path('workspaces/tauri_minimal_failure_20261109')
w.mkdir(exist_ok=False)
b = gzip.decompress((d/'minimal_start0_A.extxyz.gz').read_bytes())
assert hashlib.sha256(b).hexdigest() == json.loads((d/'MINIMAL_CASE.json').read_text())['uncompressed_sha256']
(w/'input.extxyz').write_bytes(b)
c = yaml.safe_load((d/'md_scout.yaml').read_text())
c.update(input_file=str(w/'input.extxyz'), output_file=str(w/'trajectory.hdf5'),
         log_file=str(w/'native.log'), restart_save_path=str(w/'restart.npz'))
(w/'settings.yaml').write_text(yaml.safe_dump(c))
PY

sbatch -N1 -n1 -c12 --partition=gpu_h200 --gres=gpu:1 --qos=standard_expand \
  --output=workspaces/tauri_minimal_failure_20261109/slurm_%j.out \
  --error=workspaces/tauri_minimal_failure_20261109/slurm_%j.err \
  --wrap='export OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JAX_PLATFORMS=cuda,cpu JAX_ENABLE_X64=true XLA_PYTHON_CLIENT_PREALLOCATE=false; timeout --kill-after=10s 3600s srun -n1 -c12 --cpu-bind=cores workspaces/so3lr_tauri_screen_20261106/software/env/bin/python -m so3lr.cli.so3lr_cli --settings workspaces/tauri_minimal_failure_20261109/settings.yaml'
```

Expected observation: the trans carboxylate group's original C6–C7, C7–O8 and
C7–O9 distances leave the declared integrity range during native FIRE, before
MD. Compare minimum-image distances in the 42 Å periodic cell. The capped
ASE-only comparison is separately archived; its stopped result is an integrity
threshold crossing, not a completed scission calculation.
