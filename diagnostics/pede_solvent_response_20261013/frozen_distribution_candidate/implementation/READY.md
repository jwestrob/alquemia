# Native component adapter — prepared, unexecuted

The adapter returns qualified native `G(qvac)` components. Root owns candidate
selection, the unchanged canonical calibration rule and grouped comparisons.
It does not create a substituted legacy water-minus-vacuum energy field.

## Frozen scope

The inventory contains150 canonical and1392 transfer native inputs. Canonical
requires300 component evaluations. Four exact PedE transfer inputs reuse the
eight qualified v3 results;1388 other transfer inputs require2776 new component
evaluations. These are saved-charge evaluations only: zero SCC, DFT, model,
gradient, optimization or folding calls. Atom counts and electronic states come
from each exact inventory entry; no192-atom assumption remains.

`PREPARE.py` freezes metadata and copies the implementation. It can run on the
login host because it hashes existing artifacts without scientific parsing.
All source-state/coordinate/charge parsing, compilation and evaluation are inside
the Slurm-only `NATIVE.py` entry point. No such test or evaluation was run during
preparation. The first actual canonical execution compiles only the thin wrapper
against the successful v3 static libraries and module files, then runs the
declared endpoint qualification/replays. No engine or environment is rebuilt.

## Admission and artifacts

The exact native settings, eight-decimal charge bound formulas,1e−9Eh component
allowance and1e−6Eh informative ceiling are inherited byte-for-byte from v3.
For each input, the archived water-charge replay must reproduce every printed
component before its vacuum-charge replay starts. The two distributions must
have byte-identical native cavity and Born-matrix files. Each vacuum-charge
component must also satisfy the unchanged informative uncertainty ceiling;
root confirmed this before any candidate evaluation. Every failure retains
its source, raw logs, bounds/residuals and explicit unavailable cell; no baseline
fallback or retry is implemented. Four reused pairs are rechecked against the
same input identities and diagonal gates without native calls.

`FINAL_COLLECTION.json` contains `protocol`, `manifest`, `inventory`, `stage`,
`build`, `terminal`, `status` and `tasks`. Each task maps to `status` (available or
unavailable), `G_qvac_hartree` (number or null), `result` pin and `reason`.
`status: complete` means all finite tasks reached receipts; it does not assert
all cells passed. The root scorer must require complete relevant pools.

Each task's `RESULT.json` preserves the full source identity, original archived
components, native water/vacuum results, qualification checks and exact source
charge/cavity/matrix pins. GNU time records per-child CPU/wall/RSS rather than
subtracting process-global resource counters during parallel calls. The process
pool has at most the actual Slurm allocation's12 workers, each native child with
one numerical thread. Signals stop owned calls and retain unavailable receipts.
Fresh output directories prevent historical overwrite; no automatic restart
or hidden recovery is implemented.

## Two-stage command contract

After metadata preparation, root supplies the manifest SHA and a fresh output:

```bash
sbatch /absolute/workspace/implementation/RUN.sbatch canonical \
  /absolute/workspace/manifest.json MANIFEST_SHA /absolute/workspace/canonical_v1
```

The transfer command uses the same manifest and canonical wrapper build:

```bash
sbatch /absolute/workspace/implementation/RUN.sbatch transfer \
  /absolute/workspace/manifest.json MANIFEST_SHA /absolute/workspace/transfer_v1 \
  --calibration /absolute/CALIBRATION.json --calibration-sha CALIBRATION_SHA
```

Root's calibration receipt must contain:

```json
{
  "protocol": "frozen_distribution_candidate_v1",
  "inventory": {"path": "...", "sha256": "..."},
  "native_manifest": {"path": "...", "sha256": "..."},
  "canonical_native_collection": {"path": "...", "sha256": "..."},
  "calibrations": {"operational": {"status": "available"}}
}
```

Transfer rechecks these identities and all150 canonical native results before
any dependent native dispatch. The full root receipt also records the other
policies and frozen bands. A failed canonical calibration closes dependent
transfer for this model; it does not authorize changing the gates.

## Exact prepared manifest

`/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/manifest.json`

SHA256: `bcd9fc9d8923139dfcb8d7b8d0191c5c29484e700ddc8ffd488c1dea09c10b41`. Exact canonical command and gated transfer template are
in `READY.json`. All390 compiler/module/library/header pins are from the finite
qualified build include/runtime paths; no archive-wide hash sweep was needed.
Curie independently reviewed the final source and manifest binding. Execution
qualification remains the actual Slurm stage, not this static review.
