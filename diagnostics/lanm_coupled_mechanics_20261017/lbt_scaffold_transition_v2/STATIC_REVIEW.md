# Independent V2 representation-replay review

Reviewed 2026-10-19T07:49:59Z by `/root/mex_health_review`.

**PASS for the narrow technical correction. Actual Slurm replay remains pending.**
Only source, saved JSON/text, diffs and hashes were inspected. No scientific
imports, numerical calculations, model calls or submission were performed.

## Frozen pins

| File | SHA256 |
|---|---|
| ANALYZE.py | `03e982cd63dc03b9fdae06c96fc077b078d9a5789b427d7b2dc2134e10cdbbe8` |
| SPEC.md | `ef19a514d7479522c3ab14001ba492dc84928178ee88a8a18ed1a67435aa1d55` |
| CONFIG.json | `5ef51013ff7a0e6ea3d9ceed84cb989672f9b2ab4bc67b0cb113010be37e6d6c` |
| READY.json at review | `c4efb1b9de7faad466cb39c09ae308dd5496904aa5829e12a63c3d295875976a` |

All embedded CONFIG file pins passed byte-hash verification. READY reuses the
unchanged V1 launcher (`81c78eb8…`) and selects a fresh V2 output directory.

## Failure diagnosis and correction

Actual job14785 result `626252ba…` failed at the whole-comparator/ALPB exact
replay gate with zero attempted contexts and queries. The V1 source, result,
SPEC and earlier static review remain unchanged. The earlier review missed
this operation-order mismatch.

The pinned producers have the same force sign and physical units but distinct
floating-point conversion expressions:

- REFERENCE.py `admit_alpb` (source `c22c9896…`, line243) stores the solvent
  gradient as `(g * HA) / BOHR`.
- MODEL_COMPLETION.py `package` (`f35e2225…`, line260) stores solvent force as
  `-g * (HA / BOHR)`.
- The actual full collector CONTINUATION.py (`06e40d37…`, lines278–286) replays
  the second expression and copies its admitted whole composite into the table.

Saved text already exposes the invalid cross-representation equality: before-La
atom0 ALPB force x is `+0.5506655383113366`, whereas the stored reference solvent
gradient x is `-0.5506655383113365`. No numerical residual was calculated here.

V2 keeps both original representations. It requires exact raw-MACE conversion,
exact grouped-factor solvent-force conversion, exact component sum and exact
whole-composite/table replay. It separately requires exact sequential reference
gradient conversion from the same admitted Hartree/Bohr vector. Original
collector/model/constant sources and the zero-query V1 failure are bound.
No tolerance is substituted for any exact comparison; admitted output arrays
are preserved.

## Unchanged scope and remaining gate

The source/config diffs leave all four physical/state joins, missing comparators,
187-interaction inventory, CMAP handling, MM closures, fixed directions,
reporting and interpretation unchanged. The ceiling remains 36 contexts and
816 queries, with original 1e-8 kcal/mol energy and 1e-5 kcal/mol/Angstrom force
closure thresholds. The one-CPU Slurm launcher, no-requeue setting and thread
limits are unchanged. No new molecular method, geometry, parameter or fit is
introduced.

No outstanding static blocker was found. Root should judge the actual V2 replay
and complete four-cell closure evidence before admitting a transition result;
this review establishes neither numerical execution success nor physical
accuracy. All inherited reference-state and force-qualification caveats remain.
