# Ca-conditioned preparation: 106 new sources admitted

All 107 missing requests were attempted using the existing frozen canonical
La-source4.3 Å union and unchanged source-graph/physical-map builders.
**106 passed.** Together with the three saved probes,109 of 125 Ca-conditioned
sources are ready; 15 original PQQ exclusions and one new preparation failure
remain visible. No energies, forces, proposals or optimizations were evaluated.

The new failure is `mmol_2048-pqq-la_model__conditioned_Ca__seed-1_sample-0`:
`unsupported disulfide partner`. The unchanged cysteine builder requires a
CYS partner with the specified hydrogen-parent inventory. Its failed admission
is preserved; neither a missing partner nor a unique physical defect is inferred
from that message. No altered fragment or protonation was substituted.

## What was preserved and checked

- Exact archived per-protein union, with its original first-member state anchor.
- Source heavy-atom geometry, core chemistry, water/proton inventory, formal
  charges and paired Ca/La coordinates under the existing builder checks.
- Exact state/protein keys, source atom/cap/bond identities, physical maps and
  successful mapping-derivative gates.
- The actual125-source denominator and original class labels; no new threshold
  or source selection based on computed preference.

The three existing Ca probes were reconstructed before attempting new sources.
All passed the established 1e-12 Å coordinate-copy rule; maximum difference was
8.88e-16 Å in copied cap arithmetic. **Their rebuilt XYZ hashes differ**, so
only the original three saved artifacts are used for exact scientific cache
reuse. Source/state/discrete mapping identity remains exact.

The first two wrapper attempts stopped during those replay checks: their
comparison initially omitted copied cap vectors/derived lengths and numerical
map-check residuals from cross-runtime roundoff handling. Version 3 compares
those coordinate/diagnostic quantities at 1e-12 while retaining unchanged builder
admission gates and exact noncoordinate metadata. Failed v1/v2 outputs and code
snapshots are preserved; no 107-source work ran in those attempts.

## Ready artifacts and next decision

Final workspace:
`workspaces/orca_parameter_audit_20261012/ca_conditioned_preparation_v3/`

- `PREPARATION.json`: all 125 statuses, 106 new admitted rows and receipt.
- `ORIGIN_INPUTS.json`: 212 paired origin inputs, exact maps/state pins; inputs
  only, not an automatic molecular-execution manifest.
- `REPLAY_CHECK.json`: all three real-probe comparisons.
- `manifest.json` and `implementation/`: frozen scope/runtime/code/plan.

The predeclared first sorted admissible scout is
`a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-0`.
Its two endpoints and mappings are pinned in `future_scout` and `scout_tasks`.
Root owns the CPU-throughput feasibility test and any later proposal execution.
At most 212 searches and 1272 GFN2 inputs per backend remain for the 106 new
sources, conditional on successful proposal/runtime qualification. No new
Ca-specific context discovery or calibration is introduced.

Successful v3 preparation/replay took 117.883 seconds wall and 219.423 process
CPU-seconds, with two one-thread workers. Peak observed child RSS was 125260 KiB;
parent RSS 159080 KiB. This receipt excludes initial imports and the earlier
failed replay attempts, whose process costs were not instrumented.

Reproduction into a new unused workspace:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /home/jwestrob/miniconda3/envs/fep/bin/python diagnostics/orca_parameter_audit_20261012/ca_conditioned_preparation/PREPARE.py --output workspaces/orca_parameter_audit_20261012/ca_conditioned_preparation_replay_v1 --workers 2
```

This command performs preparation only. Do not repeat completed work merely to
start molecular evaluation; use the already pinned v3 artifacts. Production and
all existing scores remain unchanged.
