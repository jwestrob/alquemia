# Reusable prepared-input descriptor interface

Implemented `scripts/mace_omol_prepared.py` with audit, prepare and report
operations, using the existing executor for dry-run, allocated execute and
collection. No production/default change. No new molecular forward was launched
for this interface work; it reuses actual job1200828 receipts explicitly.

## What runs

- Exact source replay succeeds for real MxaF1H4I and GGR1GLG preparations,
  including protein topology, hydrogen handling, cofactor/waters, formal charge,
  atom mapping and paired coordinates. No force or energy is calculated by this
  preparation check.
- The PQQ interface manifest has zero new tasks and four explicitly verified
  reuses. Its complete report exactly reproduces17.579955566129197 model kcal.
  Classification stays unavailable because no compatible calibration is supplied.
- A GGR manifest with four new tasks and zero reuses is prepared as a dry-run
  example only. It has not been submitted; do not duplicate the existing GGR
  calculations merely to test the interface.
- Four real-fixture tests pass28.325s. They include rejection of the actual1KB0
  false peptide connections and an explicitly corrupted copy of real coordinates.
  No fabricated successful scientific fixture is used.

## Environment and cost

Preparation replay needs OpenMM, already installed in the existing `lanm_qmmm`
preparation environment. The MACE worker environment intentionally has no OpenMM.
The first audit attempt used the worker Python and failed at import, before any
preparation or model forward; its empty output/timing record is preserved asv1.
The existing preparation Python then succeeded; no environment was changed.

The manifest records its preparation driver executable and OpenMM/NumPy
versions. Execute the driver with that recorded Python. The established executor
launches actual GPU workers with the separately pinned MACE Python stored in
`software`, so no MACE package or checkpoint is substituted.

Measured1H4I audit4.50wall/4.29CPU seconds; prepare with source verification
and four exact reuses21.22wall/19.28CPU seconds; report16.83wall/15.00CPU
seconds. Peak processRSS <=591844KiB for these operations. Local work is
additional to earlier GPU receipts; no inference time is implied here.

## Scope

Input is an existing source-backed whole-chain preparation.json. This does not
silently convert an arbitrary structure/XYZ or supply absent protonation.
The implemented scope is one selected metal-bearing chainA with the recorded
protein/cofactor/site-water inventory. Source exclusions remain visible in the
audit; replay verifies the declared preparation, not every chemical assumption.
Multisite source chains and unsupported assemblies fail explicitly. A PQQ
manifest cannot be relabeled as non-PQQ to drop its cofactor.

Exact1024edge/product execution and the charge-feature mask must match the
completed development qualification. Native energies cannot satisfy descriptor
caches. A new manifest contains four states or explicit matching reuses, keeps
missing results unavailable, and never applies old ORCA bands or a universalzero.
Outputs are modified learned descriptors, not quantum energies or free energies.

The canonical calibration job1200830 is still separate and running. A future
passing compatible report is required before implementing research classification
for these PQQ inputs. No generic/non-PQQ affinity band is provided.

Artifacts under `workspaces/mace_omol_20260917/`:
`prepared_interface_pqq_v1`, `prepared_interface_pqq_report_v1`,
`prepared_interface_ggr_fresh_v1`, plus audit/test/preparation/report receipts.
