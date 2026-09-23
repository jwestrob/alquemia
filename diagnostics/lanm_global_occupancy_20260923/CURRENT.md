# Small LanM occupancy pilot — queued, not yet a scientific result

Jacob explicitly requested a few two-ion and four-ion structures, beginning with
Hans-LanM, before extending across the library. He then requested email when a
meaningful result becomes available. The [frozen plan](PLAN.md) is the finite
execution scope. Existing PQQ production and manuscript delivery are unchanged.

## Prepared and checked

- Three complete chain-A structures: Hans8DQ2, Hans8FNR, Mex8FNS.
- EF1+EF2, EF2+EF3, EF1–EF4 for each source: nine conditional systems,18 La/Dy endpoints.
- No synthetic carve atoms, changed proton inventory, invented experimental labels,
  extra folding, new DFT or reserved SpyCI-LAMBS scoring.
- Five preparation tests and seven scorer tests pass on real pinned fixtures.
  These are preparation/parser/algebra checks, not successful molecular evaluations.
- One warm native MACE proposal per metal, all physical coordinates allowed within
  the declared small box. Every admitted geometry is cross-scored for both metals
  using native MACE + nativeGFN2(ALPB−vacuum).
- Dy physical multiplicities11/21 and native effective singlet are distinct,
  explicitly checked states. Multiple-spin coupling/ground-state validity remains
  untested. This is complete-monomer motion, not folding or dimerization equilibrium.
- Different source water inventories are retained. No raw cross-source energies
  or two/four-ion populations are compared.

## Jobs and exact continuation

At this checkpoint both jobs are pending; no new molecular evaluation has run.

1. **1213018**: first Hans8DQ2 EF12 proposal/common pool plus four native origin
   scalar calls. One H200,32 CPUs,200000MiB; at most24 native MPI ranks +8 MACE host
   threads simultaneously. Native/adapter full-source analytic forces are compared.
2. **1213040**, `afterok:1213018`, kill invalid dependency: verifies actual first
   complete MACE matrix, adapter equality, and four normal nativeSCF/state receipts.
   Only then runs the other eight systems on one H200/32 CPUs/200000MiB. It creates
   the exact remaining native scalar manifest and submits `run_native_matrix.sbatch`.
3. That future CPU-only stage uses32 CPUs/200000MiB, four concurrent eight-rank
   native endpoints, with four exact successful origin reuses. Maximum original
   matrix is54 composite endpoints/108 native scalar calls before deduplication;
   MACE optimization/qualification calls are counted separately. No retries or
   additional optimizations are automatically added.
4. CPU completion collects valid/failed cells, writes the conditional comparison
   and vault note, and emails the actual report to Jacob. Failed executor return
   codes do not suppress partial-result collection. Delivery is recorded as local
   relay acceptance, not assumed inbox delivery.

Do not duplicate these jobs or restart the old completed LanM pocket study.
Inspect actual receipts before interpreting results. A failed required cell/search
leaves that accommodated result unavailable; any valid static contrast remains a
separate field. No favorable-label gate controls the continuation or email.

## Artifacts

All paths below are beneath the repository root:

- `workspaces/lanm_global_occupancy_20260923/prepared_v1/manifest.json`
- `workspaces/lanm_global_occupancy_20260923/scoring_v2/manifest.json`
- `workspaces/lanm_global_occupancy_20260923/native_feasibility_v2/manifest.json`
- `workspaces/lanm_global_occupancy_20260923/AUTOMATION.json`

`scoring_v1`/`native_feasibility_v1` were preparation-only drafts, never submitted;
the v2 snapshot fixes missing-cell/search handling before any molecular call.
The native Hamiltonian and declared scientific search did not change.

Future result paths (not evidence of completed calculations yet):

- `workspaces/lanm_global_occupancy_20260923/RESULT_v2.json`
- `diagnostics/lanm_global_occupancy_20260923/PILOT_REPORT.md`
- `diagnostics/lanm_global_occupancy_20260923/pilot_result_email_receipt.json`
- Vault `agent-captures/2026-09-23_Nikasha-whole-chain-LanM-occupancy-pilot.md`

Root owns scoring/integration/jobs; Khoury's preparation and independent code
review are complete. Preparation committed73e9c7f. No remote push or baseline
promotion. Scheduler walltime defaults are not an agent-imposed compute budget.
