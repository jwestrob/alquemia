# Ca-conditioned reference recovery — completed 12 October 2026

**All eighteen unavailable accommodated results now have separately recorded
recovery results. The full 125-source ledger gives 119 correct, 3 wrong and
3 inconclusive calls.** Coverage recovery succeeded; accommodation introduces
two wrong PedE calls, so this is not an unqualified accuracy improvement.

The original primary result remains **105 correct / 1 wrong / 1 inconclusive /
18 unavailable**. The recovered ledger is an annotated alternative-source
comparison, not a replacement for that experiment. These are correlated model
sources across 25 reference proteins, not 125 independent biological observations.
Calibration bands remain frozen. Production was not changed.

## Repairs and outcomes

| Recovery | Sources | Accommodated outcome | Static-origin outcome |
|---|---:|---|---|
| MxaF P12293 interrupted La searches | 2 | 2 correct Ca | 2 correct Ca |
| MMOL_2048 XoxF4-1 hydrogen mapping | 1 | 1 inconclusive | 1 correct La |
| ExaA Q9Z4J7 two-Ca sources | 5 | 5 correct Ca | Same |
| ExaF C5AXV8 two-Ca sources | 5 | 4 correct La, 1 inconclusive | Same |
| PedE Q88JH5 two-Ca sources | 5 | 3 correct Ca, 2 wrong La | 3 correct Ca, 2 inconclusive |

MxaF recovery continued only interrupted searches under the original policy,
reusing compatible cells. MMOL_2048 had an intact disulfide: one existing Cys128
HB2 had an overlong C–H bond and no source-graph parent. Moving only that H to
the force-field bond length repaired preparation. Heavy atoms, proton count,
water inventory and the canonical core were retained.

The fifteen original ExaA/ExaF/PedE models placed their only Ca at a distant
structural site, leaving PQQ empty. The finite two-Ca AF3 test generated five
models per protein and filled both sites in all fifteen. Actual geometry selects
the PQQ ion; both ions remain in the full structures. Preparation preserves
all source heavy coordinates and exactly reproduces the existing single-Ca
regression. The distant spectator is explicitly omitted from the local
Hamiltonian after passing the fixed distance rule. This is not a global
spectator electrostatic or scaffold-response model.

ExaA 1FLG supports two Ca sites; PedE/ExaF occupancy remains a homology-based
structural hypothesis. These are modeled alternatives, not new experimental
structures. AF3 3.0.4, seed 1, five samples and ten recycles are recorded; the
original fold runtime is unverified, so this is not a controlled one-variable
folding comparison. No source was selected by score.

## Accommodation: benefit and remaining failure

On the same 125 available sources, under each policy's frozen canonical bands:

| Policy | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Static origins | 117 | 1 | 7 | 0 |
| Accommodated, operational | 119 | 3 | 3 | 0 |
| Accommodated, mathematical minima | 119 | 3 | 3 | 0 |

Accommodation adds two correct decisions overall, but also two wrong decisions.
PedE samples 1/2 start inconclusive and become incorrectly La-like. La gains more
energy than Ca in all five PedE models, including the three correct ones. The
wrong samples do not uniquely exhibit a boundary contact or unusually large
relaxation work. Saved solver, state, paired-coordinate and selection checks pass.
No demonstrated numerical failure explains these two errors.

[The PedE inspection](two_ca_scoring/PED_E_RECOVERY_FINDINGS.md) retains each
metal's work, selected geometry, residual loads and component checks. The physical
cause remains unresolved: these outputs cannot distinguish limited search,
missing scaffold/spectator response, or electronic/composite-model error.
No new optimization, threshold change or five-fold voting rule was introduced.

All six non-correct accommodated sources remain visible:

- Wrong: A0ACD6B9F2 Ca sample 4, the previously known error; PedE two-Ca samples 1/2.
- Inconclusive: ExaF two-Ca sample 4, MMOL_1770 Ca sample 0, repaired MMOL_2048 Ca sample 0.

## Execution and checks

The existing float64 MACE-OMOL model proposes geometries with the original four
force-selected physical modes, bounded SLSQP and shared origin/Ca/La pool.
Exact solvent components select endpoints with the existing 0.1 kcal/mol origin
retention rule. Published-parameter standalone xTB 6.7.1 is independently checked
with tblite 0.7.0 under component-v2 acceptance. These are not uncorrected native
ORCA La-parameter results. MACE weights and calibration remain unchanged.
Historical absolute-v1 parity remains failed; component-v2 does not supersede it.

New preparation identity:
`Nikasha_two_Ca_source_explicit_spectator_frozen_envelope4p3_v1`.
[RESULT.json](RESULT.json) pins the other versioned adapters, canonical calibration
and full workspace ledger. [SOURCES.tsv](SOURCES.tsv) preserves all 125 decisions
and unrounded scores alongside original/recovered identity.

Actual Slurm execution includes all 15 AF3 models, their geometry/paired-preparation
checks, the old single-Ca regression, both MxaF restarts, the MMOL hydrogen repair
and scoring, all 15 new shared-pool evaluations, and final collector 3049. Every
required new scoring cell is available. Each native lane passed two archived
energy/force fixtures; required corrected-solvent component gates passed.
Independent metadata review confirmed exact original 125 rows, prior 3 repairs,
disjoint 5+10 additions, all decisions and costs. Syntax checks are distinct from
these molecular evaluations. Failed startup jobs 2877/2878 remain recorded:
node-local launcher hash mismatch, then older Apptainer namespace failure; neither
produced molecular output. Compatible H200 execution resolved startup without
changing scientific admission rules.

## Measured cost

[Scheduler receipts](RECOVERY_SCHEDULER_RECEIPTS.json) retain all 28 jobs in this
recovery: **74,388 allocated CPU-seconds and 387 allocated GPU-seconds**, including
preparation, failed starts, folding, scoring and collectors. GPU allocation comprises
372 seconds folding and15 seconds input/schema support. These measure allocation,
not utilization, and exclude the original benchmark and unrelated LanM jobs.

For the 15 new sources, five native lanes took 1,075–1,660 seconds on 7 CPUs each,
with 470 actual MACE forwards including qualification fixtures. The two solvent
batches took 394/724 seconds on 12 CPUs, with 180 primary xTB and 180 independent tblite
calls total. Maximum recorded individual backend RSS was 113,272 KiB, not whole-job
memory. No new DFT ran. Detailed process times, reuse and earlier three-source
costs remain in the pinned receipts. No GPU inference speedup is inferred.

## Next decision

**Retain the recovered preparations and frozen-calibration results; do not
promote accommodation as uniformly better from this test.** The next question
is why accommodation converts two PedE abstentions into errors. A follow-up
needs a declared physical hypothesis; more optimizer freedom alone is not
justified by these outputs.

No recovery chemistry remains queued. Existing LanM jobs/watchers are preserved;
no new La/Dy run follows automatically. Full outputs remain under
`workspaces/ca_coverage_recovery_20261012/two_ca_scoring_final_v1/`.
From the repository root, the next read-only command is:

```bash
cat diagnostics/ca_coverage_recovery_20261012/DELIVERY.md
```

Exact executed manifests, wrappers and submit receipts remain in each workspace.
Do not resubmit completed chemistry to regenerate this report.
