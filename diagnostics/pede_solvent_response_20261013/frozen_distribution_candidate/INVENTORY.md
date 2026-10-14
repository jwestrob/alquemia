# Frozen-distribution candidate: complete saved-input inventory

Prepared 14 October 2026 UTC. This is metadata preparation, not a new scientific
result or an execution manifest. No candidate energies, scores, classifications,
calibrations or numerical tests were computed.

The finite existing population is ready for the separate candidate
`E_OMOL,vac + G_nativeALPB(q_GFN2,vac)`. Here G contains the native electrostatic,
hydrogen-bond, surface and reference terms. The extra self-consistent response
and intrinsic remainder are not terms of this candidate. This is a changed
approximation, not a revised decomposition of the old energy.

## Fixed population and reusable inputs

| Existing record | Declared pools | Complete saved pools | Interpretation |
|---|---:|---:|---|
| Corrected expanded-envelope population | 135 | 135 | 25 designated canonical, 3 crystals, 6 probes, 101 further source/context pools |
| Original Ca-conditioned population | 125 | 107 | Preserve 16 preparation exclusions and 2 interrupted MxaF searches |
| Separately named MxaF technical recovery | 2 | 2 | Exact original source/state/policy; complete restarted pools |
| Separately named disulfide-H preparation repair | 1 | 1 | Existing source with corrected Cys hydrogen mapping |
| Q9 two-Ca alternative sources | 5 | 5 | Every declared sample retained |
| C5/Q88 two-Ca alternative sources | 10 | 10 | Every declared sample retained |

The existing annotated current125 overlay maps each original source to its
unchanged result or its already-declared recovery. Its complete125 pools include
all15 alternative two-Ca folds. The primary125 and its failures remain a separate
cohort. Alternative folds are distinct physical hypotheses; they do not overwrite
their one-Ca predecessors or increase the number of independent proteins.

The complete inventory identifies **1,542 distinct native component inputs**:
150 canonical and 1,392 conditional transfer inputs after exact artifact-key
deduplication. Each has the full original origin/Ca-proposal/La-proposal matrix,
both-metal MACE energy receipts, qualified standalone xTB6.7.1 vacuum/water point
records, raw eight-decimal charge files, output/receipt/XYZ hashes, source maps,
state identity and existing admission status. Every cited raw artifact hash was
checked. No MACE, SCC, DFT, geometry generation or optimization is needed.

Two native evaluations per input give at most300 canonical plus2,784 transfer
component replays. The four already-qualified PedE-origin inputs from replayv3
are identified separately in `NATIVE_REUSE.json`; exact reuse can remove eight
transfer evaluations. New generic-wrapper qualification and all actual costs
belong to the execution plan. These counts are not runtime estimates.

The inventory spans169–218 atoms, including twelve actual sizes. All states are
physical singlets with Ca charges−4/−3/−2 and paired La charges−3/−2/−1.
The previously qualified wrapper's192-atom assertion therefore needs a generic
atom-count interface before this finite replay. No element/state substitution is
authorized. The old backend absolute-energy gate remains failed; the existing
accepted solvent components used their separately defined component-v2 checks.

## Fixed scoring and staged gate

1. Use precisely the25 records marked `canonical_calibration` in the corrected
   expanded manifest, preserving their14Ca/11La known classes.
2. Recalculate all six candidate cells using the new energy expression. For each
   metal use the ordered common pool `origin, adaptive_Ca, adaptive_La`; retain
   origin unless the mathematical minimum lowers energy by strictly more than
   0.1kcal/mol. Report the mathematical minimum separately.
3. Apply the existing calibration **rule**, separately for operational,
   mathematical and origin policies: all25 available; `Ca_max=max(Ca)` and
   `La_min=min(La)`; positive separation must be strictly greater than0.02kcal/mol.
   The old numerical bands are baseline evidence only.
4. Only an available operational canonical calibration releases dependent
   transfer replay. Freeze its bands before transfer outcomes. A failed canonical
   gate is a candidate result, not permission to change its members or thresholds.
5. Retain the full135 population, primary125, annotated current125 and all15
   two-Ca sources. Reuse the recorded100-triple memberships and six preparation
   exclusions, with each member's exact source/context identity. No new Ca-fold
   voting or triple rule is introduced. A valid origin remains an independent
   diagnostic when a required accommodated cell is unavailable.

All populations are consumed development references. Existing source-conditioning,
water/state and cap policies persist; no new biological evidence or unknown PLM
label is inferred. The historical corrected strict208/tenfold calibration failed
and uses a different context protocol; it is not pooled into this4.3Å candidate.

## Exact artifacts and implementation handoff

`READY.json` pins the full workspace `INVENTORY.json`, its metadata-only generator,
baseline manifests/results, calibration/selection implementation and optional
existing replay reuse. The workspace inventory retains all278 source-record
occurrences, including unavailable historical rows, with exact cohort and overlay
joins. Duplicate sources under different contexts remain distinct records.

`native_inputs[]` is the replay interface: state/XYZ, raw vacuum/water artifact
pins, and `uses[]` links to each source/candidate/metal MACE receipt.
`records[].matrix` supplies all six required cell links, including missing status.
The four archived water components are to be parsed from each pinned stdout
during Slurm preflight; they were not newly extracted on the login node.

Root owns the candidate's final plan, scoring, execution and promotion decision.
The native-replay author owns the narrowly generalized component wrapper. All
scientific execution and validation must run through Slurm. Production is unchanged.
