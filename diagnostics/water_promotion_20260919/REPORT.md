# Contextual water preparation promoted

**Jacob authorized the versioned water-preparation default.** The new
`affordable_workflow.py baseline` operations use `contextual_if_supported` and
retain original/prepared results separately. Historical inputs, references and
original-scoring access remain intact. The electronic scorer stays native
r2SCAN-3c/CPCM(Water)/DefGrid3, ORCA6.1.1.

## What is promoted and why

MACE proposes existing water H coordinates using the established source-defined
3.5Å polar context, reference water internal geometry, and source/radial-away
orientation starts. This includes normalization of water O–H lengths/angle before
rigid orientation; it is not solely rotation of the source geometry. Only inner
water H coordinates transfer back to the original DFT core. All water oxygens,
nonwater atoms, caps, charges, multiplicities and water inventory stay fixed.

This earned limited promotion because the demonstrated alpha-lactalbumin/GGR
ordering changed0/6→6/6 under native DFT and2/6→6/6 under whole-chain masked MACE,
without altering dry PQQ inputs. These are two alpha structures and three GGR
structures, **one consumed biological comparison**, not six independent labels.
[Scientific evidence](../hydration_network_20260918/REPORT.md) and
[MACE-scoring evidence](../hydration_scanner_20260919/REPORT.md) remain preserved.

Independent [parvalbumin preparation](../water_reference_validation_20260919/REPORT.md)
now also completes: CD remains exact dry identity; preparing actual EF waterA166
changes R by−3.3215921675kcal/mol. Both DFT endpoint energies decrease. Its old
water geometry had stretched O–H bonds, so energy gains combine geometry repair
and orientation. The ledger marks both site affinity directions unresolved;
this is preparation transfer, not a new accuracy win. It cost4864allocatedCPU-s
and24GPU-s, two MACE proposals plus two DFT single points.

## Scope and compatibility

| Operation | Policy |
|---|---|
| Supported wet source-backed amide-v3 site | Select water-prepared score; retain original alongside |
| Dry canonical PQQ | Exact original coordinates, recipe, protocol and released decision bands |
| Explicit `--water-policy original` | Original preparation only |
| Unsupported wet cofactor/legacy core | Explicit error; no silent substitution |
| Missing or failed prepared endpoint | Prepared score unavailable even when original is complete |
| Original successful chemistry | Reuse only compatible actual input/recipe/charge/receipt |

New workflow ID: `baseline_contextual_water_v1`.
Prepared DFT protocol: `amide_v3_native_r2scan3c_context_prepared_water_H_v1`.
Canonical unchanged protocol: `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`.

The promoted interface takes explicit whole-protein prepared-site manifests with
source mappings. The older raw inbox uses legacy core representations and has
not been silently migrated; its watchers/jobs were not restarted. Wet PQQ/cofactor
preparation remains unsupported, and no retained-water PQQ benchmark is available.
No new water, occupancy probability, proton-transfer state, entropy, environmental
score or altered electronic method is part of this release.

Scores use R=(E_Ca−E_La) with hartree converted once. Changed wet scores have no
inherited absolute bands or compatible aquo gauge; S/decision remain unavailable
or explicitly uncalibrated. Unchanged canonical PQQ retains its released
S=R−A and bands. Do not interpret the large raw elemental offset as an affinity.

## What actually ran

Promotion replay used the actual two alpha preparations, four archived successful
MACE proposals, eight original/prepared native DFT endpoints, and four archived
1H4I/4MAE endpoints. **Zero new MACE or DFT evaluations were needed for integration.**
The independent parvalbumin calculation above is separate and fully receipted.

Final replay products are under `workspaces/water_promotion_20260919/`:

- `prepared_v2/plan.json`, `dft_v2/manifest.json`, `result_v2.json`, `REPORT_v2.md`.
- `dry_prepared_v1/plan.json`, `dry_dft_v3/manifest.json`, `dry_result_v3.json`.
- Release hashes and source pins: `params/baseline_water_v1.json`, `VERIFICATION.json`.

Earlier replay versions remain as development records. Two dry-preparation
preflights rejected harmless archived comments/resource directives; parser support
was corrected without changing the Hamiltonian or launching a calculation.
The first two test errors were fixture identity/path mistakes, corrected without
changing expected energies or weakening the runner's directory boundary.

**50 real-fixture tests pass, zero skips:**10integration tests (28.311s) plus40
existing/component tests (57.215s). Checks cover real energy extraction, signs,
units, released PQQ decisions, unchanged dry coordinates, source water mapping,
explicit original access, unsupported wet chemistry, charge/input mismatch,
missing prepared endpoints, cache reuse and fresh retry workspaces. These are
parser/preparation/integration tests over real artifacts; they are distinct from
actually executed molecular calculations. Local preparation/report CPU time was
not fully metered and is not claimed to be zero.

## Operations and cost

[Commands](COMMANDS.md) provide request, prepare, water-execute, prepare-dft,
dry-run, execute, collect, report and prepare-retry without source edits.
Executors reuse the existing MACE and ORCA runners. Partial attempts are preserved;
retry uses a fresh workspace and reuses compatible successes.

A wet case with archived original scores needs two new prepared DFT endpoints
plus two MACE water proposals. Without original scores, retaining both arms means
four DFT endpoints. All tasks are visible in the manifest. This is an explicit
cost tradeoff, not a claim of speedup. Original alpha practical preparation used
96GPU-s plus118s×64CPUs for four new prepared endpoints across two structures;
parvalbumin's separate measured transfer is above. No production-wide rescore ran.

**Recommendation:** use the promoted water-H preparation on supported wet sites;
retain the unchanged canonical PQQ reference path. Expand accuracy evidence through
the separate second-shell/electronic experiments. Successful preparation and lower
energy alone do not establish broad La/Ca discrimination.
