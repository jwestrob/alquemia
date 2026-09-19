# Accuracy tracks — 2026-09-19 latest checkpoint

Jacob authorized contextual water preparation promotion and discretionary next
experiments. Accuracy and PQQ fidelity remain the objective; neither catalytic
competence nor imitating unchanged DFT is a substitute.

## Promoted water preparation

The versioned `affordable_workflow.py baseline` interface defaults to supported
contextual water-H preparation. Both original and prepared scores are retained;
dry canonical PQQ is exact identity with released bands. Unsupported wet
cofactors/legacy source maps fail explicitly. Raw legacy inbox watchers remain
untouched. [Release report](../water_promotion_20260919/REPORT.md),
[commands](../water_promotion_20260919/COMMANDS.md).

Demonstrated alpha/GGR discrimination: native DFT0/6→6/6, whole-chain masked
MACE2/6→6/6, one biological comparison across structural replicas. This is useful
preparation improvement, not broad affinity validation. Reference water internal
geometry normalization is part of preparation, alongside contextual orientation.
No new occupancy, protonation, entropy or absolute wet-site bands are promoted.

Independent parvalbumin transfer completed (jobs1202475/1202476): dryCD unchanged;
wetEF R shifts−3.3215921675kcal/mol. Original O–H bonds were stretched, so absolute
energy drops combine normalization and orientation. Its two site directions remain
unresolved in the evidence ledger. Cost4864allocatedCPU-s/24GPU-s. Six real-fixture
tests pass. [Report](../water_reference_validation_20260919/REPORT.md).
No qualified existing wet-PQQ preparation was found in the bounded source inventory;
do not invent waters or advertise wet-PQQ support.

## Environmental accuracy track

| Experiment | Result | Actual new allocation |
|---|---|---|
| Context-supported coordinate proposals | No improvement; PQQ gap30.0719→26.2757, alpha/GGR margins weaken, alpha replica spread grows3.649→10.155 | 30560CPU-s/54GPU-s including one import failure |
| Static complete polar context on GGR replicas | All6 native DFT margins improve; weakest0.654840→4.564286kcal/mol. Native OMOL2/6→6/6 | 35840CPU-s/20GPU-s |
| Complete compact-context OMOL PQQ panel | 25/25canonical and3/3consumed crystal transfers; gap79.030608→2.315459model-kcal | 1120CPU-s/70GPU-s |
| Native CPCM DFT on all4charge−1 PQQ expansions | Job1202478 running;8new endpoints, old core and neutral controls reused | Final cost pending |

These costs exclude unmetered local preparation/testing. They describe development,
not a matched end-to-end production speed benchmark. The static scorer remains
research-only. Expanded PQQ no longer has perfect label/total-charge separation:
three Ca-family contexts share charge with ten La-family contexts and still order
correctly. The much narrower gap is a robustness concern; this harder equal-charge
comparison is not itself a new biological validation set. Actual charged residues,
solvation differences and model charge-conditioning can all contribute. The native
CPCM check tests the effect without tuning charges or thresholds.

Completed reports: [geometry](../coordination_preparation_20260919/REPORT.md),
[replicas](../environment_replicas_20260919/REPORT.md),
[PQQ panel](../environment_pqq_20260919/REPORT.md).
Running scope: `diagnostics/environment_pqq_dft_20260919/`.

## Chemical-state accuracy track

All24native analytic-gradient points on three source-water/carboxylate proton
paths completed normally. Every sampled point is uphill; transferred-state costs
range20.8–90.0kcal/mol. Tiny near-origin response is possible, but no competitive
transferred basin, pH population or useful correction was established.
Both1F6Sand6IP9are bovine alpha-lactalbumin, one biological group.
Cost154752allocatedCPU-s, zeroGPU-s. Batch1202444 failed only in postprocessing;
all actual quantum results were recovered without scientific retries. Eight tests
pass. [Report](../chemical_states_20260919/REPORT.md).
Do not add an unfavorable state's attractive relative shift to the classifier.

## Independent electronic diagnostic — unfinished

Job1202429 continues the declared six DLPNO-CCSD(T1)/TightPNO/CPCM-PTES endpoints
on exact real GGR1GLG and original/water-prepared alpha1F6S cores. No final pair
was available at the latest checkpoint; four initial HF references converged.
Four16-rank tasks share64CPUs/256GiB; diffuse-basis solvent/integral transformations
are expensive (hundreds of seconds per batch,25–30batches on alpha).
Transient scratch exceeded513GiB; active files must not be removed.
This is a one-time reference diagnostic, not an affordable routine scanner method.

Actual Hamiltonian, La46ECP, Ca3s/3p and La5s/5p correlation, electron accounting,
PTES and TightPNO were checked. No new PQQ/high-level expansion was submitted.
Five source/parser tests pass; completed new CC scientific regression is unavailable.
Recover [checkpoint](../electronic_accuracy_20260919/CURRENT.md) and live receipts.
Automatic final collection is installed: detached monitor PID1633721 waits for
terminal job/accounting and free execution lock, then writes actual results/cost,
report/vault and a final email. It never launches chemistry. Root reran6real
electronic-source/incomplete-result tests successfully (0.347s). Charged-PQQ has
its existing successful-job collector plus live completion observer PID1638834
(session22209) to retain failures and allocation cost. Root reran its2source tests
(1.044s); actual new scientific outcomes remain pending. Installation receipts
and recovery commands are in each diagnostic checkpoint. Do not duplicate them.

## Preservation and validation

The previously completed water-basin entropy approximation remains unavailable;
it failed its independent DFT energy checks. Historical identity erratum is
preserved. Existing Hans-LanM EF1/EF2/EF3 scores remain ordered site vectors within
one protein-level measurement, not three independent labels.

The root reran integration and completed-component tests on actual artifacts.
See the promotion VERIFICATION.json for current totals and exact logs. No fake
scientific results, source resets, shared-job interruption, threshold refit,
production rescore, push or environmental-scorer promotion occurred.
