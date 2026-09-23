# Fixed context membership: no net discrimination improvement

**Retain the released static-context scorer.** Holding the same chemically
complete surrounding fragments across each protein's saved folds corrects one
wrong prediction, but produces six new inconclusive calls on the same 207
scorable structures. No default, historical score, input or reference was changed.

## Full result and denominator

All 225 declared noncanonical sources are retained (25 proteins, 100 La- and
125 Ca-conditioned folds). These are consumed structural replicas of the
calibration proteins, not 225 independent biological observations or new blind
validation. The original 25 canonical sources calibrate the new protocol;
three consumed crystal controls are separate identity transfers.

| Method | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Released static context | 203 | 2 | 2 | 18 |
| Fixed membership | 199 | 1 | 8 | 17 |
| Preserved DFT | 199 | 3 | 6 | 17 |
| Native MACE core | 199 | 0 | 9 | 17 |

On exactly the **same 207 scorable sources**, released static context gives
203 correct / 2 wrong / 2 inconclusive, compared with **198 / 1 / 8** for fixed
membership. The new result's additional correct source is a technical SCF
recovery, described separately below. Neither abstentions nor the wider training
gap are counted as accuracy gains.

The only remaining wrong union prediction is A0ACD6B9F2, Ca-conditioned sample4.
The simpler native-MACE core already has the same 199 correct calls, zero wrong
calls and nine inconclusives with the same preparation coverage. Fixed membership
has not demonstrated an advantage over that inexpensive comparator either.

## Strict fold summaries

- La4: unchanged at **23 correct / 2 unavailable** per 25 proteins.
- All 100 La three-fold subsets: static context, DFT and native core each give
  **94 correct / 6 unavailable**. Fixed membership gives **92 correct /
  2 inconclusive / 6 unavailable**, with no wrong calls.
- Ca5: static 21 correct / 4 unavailable → union 22 / 3.
- Equal mean of complete La4 and Ca5 medians: static 20 correct / 5 unavailable
  → union 21 / 4.

The last two coverage gains arise solely from the unchanged-context A8R3S4 SCF
recovery. No previously complete group changes classification. Every required
member remains mandatory; a failed member is never replaced by a favorable fold.
These medians are structural descriptors, not thermal populations.

## What changed and what the experiment learned

The new protocol is `native_OMOL_GFN2_fixed_group_source_fragment_union_v1`. Fragment membership is the
label-independent union of the existing source-backed selections across the
supported saved folds. Complete graph closure, physical atom/H identities,
cap boundaries, formal charges and paired states must agree within each protein.
All 233 previously supported sources remain supported; all 17 prior failures stay
unavailable. Source coordinates, original core, water inventory, protonation,
MACE checkpoint and native GFN2 Hamiltonian remain fixed. This changes context
composition and cavity together; it is not an isolated geometry intervention.

All 25 canonical and all three crystal controls classify correctly. The new
canonical-only gap is **11.066622 model kcal/mol** versus the
static gap of 5.076366. The exact class-extrema reference was frozen before the
pilot and full-transfer calculations. It was never adjusted to their outcomes.
The larger calibration gap did not improve transfer utility.

| Source | Classification change | Raw R shift, model kcal/mol |
|---|---|---:|
| a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-0 | correct → inconclusive | -1.004855 |
| a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1 | wrong → correct | 8.740205 |
| a0a3f2yly8-pqq-la_model__conditioned_La__seed-1_sample-0 | correct → inconclusive | -1.309832 |
| a0a3f2yly8-pqq-la_model__conditioned_La__seed-1_sample-3 | correct → inconclusive | -1.037021 |
| c5axv8-pqq-la_model__conditioned_La__seed-1_sample-3 | correct → inconclusive | 0.000000 |
| mmol_1770-pqq-la_model__conditioned_Ca__seed-1_sample-1 | correct → inconclusive | -2.117111 |
| q88jh0-pqq-la_model__conditioned_Ca__seed-1_sample-0 | correct → inconclusive | -1.685363 |

For the corrected A0A3F2YLY8 Ca-sample1, the +8.740205 total shift contains
−5.731701 native-OMOL and +14.471905 solvent-transfer contributions. That cannot
be assigned uniquely to a hydrogen bond or mechanical restraint. Fold score
spreads improve for some proteins and worsen for others. C5AXV8 La-sample3 has
an exactly unchanged raw score; its new inconclusive call comes from the new
independently calibrated reference interval.

The score remains `R = E_Ca − E_La`, with
`E_M = E_OMOL,vac,M + E_GFN2,ALPBwater,M − E_GFN2,vac,M`.
Native eV and quantum hartree components are converted once. Absolute aquo and
binding-free-energy claims remain unavailable. This panel tests PQQ functional
class and structural robustness, not direct La/Ca affinity across protein families.

## Useful independent SCF recovery

A8R3S4 Ca-conditioned sample3 has the same actual selected fragment identities,
atom ordering, charge/multiplicity and both original native-MACE receipts. The
maximum coordinate difference is only **1.78e−15 Å**, within the already frozen
1e−12 Å arithmetic allowance; this is not byte-identical XYZ text. Its previously
failed **La vacuum** endpoint now converges under the already qualified
MaxIter500 policy. The repeated La ALPB energy exactly matches the old printed
energy; both Ca endpoint energies were reused.

The recovered raw contrast is **-405474.2087136921 model kcal/mol**,
which is **Ca-supported under the unchanged released static bands**. This record
can be reused independently of adopting union membership. Historical failed
collections and baseline tables remain untouched. See SCF_RECOVERY_v1.json for
source pins, exact states, original and recovered endpoint values.

## Execution, costs and practical limits

All **258 new native-MACE and 518 new native-GFN2 calls** succeeded. No new DFT,
folding, optimization, waters, protonation states or PLM scoring ran. Reuse was
explicit: the full225 stage alone reused 228 MACE and 454 GFN2 endpoint values,
including all successful pilot values. It performed 188 and 378 new calls.

Measured scheduler allocation: **83,513 core-seconds
and 72 requested GPU allocation-seconds**,
including the full collector and solver preflight. Pending cancellations had
zero elapsed time and no molecular work. Nested ORCA receipts are not added a
second time. Complete receipts are in COST_v1.json and SACCT_FINAL.txt.

The largest stage used a single H200 warm worker with 16 CPUs / 200000 MiB:
188 fresh native-MACE endpoints took 23.248 s worker wall including 1.281 s model
load, with peak CUDA allocation 5,676,625,408 bytes. The original 32-CPU pending
allocation was replaced using otherwise available CPUs; no other job was
interrupted. The 64-CPU / 128-GiB solvent allocation took 919 s wall including
preflight; the read-only collector took 185 s on one CPU. These are development
batch measurements with mixed reuse, not preparation-inclusive production costs
or matched-hardware speedups over DFT.

Initial source-graph preparation took 56.415 s wall with four
local workers. Local manifest construction and test CPU time were not fully
metered and remain additional one-time costs, not zero. No broad speed claim is
made. Memory and incomplete scheduler fields remain as actually reported.

The 12 real-fixture component/preparation/transfer checks and two final-record
checks pass. They cover graph
union, unchanged core and paired coordinates, original selection/crystal replay,
actual energy algebra, canonical-only reference, missing-member summaries,
corrupted-reference rejection, exact reuse of pilot receipts, full matched
denominators and independent recovery under released bands. Final-record checks
are recorded in FINAL_RECORD_TESTS_v1.txt. No fabricated scientific output was
used. The preserved first calibration status-merge error was repaired by parsing
the same original energies; it caused no molecular rerun.

## Disposition and runnable operation

**Do not promote fixed membership.** Keep the preparation adapter and records for
reproducibility; reuse the independently compatible A8R3S4 SCF recovery if useful.
No further union sweep is needed to establish this result. Other approved physics
tracks remain independent.

Read [RESULT_v1.json](RESULT_v1.json) for machine-readable counts and source pins,
[FROZEN_REFERENCE_v1.json](FROZEN_REFERENCE_v1.json) for the untouched new
canonical reference, and [COMMANDS.md](COMMANDS.md) for read-only inspection,
explicit preparation, validation, collection and comparison commands. All jobs
for this branch are complete; no production rollout or default change occurred.
