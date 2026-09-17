# Fixed-field DFT/MACE screen: useful reduction, failed gate

**Adding the learned local scaffold term reduces the partition jump but does
not meet the frozen 2 kcal/mol criterion.** This saved-output screen used no new
DFT, MACE, charge-fit or solvent calculations. Baseline remains unchanged.

| Contribution to qm36 minus qm33 contrast | Medium | Large |
|---|---:|---:|
| Vacuum DFT core | +61.738603357 | +61.738603357 |
| CHELPG coupling to permanent protein charges | -51.498665785 | -51.498665785 |
| Full-minus-core MACE short correction | -3.637970036 | -4.680450290 |
| Result, kcal/mol | **+6.601967536** | **+5.559487282** |

Without MACE, this same CHELPG vacuum expression has a 10.239937573 kcal/mol
jump. Replacing fitted coupling with the saved exact density diagnostic yields
7.009438579 / 5.966958325 kcal/mol after adding MACE. That comparison does not
select a charge scheme per case; medium/CHELPG remains the declared primary.
The residual cannot be attributed uniquely to the fitted monopoles, caps,
classical residue approximation or missing field response from this arithmetic.

The expression is native vacuum r2SCAN-3c core + direct CHELPG/environment
coupling + MACE interaction_energy(full) minus interaction_energy(core).
It includes neither MACE's learned global electrostatic nor its field-dependent
electron-energy scalar. The short/field split was learned jointly and is not a
unique physical decomposition; transferring it to external charges is itself an
approximation. This failure rejects this frozen screen for expansion. It does
not establish that every hybrid model or environmental contribution is useless.

No solvent or affinity score was calculated. Those fields are null, not zero.
No existing CPCM, aquo offset, threshold, entropy or relaxation correction was
attached. All four consumed states retain original geometry, caps, charge,
protonation, environment, assembly and empty water inventory. Both full inputs
are the same across partitions; the full short contrast cancels exactly from
the partition difference. All 38 recorded source/accounting/rigid checks pass.

The old 1H4I boundary trial moves the complete capped Asp303 sidechain, not an
entire peptide residue. This screen retains that inherited limitation and
cannot claim a stronger complete-residue partition test.

## Verification and measured cost

Four real-artifact tests pass: archived component algebra, real charge coupling
and joint rigid transformation, corrupted actual atom mapping rejection, and
roundoff-only replay versus altered actual output rejection. These are software
and saved-output checks; no unavailable solver was replaced by invented output.

Successful local analysis: 10.084510429 wall seconds, 9.741181755 CPU seconds.
No new cluster allocation. Three initial preflight invocations rejected exact
floating comparisons before producing the candidate result. They exposed only
summation roundoff (about 2e-12 kcal/mol in replay and 5e-14 e in formal charge
closure); the adapter now allows 1e-9 absolute floating replay error with exact
metadata, while retaining the frozen scientific tolerances. Failed preflight
and local testing costs were not fully profiled and are unavailable, not zero.

Cached costs remain real: medium pilot1200470 used329GPU-s/5264core-s;
large1200525 used586GPU-s/9376core-s. Original vacuum DFT/MBIS1199949 used59648allocatedcore-s; native CHELPG
utility1199984 used296allocatedcore-s. These are prior development allocations,
not costs incurred again by this screen or a matched production speed estimate.
Read the original source reports for all initial utility attempts and timing.

**Recommendation:** retain baseline; do not expand this frozen vacuum mixture
into an accuracy trial. MACE's local term has a measurable, consistent effect,
but the residual representation dependence still requires a coherent fix.
The broader discriminator goal remains active. No fitting, parameter rescue,
production promotion or new biological comparison was performed.

## Replay from the repository root

```bash
python scripts/mace_fixed_field.py --medium workspaces/mace_analytic_20260916/pilot_v1/collection_job_1200470.json --large workspaces/mace_large_20260916/pilot_v2/collection_job_1200525.json --charge-fit workspaces/density_embedding_20260916/chelpg_v1/comparison_v1.json --agreement diagnostics/mace_fixed_field_20260916/PLAN.md --output workspaces/mace_fixed_field_20260916/screen_replay_v1
python -m unittest discover -s tests -p test_mace_fixed_field.py -v
```

The output directory must be fresh; original scientific outputs and screen_v1
are immutable. The report copies its executed implementation into that directory.
[result.json](result.json) retains endpoint components, source pins, all38checks,
null unsupported fields and measured analysis time.
