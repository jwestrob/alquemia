# Nikasha: adaptive accommodation earns further development

**Adaptive accommodation improves individual-fold discrimination on matched
coverage: 200→202 correct out of 204.** It repairs one wrong prediction and two
inconclusives, at the cost of one new inconclusive. The strongest example is
A0ACD6B9F2, where fold disagreement also decreases substantially. All 111 available
Ca-class structures remain correct. This is a modest, concrete gain on consumed
structural-development data, not broad independent biological validation.

Recommendation: **pursue adaptive accommodation further; retain the released
fast MACE/GFN2 scorer as the production default.** The adaptive version loses
coverage and requires more work. Consistent pocket membership alone does not
earn promotion; the tested CPCM route is numerically unqualified. All three
approved branches are complete, with no new DFT, folds or PLM cohort rescore.

## Full comparison

The 225 noncanonical structures come from the same 25 reference proteins used in
development, including both Ca- and La-conditioned folds. They are structural
replicas, not 225 independent proteins. Each changed method uses its own reference,
calibrated on exactly the original 25 designated canonical structures and frozen
before its fold-transfer results. No fold outcome changed a threshold.

| Method | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Released MACE/GFN2-ALPB | 203 | 2 | 2 | 18 |
| Preserved DFT | 199 | 3 | 6 | 17 |
| Consistent pocket membership | 199 | 1 | 8 | 17 |
| Adaptive common geometry pool | 202 | 1 | 1 | 21 |

On the **same 204 scored structures**, released MACE/GFN2 gives200 correct,
2 wrong and2 inconclusive; adaptive gives202 correct,1 wrong and1 inconclusive.
Of the 200 previously correct calls, 199 remain correct and one becomes
inconclusive. The other three adaptive losses relative to the full released
coverage are missing calculations, not incorrect class predictions.

All 25 canonical members and all 3 already-consumed crystal controls classify
correctly under the new adaptive reference. Its calibration gap is 2.5466901856
model kcal/mol; this calibration success is separate from the fold-transfer gain.
Adaptive old-band transfer gives 199 correct / 0 wrong / 5 inconclusive / 21 unavailable.
The frozen adaptive bands give 202/1/1/21; both results remain in the ledger.
Mathematical row minima and the operational 0.1 kcal origin-retention rule agree
on all 204 available transfer sources.

## What accommodation actually improves

| Source | Released → adaptive | Change in raw R, model kcal/mol |
|---|---|---:|
| A0ACD6B9F2, Ca-conditioned sample4 | wrong → correct | +17.828009 |
| A0ACD6B9F2, La-conditioned sample4 | inconclusive → correct | +21.662208 |
| A0A3F2YLY8, Ca-conditioned sample3 | inconclusive → correct | +6.889575 |
| C5AXV8, La-conditioned sample3 | correct → inconclusive | −0.165259 |

All three repaired calls also become correct under the old bands: their benefit
is not created solely by recalibration. The C5AXV8 raw contrast changes little;
its new abstention arises primarily from the changed, independently frozen
reference interval. The A0A3F2YLY8 Ca-conditioned sample1 error remains wrong
under the adaptive reference.

For the repaired A0ACD6B9F2 Ca-conditioned structure, selected composite energies
fall by 11.013156 for Ca and 28.841165 for La. Their difference shifts R by +17.828009,
containing +10.174326 from native MACE and +7.653683 from solvent transfer. Both
metals select their own adaptive geometry from the common pool. This supports a
local accommodation explanation within this model; it does not establish a
binding free energy or supply new DFT validation of these particular motions.

Across that protein's complete source arms, the La-conditioned score range
shrinks 22.7467→5.9950 and the Ca-conditioned range 31.8054→13.0035 model kcal/mol.
This is a real structure-consistency benefit on the named problem case. Across
all proteins, spread does not generally improve: 11/22 comparable La4 groups
narrow and 11 widen; 6/19 Ca5 groups narrow and 13 widen. Median range changes are
−0.0880 and+0.9325 respectively. These spreads describe structures, not uncertainty
estimates or equilibrium populations.

The existing strict aggregate decisions were already correct wherever available.
Adaptive accommodation adds no correct aggregate decision on common coverage.
Strict La4 coverage becomes 22/25 versus 23/25 released; Ca5 becomes 19/25 versus
21/25; balanced summaries 18/25 versus 20/25. Correlated La triples become 91 correct /
9 unavailable versus 94/6. All available aggregates remain correct; the loss is
coverage. The four-member and three-member summaries require every declared
member, with no favorable subset selection.

## What was implemented and run

Existing source graphs, kinematics, MACE proposals, native GFN2 runner, archived
base pools and exact compatible endpoint energies were reused. The generic
SLSQP completion policy expresses the same native MACE objective in hartree
units, avoiding the old unsafe-trial failure without changing the four selected
physical angular coordinates, ±0.8 rad bounds or final 0.8 Å displacement domain.
Its source protocol is `common_four_angular_native_OMOL_SLSQP_hartree_units_v2`.
The parameter scaling and SciPy tolerance implications are documented in
[the completion report](../adaptive_completion_20260922/REPORT.md).

Both metals receive the same admitted geometry pool:

`{origin, old_Ca_proposal, old_La_proposal, adaptive_Ca, adaptive_La}`.
Each matrix cell uses the appropriate metal/charge and unchanged physical atoms,
water/proton/cofactor state and cap mapping. Selection uses
`E_M = E_OMOL,vac,M + E_GFN2,ALPB,M − E_GFN2,vac,M`;
`R = min(E_Ca) − min(E_La)`. These are protocol-specific raw electronic
contrasts; no compatible aquo reference or affinity interpretation was invented.

The original 30 completed all 60 searches, 60 fresh cross-MACE and 240 native GFN2
calls. The 225 transfer ledger retained20 inherited unavailable sources, yielding
410 searches: 409 produced accepted candidates and one hit the 200-iteration
limit. The 204 complete pairs needed 408 new cross-MACE and 1632 native GFN2 calls;
**every one of these new singlepoints succeeded**. No rescue or fallback was used.

The failed Q60AR6 La-conditioned sample0 La search remains unavailable. Additional
losses relative to released coverage are inherited old-proposal failures for
P12293 Ca-sample4 and Q60AR6 Ca-sample1. The baseline's own unavailable A8R3S4
source remains separate. The full 17 preparation failures also remain visible.

172/409 accepted transfer candidates reach a search boundary. They are bounded
proposals, not proven stationary minima. 393/408 selected metal endpoints use a
new adaptive geometry; 109 select the other metal's adaptive proposal. Median
composite work is−5.536 for Ca and−9.776 for La. More energetic relief is not by
itself predictive evidence; the classification and coverage comparison above
determines the practical benefit.

## The other branches

**Consistent membership:** the same source-backed complete context fragments
were held across each protein's supported folds. On 207 matched scored sources,
203 correct / 2 wrong / 2 inconclusive becomes 198/1/8. It fixes the A0A3F2YLY8 error,
but six previously correct calls become inconclusive. A wider canonical gap
(11.066622) does not compensate for that transfer result. An independent,
unchanged-context SCF recovery for A8R3S4 is correctly Ca-like under the released
bands and is saved separately. [Full report](../consistent_context_20260922/REPORT.md).

**CPCM:** all 16 declared calculations were executed, including matched ordinary
vacuum controls. ORCA's CPCM path forces ordinary SCF; only 4/8 CPCM and 1/8 vacuum
endpoints converge, giving 0/4 complete metal contrasts. All 11 failures explicitly
exhaust SCF iterations. The cavity/reaction field exists, but the numerical solver
is unqualified and this attempt yields no discrimination result. The wider
calibration and 225-source expansion were not run. Close this implementation
branch. [Full report](../solvent_cpcm_20260922/REPORT.md).

## Cost and qualification

| Branch | Allocated core-seconds | Requested GPU allocation-seconds |
|---|---:|---:|
| Consistent context, including recovery | 83,513 | 72 |
| CPCM and matched vacuum, including failures | 290,048 | 0 |
| Adaptive original30 and calibration pool | 60,864 | 548 |
| Adaptive225 proposal searches | 110,688 | 3,459 |
| Adaptive225 cross-scoring and in-job collection | 367,680 | 286 |
| **Total** | **912,793** | **4,365** |

This is 253.554 allocated CPU-hours and 1.2125 requested GPU-hours, separate resource
measures. It includes failed attempts, reserved idle capacity and in-job collection;
unmetered local setup/tests/reporting are additional, not zero. Historical reused
calculations are excluded. Total fresh work is 13,742 MACE evaluations and 2,406 GFN2
attempts, with zero new DFT. The finite declared jobs are all terminal.

Full 225 cross-scoring took 38 min 43 s from first allocation start to last end,
including staggered scheduling and collection. Its warm cross-MACE worker used
68.393 s for 408 calls, loaded the model in 1.321 s, and peaked at 5,645,979,136 bytes
of CUDA allocation. Proposal-search median was 4.751 s/endpoint; one failed search
took 579.392 s. None of these figures is a measured cold full-scanner latency.

A cold five-geometry score needs up to 20 native GFN2 singlepoints per source,
versus four for released static scoring, plus proposal work. The current gain
does not justify blanket production use at that overhead. Numerical credibility
is supported for the admitted finite pool; scientific benefit is modest and
localized; routine affordability of the full candidate workflow remains unproven.

## Validation, delivery and next step

Real-artifact checks cover matrix algebra/signs, paired coordinates and charges,
exact reuse, source/state mismatches, missing candidates, canonical-only bands
and full denominator preservation. Existing shared-pool regression: 27 passed;
new population-join checks: 3 passed; new comparison-join checks: 2 passed.
The context branch reports 14 checks; CPCM reports 7 passed and one already-covered
prelaunch-state skip. Adaptive completion and sharding have their separate passing
real-fixture checks. These parser/algebra tests are distinct from the actual
molecular execution above. No fabricated scientific fixture was substituted.

The final three reporting commands all exited 0; direct component checks reproduce
`delta_R = work_Ca − work_La` and its native/solvent sum for every available pair.
The final collector snapshot was reused byte-for-byte, without duplicate parsing
or molecular work. Detailed ledgers, costs, report-code pins and commands live in
`workspaces/nikasha_next_phase_20260922/`; all 225 individual outcomes, 75 strict
groups and 100 correlated triples are retained in `COMPARISON.json` and `.md`.
[`COMMANDS.md`](COMMANDS.md) gives executable reproduction operations.

Current protocol IDs remain distinct:

- Production: `fast_PQQ_OMOL_GFN2_ALPB_v1`, unchanged.
- Context challenger: `native_OMOL_GFN2_fixed_group_source_fragment_union_v1`.
- Adaptive challenger: `nikasha_scaled_angular_common_geometry_native_OMOL_GFN2_ALPB_v1`.
- Unqualified solvent route: `native_OMOL_plus_native_GFN2_CPCM_water_transfer_v1`.

**Next proposed scientific test:** combine consistent source-backed pocket
membership with adaptive accommodation under one newly calibrated protocol.
The two isolated branches corrected different errors, providing a concrete reason
to test their interaction. Their gains cannot be assumed additive, and that
combined analysis has not been executed or promoted. Its scope should be agreed
before launch. The broader accuracy goal remains open; this phase provides a
specific positive development result and closes two unsuccessful standalone routes.
