# PedE failures: solvent response is the leading discrepancy

13 October 2026. Saved-artifact investigation requested by Jacob before more
molecular experiments. Production, calibration and original results are unchanged.

## Scientific judgment

The remaining two PedE errors are **not explained by a demonstrated bad fold or
preparation defect**. The gross missing-site-metal problem was repaired by the
separately recorded two-Ca source hypotheses. Their finer geometry and confidence
do not distinguish the two wrong calls from the three correct ones. Instead,
the largest starting-score difference is almost entirely in the **GFN2
ALPB-minus-vacuum solvent contribution**, with very little difference in the
MACE Ca-minus-La contribution. Investigate that response before commissioning
more folds or a larger accommodation search.

This identifies where the model discrepancy enters. It does not establish that
the solvent term is physically wrong, that AF3 coordinates are experimentally
correct, or that removing solvent would improve a valid classifier. Small
coordinate changes, electronic charge redistribution and the representation of
a capped pocket in bulk-water continuum remain distinguishable hypotheses.

## 1. No measured geometry defect specific to the wrong sources

All five recovered sources contain Ca at PQQ and retain the remote second Ca.
Their same prepared local state has 192 atoms, 14 caps, Ca/La charges −3/−2,
and physical singlets. Both metals access the same saved candidate pool.

| Two-Ca PedE sample | Accommodated outcome | Typed coordination number | Ca–Asp358 OD1, Å | PQQ RMSD to canonical source, Å | Metal pLDDT |
|---|---|---:|---:|---:|---:|
| 0 | Correct Ca | 7 | 3.051 | 0.162 | 91.24 |
| 1 | Wrong La | 7 | 3.053 | 0.147 | 91.30 |
| 2 | Wrong La | 6 | 3.107 | 0.210 | 91.34 |
| 3 | Correct Ca | 6 | 3.179 | 0.259 | 91.46 |
| 4 | Correct Ca | 6 | 3.108 | 0.189 | 91.74 |

RMSD uses the same five-residue backbone fit and chemically equivalent PQQ
carboxyl-oxygen swaps. The canonical source is another prediction, not an
experimental truth structure. All five have high PQQ confidence, similar ring
planarity and the same PQQ/Asn/Glu donor pattern. Correct sample3 is more displaced
from the reference than either wrong sample. The 3.1 Å contact definition makes
the Asp contact marginal in several models of both outcomes.

The admitted La proposals also resemble one another: wrong2 and correct3 rotate
Glu221 chi3 by approximately +20.9°/+20.7° and reach similar oxygen distances.
PQQ and metal positions are fixed during these searches; mapped bonds are
preserved. A boundary stop or large La relaxation is not exclusive to an error.
See [geometry report](geometry/REPORT.md), which retains all ten raw sources and
nineteen saved local geometries. Old La sample1 remains raw-only with its metal
outside the PQQ site; it is not assigned a replacement score.

High confidence is not proof of energetic accuracy. These observations reject
simple explanations such as “only the wrong models have the metal misplaced,
poor confidence, distorted PQQ, or one lost donor.” They do not exclude errors
in sub-ångström coordinates or source chemical state.

There is also a shared preparation/QA limitation: PQQ C5–O1 is approximately
1.34 Å in both these sources and the canonical model. The executed preparation
declares oxidized quinone chemistry but applies broad connectivity checks,
not bond-order-specific geometry validation. No completed earlier audit found
here validates that exact fine geometry. The [state-hypothesis record](../pqq_state_hypothesis_20260922/PROPOSED.md)
keeps reduced-PQQ testing proposed-only. A bond length alone is not a redox
assignment or an explanation specific to the two errors; it is a reason not
to equate source admission with fully validated cofactor chemistry.

## 2. The solvent contribution already separates the origins

The actual corrected composite is

`E_M = E_OMOL,vac,M + E_GFN2,ALPB-water,M − E_GFN2,vac,M`;

`R = E_Ca − E_La`, with larger R more La-like on this calibrated scale.
Current solvent components use standalone xTB6.7.1, checked against tblite0.7.0.
This is not an unchanged run of the defective native-ORCA La parameters.

Origin differences relative to the **same canonical PedE source**, kcal/mol:

| Sample | MACE contribution to ΔR | Solvent contribution to ΔR | Total ΔR |
|---|---:|---:|---:|
| 0 | −2.697 | −6.673 | −9.370 |
| 1 | −2.514 | +11.642 | +9.128 |
| 2 | −0.722 | +5.391 | +4.669 |
| 3 | −2.363 | −9.615 | −11.978 |
| 4 | +1.828 | −3.549 | −1.720 |

Using the sources defining the total-score extrema, correct3→wrong1 changes R
by **+21.106215 kcal/mol = −0.150413 MACE +21.256628 solvent**. This is a
decomposition of the already-observed spread, not a claim that total MACE
energies or the structures are identical. Four compatible La-conditioned PedE
controls are also retained; no incompatible historical score fills missing La1.

All54 saved metal/candidate component checks passed the existing gates. Maximum
component residual is 1.29004e−7 Eh and population residual 8.33357e−8 e.
Agreement between implementations does not independently establish the correct
electronic solution or a physically accurate solvent response. The previous
absolute-v1 parity failure remains a separate failed test. See
[energy decomposition](energetics/REPORT.md) and its pinned tables.

### Saved charges narrow the question further

All36 saved origin calculations, each with both implementation records, have
consistent atom maps, charge sums and converged SCF evidence. Correct3→wrong1
changes the solvent contribution by −56.7293 kcal/mol for Ca and −77.9859 for
La. The much greater La stabilization produces the +21.2566 difference.

The accompanying charge redistribution is mostly on PQQ and surrounding
fragments, not a large jump in metal charge. Selected xTB population sums in
ALPB water (electrons of net charge; these are not formal oxidation states):

| Source / endpoint | Metal | PQQ, including its H atoms | Arg386 fragment | Pro309 fragment |
|---|---:|---:|---:|---:|
| Correct3 / Ca | +0.613 | −2.923 | +0.882 | +0.106 |
| Wrong1 / Ca | +0.603 | −3.752 | +1.406 | +0.888 |
| Correct3 / La | +0.751 | −2.781 | +0.922 | +0.145 |
| Wrong1 / La | +0.724 | −3.639 | +1.532 | +1.089 |

Wrong2 shows the same direction with a less extreme PQQ population. Vacuum
PQQ populations differ much less between1/3. **Excessive solvent-enabled charge
redistribution across the capped pocket is therefore a concrete candidate
failure mechanism.** These observations do not prove spurious charge transfer,
a redox-state change, or a discrete SCF-branch switch. Small electronic gaps
and fractional occupations also occur in correct controls.
The correct canonical La-conditioned PedE source also has strongly negative
PQQ populations (−3.447 for Ca, −3.316 for La); this is not a universal
wrong-case signature or a basis for a new charge-based rejection rule.

One command-line warning is retained: xTB rejects the literal reference token
`gsolv`, but the actual calculation setup still prints the effective reference
`gsolv [1 M gas/solution]`. A shared coordinate-independent reference shift
cancels in the within-metal source difference; this warning is not evidence
for the observed structural sensitivity. Correcting future syntax should be
separate from changing or reinterpreting these completed results.

## 3. Accommodation and calibration also matter

The wrong sources were already inconclusive at their origins. Their raw score
changes with accommodation are +9.419/+10.203 kcal/mol, but the canonical Ca
band edge also shifts +8.189. Thus their Ca-band margins deteriorate only
1.230/2.014 kcal/mol; correct sample3 deteriorates more, by 2.394.

The canonical abstention interval contracts from **10.298 to 3.472 kcal/mol**.
One wrong source ends only **0.277 kcal/mol** inside the La band. These facts
explain why ordinary-sized response changes can convert pre-existing ambiguity
into a wrong decisive call. The bands are descriptive reference bounds, not
validated probabilities or uncertainty intervals. All25 calibration structures
were La-conditioned; two-Ca PedE tests transfer across source geometry and an
occupancy hypothesis absent from that calibration.

Across the recovered125 sources, static scoring gives117 correct/1 wrong/7
inconclusive; accommodation gives119/3/3. This is a coverage–error tradeoff on
consumed, correlated structures, not an unqualified accuracy gain. No threshold,
vote, source subset or classification was changed during this investigation.
See [frozen-band analysis](CALIBRATION_RESULT.json).

## 4. Next discriminating test

Use the existing correct3/wrong1 origins to check **each metal's solvent response
and ligand charge redistribution with an independent electronic treatment**.
The saved-charge audit above found no missing convergence or mapping mismatch.
A matched two-origin × two-metal × vacuum/continuum matrix is
at most eight endpoints before valid reuse. Freeze method, state and numerical
settings before execution. Compare origin-to-origin response double differences,
not absolute cross-method totals or a desired PedE sign.

This can distinguish a solvent/electronic-model discrepancy from a response also
present in the independent treatment. It cannot by itself validate the absent
protein environment. The present local model includes strain within retained
fragments but has no omitted-scaffold deformation cost, freezes metal/PQQ pose,
and proposes motions in vacuum before solvent selection. Those are real limits,
but none uniquely explains the large **origin** spread. Prior DFT qualification
used different donor motions and does not qualify these PedE candidates; see
[model review](model_review/REPORT.md).

No independent-reference endpoints, new folds, optimizations, production changes
or recalibration were launched in this investigation. The next molecular test
is a recommendation for discussion, not an automatic continuation.

## Execution, failures and restart

All work reused existing scientific outputs. Seven one-CPU Slurm jobs accounted
for 12 allocated CPU-seconds at scheduler integer-second resolution; GPU use
and new molecular calls were zero. Subsecond process times, memory and exact
commands are retained in [EXECUTION.json](EXECUTION.json) and its workspace
submission records. These costs exclude the already-completed molecular work.

- 3338: saved-score reader failed on a mixed-role collection. Immutable failure
  retained;3339 applied the predeclared canonical-role filter and completed.
- 3340: all raw geometry/confidence results completed; a 192-scored-atom versus
  189-physical-map assumption blocked local metrics.3344 reused the completed
  raw results and resolved the explicit links for all19 local geometries.
- 3341: all54 saved metal/candidate cells decomposed; no energy calls.
- 3345: all36 origins/both implementations audited for charges, maps and SCF.
  Its flat charge table had a field-name collision; JSON was correct.3348
  re-exported the unchanged JSON with72 preserved metal identities and prefixed
  charge fields. The original table remains archived and superseded.

All seven jobs are terminal and root has read the results. No owned diagnosis
job or molecular continuation remains. Read [CURRENT.md](CURRENT.md) before
resuming; do not rerun completed chemistry on a delayed notification.
