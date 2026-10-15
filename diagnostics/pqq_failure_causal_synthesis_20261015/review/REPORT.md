# Why nine PQQ sources remain inconclusive

**There is no demonstrated single failure shared by all nine.** Six new
abstentions are an established solvent-model/calibration tradeoff. Three other
sources remain unresolved for different reasons. The subsequent ExaA experiments
rule out several proposed remedies on one matched source pair; they do not explain
every inconclusive protein or establish a universally correct chemical state.

Read-only synthesis; no analysis, molecular calls or submissions. These are five
consumed proteins, not nine independent biological tests.

## Latest mechanistic conclusion — 15 October

The XoxF1 case now has a concrete **local-search limitation**: its original Asn290
terminal amide presents nitrogen toward the metal while oxygen points away. An
intact 180° alternative presents oxygen instead. The old bounded torsion search
could not reach this separate orientation. This is a missed candidate geometry,
not evidence that the whole fold is defective or that atoms were mislabeled.

The same rule was tested on all ten XoxF1 sources. Adding the fixed alternative
changes the ambiguous Ca-conditioned sample4 paired contrast by +10.975 kcal/mol,
the correctly called La-conditioned sample4 by +0.673, and the other eight by
zero. Both sample4 sources have almost the same favorable La solvent change;
their different vacuum competition explains most of their different net gains.
A crossed geometry/charge test shows the solvent benefit survives either saved
charge distribution, while changing charges opposes it. The large mixed term
prevents unique additive attribution. This constrains the mechanism within the
model; it does not independently validate its solvent physics.

The fixed alternative leaves a short represented amide-H/backbone-H contact.
Consequently, its favorable score is a lead, not an accepted relaxed structure
or a new class decision. The first bounded-relaxation control passes physical
and optimizer checks but keeps its original winners and score. The unchanged
remaining nine-source test is now complete. Allten common6 pools are available;
Ca4 and La4 gain +16.527 and +11.464kcal/mol in paired vacuum contrast relative
to the fixed-flip pool, with the other eight unchanged. Both relevant sources
retain oxygen-facing Asn and lengthen the named H contact to2.17–2.28A; other
short included H–H contacts remain. The selected new proposals are interior.
Matched solvent evaluation is running after the frozen gate and actual preflight
passed. This supports a missed local alternative in the old search, but does not
yet establish a corrected composite classification or the true equilibrium state.

For MMOL1770, all six native DFT endpoints are now complete as energy diagnostics.
Both methods choose the same metal-specific geometries. DFT reproduces a stronger
Ca-directed pool response (−5.222 versus MACE −4.005 kcal/mol), so replacing MACE
with this DFT treatment would not rescue that response in the same representation.
This is directional support, not exact quantitative agreement: the pool discrepancy
is −1.217 and the Ca-generated-path discrepancy is −1.687 kcal/mol, exceeding the
declared 1 kcal/mol interpretive resolution. Density residuals and nonuniform
frontier-ordering warnings leave electronic stability and same-state continuity
unqualified. A MACE-only explanation is unsupported; source/state/environment and
the adequacy of the shared local representation remain unresolved.
The small PQQ rotation probe is complete and closed: neither source reaches its
predeclared solvent-continuation threshold. These results strengthen a
case-specific explanation; they do not supply one universal solvent or folding
repair for all nine abstentions.

Evidence: [Asn orientation and contacts](../a0ac_asn_geometry_response/REPORT.md),
[full fixed-alternative comparison](../a0ac_asn_solvent/REPORT.md),
[geometry/charge factorial](../a0ac_asn_solvent_factorial/REPORT.md),
[first relaxation control](../a0ac_asn_relaxation/SCOUT_REPORT.md),
[full alternative-search result](../a0ac_asn_relaxation/FINAL_MACE_REPORT.md),
[final named-contact evidence](../a0ac_asn_relaxation/contact_review/FINAL_REPORT.md),
[complete MMOL reference](../mmol_reference/FINAL_REPORT.md),
[closed PQQ pose probe](../pqq_pose/REPORT.md).

## 1. What changed the classifications

Changing from self-consistent GFN2 solvent subtraction to frozen-distribution ALPB
also changed the canonical bands under the unchanged calibration rule. The Ca edge
moved −4.398647 kcal/mol and the La edge +0.579505. Eight of the nine targets retain
their previous selected geometry IDs. ExaA La2 alone reselects, and that reduces
its contrast by 1.007315 kcal/mol, helping its Ca margin. Thus a newly broken
optimization cannot explain the six new abstentions.

The candidate fixes both wrong two-Ca PedE calls but loses Ca-control decisions.
Independent PedE DFT gives +10.0342 kcal/mol solvent contrast versus corrected
GFN2's +21.2566. Replay assigns +18.7467 of the latter to additional redistribution
plus intrinsic-energy compensation. That supports investigating electronic response,
not removing all response or fitting a response fraction to these labels.

## 2. Exact cases and strength of the diagnosis

Distances below are shortfalls to each expected-class band, copied from the saved
review; they are not physical error bars.

| Exact source | Established result | Still unresolved |
|---|---|---|
| XoxF1 A0ACD6B9F2 Ca4, La class; 5.815 kcal/mol short | Wrong Ca becomes inconclusive; original accommodation helps by +12.176. A separately proposed intact Asn orientation, unreachable by the old domain, adds +10.975 to the current model's paired contrast. Full relaxed competition adds a further +16.527 in vacuum relative to the fixed-flip pool; the named short H contact lengthens while oxygen-facing Asn persists. | Relaxed composite competition is running; calibration remains untested and other short H contacts persist. Correctly called La4 shares the starting orientation and also benefits from its alternative, so orientation is not a sufficient failure/exclusion rule. |
| XoxF4-2 MMOL1770 Ca0, La class; 0.812 short | Correct static origin becomes inconclusive after −4.186 accommodation, mostly MACE (−4.005). Ca proposal compresses Glu200 to 2.151 Å and reaches the movement bound; La does not. Complete same-geometry DFT gives a stronger Ca-directed pool response (−5.222), with the same winners. | Quantitative differences exceed the declared resolution on the pool and Ca path; electronic stability remains unqualified. Shared local representation, source/state/environment or calibration may still be inadequate. Ca1's extra-Asp experiment and La4's hydrogen repair concern different sources. |
| ExaF C5AXV8 two-Ca4, La class; 0.196 short | Already inconclusive before the solvent change. Both metals select the same La proposal; no bound hit or distinctive gross geometry defect. | Fine energetics, occupancy/state or calibration coverage; no case-specific physical defect is established. Preserve abstention. |
| ExaA Q9Z4J7 two-Ca3, Ca class; 1.191 short | Newly inconclusive with unchanged selections. Large donor response is reproduced by independent vacuum DFT. Metal/PQQ placement differs from same-condition peers. | Whether that source geometry or the shared dry/cofactor/environment representation is physically inappropriate. Tested simple repairs do not resolve it. |
| ExaA La2, Ca class; 1.182 short | Newly inconclusive; reselection helps. Both proposals hit bounds, but origin is correct under static bands. | Its own electronic motion response and chemical environment; the two-Ca DFT result cannot be silently transferred to this source. |
| ExaA La3, Ca class; 4.005 short | Newly inconclusive; origin short-Glu warning is real and proposals rebalance its donor contacts. | A short contact or a bound is not sufficient: successful peers share those features. Its own motion reference remains missing. |
| ExaA La4, Ca class; 1.583 short | Newly inconclusive; static origin was already inconclusive. Both rows select La proposal, and both searches hit bounds. | No unique malformed source or failed optimizer has been established; bound widening is unsupported. |
| PedE Q88JH5 La3, Ca class; 1.835 short | Newly inconclusive; static origin was already inconclusive. Complete CN7 pool, high confidence, neither proposal hits the bound. | Fine solvent/density/state response; no shared gross fold defect or demonstrated need for another start. |
| PedE La4, Ca class; 1.425 short | Newly inconclusive; static origin already inconclusive. One bound hit also occurs in correct peers. | Same limitation as La3. The completed PedE DFT test used two-Ca1/3, not these La-conditioned motions. |

## 3. What the later ExaA tests establish

- **Electronic surface:** DFT gives target/control differential works
  +17.13547/+5.37792 versus MACE +16.40348/+5.26150 kcal/mol on the same
  La-generated path. A MACE-only exaggeration is unsupported. Two energy-only
  SCF endings and state-stability warnings remain; both methods share the same
  conditional local representation.
- **Packing:** the tested complete added fragments change the target differential
  by only −0.374789 kcal/mol without reselection. That particular omitted-packing
  explanation is weak; all scaffold mechanics were not tested.
- **Metal coordinate:** reciprocal placement changes the donor-path response but
  adds zero gain against the full original pool. The midpoint's target paired
  change is only −0.176971 kcal/mol. The tested reciprocal and midpoint probes did not provide a useful pool improvement.
- **Cavity exclusion:** origin shifts are mostly common. The matched follow-up
  changes source-dependent accommodation by only −0.076458 kcal/mol and leaves
  all selections unchanged. Aggregate precision flags remain false. This closes
  the tested fixed-charge Gelec exclusion remedy, not HB/nonpolar burial,
  permanent protein electrostatics, density response or hydration.
- **Crystal:** conditional 1FLG A has smaller accommodation and a Ca call.
  Its static→accommodated call transition uses different bands; raw motion is
  La-directed. Crystal B remains unsupported. Neither establishes the physical state.

All nine current states are dry and use oxidized PQQ3− while folding supplied a
neutral-acid PQQ graph. The PedE audit found complete PQQ, correct mapped H and
no cofactor cap or source-to-core distortion; unusual quinone-region bonds also
occur in successful predictions. Fixed PQQ geometry could still couple to donor
strain. No result proves a wrong protonation/redox state, missing-water cause or
unique defective fold. Crystal bond lengths must not become hand-fitted targets.

## 4. Declared reference tests and current disposition

1. **MMOL1770 Ca0 reference is complete, six of six.** Exact admitted 193-atom
   state and saved origin/adaptive_Ca/adaptive_La under both metals, native
   r2SCAN-3c/DefGrid3/TightSCF vacuum. Both methods choose Ca/adaptive_Ca and
   La/adaptive_La; DFT strengthens the Ca-directed response rather than reversing
   it. Do not describe the quantitative discrepancies or electronic-state warnings
   as a perfect reference match. This closes the finite MACE-versus-DFT motion
   test, without proving a bad fold or identifying a unique missing interaction.
   No automatic rerun or method sweep follows.
2. **The finite rigid-PQQ response probe is complete and closed.** The saved-force
   audit identified a coherent rigid mode in ExaA target and control. Its tested
   ±1° movement did not meet the frozen continuation criterion; do not rerun or
   widen it from this old proposal. The broader initial question was cofactor strain,
   including Ca/La coupling. Internal cofactor flexibility remains untested, not
   established as the missing cause. No guessed carbonyl
   shortening, crystal transplant or redox switch. Existing ExaA DFT inputs are
   SP, not EnGrad; cached DFT analytic forces are not established. A shared force
   alone does not explain selective failures or justify recalibration.

ExaF's near-boundary abstention and XoxF1's distinct Asn-basin hypothesis should
remain explicit rather than forcing these two tests to serve as universal fixes.
Another global solvent factor, blanket refolding or automatic search enlargement
is not supported by the completed evidence.

## Main evidence

- [Nine-source scores and identities](../../pqq_inconclusive_review_20261014/REPORT.md),
  [geometry and corrected-contact findings](../../pqq_inconclusive_review_20261014/geometry/REPORT.md).
- [Exact model/selection/band decomposition](../../pede_solvent_response_20261013/frozen_distribution_candidate/selection_diagnosis/REPORT.md),
  [full candidate denominator](../../pede_solvent_response_20261013/frozen_distribution_candidate/REPORT.md).
- [PedE independent solvent reference](../../pede_solvent_response_20261013/REFERENCE_REPORT.md),
  [native replay](../../pede_solvent_response_20261013/fixed_charge_replay/REPORT.md),
  [cofactor-state audit](../../pede_solvent_response_20261013/cofactor_review/REPORT.md).
- [ExaA vacuum reference](../../exaa_donor_response_20261014/FINAL_REPORT.md),
  [packing/crystal](../../exaa_source_packing_20261014/REPORT.md),
  [full-pool/force test](../../exaa_metal_position_20261014/FULL_POOL_FORCE_REPORT.md),
  [midpoint](../../exaa_metal_midpoint_20261014/REPORT.md),
  [cavity-dependent accommodation](../../native_cavity_accommodation_20261014/REPORT.md).

These links point to completed evidence and its source pins; historical proposed
next steps inside them are superseded by their later actual results above.
