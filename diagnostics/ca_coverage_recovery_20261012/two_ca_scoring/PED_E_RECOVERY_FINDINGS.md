# PedE recovery: two source-dependent wrong calls

**Samples 1 and 2 are genuine wrong La calls under the frozen accommodated
policies.** Their origin calls are inconclusive; samples 0, 3 and 4 are correct
Ca at origin and after accommodation. All five are scored and retained. The
original exclusions remain in the primary ledger, alongside these alternative
source results. No new five-fold vote or recalibration is introduced.

## What changed in the score

`R = E_Ca - E_La`, kcal/mol. Work is each metal’s selected composite energy
minus its origin energy. The identity is `R_after = R_origin + W_Ca - W_La`.
Operational and mathematical selections agree for every sample.

| PedE sample | Origin R / call | Accommodated R / call | Ca selected / work | La selected / work |
|---|---|---|---|---|
| 0 | -405478.057649 / Ca | -405468.388555 / Ca | Ca proposal / -2.122885 | La proposal / -11.791978 |
| 1 | -405459.559618 / inconclusive | -405450.140587 / La | Ca proposal / -3.417673 | La proposal / -12.836704 |
| 2 | -405464.018944 / inconclusive | -405453.815599 / La | La proposal / -2.101846 | La proposal / -12.305191 |
| 3 | -405480.665833 / Ca | -405470.082612 / Ca | La proposal / -1.079651 | La proposal / -11.662872 |
| 4 | -405470.408100 / Ca | -405464.435887 / Ca | Ca proposal / -1.291723 | La proposal / -7.263936 |

Both wrong sources start more La-like in raw contrast than the three correct
sources. La gains more energy than Ca in **all five** accommodation comparisons.
The wrong samples’ La works are similar to correct samples 0 and 3; these records
do not identify an exceptional relaxation magnitude exclusive to the failures.
The supported description is pre-existing source sensitivity followed by a
La-favoring accommodation shift. Its physical cause remains unresolved.

Origin uses its own frozen bands (Ca ≤ −405465.7535534819; La ≥
−405455.45587457). Operational/mathematical bands are Ca ≤ −405457.564455343
and La ≥ −405454.092284296. These are protocol-specific raw scales, not affinity.
Sample 2 lies just on the La side of the accommodated boundary; no boundary is
moved in response. Unrounded values and the exact calibration pin are in JSON.

## Optimizer and boundary evidence

All ten bounded proposal optimizations report success and physically feasible
final geometry. They minimize vacuum OMOL energy; the composite is subsequently
used to select a finite geometry. No composite gradient or unconstrained minimum
is claimed. Each source uses the same four selected modes for Ca and La.

| Sample | Ca / La iterations | Ca proposal at bound | Ca / La maximum heavy displacement (Å) | Selected geometry at bound? |
|---|---:|---|---|---|
| 0 | 12 / 14 | False | 0.506930621 / 0.509457568 | No |
| 1 | 13 / 16 | False | 0.628989272 / 0.610311300 | No |
| 2 | 13 / 14 | True | 0.800000000 / 0.783920834 | No |
| 3 | 14 / 13 | True | 0.800000005 / 0.720464854 | No |
| 4 | 13 / 15 | False | 0.279974435 / 0.432665944 | No |

- Wrong sample 1 selects its own interior proposal for each metal. Its selected
  modes are A/221/chi1, A/358/chi2, A/221/chi3 and A/358/chi1.
- Wrong sample 2 reaches the 0.8 Å heavy-displacement boundary in its Ca search,
  but **both metals select the interior La proposal**. Its selected modes are
  A/221/chi1, A/358/chi1, A/308/chi1 and A/221/chi3.
- Correct sample 3 also reaches that boundary during its Ca search and selects
  the interior La proposal for both metals. The boundary flag alone therefore
  does not separate correct and wrong calls or establish the failure mechanism.
- Rejected infeasible trials remain recorded. The final selected geometries pass
  existing admission checks; failures are not being converted into successful
  native energies or origin fallbacks.

The vacuum-force diagnostics retain omitted angular loads. For wrong sample 1
Ca/La omitted maximum loads are 9.700574476676383 / 18.17421288098789 kcal/mol/Å;
for wrong sample 2 they are 6.647920215440698 / 7.43385728831857. Correct sample 3
also has a La omitted load of 7.422492104369335. These indicate limited search
freedom under the proposal Hamiltonian; they do not establish how expanding
freedom would change the composite contrast. No expansion was run here.

## State and numerical gates

Every source retains the same 192 atoms, 14 caps and frozen state key
`5258d9144f78a32f95917b3637f62adc8063ef9bd1774d72bfd25fed84033360`.
The exact state-anchor comparison, physical maps and paired maps pass. Tasks use
Ca charge −3 / La charge −2, both multiplicity 1. Original-selection replay and
all paired origin/proposal coordinate errors are exactly zero. No proton, water,
cap or local membership change is recorded between the compared preparations.

All thirty required metal/candidate solvent-component records are accepted.
The saved xTB–tblite component residuals and population residuals pass the frozen
1e-6 Eh / 1e-5 e gates, including both wrong samples. All three native lanes passed
their two real force/energy fixtures. These checks support numerical consistency
of this execution; they do not validate biological accuracy or the Hamiltonian.
Historical absolute-v1 parity remains failed and is not repaired by component-v2
acceptance. Exact individual residuals and raw-result pins are retained in JSON.

The spectator-admission records pass for every source. Their saved nearest
retained-heavy distances are about 23.7–23.8 Å. Both Ca ions remain in full
sources; the distant spectator has no local energy term. No global electrostatic
or mechanical response was evaluated, and these errors do not prove that missing
spectator physics is their cause.

## Conclusion and limits

The two errors survive this bounded saved-output inspection. No demonstrated
solver, paired-state or endpoint-selection failure explains them. Starting
geometry and the finite accommodation model are live hypotheses, alongside
electronic/composite model error; this inspection cannot identify one uniquely.
It contains no structural-distance recomputation, energy decomposition beyond
saved components, new curves or molecular calculations. Do not retune bounds or
bands, choose favorable sources or invent a five-fold aggregate around these
outcomes.

For context, the actual final 125-source annotated ledger records
**119 correct / 3 wrong / 3 inconclusive / 0 unavailable** under operational and
mathematical scoring; origin records **117 / 1 / 7 / 0**. Parent delivery owns the
full comparison and any next scientific decision. Five PedE samples remain one
consumed protein, with alternative fold settings and occupancy explicitly noted.

## Evidence

- [Pinned copied-field record](PED_E_RECOVERY_FINDINGS.json): all five source
  decisions, shared pools, works, state admissions, exact optimizer reports,
  residual loads, coordinate pins and component gates.
- Native source: `workspaces/ca_coverage_recovery_20261012/two_ca_scoring_c5q88_v1/NATIVE_COLLECTION.json`.
- Corrected solvent source: `workspaces/ca_coverage_recovery_20261012/two_ca_solvent_c5q88_v1/FINAL_COLLECTION.json`.
- Final denominator: `workspaces/ca_coverage_recovery_20261012/two_ca_scoring_final_v1/RESULT.json`.
- No scientific computation, reclassification, submission, rerun or production
  mutation was performed for this report. All displayed values were copied from
  the actual completed Slurm outputs; proposal-file hashes were checked as metadata.
