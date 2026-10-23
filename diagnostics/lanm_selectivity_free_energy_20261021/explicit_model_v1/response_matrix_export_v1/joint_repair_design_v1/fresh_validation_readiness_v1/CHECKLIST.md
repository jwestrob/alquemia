# Fresh validation: existing gates and remaining binding work

Source-only review while fit18652 runs. No fit trajectory, held validation
outcome, molecular calculation or new threshold was used. This is a delta to
the existing declarations, not another validation protocol or an execution release.
Full paths and SHA-256 pins are in `SOURCES.json`.

## Decision

**The 25 water-development gates plus consumed trimer regression do not by
themselves release any fresh calculation.** If all25 pass for the one terminal,
frozen candidate, its actual numerical/input/raw review passes, and root accepts
the complete consumed response/regression report, the existing four-cell W9
radial test is the next compatible fresh test. Its old adapters still require
job17676 and six parameters; they cannot bind this twelve-parameter candidate.
Version only those candidate/source interfaces and retain the frozen geometry,
reference recipe, scorer, limits and full denominator.

The six additional trimer evaluations complete eight **consumed** E3/mu3 cycles.
There is no declared physical pass threshold for E3, dipoles or susceptibility.
Their descriptive improvement, degradation or small joint loss cannot replace
all25 E/F gates or be relabeled an independent validation pass. The regression
adapter intentionally runs without an E/F-pass prerequisite, so it can diagnose
a failed candidate. Root must make the scientific progression decision; this
note does not invent a new response-error threshold.

Donor/protein validation has a separate, concrete blocker: donor17621 fails
21/33 energy checks and27/36 force fixtures. A water-only pass does not satisfy
the consumed donor gate. The current nonzero local-response implementation also
rejects all atom types except36/37/500/501 and systems beyond its vacuum,
single-metal, complete-water, at-most28-atom scope.

## 1. Four fresh W9 radial cells: smallest compatible release

Reuse `water_coordination_radial_validation_v1/ROOT_DECLARATION.md`48938d14,
geometry manifest a34e9d05, reference panel04d99bec and scorer92678201.
These freeze inward/outward0.10Å rigid translations of all nine waters along
their original metal–O directions, identical La/Dy coordinates, four28-atom cells.
Do not optimize, reselect a geometry or replace an unfavorable sign.

Before root releases either molecular arm:

1. Freeze exactly the terminal selected twelve-parameter vector and its full
   selected46 artifact. Admit source/build/analyze/testgrad pins, exact effective
   inputs, native outputs, current531-component score and all25 original gates.
   Bind the passed57-program numerical qualification to the same native source
   and executable. Preserve actual selected-cell stability/induction/readback
   evidence; the fixed numerical probes are not physical acceptance tests.
2. Complete/report all eight consumed trimer cycles with the same candidate,
   retaining each outcome and any unavailable cell. Freeze the candidate before
   opening radial outcomes. No alternate iterate chosen using regression or
   validation outcomes is implied by the original strict-J fit selection.
3. Add a narrow versioned model binder using the current evaluator's rendering,
   input admission and E/F parsing. Reuse exact prepared displaced XYZ/graphs;
   copy/reproduce each selected central-W9 effective key/parameter bytes. Use
   the **new selected candidate's** W9 central E/F as model anchors. Reuse the
   original17670 central reference E/F and unchanged ECP-PBE0-D4 recipe.
4. Version the old reference release check: `reference/RUN.py`121–124 hardcodes
   job17676, six parameters and `selected_ten_raw_review_passed`. Update the
   model/reference/collector candidate joins together to the actual new terminal
   review. Do not bypass the old failed gate by changing only a job label.
   Retain existing decks, parser/collector missing-slot behavior and fixed scorer.
5. Root binds the reviewed source/config/result inventory and exact finite scope:
   four reference E/gradient calls, four MPI each; eight model programs
   (analyze plus analytic testgrad), four single-thread workers; zero new Cartesian
   finite-difference energies. Existing default scheduler resource policy applies.
   Preserve first-failure stops, all four slots, raw warnings and actual costs.

For model/reference X, define W_X,m,sign = E_X,m,sign − E_X,m,central.
Require all four individual work errors and both
`(W_Dy−W_La)_model − (W_Dy−W_La)_reference` errors to have absolute value≤1kcal/mol.
For each sign retain La, Dy and Dy−La rigid responses: all six require category-
balanced RMS≤2 and maximum component≤5kcal/mol/Å. Use the unchanged net-force/
centroid-torque helper (torques divided by1Å). **All12 comparisons are required.**
Missing/reference-limited rows are unavailable, never passing zeros.

These are model-to-provisional-reference agreement criteria, not uncertainty
bounds or proof of bulk hydration, electronic-state accuracy or exchange affinity.

## 2. Before the twelve held donor orientations

Reuse the consumed28-row training authority99ea27f8 and unchanged donor
`OBJECTIVE.py`86e85fcb with body helperbf65a20a:

- Four water/aquo centers,18 hydrated donor cells (mono-acetate, bi-acetate,
  acetamide at minus/central/plus for both metals), and six matched-dry cells.
- All33 energy gates:12 aquo-referenced contrasts,3 wet-minus-dry contrasts,
  12 individual radial works and6 paired radial works; every absolute error≤1.
- All36 La/Dy/paired donor-force fixtures: each RMS≤2 and maximum≤5. Preserve
  full water qualification separately; the old four-center objective does not
  stand in for all eight original water rows or the current ten-cell/25-gate panel.

The pure objective already requires an explicitly supplied **future frozen
candidate** water baseline; it can be reused without changing its arithmetic.
The old evaluator and validation binder hardcode17605 water/17621 donor results
and cannot be reused as candidate authorities. A new consumed replay must bind
the current water candidate and explicitly declared fixed donor pair terms,
epsilons and La amide-carbon override. Do not silently combine failed17621 donor
values with the new water fit, or substitute old aquo anchors. Exact dry-cell
reuse needs matching effective inputs and qualified implementation behavior.

First implement and numerically qualify the mixed-donor path while preserving
the same h/H definition, native donor tensors/exclusions/local frames and
unrelated interactions. Merely removing the native type/size guard is not a
qualification. Retain the existing exact-candidate donor numerical stage:
twelve consumed wet-central/dry cells, at most12 testgrad programs/1,344 internal
FD energies, total step1e−5Å, energy replay1e−4kcal/mol and Cartesian/force replay
1e−3kcal/mol/Å. All consumed donor controls and water checks must pass before any
held donor outcomes are opened. A fixed-parameter replay is sufficient to test
this prerequisite; no donor refit is authorized by this note.

Then reuse the original24-cell reserved manifest eed99b0f and native panel
75d3fc69, with the reviewed StageA scorer. Keep the full16-cell/44-check StageA
report: eight aquo-referenced contrasts, twelve plus-minus La/Dy/paired works,
and24 rigid fixtures at the same1/2/5 limits. Its four aquo-orientation cells and
11 checks are now **historical regression**. The remaining twelve donor cells
provide15 energy and18 rigid checks that are still fresh. This is a provenance
partition, not permission to drop the water checks or choose a donor subset.
Use candidate-specific aquo anchors and water-model replay; old model scores
cannot substitute when the water candidate changes.

## 3. Eight held protein slots and periodic/GPU separation

The existing donor outline requires complete StageA success before StageB. Keep
all eight V3 ACE–Asp/Glu/Asn/Trp cells, exact state/graph/cap/contact evidence;
V2 is withdrawn. Prospective existing limits are four aquo-referenced protein
energy errors≤1 and twelve La/Dy/paired rigid fixtures at2/5. No protein radial
work target is present. A concrete StageB executor/scorer and mapped native
qualification remain to be supplied; the old outline is not executable support.

Required mapping is substantive: BIO900/901 metals; water O402/class73 and
H403/class74; carboxylate O143/161/class31; represented amide O11/153/**227**/class5;
La amide-carbon class3 includes9/152/**226**. Preserve BIO multipoles, local frames,
intrinsic polarizabilities, bonded terms and the actual represented inventory.
Small-system acetate403 is not BIO water-H403. Match a BIO-water aquo anchor;
SI/BIO Urey distances1.5337/1.5537Å remain distinct. The new response selector
and effective metal response require an explicit mapped definition and numerical
qualification, not only VDWPR/POLPAIR remapping. No protein outcome read is needed
to prepare that source support.

Arendt's separate periodic/GPU extension is not released by any finite vacuum
pass. Native periodic/GPU implementation/parity, mapped water/salt physics and
aqueous stability/sampling remain separate requirements for exchange free energy.

## Failure decisions

- Any failed current development gate: keep fresh W9/donor/protein outcomes held;
  use the admitted consumed residual/regression evidence for a separately declared
  mechanics decision. A local/budget stop does not prove family impossibility.
- Technical/native/reference admission failure: preserve its unavailable slots
  and costs; no automatic molecular retry, threshold change or physical pass.
- Fresh radial or later donor/protein scientific failure: retain the whole
  affected fixed panel and stop progression under its declaration. Using it for
  a later repair consumes its independence and requires a new declared test.
- Passing W9 transfer does not resolve failed donor development. Passing donor
  StageA does not supply missing BIO implementation, periodic/GPU support or an
  exchange-free-energy result.
