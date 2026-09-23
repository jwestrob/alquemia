# Collective scaffold response: a matched-target proposal experiment

**Proposed only — no new molecular energy, force, optimization or job.**
September 22, 2026. This implements the *design* part of Jacob's request for
higher-upside scaffold accommodation. It does not authorize a new production
score. The released scorer and completed experiments remain unchanged.

## Recommendation and discriminating question

First test **whether the same donor movement becomes more useful when the
surrounding protein can respond**, using a standard protein force field to
*generate* collective geometries and the unchanged compact composite to assess
them. Do not begin by adding an unqualified scaffold penalty to the score.

For each source, compare its archived adaptive donor geometry against a new
geometry having **the same donor positions**, but allowing nearby complete
residues, peptide backbone and second-shell contacts to move together. Also run
an otherwise identical response with donors held at their original positions.
This q0-target control distinguishes ordinary source-to-force-field relaxation
from response specifically associated with donor accommodation. It does not
perfectly remove force-field bias, and its result must be shown separately.

This changes the mechanical space: multiple peptide units and neighboring
residues can respond simultaneously, with real source connectivity across their
boundaries. It is not another independent donor torsion or the existing eligible
single peptide-crankshaft mode. It deliberately keeps donor displacement size
fixed; merely enlarging the old angular box would not answer this question.

The completed adaptive experiment motivates this test: matched individual folds
improved 200/204 to 202/204 correct, but one error remains, many structural spreads
widen, and 172/409 accepted proposals touch a bound. The bounds do **not** establish
that omitted scaffold compliance is the cause. See the
[actual comparison](../nikasha_next_phase_20260922/REPORT.md).

## Exact proposed construction

Let x be source protein coordinates, x0 their original values, D the actual
protein donor atoms identified by the existing core role/source map, and L(x)
the **unchanged** local context reconstructed from physical source atoms.
For t in {0, Ca-adaptive, La-adaptive}, let d_t be the donor coordinates in the
original or the corresponding archived, accepted adaptive candidate. Generate:

    x*_t = constrained_local_minimize U_ff19SB,protein(x)
    subject to x_D = d_t and the source-geometry constraints below.

The parent is the already inventoried protein-only ff19SB function: bonds,
angles, proper/improper torsions, CMAP, and native NoCutoff Coulomb/LJ with exact
1–4 exceptions. It has no invented La/PQQ parameters and no added solvent.
Initialize each search from that target's actual source-connected candidate;
use exactly one start, with no outcome-directed restart. Full nonlinear forces
retain the source's nonzero affine force; there is no zero-gradient assumption,
Hessian inversion, negative-mode clipping, spring fitting or entropy term.

Proposed mobile set, chosen from x0 before energies: complete protein residues
having any heavy atom within 8 Å of a canonical-core heavy atom, plus complete
core-contributing residues and their directly bonded peptide neighbors. Use
actual bonds, not residue-number arithmetic. Freeze all other protein atoms.
This radius is a **new proposed mechanical-domain choice**, not an established
physical optimum. Keep the former 0.8 Å final source-heavy displacement limit;
report boundary outcomes rather than enlarge it after inspecting scores.

Preserve source covalent bond lengths, including X–H, via explicit sparse
geometric constraints, and preserve chirality, peptide cis/trans identity and
independent overlap guards. Fix PQQ and metal at the original coordinates;
fix carver-generated PQQ H. The existing angular candidates do not translate the
metal or change the cofactor. Donor targets therefore differ only by the same
already admitted protein accommodation. Full source bond connectivity must be
validated at the initial target before any force call. Do not normalize H lengths
or interpret the known source/FF H-length convention difference as scaffold strain.

Implementation needs a small constrained Cartesian adapter around the installed
OpenMM energy/force function, with analytic sparse distance/target Jacobians
and fixed coordinates eliminated. Existing local Kinematics can supply the
actual target and cap maps, but does **not** already implement this general
collective constrained search. Proposed one-search algorithm limit: 200 accepted
optimizer iterations, explicit stationarity/constraint residuals, no claim that
a boundary or nonconverged candidate is a protein minimum. The exact solver and
residual acceptance values must be fixed in an implementation plan before forces;
this note does not pretend that adapter has been validated.

Map every moving physical atom back by chain/residue/insertion-code/atom identity.
Caps are reconstructed from their original retained/omitted source bond and fixed
cap length, never independently optimized. Context membership, charge,
multiplicity, proton inventory, water inventory and original cap definition stay
unchanged within each source. Failed parent templates or missing mapped atoms
make that source unavailable; do not add new terminal/protonation repairs.

## Scoring and the essential limitation

For every admitted new geometry g = L(x*_t), score both metals with:

    E_M(g) = E_OMOL,vac,M(g) + E_GFN2,ALPBwater,M(g) - E_GFN2,vac,M(g).
    R_pool = min_g E_Ca(g) - min_g E_La(g).

The shared pool retains the existing origin/old/adaptive candidates and appends
all three new geometries. There is no donor-target-specific winning rule. Show
R for each matched target and the q0-only response control **before** interpreting
the augmented pool. Keep missing required cells explicit. This is an electronic
accommodation descriptor, not an equilibrium mixture or a binding free energy.

**Do not add U_ff19SB to this score.** The outer protein does mechanical work in
proposal generation, but that deformation work is absent from the reported
composite. This first version can reveal a useful geometry mechanism; it cannot
establish complete scaffold-coupled energetics. Fixed PQQ/metal also exert no
force on exterior FF atoms in the generator. Independent geometrical guards
prevent overlap but do not replace those missing interactions. Any observed
benefit therefore requires later scrutiny, not an immediate physics claim.

## Why an additive score is not ready

The [existing real inventory](../scaffold_environment_20260922/REPORT.md)
parameterizes 1H4I/4MAE/PLM8344 protein parents (9,060/8,826/8,678 atoms). Its
source/context/cap maps are usable; 1H4I and 4MAE reuse exact archived distant OXT
completions. No parent force or energy has yet been measured.

1. Subtracting terms wholly inside the local protein is not a matched reference:
   the capped MACE context already represents part of boundary angular/many-body
   response. Actual chopped local fragments fail ff19SB templates for 10/9/11
   residues. Adding every mixed parent term could count boundary stiffness twice.
2. Installed ff19SB XML contains genuine ACE and NME templates. A **different**
   complete-peptide local context could in principle support the expression
   E_composite(new local) + U_FF(full protein) - U_FF(same capped local protein).
   That would require a new source closure/cap/proton map and changes both the
   learned system and solvent cavity. Template presence alone does not validate
   the matched subtraction or its cost.
3. Even that expression omits exterior protein–PQQ/metal coupling. When exterior
   atoms move, it is not a canceling constant. Missing cofactor/metal parameters
   and a consistent electrostatic/solvent reference remain real blockers.
   Whole-protein native OMOL/masked subtraction has already been tested; it is
   not a new easy bypass.

Thus the proposal-only route is the smaller first experiment. A broader hybrid
energy is deferred explicitly, rather than smuggled into a force-field penalty.

## First proposed panel, tasks and decision

Use the parent's common eight consumed sources, with no new folds or labels:

- 1H4I (Ca control) and 4MAE (La control).
- Canonical Q9Z4J7 and Q88JH5 (Ca guardrails).
- A0A3F2YLY8 Ca-conditioned samples 1 and 3 (La class).
- A0ACD6B9F2 Ca-conditioned sample 4 and La-conditioned sample 4 (La class).

These are four biological reference groups plus two crystal control identities,
not eight independent biological tests. Selection is explicitly developmental.
Only the two crystal protein parents above are already demonstrated available;
parameter availability for the other six exact sources remains a preparation
gate. PLM8344 supplies an existing mapping example, not a labeled success case.

Maximum: **24 constrained protein searches, 48 fresh compact MACE singlepoints,
96 native GFN2 singlepoints**; archived original/adaptive pool cells reused only
with exact state/recipe/coordinate compatibility. Every search's actual FF calls,
failures and allocated CPU/GPU time are counted. This adds three geometries to
an existing five-geometry pool: compact cross-scoring grows by at most 1.6x that
research pool's total, not 1.6x the released two-endpoint scanner. No parent-FF
runtime has been measured, so overall constant-factor affordability is a
hypothesis. Use the existing GPU MACE worker and CPU native-GFN2 runner; benchmark
protein-force throughput within the first declared search instead of creating
a separate static-energy campaign.

Primary readout: matched-target changes in La/Ca contrast and known-class margins,
with the q0 response and both donor targets reported individually; then common
pool benefit, failures, structural spread for the two paired source groups, and
Ca guardrail preservation. Report old-band transfer as developmental only;
no calibration fitting on these eight or claim of full-panel classifier fidelity.
Energy lowering alone, a moved unknown prediction, or movement away from a bound
is not success. A useful signal would be repeatable benefit specific to donor
accommodation beyond the q0 control, with no Ca guardrail regression. A null
result limits this proposal mechanism; it does not disprove all scaffold physics.

**Remaining decisions before execution:** agree the 8 Å mobile-domain rule and
matched-target design; then freeze the constrained solver, feasibility tolerances,
stationarity test and near-PQQ overlap checks. If source-exact constraints cannot
be satisfied by the existing candidates, report the precise mapping failure;
do not silently loosen geometry or introduce a new chemical state. No molecular
execution is authorized by this design note alone.
