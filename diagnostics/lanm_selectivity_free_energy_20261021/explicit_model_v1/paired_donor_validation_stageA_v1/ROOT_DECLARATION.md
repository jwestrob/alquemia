# Conditional orientation validation: source preparation

## Question and entry condition

Does one frozen donor candidate, with the passing six water parameters fixed,
reproduce untrained orientations of hydrated ions and donor ligands? This is a
transfer test against the existing provisional PBE0-D4 reference model, not a
protein affinity calculation or proof that this reference is physically exact.

Prepare adapters while donor job17621 runs. Release numerical work only after
its actual terminal selected candidate passes every original donor development
control and all eight water invariants, with independent raw/input/objective
replay admitted. A local or budget stop with failed controls leaves this branch
unexecuted. Preparation must not read reserved energy/force outcomes, bind an
unpassed candidate, select a different training iterate or change the fit.

## Exact reserved inputs

Preserve the original24-cell inventory:

- `RESERVED_PANEL.json`: `eed99b0fb6dbe7b49857d68853fe809481f6b657c928a7e75197d77c61f3432b`.
- `NATIVE_PANEL.json`: `75d3fc69703c8a5acca61f547e3fcf5d614cfe6ac653d2f3ecd3f34a93bcaf5a`.
- Both reside in workspace `explicit_model_v1/matched_dry_v1/reserved_validation_v1`.

Stage A contains exactly16 cells: aquo8 water-H rotations ±15°, and hydrated
mono-acetate, bi-acetate and acetamide ligand rotations ±15°, each for La and Dy.
Preserve every frozen coordinate, chemical state, atom order and geometric gate.
The other eight V3 protein cells remain separate, unexecuted Stage B slots.

Use unchanged small-system typing, intrinsic metal terms and the selected six
water/eight donor values. Preserve fixed donor epsilons, La amide-C override,
all other pair terms, permanent multipoles, polarization and bonded parameters.
Reference cells retain the frozen native REFERENCE.py recipe, basis/ECP assets,
La46/Dy55 cores, effective singlets and separate physical La1/Dy6 labels. Do not
splice counterpoise energies into the uncorrected energy/gradient surface.

## First: consumed numerical qualification

On the exact selected candidate, check both metals for the three hydrated central
donor geometries and three matched-dry donor geometries:12 consumed cells.
Reuse their exact selected analyze/analytic-force records and inputs. Run only
the established Cartesian `testgrad system.xyz Y Y 0.00001` check for each cell.
Maximum12 new native programs and1,344 internal finite-difference energies:
six wet cells total172 atoms, six dry cells total52 atoms, six differences/atom.

Require all12 checks, unchanged actual inputs, saved energy replay within1e-4
kcal/mol and force-component/Cartesian consistency within1e-3 kcal/mol/Å.
Keep existing solver/output admission rules. No optimizer or retry is included.
These numerical checks must pass before reserved evaluations are released.

## Then:16 reserved cells and fixed scoring

Evaluate every Stage A cell with the frozen candidate: analyze plus analytic
testgrad, at most32 native programs. Evaluate the same16 fixed cells with the
unchanged ECP-PBE0-D4 energy/gradient recipe. No free fragments or new geometries.
Retain missing, failed or reference-limited cells in the declared denominator;
they do not count as passes and cannot be replaced by a different orientation.

Freeze the scoring from `paired_donor_validation_outline_v1/OUTLINE.md`:

1. Eight paired aquo-referenced contrasts, using exact admitted central-aquo
   model and reference anchors, with model-minus-reference error≤1 kcal/mol.
2. For each of four orientation families, plus-minus work for La, Dy and their
   difference:12 errors, each≤1 kcal/mol. These require no extra molecular cells.
3. All eight configurations' La, Dy and paired rigid-body responses:24 fixtures.
   Preserve existing category/body weights, net forces, and centroid torques
   divided by1Å. Every fixture requires RMS≤2 and max component≤5 kcal/mol/Å.

Use absolute energy errors. Atomic/internal force residuals remain diagnostics.
Do not fit on reserved outcomes, average away failing fixtures or change limits
after readout. Report each family and metal. Original reference, electronic-state,
counterpoise and derivative-qualification limits remain explicit.

## Execution and ownership

This declaration authorizes source/metadata preparation only until entry gates
are met and root binds the actual candidate and reviewed executable manifests.
Root owns all submissions. Model work uses up to eight single-thread workers in
one eight-CPU standard_expand allocation; reference cells use independent four-CPU
standard_expand allocations, default memory/time and allocation-derived memory
settings. Respect actual QoS concurrency limits. Arm completion watchers and
record actual costs; no polling collector allocation or automatic scientific retry.

Singer owns the minimal reference adapter, Arendt model binding/numerical adapter,
Planck scoring and independent adapter review. Root integrates and reviews the
scoring adapter. Local numerical preparation shares: root≤2, Arendt≤2, Singer1,
Planck1, within the authorized16-thread total. No native calls during preparation.

Stage B requires Stage A success and a separate concrete protein-parameter mapping
and matching BIO aqueous anchor. Neither Stage B nor periodic dynamics, parameter
refitting, production promotion or an exchange-affinity claim is released here.
