# Frozen TABI-PB numerical protocol — 2026-09-16

Frozen before scientific surface-solver execution, under Jacob's approval of
`diagnostics/accuracy_strategy_20260915/GLOBAL_ELECTROSTATIC_PROPOSAL_20260916.md`.
This records numerical implementation choices, not a new physical model.
The numerical schedule was confirmed by the parent executor before execution.

## Pinned implementation and actual controls

Use the unmodified archived APBS 3.4.1 distribution's standalone TABI-PB and
NanoShaper binaries, copied into this campaign's software directory. Leave
archived files and permissions untouched. SHA256 values:

- `tabipb`: `2e3b2ef3016982438401896ad0817fda5095b3042dc23a601e4e9e96406e0a9b`
- `NanoShaper`: `6d41b6f0ddd400aa5713d08b3a178884bfff4c7ae826530ab73549aa93584baf`

APBS v3.4.1 pins TABI source commit
`fe1c237b057418fed48535db125394607040d9de`. The inspected source is preserved
under `workspaces/global_electrostatic_20260916/software/source_inspection/`.
The installed headers and native output strings agree with that interface.

The standalone input `sdens` is passed directly to NanoShaper `Grid_scale`.
This is a linear sampling scale in inverse Angstrom, with spacing its inverse;
it is **not** the actual number of mesh triangles per square Angstrom. Preserve
the actual vertex/face counts and mesh files. The NanoShaper SES triangulation
uses its grid only for surface generation; this is not the old volume-grid PB
solver.

| Level | Grid_scale, A^-1 | Nominal surface-generation spacing, A | Tree degree | Tree theta |
|---|---:|---:|---:|---:|
| primary | 2 | 0.5 | 5 | 0.5 |
| refined | 3 | 1/3 | 5 | 0.5 |
| tree | 2 | 0.5 | 7 | 0.3 |

All use tree leaf size 500, preconditioner off, and native GMRES relative
tolerance 1e-4, restart 10, maximum 100 iterations. These GMRES controls are
hard-coded in this pinned executable and are not user-adjustable input
keywords. The iteration maximum is the solver's convergence criterion, not
a project compute/time budget. Nonconvergence remains failure; no substitute
energy or extra numerical sweep is authorized.

Physical settings remain: interior dielectric 1, exterior 78.54, zero salt,
298.15 K, SES probe radius 1.4 A, previously recorded Bondi CHNOS radii, and
common Ca/La radius 1.80 A. Native SES settings include smoothing, accurate
triangulation, Grid_perfil 90, and no cavity filling; retain the generated
configuration verbatim. Set NanoShaper `Number_thread=1` and BLAS/OpenMP thread
environment to one for resource-aware independent task concurrency.

## Surface and charge separation

The source physical atoms, sorted by stable source ID, produce a separate
XYZR file with coordinates at ten-decimal Angstrom precision and radii at six
decimals. Production meshes never include synthetic cap spheres. Charged
QM caps remain independent charge sites in the solver PQR.

Native TABI generates its own XYZR and invokes `NanoShaper` through PATH.
A pinned adapter bridge verifies its meshing parameters, replaces that XYZR
with the source-only file, then runs the unchanged pinned NanoShaper. It saves
the input/configuration, output, vertex and face files before TABI's native
cleanup. Preserve a geometry hash excluding only descriptive MSMS headers;
the actual mesh must match across endpoints/partitions for each identical
source geometry and numerical level. The tree check must share the primary
mesh. No charge-dependent meshing or cap-dependent cavity is permitted.

The two explicitly labelled isolated-reduction controls are the exception:
they use the same QM charge sites with no environmental charges and the
declared isolated core-plus-cap surface. They are numerical reduction checks,
not terms subtracted from the proposed production descriptor.

## Frozen solver task schedule: 25 tasks

- Four 1H4I endpoints (qm33/qm36 by La/Ca), each primary/refined/tree: 12.
- qm33 La/Ca translated at primary, and rotated at primary: 4.
- Independent repeat of primary qm33 La: 1.
- Isolated reduction of qm33 La/Ca, primary: 2.
- Core-only charges on the unchanged full physical surface for all four: 4.
- Environment-only charges on that full surface once per partition: 2.

Rigid transformations move **every** physical and charge coordinate together.
Translation is (0.173, 0.271, 0.389) A. Rotation is 37 degrees about the
normalized (1,2,3) axis through the coordinate origin. No re-centring of only
one inventory is allowed. The orchestrator may materialize these exact
transformed states and request the adapter's identity transform, avoiding
double transformation.

## Acceptance and energy accounting

- Absolute change in paired Ca-minus-La contrast under mesh refinement,
  tree refinement, or either rigid transformation: at most 0.5 kcal/mol.
- Difference between qm33 and qm36 global contrasts: at most 2 kcal/mol
  separately at **each** of primary, refined, and tree levels.
- Repeat and isolated-reduction/algebra consistency: at most 0.01 kcal/mol.
- Preserve all individual endpoint changes even where pairing cancels them.
- Require finite outputs, successful mesher/solver receipts, charge/state
  closure, and actual identical physical surface where required.

Use native CSV solvation energy (12-digit scientific output) in kJ/mol;
convert once by dividing by 4.184. Preserve native Coulomb and total/free
energies for auditing; **never add them** to the QM-plus-direct-coupling score.
Require solvation + Coulomb = native free energy within output precision.
For component accounting, RF(total)-RF(core)-RF(environment) uses the same
full surface and exposes the reaction-field cross term. It does not revive
the old isolated-core CPCM/PB transfer subtraction.

At primary level, also require
`[Coulomb(total)-Coulomb(core-only)-Coulomb(environment-only)]/4.184`
to agree with the independently evaluated full core/environment Coulomb
interaction within 0.01 kcal/mol. Environment-only output is shared within
each metal pair because its charges and physical surface are identical.

Before any solver result, the pinned installed `include/tabi/constants.h`
was checked: `UNITS_COEFF = 1389.3875744` kJ A/(mol e^2), equivalent to
**332.0716 kcal A/(mol e^2)**. Native `UNITS_PARA = 8729.779593448` is the
reaction-field prefactor (approximately `2*pi*UNITS_COEFF`). The new
descriptor's direct Coulomb and its independent accounting check both use
332.0716 for consistency with this backend. The earlier module's constant
332.063713299 remains unchanged: the absolute difference is 0.007886701,
approximately 0.002375 percent. This pre-result convention reconciliation
does not adjust a scientific acceptance tolerance or rewrite old results.

Zero/nonfinite results are not inferred from failed runs. Tests of archived
parsers and input preparation are distinguished from actual scientific
integration. No solver execution had occurred when this document was frozen.

## Primary implementation references

- [Pinned TABI surface construction](https://github.com/Treecodes/TABI-PB/blob/fe1c237b057418fed48535db125394607040d9de/src/particles.cpp)
- [Pinned native GMRES settings](https://github.com/Treecodes/TABI-PB/blob/fe1c237b057418fed48535db125394607040d9de/src/boundary_element.cpp)
- [Pinned native CSV/energy output](https://github.com/Treecodes/TABI-PB/blob/fe1c237b057418fed48535db125394607040d9de/src/output.cpp)
- [APBS TABI-PB documentation](https://apbs.readthedocs.io/en/latest/using/input/old/elec/tabi.html)
