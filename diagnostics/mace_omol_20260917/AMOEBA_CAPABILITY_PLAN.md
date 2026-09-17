# Next contained step: explicit protein polarization with the existing MACE context

Declared after the responsive-field result. This is a new engineering and
physical-accounting investigation, not a backend switch inside that failed
candidate. Active-goal autonomy applies; baseline/default remains unchanged.

## Why this direction

The fitted distributions now pass actual near-boundary coupling checks, and
core electronic response plus a consistent solvent update changes the ordering
from1/4to3/4. The remaining4.44kcal boundary sensitivity still disqualifies that
candidate. Neither observation identifies its unique cause. Permanent protein
charges and the one-way reaction field remain substantive approximations.

AMOEBA was explicitly raised by Jacob and is already listed in the original
MACE plan. It is not a new invention here. The read-only installed inventory
workspaces/mace_omol_20260917/amoeba_installed_inventory_v1.json confirms
OpenMM8.5.1, AMOEBA2018/GK files, and multipole/GK APIs. It does not prove that
this hybrid can be implemented coherently, that La parameters are appropriate,
or that a model will improve accuracy. The original plan references
https://doi.org/10.1021/acs.jpcb.2c07237; verify the actual primary paper and
parameter definition before relying on it.

## Question and exact first-stage scope

Can the installed maintained solver represent the real protein's permanent
multipoles and induced electronic response around a QM charge distribution,
with an auditable core subtraction and compatible reaction field, at the
existing affordable inference scale?

Read installed source/API documentation and the primary parameter evidence.
Use precisely the3existing full normalized physical systems: GGR1GLG,
alpha1F6S andalpha6IP9, covering the4QM representations already tested.
Perform at most3real protein topology/parameterization preparations, one per
full structure, using AMOEBA2018. This stage evaluates **no energies or forces**,
performs no minimization, and makes no scientific score or new biological
comparison. No newDFT, charge-fit, native-potential, MACE or GPU calls.
No installation over the working environments. Record actual preparation cost.

Keep all source heavy/H coordinates, protonation, disulfides, water inventories,
assembly and physical IDs. A documented atom-name mapping is permitted when
identity and connectivity are unchanged. Do not add/delete atoms, generate
missing geometry, alter protonation, invent a La template or omit a cofactor
because parameter matching fails. Report unsupported chemistry explicitly.
The quantum ion/selected sites may need an external-source representation,
not a classical metal force-field template; determine this from the actual
Hamiltonian rather than assigning defaults. Published La parameters, if used
for damping/response, require their exact definitions and units.

## Required deliverable before any new energy pilot

A compact capability result and a concrete additive/subtractive expression:

-Which native force components contain permanent Coulomb, induction, and GK
 reaction-field energies; whether induction and solvent response are coupled,
 and which self terms are included. Verify the installed implementation.
-How a frozen QM-derived distribution enters the field without retaining a
 second force-field charge/multipole on the same site. Preserve the existing
 cap/source mapping and distinguish actual source atoms from synthetic caps.
-Which classical core interactions are subtracted because DFT already contains
 them. A generic ONIOM expression is insufficient: specify actual source sets,
 covalent exclusion/scaling groups, damping and core/environment cross terms.
-How protein-only permanent multipoles, induced dipoles and environment-only
 constants cancel or remain. Zero QM polarizability alone may not define the
 correct damping or coupling; verify independently.
-A single traceable covalent-boundary rule with charge closure and no global
 neutralization or score-dependent redistribution. Axis atoms used to orient
 multipoles remain physical even if their charge representation changes.
-A common physical full boundary for GGR extended/connected. Replacing OBC2
 by GK changes the model and cannot inherit an old absolute reference/band.
-What electronic response is already present in the selected QM density and
 MACE short readout, and what would be counted twice. The trained local readout
 is not a uniquely separable non-electrostatic energy.
-Which charge, electric-field and partition checks can actually qualify the
 representation, and the expected solver/cost scale from documented algorithms
 and existing receipts. Potential-only agreement is not an induced-field test.

Keep the current MACE full/core short terms available as matched context, but
do not append a new environmental scalar until this accounting is coherent.
Direct AMOEBA scoring or relaxation is not automatically part of this stage.
No fitted labels/weights, flexible classifier, assumed curvature or entropy.

If a coherent affordable route exists, declare its exact parameters, finite
energy/force task inventory, independent numerical tolerances and comparison
rules BEFORE inspecting its outputs, then implement/run it autonomously under
the standing goal authorization. If it does not, record the specific failed
capability and pursue another defensible approach. No invented successful
scalar or automatic replacement of the baseline. The goal remains broader
than this backend and cannot be marked complete by an engineering pass.
