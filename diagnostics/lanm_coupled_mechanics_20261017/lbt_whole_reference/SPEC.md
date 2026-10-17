# Whole-LBT stiffness reference: finite feasibility proposal

**Proposal only. No runner, preflight, molecular call or submission is authorized
by this document.** Root must review the actual ongoing Mex DFT cost and this
scientific scope before binding a first scout. The frozen Mex eight-cell
inventory is unchanged.

## Question and judgment

Does whole-peptide electronic structure support the large restoring forces and
deformation work introduced by the matched hybrid's MM scaffold?

Saved audit 13032 compares the same coordinates. At the primary 7CCO/La start,
Phe1 CA is wholly MM: whole-MACE composite force magnitude is 1.25991 versus
51.06301 kcal/mol/Å under the hybrid. The saved components identify MM as the
dominant contribution. Root's separate MM ledger review has not found an obvious
implementation error. This motivates an independent reference; it does not prove
that either force is physically correct.

A **whole-256-atom vacuum PBE0-D4 calculation plus the same fixed-charge ALPB**
is a defensible bounded test of this modeled deformation. It replaces the whole
internal energy, so the reference does not inherit the MM stiffness being tested.
At identical coordinates the same ALPB term cancels from reference-minus-model
force and energy differences. Report vacuum and composite quantities separately.

This is not a clean experimental measurement of an MM force constant. Replacing
the hybrid also changes electronic extent, polarization and the partition, so
agreement or disagreement cannot uniquely identify one MM parameter. Four
endpoints provide forces and finite deformation work, not a Hessian or a fitted
stiffness. PBE0 is an independent vacuum reference, not a claim of reproducing
the exact MACE training functional or spin representation.

## Fixed four-cell selection

Use only the 7CCO primary-region **La-proposed** deformation, regardless of its
energy sign. Both metals see the same before coordinates and the same after
coordinates. Do not substitute the Dy-optimized endpoint.

| Order | Cell | Fixed geometry | Role |
| --- | --- | --- | --- |
| 1 | `before_La` | Exact hybrid start `3fac44a7...`; also whole-MACE-selected `cartesian_01` | First cost/analytic-force scout |
| 2 | `after_La` | Actual 13017 La terminal `2db7d738...` | Same-metal deformation work |
| 3 | `before_Dy` | Same before positions, target substituted as Dy | Paired initial state |
| 4 | `after_Dy` | Same after positions, target substituted as Dy | Paired deformation work |

The before and after POINT/physical/coordinate pins are in `INVENTORY.json`.
The hybrid endpoint stopped at 80 accepted steps, with max force 3.31740 and
atom-vector RMS 0.901487 kcal/mol/Å. It is admitted but **not converged**. The
initial whole-MACE state was also not converged. Neither is an equilibrium
structure, and the proposal does not optimize or repair either.

The finite ceiling is four new whole-peptide `EnGrad` endpoints: one scout, then
at most three further cells only after actual native admission and measured
tractability. Preserve failed/unavailable slots. There is no automatic retry,
longer iteration limit, altered state, different basis/functional, additional
source, optimizer or library campaign.

## Physical and electronic state

Preserve 255 peptide atoms plus one target, in original physical order, with
target index 255. Preserve all represented hydrogens and the original covalent
graph. There are no represented waters, spectators, external point charges or
synthetic partition caps in this whole-system calculation. Do not send the
172-atom hybrid core XYZ, its metal index 0, or its shifted field to the runner.

The existing model declares holo charge −2, La physical multiplicity 1 and Dy
physical multiplicity 6. Preserve the declared free NH3+/COO− termini and five
deprotonated acidic side chains. The experimental terminal chemistry remains
unresolved in `LBT_TERMINAL_CHEMISTRY.md`; this calculation cannot resolve that
ambiguity. No proton transfer, neutralization or terminus replacement is proposed.

Use the established spin-free frozen-f reference convention: La ECP46 and Dy
ECP55, each with effective singlet bookkeeping; physical La1/Dy6 remain separate.
No spin-orbit treatment or stable electronic-state claim follows. The **actual
whole-system** all-electron/explicit-electron counts, alpha/beta populations,
even parity, both species' charge/state and basis/ECP coverage must be derived
and checked in a future zero-call Slurm preflight. They have not been calculated
or qualified by this inventory. Never copy the LBT48 count of 200 explicit
electrons or the 172-atom hybrid's prospective DFT state into the whole model.

## Recipe and solvent separation

Proposed native recipe is the pinned existing vacuum route:

`PBE0 D4 def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF EnGrad`,
with `Guess PModel`, the exact same metal lcecp1-TZVP and AutoAux-J bytes, and
ORCA 6.1.1. There is no CPCM, SMD, external field, MM energy, GFN embedding
correction, extra dispersion, geometry optimization or self-consistent solvent.

Define `E_ref = E_whole_PBE0_vacuum + E_saved_fixedFF_ALPB`. The solvent term
uses the same physical coordinates, unchanged ff19SB charges and formal +3
target, native metal identity and existing fixed-charge binary. Do not replace
its charges with DFT populations. Reuse a cell only after exact source/order,
metal, coordinates, charges and raw gradient provenance are admitted.

This vacuum dianion may have electronic-state and charge-localization limitations;
post-added fixed-charge solvent does not polarize or stabilize the vacuum SCF
density self-consistently. Preserve density/orbital residuals, small-gap and other
native warnings. Finite-basis SCF convergence does not establish electron binding
or a physical ground state. No silent CPCM change is a recovery path. These
limitations bound the stiffness interpretation even if the calculation finishes.

## Existing components and explicit gaps

- Both before-metal hybrid cells and their same-coordinate whole-MACE/ALPB
  evidence exist. The after-La hybrid endpoint and full ALPB energy/gradient exist.
- The same **after-La coordinates evaluated as Dy** are not yet supplied by the
  owned response results. Do not substitute the separate own-Dy terminal or reuse
  La ALPB as Dy. Bind a future admitted paired-pool cross cell if it appears.
  Otherwise a separately declared single fixed-charge ALPB evaluation is needed
  for `after_Dy` before a complete composite reference; this proposal launches none.
- Whole-MACE values at the new after geometry are not supplied by the start audit.
  Report only existing matched observations. Any additional whole-MACE calls need
  their own explicit finite declaration; they are not part of the four DFT cells.
- Existing LBT48 PBE0 results have different atom inventory, caps and charge.
  They cannot replace any whole-256 reference endpoint.

## Thin adapter requirements, before any execution

Reuse the native runner, allocation-aware renderer, raw `engrad` parser and
existing state/SCF/normal-termination admission. The current
`metal_environment_frozen_f.state/input_text` hardcodes charge −1, while
`LBT_LOCAL_REFERENCE_V2.inputs` enforces 48 atoms and 200 explicit electrons.
Those preparation functions cannot be called unchanged. A narrow whole-system
state/deck wrapper must use the actual full physical XYZ and charge −2; the
generic analytic-gradient parser can consume the newly derived state.

Future Slurm preflight must freeze all four cells, reproduce both actual physical
arrays and the target substitution, check native assets/atomic-order/electron
bookkeeping, and bind a single-cell scout manifest. It must record no external
field, no caps and no MM addition. Use 12 actual MPI ranks and memory derived
from the actual allocation through the existing renderer. Keep ordinary shared
Slurm limits; the scout measures cost rather than predicting it from atom count.

The first actual native result must demonstrate the requested charge, explicit
electron count, effective multiplicity, ECP and method, complete finite 256-atom
analytic gradients including ECP/dispersion, energy/geometry correspondence and
normal termination under existing native convergence semantics. Retain printed
residual flags. No new numerical derivative qualification is included in the
four-call budget; force comparisons remain descriptive unless independently
qualified. Root records admission and actual CPU/memory/time before scheduling
at most the remaining three endpoints.

## Predeclared comparisons

1. At `before_La`, compare the full Phe1 CA force vector, the full Phe1 residue
   force table and all-atom residual distributions against the two saved models.
   Retain the whole-peptide table so a local anecdote cannot conceal disagreement
   elsewhere. Phe1/CA selection comes from the saved audit, before DFT outcomes.
2. Repeat available same-geometry comparisons for the remaining declared cells.
   Report raw vacuum and composite forces, units/signs and SCF limitations.
3. Define `W_M = E_M(after) − E_M(before)` using native energy differences before
   conversion, and `delta_W = W_Dy − W_La`. Show reference and available hybrid
   electronic/MM/ALPB contributions separately. Retain the plan's 1 kcal/mol
   response resolution; it is not a statistical error bar.
4. In zero-call Slurm preparation, freeze the normalized full before-to-after
   displacement chord and the saved initial Phe1 CA hybrid-force direction.
   Report endpoint force projections and raw vectors without fitting a force
   constant or treating those directions as independently qualified derivatives.

Interpret only the same physical deformation and state. Favorable metal ordering,
extra energy lowering or force reduction is not an admission criterion. Strong
reference support for the restoring response would favor the hybrid's mechanics
over the relaxed whole-MACE surface on this conditional test; disagreement would
motivate locating the model error, not automatic parameter fitting. One source,
one deformation and the frozen-f/vacuum assumptions cannot establish an accurate
LanM preference profiler.
