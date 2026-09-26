# Real 1H4I embedding preparation: ready for a component-response scout

**A traceable 54-atom electronic core and two matched 9,087-charge protein fields
are ready. A complete additive QM/MM model is not ready.** No molecular energy,
gradient, optimizer or scheduler job was executed by this preparation task.

Input manifest: `workspaces/metal_environment_response_20260926/preparation/scout_v3/INPUTS.json`.
`INVENTORY.json` pins the exact inputs, inspected implementations and parameter
inventories. Root owns the energy ledger, reference execution and interpretation.

## Region, state and source

Reuse the exact archived `1h4i_qm36` pair from the density-embedding diagnostic:
metal + complete PQQ(3−), Glu177, Asn261 and Asp303 side chains, with three 1.09 Å
CB–CA link hydrogens. This is 54 atoms: 27 PQQ, 23 physical protein atoms, three
links and one metal. Ca charge −3 and La charge −2, both multiplicity 1. No water
is present in this conditional prepared source. Each endpoint retains its actual
coordinates, rather than reconstructing an older input from prose. These are
consumed development fixtures, not new validation observations.

The actual ORCA 6.1.1 archived La output documents Def2-ECP replacing 46 core
electrons and 280 explicit electrons at total charge −2 for this region. The
Ca output reports 290 explicit electrons at total charge −3. The formal core
charges follow metal + PQQ(−3) + Glu177(−1) + Asn261(0) + Asp303(−1):
La(+3) gives −2 and Ca(+2) gives −3. With the common −6 environment, total
conditional-system charges are −8 and −9, respectively. The smaller qm33 region
omits Asp303 and has the different −1/−2 core charges; those must not be copied
into this qm36 preparation. Root must qualify the
chosen current gradient Hamiltonian independently; earlier converged energies
are not a gradient qualification. Native runner and runtime-renderer paths are
in `INPUTS.json`, along with the executable SHA-256. Original r2SCAN-3c inputs use
DefGrid3, MBIS, `%pointcharges`, `DoEQ false`, and no continuum solvent in the
completed density diagnostic. The older balanced-embedding experiment used CPCM
and is a different Hamiltonian, not an interchangeable reference.

The protein environment is chain A of `workspaces/mxaf_qm/1H4I_protonated.pdb`.
Its selected PQQ and Ca source records are **heterogen chain B**, identified by
actual coordinate matches. The core map retains those identities; it does not
mistake their chain labels for an extra protein environment. All 24 PQQ heavy
atoms and the metal have unique source matches. Three PQQ hydrogens are
carver-local atoms. Other deposited protein chains are outside this explicitly
conditional monomer model.

## Boundary accounting already available

The archived charge map has 9,113 full protein atoms and ff19SB chain charge −8.
QM side-chain atoms and the adjoining MM1 C-alpha charges are removed from the
field. The residue-local correction is redistributed equally onto that residue's
retained backbone N and C. The full per-fragment ledger lists removed identities,
formal charges, recipients and increments; environment charge is −6 e. No charge
is added to a cap, no global neutralization is applied, and both metals consume
the same charge file. Boundary validity beyond this bookkeeping is still a
scientific qualification question.

Three link maps explicitly record retained CB and omitted CA identities,
coordinates, link length, and both physical-coordinate Jacobians. Analytic link
Jacobians agree with coordinate finite differences to 4.59e−10, below the
predefined 1e−8 mapping tolerance. Stored six-decimal core coordinates have
sub-microangstrom serialization differences from exact reconstructed link
positions; those differences remain recorded rather than silently modifying the
historical core.

**Do not splice in the newer scaffold inventory.** Its source has 9,060 protein
atoms and charge −9, while this older exact embedding source has 9,113 and −8.
The newer inventory provides useful general parent mechanics and source mapping,
but its charges/coordinates are not interchangeable with this scout.

## Deterministic environmental perturbation

Before any new energy was calculated, select the closest complete exterior
Ser/Thr/Tyr hydroxyl oxygen to any QM heavy atom, with lexicographic identity
breaking ties. This selects **Thr159 OG1**, 2.6791 Å from the core. Environment B
rotates only its actual HG1 by **+10°** around the CB→OG1 axis, preserving the
O–H distance and C–O–H angle. No core coordinate, atom, proton inventory, force-field
charge, donor identity or other environment coordinate changes. This is a modeled
orientation perturbation, not an equilibrium sample or experimental structure.

The hydrogen moves 0.189844 Å. Its closest QM-heavy distance changes from 2.1930
to 2.0139 Å; its closest nonbonded exterior heavy distance remains 2.1876 Å.
Both exceed the declared 1.2 Å H-heavy clash screen. The inherited source O–H
length is **1.18537 Å**, longer than a typical equilibrium hydroxyl. This scout
preserves the old source hydrogen convention exactly; it does not claim to
have repaired or equilibrated the protein's H geometry. No archived whole-protein
radial H normalization was applied. That limitation belongs beside any response
result, particularly if the magnitude depends strongly on this hydrogen.

A is byte-identical to the archived point-charge file. A and B have identical
charge and identity lists; only point-charge index 2412 changes coordinates.

## What is still missing

This preparation supports the explicitly declared **embedded electronic energy
component** and its coordinate derivatives. It is not a complete hybrid energy.
Protein-only ff19SB parameters exist in the newer parent inventory, including
bond, angle, torsion, CMAP and nonbonded exception terms. They do not supply a
validated metal/PQQ–exterior repulsion/dispersion model, a matched capped-local
subtraction, or an established treatment of all boundary mechanical terms.
Neither zero-filling those interactions nor assuming their metal differences
cancel is justified. Do not optimize the protein on this component or report its
contrast as affinity.

A chemically useful expanded-region qualification would promote a complete
hydrogen-bonding group such as the selected Thr159 side chain with its actual
source CB–CA boundary, retain the exact physical A/B states, and remove/redistribute
its MM charges through a declared corresponding map. It has **not** been prepared
or scored here. Under that changed partition the moving HG1 becomes a QM atom:
it is a physical-system boundary test, not another fixed-core A/B test. The
already prepared 154-atom context includes other second-shell residues but uses
the different newer source and cannot substitute without a new preparation.

## Reproduce preparation and checks

From the repository root, use a new output directory (writers refuse overwrite):

```bash
python diagnostics/metal_environment_response_20260926/preparation/PREPARE.py \
  --repository "$PWD" \
  --output "$PWD/workspaces/metal_environment_response_20260926/preparation/scout_reproduction"
python diagnostics/metal_environment_response_20260926/preparation/CHECK.py \
  --repository "$PWD" \
  --inputs workspaces/metal_environment_response_20260926/preparation/scout_v3/INPUTS.json \
  --output workspaces/metal_environment_response_20260926/preparation/check_reproduction.json
```

Five real-fixture checks pass: paired coordinates, unchanged inventory/charge,
complete source/synthetic mapping, physical cap Jacobians and hydrogen-contact
screen. These are parser/geometry/algebra checks, not scientific integration tests.
Two earlier preparation-only attempts remain visible: v1 stopped on a file-path
comparison for identical metal-specific boundary records; v2 retained incorrect
chain-A-only heterogen identity audit labels. v3 fixes those bookkeeping issues,
with the same intended coordinates/perturbation and no molecular calls.

Preparation was a few seconds of local CPU work; no formal allocation receipt
was produced, and process CPU/memory were not measured. Scientific compute cost
is exactly zero for this task. No jobs, production changes, push or email occurred.

### Documentation correction

The initial report prose mistakenly gave the smaller qm33 core charges. Corrected
to the actual qm36 La −2 / Ca −3 after checking both archived output electron
counts and the unchanged prepared `INPUTS.json`. No prepared inputs, point
charges, maps, queued calculations or archived results changed.
