# Saved Mex electronic-response comparison

Declared 16 October 2026 before calculating these projections. Source-only
preparation; root supplies the admitted static RESULT pin and submits one saved-data
Slurm job. No MACE, SCC, DFT, MM or ALPB evaluation, displaced geometry, optimization,
Hessian, new charge assignment or empirical fitting occurs.

## Question and exact comparison

At the identical repaired Mex8FNS EF3 origin, does the new electronic approximation
predict the same local donor and target response as the archived embedded PBE0
reference, including the La–Dy difference?

Compare `g_reference = g[PBE0-D4_embedded]` against
`g_approximation = g[MACE_vacuum] + g[GFN2_embedded − GFN2_vacuum]`.
Use the exact existing 209-atom capped region, 1840 shifted external sites and
2049 physical atoms. Local charge is 0; the three spectator ions remain Nd(III).
Map real atoms, both cap-anchor derivatives and all external gradients onto the
same physical source-ID order. No unlike-metal total-energy offset is interpreted.

Source contract:
`workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/Mex8FNS/INPUTS.json`,
SHA256 `4b3b98983af0cb5fcdbf7b8f52e35063308201ce9df777a30484acf5146b7930`.
Archived collection:
`workspaces/metal_environment_response_20260926/lady_repaired_origins_Mex8FNS_v1/FINAL_COLLECTION.json`,
rows `La_A` and `Dy_A`. Both raw `endpoint.engrad` and
`endpoint.runtime.pcgrad` files and execution receipts are available. They have
the same XYZ and field hashes as preparation 12569 and local-field result 12570.
Require admitted static La/Dy rows before any projection; the earlier failed
static replay is not an eligible replacement.

## Six directions fixed independently of force outcomes

Use this order, without ranking, replacement or additional directions:

1. Target translation along source Cartesian +x.
2. Target translation along source Cartesian +y.
3. Target translation along source Cartesian +z.
4. Asp84 rotation about `A/84//CA → A/84//CB`.
5. Asp86 rotation about `A/86//CA → A/86//CB`.
6. Glu95 rotation about `A/95//CA → A/95//CB`.

For each donor, cut that actual native graph bond and take the complete component
on the CB side, including hydrogens and the intact carboxylate. Require a graph
bridge, exactly the expected sidechain atoms and residue identity. The positive
tangent is `unit(CB−CA) × (R_atom−R_CB)` in Å/radian; the target columns are unit
translations. These origin tangents are not finite displacements or candidates.
No synthetic cap receives its own degree of freedom.

Normalize each column using its Euclidean norm over **all physical atoms and
Cartesian components**, including H. This gives a unit physical-displacement
direction and comparable loads in kcal/mol/Å. Save the native Jacobians, column
norms, normalized Gram matrix, rank and conditioning. Require all six nonzero
columns and rank six at relative eigenvalue cutoff 1e−10. The same directions
and normalization must apply to both metals and methods; no heavy-only metric,
source alignment, force-selected basis or direction rescaling is introduced.

## Admission and outputs

Verify the comparison implementation/plan/dependency pins; the static RESULT,
config, source implementation and both row receipts; the exact archived collection,
manifest, raw inputs/gradients/outputs and successful execution receipts; and the
shared source graph, maps and field. Reconstruct native PBE0 physical gradients
with the existing `metal_environment_force_assembly.assemble/scatter` helpers.
Replay the static MACE mapping from its saved negative forces and the GFN2 field
difference from its saved two physical gradients. Compare saved component arrays
with absolute tolerance 1e−10 plus relative tolerance 1e−12, rather than requiring
bitwise identity after floating-point matrix multiplication. Require exact input
hashes, state metadata and raw file hashes. Report replay residuals separately;
these tolerances do not set a physical agreement threshold.

For each component/method report the signed energy-gradient projection `g·v` and
the opposite force load `−g·v`, separately for La and Dy. Also report common
`(La+Dy)/2` and differential `Dy−La` projections. Include approximation minus
reference, the separately saved MM and global fixed-FF-charge ALPB loads, and
each electronic method plus those identical shared terms. Their addition exposes
the represented common mechanical/solvent load; it cannot change the electronic
method discrepancy. Do not sum the six locally normalized loads into an energy.

Write RESULT.json, a short Markdown table and a pinned Jacobian artifact. Missing
or incompatible inputs make the comparison unavailable, with no partial-method
substitution. Report signs and continuous discrepancies without inventing an
accuracy pass threshold, fitted scaling, favorable-site selection or biological
classification. Zero new molecular calls and actual analysis cost are explicit.

## What the reference can and cannot establish

The archived ORCA 6.1.1 calculation is PBE0-D4 with light-atom def2-TZVP,
explicit lcecp1-TZ metal bases and AuxJ, RIJCOSX, DefGrid3, VeryTightSCF,
EnGrad and DoEQ false. La uses ECP46; Dy uses ECP55. Both have 832 explicit
electrons and effective restricted-singlet valence. Physical La(III) singlet and
Dy(III) 4f9/sextet metadata are distinct from that spin-free frozen-f model; no
spin–orbit treatment or stable-state proof is present. Both archived
`physical_state_qualified` and `force_consistency_qualified` flags remain false,
and their early small-gap warnings remain preserved.

MACE receives physical multiplicities 1/6. GFN2 uses its own spin-independent,
f-in-core effective singlet parameterization. The methods therefore differ in
their electronic representation as well as their approximation to the energy.
Agreement is independent local response evidence conditional on these models;
disagreement is a useful warning, not a unique attribution to MACE or the GFN2
field correction. No matching isolated PBE0 counterpart is available in this pair.
Neither result proves affinity, source robustness, a stationary structure, solvent
accuracy, full-gradient qualification or transferable within-series selectivity.
This saved comparison may run before the separately declared finite derivative
campaign; it does not replace that consistency check or authorize new DFT.

## Reproducible interface

`python MEX_REFERENCE_RESPONSE.py --config CONFIG.json --output FRESH_DIRECTORY`

Config keys: `repository`, `implementation`, `plan`, `reference_collection`,
`static_result`, and `dependencies` mapping the nine named script files listed
in the adapter to their exact path/SHA256 pins. Do not populate `static_result`
until root supplies the admitted receipt. Runtime needs Python and NumPy only;
imported project helpers parse and map saved data and instantiate no energy engine.
One shared CPU, no GPU and no memory/time override. No new test framework.
