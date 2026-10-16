# Mex physical directional derivative check

Declared 16 October 2026 from source mappings and implementation, before reading
the new assembled force values or any displaced energies. This is a design only.
Root owns execution; implementation waits for the admitted `MEX_STATIC_PLAN.md`
result and its frozen source. No molecular work or geometry calculation was run
to choose these directions.

## Question and fixed model

Does the assembled analytic gradient differentiate the actual whole-system energy
when an external charge site moves, and when an omitted physical cap anchor moves?

`E = MACE_vac(Q) + GFN2_emb(Q;MM) - GFN2_vac(Q) + MM_ledger + ALPB(qFF,R)`.

Use exactly the static model: 2049 physical atoms, 209 local atoms including four
synthetic caps, 1840 external charges, local charge 0, physical La multiplicity 1
or Dy multiplicity 6, GFN2 effective singlet/UHF 0, three fixed Nd(III) spectators,
and the original global FF-charge vector of total +2. No charge redistribution,
occupancy change, new solvent recipe, additional classical Coulomb or fitted term.

## Two source-selected directions

Indices below are zero-based checks on the actual 12569 mapping; source IDs are
authoritative. Both metals use the same physical source coordinates and directions.
Normalize the following vectors once at the original source in the Slurm fixture
preparation. Do not recompute a direction at a displaced endpoint.

| Fixture | Moving physical atom | Unit direction and reason |
|---|---|---|
| `external_Leu97_O` | `A/97//O`, physical 1019, field row 832, charge −0.5679 e | `(R[O97]−R[C97])/norm(R[O97]−R[C97])`. A real nonboundary MM carbonyl oxygen; no local QM atom or cap moves. |
| `boundary_Leu83_CA` | `A/83//CA`, physical 801, omitted anchor of cap QM index 205 | Let `u=unit(R[CA83]−R[C83])`, `w=R[N83]−R[CA83]`; use `unit(w−u*(u·w))`. Transverse motion changes the cap direction, avoiding a radial anchor displacement whose fixed-length cap would remain stationary. |

For the boundary fixture, retained `A/83//C` is physical 803; frame atom
`A/83//N` is physical 799. Cap 205 has length 1.09 Å. The omitted CA is present
in the physical/MM/cavity inventory but absent from the external-charge list.
Require these identities and graph bonds, a nonzero projected frame vector
(norm greater than 1e−8 Å), and the declared cap relationship. Any mismatch makes
the fixture unavailable; no alternative atom or direction is selected.

For each direction and metal use only `h = ±0.0001, ±0.00005 Å`:
`R(h)=R0+h*v`, where `v` is nonzero on exactly the named physical atom.
There are four directional tests and sixteen displaced combined cells. These
single-atom perturbations are derivative fixtures, not proposed finite protein
motions: they deliberately permit the corresponding tiny bond/angle changes.
All other physical coordinates, including waters, target and spectators, remain
exactly equal to the source. Preserve the full denominator if any cell fails.

## One coherent coordinate update

Start every endpoint from the original physical source, never from another
perturbation. Update the physical/export coordinate records by source ID. Update
every external-site coordinate from that same physical array while keeping the
existing shifted charges and charge-shift membership fixed. Do not rerun charge
assignment or constrain the old dipole diagnostics as the geometry changes.

Rebuild every synthetic cap from its physical anchors, using the existing rule
`H=Rr+l*(Ro−Rr)/|Ro−Rr|`. Recompute both cap Jacobians:
`Jo=l/|Ro−Rr|*(I−u*u^T)` and `Jr=I−Jo`. Source-derived real QM coordinates follow
their physical IDs. Caps are neither independent particles nor solvent sites.
Update cap coordinates, anchors, Jacobians, boundary-coordinate bookkeeping,
core XYZ, field file and mappings together. Serialize coordinates with 17 digits
and verify round trips and exact unchanged-slot identity.

For `external_Leu97_O`, only field row 832 moves; the complete local XYZ and all
caps must remain byte-identical to the appropriate origin. Reuse that metal's
origin MACE and isolated GFN2 energies. Recalculate embedded GFN2, the full MM
ledger and global ALPB. For `boundary_Leu83_CA`, field coordinates/charges must
remain byte-identical; cap 205 changes and the other three caps do not. Recalculate
MACE, embedded and isolated GFN2, MM and ALPB. Do not reuse a local component merely
because all *real* QM atoms stayed fixed: its cap geometry changed.

Evaluate the same full classical ledger in its native particle order, joined by
source ID with explicit target-alias handling. Reuse its topology, parameters,
exceptions and retained-term policy. Set all 2049 trial coordinates, including
zero-field anchors. Preserve retained radial cut bonds and existing omissions of
cap-represented angular terms. Record all four ledger energy groups separately.
Use the existing Reference OpenMM platform; no constraints, integration or added
force-field energy beyond this ledger.

ALPB receives the entire physical 2049-atom trial geometry, the original qFF in
the same source-ID order, and the actual target element. Hold every qFF entry
fixed, including the boundary atom's original FF charge; external charge shifting
does not alter this vector. No cap enters the cavity. Reuse the static native
wrapper/library/settings and charge-closure gate.

## Provenance, reuse and finite cost

Bind the successful static `RESULT`, config, implementation, both assembled
origin rows and all raw component receipts. Re-admit the 12570 local GFN2 rows
with their existing parser; preserve accuracy 0.0001, maximum 500 iterations,
300 K electronic temperature, author parameters, runtime and effective states.
Reuse the exact static MACE model/head/float64/runtime/physical-state checks and
native ALPB build. No origin molecular calls are needed.

Required source anchors, relative to the repository:

- `workspaces/lanm_coupled_response_20261016/mex_local_field_prepare_v1/execution/MAPPING_{La,Dy}.json` (12569).
- `workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/Mex8FNS/INPUTS.json`, SHA256 `4b3b98983af0cb5fcdbf7b8f52e35063308201ce9df777a30484acf5146b7930`.
- `workspaces/metal_environment_response_20260926/scaffold_origin_transfer_v1/Mex8FNS/LEDGER.json`, SHA256 `f37fb9d6185b987acdb4b17d824cef114fa1acbad6bc7e34240fd4588c53eb9e`, and its exact `{La,Dy}_A.json`, artifacts and native XML pins.
- The original export `atoms_charges.json` and `bonds.json` reached through the
  source INPUTS; use the *repaired physical coordinates* rather than replacing
  them with the unmodified export's coordinate array.

Reuse the coordinate/cap/field update pattern in `PREPARE.py::write_candidate`,
`metal_environment_force_assembly.prepare_mapping/scatter`,
`metal_environment_component_checks.build_system`, `MEX_LOCAL_FIELD` execution
and admission helpers, and the frozen `MEX_STATIC` MACE/ALPB adapters. The Hans
source-specific mode manifest and its 1665-protein-atom assertions do not transfer.
No new molecular engine or full Cartesian numerical gradient is required.

| New calls, maximum | External fixture | Boundary fixture | Total |
|---|---:|---:|---:|
| MACE local energy | 0 | 8 | 8 |
| Embedded local GFN2 | 8 | 8 | 16 |
| Isolated local GFN2 | 0 | 8 | 8 |
| Whole MM ledger configurations | 8 | 8 | 16 |
| Global fixed-charge ALPB | 8 | 8 | 16 |

Each MM configuration exposes four component energy queries (64 total).
Displaced energies suffice; native analytic outputs produced by the existing
GFN2/ALPB admission route may be saved but do not constitute additional finite
directions. No full-protein SCC, DFT, Hessian, restart or repeated origin call.
Cache keys must include exact component geometry/field, metal/state, parameters,
runtime and implementation; local unchanged reuse requires equality checks.
Record actual attempts, failures, wall/CPU/GPU costs and reuse separately.

## Derivative definition and frozen gates

At the origin use `a = g(R0)·v`, with `g=dE/dR` in kcal/mol/Å. Native MACE force
has the opposite sign. For each component and the total compute
`D(h)=[E(R0+h*v)−E(R0−h*v)]/(2*h)` at both positive step magnitudes. Difference
each component's endpoint energies before summation to avoid unnecessary loss of
precision. The field component is embedded minus isolated, not embedded alone.
Convert eV/Hartree/bohr once using the static constants. Preserve individual
embedded and isolated diagnostics, plus common `(La+Dy)/2` and `Dy−La` results.

1. All source/state/coordinate/component availability and native convergence
   gates from the static/local routes must pass. No fallback or silent repair.
2. At **both** step sizes require `|D(h)−a| ≤ 0.2 + 0.05*|a|` kcal/mol/Å for
   each metal's total and each of the four physical components (MACE, field
   difference, MM ledger and ALPB), preventing cancellation from hiding a failed
   component. Individual embedded/isolated residuals are also reported.
3. Require step agreement `|D(0.0001)−D(0.00005)| ≤ 0.2 + 0.05*|a|` for the
   same quantities. This is a fixed numerical consistency gate, not a fitted
   energy resolution. Apply the same additive gate to the derived common and
   Dy−La total projections; retain them even if the individual-metal tests pass.
4. For each classical group additionally use the existing component-check
   tolerance `0.001 + 0.0001*|a_group|` kcal/mol/Å at both steps. Unchanged local
   components in the external fixture must have exactly zero finite differences
   through explicit cache reuse, with zero physical projection within 1e−10
   kcal/mol/Å. Report the known native pcgrad output-rounding bound separately.

Admission requires all four directional tests and their sixteen cells. A failure
closes this finite check with its actual residual or unavailability; no smaller
step, changed tolerance, selected direction or chemistry retry follows. A pass
supports these two assembly interfaces locally, not every protein mode, electronic
state stability, solvent accuracy or binding preference. Old failed checks remain
preserved.

## Immediate bounded response pilot after a pass

Next prepare one source-defined seven-coordinate Mex basis: target xyz, whole
Asp84 CA–CB, Asp86 CA–CB and Glu95 CA–CB donor rotations, and the Leu83 N–CA
backbone hinge with its complete covalently connected downstream component.
Confirm these graph bridges and protected groups from this Mex export before
any call; no Hans axes or independently moving fragments. Reuse exact kinematics,
the full ledger/cavity and original-source contact/chirality/peptide gates,
including the current moving target against all physical atoms. From the same
origin make one projected own-total-gradient proposal per metal, with normalized
physical-coordinate Gram pseudoinverse cutoff 1e−10, maximum displacement of any
physical atom 0.02 Å and maximum torsion 0.5 degrees; no retries or difference
optimization. Cross-score both admitted proposals for both metals, retaining the
origin in the common three-geometry menu (four new combined energy cells at most).
Report predicted versus actual own work, component changes, common strain and
conditional Dy−La accommodation at fixed three-Nd occupancy. Failure of geometry,
execution or own descent stops expansion; successful finite motion is not a
stationary minimum or affinity claim. This pilot needs its own frozen execution
contract after the present check, rather than an open-ended optimization campaign.
