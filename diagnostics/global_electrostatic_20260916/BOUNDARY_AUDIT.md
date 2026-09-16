# Stage 1 boundary and energy-accounting audit

2026-09-16. Scope: the approved 1H4I qm33/qm36 partition pair, using archived
preparations unchanged. This audit inspected files and ran nine real-artifact
integrity/arithmetic tests. It did not prepare structures, launch electronic or
continuum calculations, or inspect new pilot energies.

**Finding:** no source, charge-ownership, or cap-containment defect blocks the
four approved gas-phase endpoints. This is not a physical validation of the
boundary approximation. New gas-phase charge quality, identical actual surface
meshes, solver convergence, and the partition gate remain required.

## Frozen physical system

All four skeletons refer to the same protonated source:
`workspaces/mxaf_qm/1H4I_protonated.pdb`, SHA256
`e91e6d51e0175a5f92e7eb39a446d65f3e869db7d475b643fc8a5eb358127bab`.

The declared assembly is **deposited catalytic chain A**, not the complete
biological assembly. Deposited chains B, C, D, E, and F are explicitly excluded.
The physical cavity contains 9,141 atoms: 9,113 protein atoms, 27 atoms of the
frozen oxidized PQQ(3−) state, and one metal. There are no explicit waters.
The source protein is parameterized with the existing ff19SB file, SHA256
`086bfdb1c05e7d5cb1d330e6c39fab5056a95e89492cdbd3044c746db13e9d09`.
Its net charge is −8 e to floating-point precision.

The physical cavity has no synthetic cap spheres. Radii remain H 1.20, C 1.70,
N 1.55, O 1.52, S 1.80, and common Ca/La 1.80 Å. The metal has common cavity ID
`metal` and element label `M`. Coordinates/radii at the existing PQR precision
(10/6 decimal places) are identical across all four states. Their common
`physical_boundary_key` is
`4c518676f0b57cd3d30b135734e226e49b98627f365908f668cd7e27c7ecb1c4`.

Five Asp303 coordinates differ between the raw JSON representations by at most
about 7.1e−15 Å because one route converted OpenMM coordinates. This is binary
roundoff, not a geometry change. Canonical serialization/hash must define mesh
identity; raw Python dictionary equality is inappropriate across partitions.

## Quantum/environment ownership

The historical names qm33/qm36 are not their total atom counts.

| Partition | QM atoms | Caps | Environment atoms | QM La/Ca charge | Environment charge | Whole La/Ca charge |
|---|---:|---:|---:|---|---:|---|
| qm33 | 47 | 2 | 9,094 | −1 / −2 | −7 | −8 / −9 |
| qm36 | 54 | 3 | 9,087 | −2 / −3 | −6 | −8 / −9 |

qm33 includes metal, PQQ, Glu177 sidechain, and Asn261 sidechain. qm36 additionally
includes the complete capped Asp303 sidechain. PQQ and the metal carry no MM
charges. Source QM atoms and environment charge atoms have disjoint IDs.
Paired La/Ca endpoints have exactly identical environment atoms/charges,
coordinates, source, assembly, microstate, water list, and cavity radii.

For each selected sidechain, all its represented source atoms **and its CA MM
charge** are removed from the environment. CA remains in the physical cavity.
Local charge closure distributes a fixed correction equally to the source
backbone N and C; it does not alter other residues or globally neutralize the
protein. The original ff19SB residue templates reproduce these ledgers:

| Fragment | QM formal charge | Increment to each backbone N/C (e) | Final N / C charges (e) |
|---|---:|---:|---|
| Glu177 | −1 | +0.07875 | −0.43755 / +0.61535 |
| Asn261 | 0 | +0.00480 | −0.41090 / +0.60210 |
| Asp303, qm36 only | −1 | +0.09000 | −0.42630 / +0.62660 |

The complete cross-partition environment change is removal of Asp303
CA/CB/CG/HB2/HB3/OD1/OD2 charges and +0.09 e on each Asp303 N/C. The QM addition
is CB/CG/HB2/HB3/OD1/OD2 plus its CB–CA link H. No other residue, donor, water,
protonation state, or source geometry changes. This is a charge-preserving
sidechain partition test with an explicitly approximate covalent boundary;
it does not move an entire uncapped peptide residue into QM.

## Caps and near contacts

Every link H lies on its recorded source CB-to-CA direction at 1.09 Å from CB.
The following numbers come directly from the frozen coordinates. Containment
margin is `r(CA) − distance(cap,CA) − r(cap)`; positive means the whole cap sphere
lies inside an existing physical CA sphere.

| Cap | Distance to cavity-only CA (Å) | Containment margin (Å) | Nearest charged environment atom | Distance (Å) |
|---|---:|---:|---|---:|
| Glu177 | 0.443888120 | 0.056111880 | Glu177 HA, +0.1105 e | 1.437421058 |
| Asn261 | 0.445570283 | 0.054429717 | Asn261 HA, +0.1048 e | 1.453939387 |
| Asp303 | 0.442435203 | 0.057564797 | Asp303 HA, +0.0880 e | 1.366865264 |

No QM source or cap charge is within 1 Å of an environment charge. That rule
does **not** mean caps are far from all physical nuclei: they intentionally
sit about 0.44 Å from the cavity-only CA whose MM charge was removed. The
synthetic link H is not a native protein atom, and charge closure does not
validate its multipole field. These short cap–HA contacts are retained in the
declared complete direct Coulomb sum; no new force-field 1–2/1–3 exclusions or
screening factors are introduced.

The TABI adapter must generate the surface from the common source-only cavity
while retaining the link-H charge locations in the charge list. The backend
owner's planned adapter replaces only NanoShaper's XYZR geometry input and
preserves the actual vertices/faces. Actual mesh equality across endpoints and
partitions must be checked after execution; containment alone is not proof of
identical numerical meshing.

## Accounting and physical checks

The new formula is `E_QM,vac + C(core,environment) + G_RF(full)`. The endpoint
energy function takes Hartree, kcal/mol, and TABI reaction energy in kJ/mol,
respectively. Convert Hartree by the recorded 627.509474 factor and kJ/mol by
division by 4.184 once. For contrasts, subtract native endpoint quantities
before conversion. Report each term and Ca-minus-La sign independently.

Before solver execution, inspection of the pinned TABI header established its
Coulomb coefficient as 1389.3875744 kJ Å/(mol e²), or 332.0716 kcal Å/(mol e²).
The new descriptor and its independent Coulomb accounting check use that same
coefficient. The older 332.063713299 coefficient and all old results remain
unchanged. The numerical protocol records this approximately 2.375e−5 relative
convention difference; no acceptance threshold was adjusted after results.

The new total must contain **neither** the old isolated-core RF subtraction
**nor** an environment RF subtraction. An environment-only RF component can
be retained for auditing. It is constant between metals within a partition;
the environment differs across partitions, so its value cannot be reused
between qm33 and qm36. Any component solves must use the same full physical
surface. TABI's separate Coulomb output is not an extra term in this model.

Checks still needed on the approved solver schedule:

1. Verify each new output is gas-phase native r2SCAN-3c with the declared
   charge/multiplicity and La ECP; obtain its own MBIS charges and ESP receipt.
   Archived CPCM charges/ESP passes are not transferable. In particular, qm36
   Ca is a −3 gas-phase core; convergence alone will not validate its density.
2. Require charge sums, exact source mapping, full cavity/mesh identity, and
   unchanged paired environment ownership before evaluating a contrast.
3. Retain total, core-only, and environment-only reaction contributions on the
   common surface if used for the planned component audit; the full RF must
   include the cross response. Do not introduce separate core-only cavities
   into the production expression.
4. In an environment-removal check, `Q=0` makes direct core–environment Coulomb
   zero but generally leaves nonzero core RF inside the retained protein
   cavity. This new descriptor has no transfer correction that must vanish.
   A uniform dielectric limit, if included in the frozen solver manifest,
   has zero RF. Do not conflate these two checks.
5. Apply the approved 0.01 kcal/mol reduction/algebra, 0.5 kcal/mol paired
   numerical/rigid-transform, and 2 kcal/mol partition tolerances. Partition
   tolerance is required at primary, refined, and tree settings individually,
   as frozen before execution. The component
   ledger can locate changes; it cannot by itself assign a unique physical
   cause to a failure. No parameter rescue follows a failed gate.
6. Require new-protocol/cache identities and explicit missing status for any
   absent term. No old aquo offset, universal zero, or released PQQ band applies.

## Immutable fixtures and tests

Skeletons are under
`workspaces/affordable_challenger_20260915/environment/<case>/skeleton.json`:

| Case | SHA256 |
|---|---|
| 1h4i_qm33_La | `b6597141c4143f854b23f2d9cfdcfce0e4c973cee08b597c81550d72c5fa9c8a` |
| 1h4i_qm33_Ca | `25f7b78bca9e000b263ce2514f6999cb248f79912f1e0b7636507e3a288df952` |
| 1h4i_qm36_La | `12e7496aece17805dc0173f714d4f9b1beb1aee535632551b3591331585d6c4b` |
| 1h4i_qm36_Ca | `f2f2d39cc9109383814a399678412ca15687f8926e983934a256b0b98ebcd971` |

Each test verifies the source, boundary-map, XYZ, and force-field hashes where
used. The source pilot manifest SHA256 is
`91ab53de327806f3e042791957f72c34ace68a9368eef0f99601efdfa4dc987a`;
the archived solver-result arithmetic fixture SHA256 is
`542d491159111114e54ba81df220b08f7bbf7d8db619e2a32b970dea3b0ee885`.

```bash
python -m unittest discover -s tests -p 'test_global_electrostatic_accounting.py' -v
```

Result: **9 tests passed**, 5.092 s unittest elapsed. Six verify actual state,
ownership, local closure, partition change, and cap containment; three exercise
the new pure arithmetic API using real archived numerical values and explicitly
corrupted missing/nonfinite copies. The old CPCM/APBS values are used only as
arithmetic inputs and are not reported as a new gas-QM/TABI scientific score.
No unavailable solver integration was replaced with a fabricated result.
