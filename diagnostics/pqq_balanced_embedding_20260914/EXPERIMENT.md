# Preregistered PQQ balanced electrostatic-embedding experiment

**Status at freeze:** design only; no embedded electronic energy has been
calculated.

**Question:** Does a charge-conserving whole-chain protein field remove the
generic La stabilization caused by moving an Asp(-) across the QM boundary,
while retaining a reproducible XoxF-specific response?

This experiment follows the completed unembedded 3.3/3.6-A diagnostic in
`diagnostics/pqq_boundary_pair_20260914`. That experiment admitted exactly one
seven-atom Asp(-) fragment at 3.6 A in every structure and found:

| structure | score at 3.3 | score at 3.6 | bare boundary response R |
|---|---:|---:|---:|
| MxaF 1H4I | +6.176 | +8.380 | +2.204 |
| C5B120 monomer | +11.193 | +15.369 | +4.176 |
| C5B120 dimer | +12.190 | +17.167 | +4.977 |

The bare selective margins over MxaF were +1.972 and +2.772 kcal/mol. Thus
3.6 A contained both a generic anionic-boundary effect and a possible modest,
geometry-dependent XoxF component. MxaF remained falsely La-positive at
3.3 A, so the radius alone did not create the broader PQQ score bias.

## Frozen panel

The experiment is a 3 x 2 x 2 matrix: **12 new ORCA single points**.

- Structures: MxaF 1H4I, C5B120 monomer, C5B120 dimer.
- QM inclusion radii: 3.3 and 3.6 A.
- Metals: La and Ca at identical nuclear coordinates within each pair.
- No apo or aquo calculations: both cancel from every primary paired endpoint.
- No geometry relaxation, alternative PQQ microstate, basis change, cavity
  tuning, or new fold is permitted in this experiment.

Frozen source structures:

| structure | donor-bearing chain | protonated source SHA-256 |
|---|---|---|
| MxaF 1H4I | A | `e91e6d51e0175a5f92e7eb39a446d65f3e869db7d475b643fc8a5eb358127bab` |
| C5B120 monomer | A | `1db4fe7872530e97fb3027ce23036096fc6fa5de49e6c9237def991397cc5531` |
| C5B120 dimer conformation | B | `579febf0b4fb1aab9198916b55336258c74dccd3d0b01ee54aade650ac8aa170` |

"Whole protein" is operationally the complete donor-bearing approximately
600-residue catalytic chain, not every crystallographic or predicted assembly
partner. This gives a matched physical object across 1H4I, the monomer fold,
and the dimer-derived conformation. It also avoids introducing the second
PQQ/metal site of the MxaF alpha2beta2 assembly or the second C5B120 site as an
unparameterized, unmatched environmental charge. Assembly electrostatics are
outside this causal test.

## Electronic model

The QM coordinates, atoms, formal charges, PQQ(3-) microstate, link hydrogens,
r2SCAN-3c Hamiltonian, native metal basis/ECP policy, DefGrid3, CPCM(Water),
and ORCA 6.1.1 settings remain exactly those of the completed boundary panel.
The only new term is a frozen protein electrostatic field supplied through
ORCA `%pointcharges`.

Protein charges will be assigned on the already protonated donor-bearing chain
with AMBER ff19SB through the installed OpenMM force-field implementation.
Coordinates and protonation will not be changed. PQQ and metal are QM species
and are never parameterized as MM point charges.

At a fixed structure and radius, La and Ca must consume a byte-identical point-
charge file. Point charges do not change the QM electron count.

## Charge-conserving QM/MM boundary

For every protein fragment in the QM carve:

1. Remove point charges for every source atom represented in the QM geometry.
2. Remove the adjacent MM1 C-alpha charge. Leaving it in place would put a
   point charge about 0.44 A from the synthetic C-beta link H.
3. Let `q_residue` be the original full force-field charge of that residue and
   `Q_fragment` its declared QM formal charge. Redistribute the required local
   correction over the retained heavy atoms N and C directly bonded to the
   omitted C-alpha, equally, so that the retained residue point charges sum to
   exactly `q_residue - Q_fragment`.
4. Record every removed atom, original charge, recipient, and increment.

This nearest-bonded-heavy-atom map is the frozen primary boundary model. It
does not neutralize or scale the protein globally. It preserves the complete
force-field chain charge when the declared QM protein-fragment charges are
added back:

`sum(q_PC at radius r) = Q_chain - sum(Q_QM_protein_fragments at r)`.

Consequently, promotion of one Asp(-) from point charges at 3.3 A to the QM
region at 3.6 A must make the PC sum exactly +1 e more positive while making
the QM scaffold exactly -1 e more negative. The total embedded-system charge
must remain invariant. Failure of any closure at 1e-6 e is a preparation
failure, not a result.

Required pre-execution records include hashes and ordered identities for all
QM atoms and PCs, chain and residue charge ledgers, every boundary shift,
minimum PC-to-QM and PC-to-link-H distances, and the point-charge potential at
the metal decomposed into <=6, 6--12, and >12 A contributions. A PC within
1.0 A of any QM center is a preparation failure.

## Primary endpoints

Use the existing score algebra and unrounded ORCA energies. For structure `s`:

`R_PC(s) = score_PC(s, 3.6 A) - score_PC(s, 3.3 A)`.

The aquo gauge cancels exactly. Define the selective difference-of-differences:

`D_PC(mono)  = R_PC(C5 monomer) - R_PC(MxaF)`

`D_PC(dimer) = R_PC(C5 dimer)   - R_PC(MxaF)`.

The balanced embedding has the requested selective-rescue capability (**YES**)
only if all three preregistered conditions hold:

1. `abs(R_PC(MxaF)) <= 1.0 kcal/mol` (generic boundary bonus suppressed);
2. `D_PC(mono) >= +2.0 kcal/mol`; and
3. `D_PC(dimer) >= +2.0 kcal/mol`.

Do not average the monomer and dimer. Thresholds will not be changed after an
energy is seen. Any failed condition is **NO** for this embedding model.

Secondary, non-decisive outputs are the absolute embedded scores, the fixed-
radius embedding shifts relative to the bare panel, and whether MxaF falls
below the current +5 prospecting threshold. The absolute MxaF sign cannot
override the primary call.

## Prespecified interpretation

- All three `R_PC` values near zero: the apparent 3.6-A rescue was chiefly a
  QM-partition artifact; answer NO.
- MxaF remains above 1 kcal/mol: the embedding failed to remove the generic
  Asp(-) boundary response; answer NO.
- C5 monomer and dimer disagree: response is conformation-dependent; answer NO.
- MxaF is suppressed and both C5 margins remain at least 2 kcal/mol: robust
  provisional selective rescue; answer YES and advance only the controls and
  strongest candidates to separately preregistered restrained QM/MM relaxation.

If a primary endpoint lies within 0.5 kcal/mol of a decision threshold, repeat
the same 12-point matrix with the only allowed boundary sensitivity: spread
the exact same correction uniformly over every retained atom of the cut
residue. A final YES then requires the same YES under both maps. This branch is
triggered by proximity to a frozen cutoff, not by a preferred result.

## Known scope limits

Moving Asp from the PC environment into QM also changes the CPCM cavity, so the
experiment tests the operational embedded pipeline rather than a mathematically
pure electrostatic term. Whole-chain net charge and distant residues can also
affect the La/Ca exchange. The distance-resolved potential is therefore
reported explicitly; no post hoc charge scaling or truncation is allowed.

This experiment is a frozen-geometry causal diagnostic. It is not yet a
relaxed QM/MM binding free energy and it does not authorize production-wide
rescoring.
