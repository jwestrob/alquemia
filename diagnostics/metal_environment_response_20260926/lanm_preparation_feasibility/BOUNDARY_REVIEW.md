# Hans8DQ2 EF3 peptide-boundary review

Status: native export, local accounting and final prepared field/geometry
independently verified. Electronic-state/energy qualification remains unrun.
No electronic energies, new force-field evaluations, coordinate repair or duplicate
preparation were performed. The preparation agent owns generated artifacts.

## Recommended local rule

Use the pinned ff19SB residue templates already identified for this exact
normalized source. The QM region comprises full residues83–94, neutral partial
peptide groups from65,66,82,95, neutral selected waters, and oneLn(III).
Its formal charge is−1. The original protein charge is−4; deposited spectators
areLa201(+3),La202(+3),Na204(+1). Full-source charge is+6, so the finite field
must sum+7.

For each partially selected residue r, define:

- Q_r: original complete-residue FF charge.
- Q_QM,r: formal charge assigned to its selected capped peptide contribution,
  zero for all four partial residues here.
- S_r: actual selected QM atoms plus the cap-omitted CA atom.
- T_r: remaining source-bonded **heavy** neighbors of S_r within that residue.
  Derive from the real graph; reject absent, duplicate, selected or unsupported
  recipients rather than introduce a fallback.
- delta_r = Q_r − Q_QM,r − sum(q_i for retained MM atoms in residue r).
- Set every S_r field charge tozero; add delta_r/|T_r| to each recipient.

This conserves each partial residue's formal partition locally. No distant
residues, metals, waters or global normalization are adjusted. Full loop residues
have exactly the expected total FF charge−4 and are wholly removed from the
field. Selected waters have zero total charge and are wholly removed.

Independent arithmetic from pinned XML templates gives:

| Residue | Removed QM atoms and MM1 CA | delta/e | Actual heavy recipients | Increment each/e |
|---|---|---:|---|---:|
|Leu65|C,O,CA|−0.0224|N,CB|−0.0112|
|Glu66|N,H,CA|−0.1830|C,CB|−0.0915|
|Ala82|C,O,CA|+0.0631|N,CB|+0.03155|
|Trp95|N,H,CA|−0.1713|C,CB|−0.08565|

These independently derived expectations now match the actual `export_v1`
parameters exactly; they are not a new parameter fit. The selected partial peptide atoms carry FF charge−0.3077;
the omitted CA atoms carry−0.0059. Before redistribution the remaining protein
field sums+0.3136, which becomeszero after the summed−0.3136local correction.
The spectator+7 then gives the required field+7 and total+6.

Simply deleting or shifting MM1 charges would not fix the separate fractional
QM-fragment versus formal-charge mismatch. Conversely, satisfying the global
sum alone would not justify an arbitrary global correction.

## Physical interpretation and limits

The heavy-only convention keeps charge increments off the nearby small boundary
H sites. All recipients are directly bonded exterior source atoms; no charge is
placed at a link cap or a nucleus not present in the source. The largest increment
is0.0915e. This is a **local residue-charge-conserving boundary convention**, not
an exact electrostatic reconstruction or a claim of preserved dipole.

ORCA documents charge-shift and charge/dipole redistribution alternatives at
covalent boundaries. Our explicit per-residue restoration is not asserted to be
its native CS implementation. Record original/shifted charges and any field-only
dipole change as representation diagnostics, not quantities adjusted to obtain
favorable metal results. [ORCA boundary documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/multiscalesimulations/qmmm-general.html#charge-alteration)

Keep the complete actual spectator identities/coordinates fixed across targetLa
and targetDy. Treating them as formal monopoles is the declared frozen electronic
scout approximation, not a spectator force field or open-shell quantum model.
No spectator charge must be moved to compensate a peptide boundary.

## Geometry and perturbation requirements

Use normalized `sources/Hans_8DQ2/atoms.json` coordinates, with source IDs matched
exactly to exported ff19SB+TIP3P topology. Do not reuse original long-H PDB
coordinates or call addHydrogens. Preserve full source occupancy and all actual
waters; the QM water subset moves representation only, not physical membership.

The +2° Asp85 CA→CB torsion must rotate CG,OD1,OD2,HB2,HB3 together. CA andCB
stay fixed, so this is a chemically intact sidechain motion: no bond lengths or
CB valence angles are intentionally strained. All other real atoms, caps and
spectators remain fixed. Verify covalent geometry, recorded atom displacement
bound and actual nearest contacts before admitting B. No atom may be selected
because of its metal-score effect.

The point-charge array must be identical for A/B and targetLa/Dy because all
moved atoms are QM. Charge/state/input hashes must still distinguish the four
endpoints and the changed core geometry. Cap source mapping and Jacobians must
remain physical, even though this particular motion leaves all cap positions
unchanged. This does not qualify other motions or a whole-protein force model.

## Separate finding for any later additive mechanical model

ORCA's native default omits link-boundary bonded terms containing only one MM
atom, controlled by DeleteLADoubleCounting/DeleteLABondDoubleCounting. The earlier
1H4I ledger retained all cross bonded terms, so it is explicitly a different
boundary approximation; passed classical derivative checks do not establish
native additive equivalence or eliminate cap double-counting concerns. Preserve
that historical test. The present EF3 **electronic-only** scout adds none of
those classical terms, so this distinction does not block field preparation.
[ORCA bonded-boundary policy](https://www.faccts.de/docs/orca/6.1/manual/contents/multiscalesimulations/qmmm-general.html#bonded-interactions-at-the-qm-mm-boundary)

## Native export and motion check completed

`BOUNDARY_EXPECTATIONS.json` pins the actual native export. All1,887 source IDs
and normalized coordinates match exactly. Its direct charge sum is−4 within
roundoff; all four boundary corrections match the independent template arithmetic.
The expected field has1,696 rows if the four zero-charge MM1 atoms are omitted:
1,887 source atoms−190 real selected nonmetal atoms−4CA+3spectators. The omitted
CA atoms must still remain in the physical source mapping.

An independent in-memory check of the proposed Asp85 motion gives maximum
physical displacement0.0829489Å, maximum covalent-length change7.8e−15Å and
maximum valence-angle change5.8e−13degrees. The closest moved-to-source nonbonded
contact (excluding1–2/1–3 pairs), OD1–backboneH85, increases2.068959→2.123500Å.
Metal distances OD1/OD2 change2.695897/2.597779→2.703276/2.596560Å. No new
coordinate file was generated for this independent geometric review. Final
prepared A/B and field products have now passed the separate review below.

## Final prepared artifact review: pass

`PREPARED_REVIEW.json` pins actual
`workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Hans8DQ2/INPUTS.json`
(SHA256 `c5886704ee541752ec2adfe7ada0b5f8b2c0df6014d5ed8676b74447d1141c3f`).

The actual field has1,696 rows and charge+7 within3.6e−15e. With QM−1 the
whole-source charge is+6. Each emitted charge agrees with independent native
parameter/local-graph accounting to1.11e−16e. ActualLa/La/Na spectators are
present at their source positions; all four endpoints use the exact same field
file. No selected real QM atom or omitted MM1 CA remains in the charge array.

All190 selected nonmetal atoms match the proposed source-ID set; the target ion
andfour caps give195QM atoms. A coordinates exactly equal the normalized source.
La/Dy coordinate arrays match within each configuration. Only the five intended
Asp85 atoms move inB; independent torsion reconstruction, fixed cap coordinates
and analytic cap Jacobians match exactly. Serialized XYZ coordinates reproduce
the maps exactly. Physical electron parity was checked, but no native Dy
Hamiltonian, spin localization or energy was evaluated by this review.

`REVIEW_PREPARED.py` reproduces this independent read-only review, accepting
explicit --inputs and a fresh --output path. It does not construct a force-field
System, prepare a new structure, or call any molecular executable. The earlier
export-only attempt is retained by its owner; its final regenerated export
retains the same native parameter identity. No production changes occurred.
