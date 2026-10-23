# Hydrated comparison: independent review

21 October 2026. Saved-data arithmetic only, one thread; no new molecular
calculation or modification of the reviewed collector.

## Result

**The saved comparison is arithmetically correct.** All 20 MACE, 20 transferred
AMOEBA and 8 frozen-f PBE0 rows are present. Independent reconstruction passed
4,582 checks and verified 188 distinct artifact pins. All 21 exchange contrasts
replay exactly; all 42 ligand derivatives agree within 1.42e-14 kcal/mol/Å and
six force-error summaries within 2.13e-14 kcal/mol/Å.

| Hydrated central donor | MACE − PBE0 contrast error | AMOEBA − PBE0 contrast error |
| --- | ---: | ---: |
| Monodentate acetate | −2.313 kcal/mol | −6.811 kcal/mol |
| Bidentate acetate | −2.592 kcal/mol | −7.838 kcal/mol |
| Acetamide | −0.569 kcal/mol | +7.388 kcal/mol |

These are discrepancies relative to the stated frozen-f PBE0 calculations,
not demonstrated errors relative to physical truth. Neither model is qualified
for selectivity prediction by these fixtures.

## Important causal limit

**The difference from the original dry panel changes radial geometry as well
as hydration.** Original central metal–donor distances are approximately:

| Donor | Original dry | Hydrated central |
| --- | ---: | ---: |
| Monodentate acetate O | 2.021 Å | 2.500 Å |
| Bidentate acetate O/O | 2.248 / 2.240 Å | 2.500 / 2.500 Å |
| Acetamide O | 2.090 Å | 2.500 Å |

The smaller MACE/reference discrepancy in the hydrated fixtures is real. It
does not isolate hydration as the cause. A later comparison would need dry
clusters extracted from the hydrated coordinates, with the same metal–ligand
placement and a consistently matched aqueous reference. No such additional
calculation was run for this review.

## What was checked

- **Reaction/sign:** `La(mixed) + Dy(aquo8) -> Dy(mixed) + La(aquo8)`.
  The collector's expression is exactly this balanced four-energy contrast.
  Both elements and total charge balance. Free ligand and released-water terms
  cancel in the matched La/Dy difference; their omission does not create an
  imbalance. A negative value favors the Dy mixed cluster relative to the
  specified aquo8 reference, not necessarily the experimental protein.
- **Coordinates and state:** all original XYZ records match. Actual Tinker
  native coordinates also match exactly, including type/bond mapping. La/Dy
  atom order and positions match. Native ORCA atom identities and printed
  coordinates agree within 4.19e-8 Å, consistent with output precision.
  Native charges, effective singlet multiplicity, electron counts and 46/55
  electron ECP conventions were checked. Physical Dy multiplicity six remains
  distinct from the effective closed-shell reference. Native admission does
  not qualify physical states, strict residuals or reference accuracy.
- **Units/sign:** actual Tinker output energies and negative analytic gradients
  replay exactly. Actual ORCA engrad energies and `−gradient ×
  627.509474/0.529177210903` forces replay exactly. MACE uses 23.060547830619
  kcal/mol per eV. Reused and new MACE execution records specify identical
  model/checkpoint metadata; the large checkpoint binary was not reread or loaded.
- **Radial derivative:** `dE/dδ = −u · sum(F_i)` over **every ligand atom**,
  including its hydrogens. The saved path moves the whole ligand by `δu`,
  with metal and waters fixed. Unit direction and actual translations were
  checked at all three positions. This is not a metal-only force or a sum
  over coordinating oxygens alone. No extra normalization by ligand size is
  appropriate for this translation coordinate.
- **Force errors:** RMS uses all Cartesian components, with denominators
  87 (29 atoms, mono acetate), 78 (26 atoms, bi acetate), and 93 (31 atoms,
  acetamide). Paired error is the model-minus-reference error in `F_Dy−F_La`.

## Preserved evidence

The matching `experimental_target_v1` workspace contains:

- `hydrated_collector_review_v1.py` and `.json`: independent implementation,
  numeric results, native-state metadata and reviewed pins.
- `hydration_causal_confound_v1.json`: pinned dry/hydrated source geometries
  and explicit donor-distance comparison.
- `hydrated_collector_review_attempt1.py` and its failure record: initial
  path-resolution error before reading scientific input; corrected in the
  executed review without changing scientific checks.

Successful saved review runtime: 0.511 s internal wall time. Zero molecular
calls, no scheduler submission, and no changes to the root collector. The
separate GPU adapter preparation and its tests likewise establish no new
physical accuracy result.
