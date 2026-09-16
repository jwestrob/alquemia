# GGR: representation sensitivity and experimental clues

## Scope and conclusion

Jacob requested: “Yep. Let's see if we can get any clues by looking into GGR.”
This follow-up audited existing preparations, the four completed GGR endpoints,
previous repair comparisons and primary literature. No new quantum calculations,
preparations, scientific variants or classifications were launched.

**The original narrow Ca-side result is sensitive to a verified chemistry repair.**
The change is reproducible and localized; we found no additional demonstrated
preparation defect. Local donor chemistry, scaffold response and hydration are
credible issues, but existing results do not isolate a missing physical term.
They do not establish that a global representation would fix the result.

## What actually changed

The original v2 fragment replaced Gln140's peptide carbonyl environment with
formaldehyde-like C/O/H/H. Existing amide v3 retains its actual bonded Ile141 N
and source hydrogen, caps N–Cα, and removes the old C–N replacement hydrogen.
This is a net addition of NH: **50 → 52 atoms**.

- All original heavy atoms retain exactly the same coordinates.
- 43 original atoms are numerically identical; another six caps match at the
  original six-decimal XYZ precision. v3 writes ten decimals; the largest
  resulting coordinate difference is 6.978887e-7 Å. Thus 49 atoms are unchanged
  at the original recorded precision, not literally bit-identical.
- Donors remain Asp134, Asn136, Asp138, Gln140 backbone, Gln142 and bidentate
  Glu205: seven oxygens. There are no explicit first-shell waters.
- Ligand charge remains −3; endpoint charges remain La 0/Ca −1, singlets.
  Source/protonation hashes and donor descriptors are unchanged; preparation pH
  remains 7. No new hydrogen placement run or heavy-atom reconstruction occurred.
- The repaired graph records 25 source heavy atoms, 19 retained heavy bonds and
  seven caps, with no duplicate source mappings. This validates bookkeeping,
  not the energetic equivalence of a capped fragment and the intact peptide.

## Existing energies: no new computation

| Quantity | kcal/mol |
|---|---:|
| Original v2 S | −1.183145669 |
| Repaired v3 S, same reporting reference | +3.057472569 |
| Repair effect on R = E_Ca − E_La | +4.240618238 |

Original GGR passed its frozen direction test. **v3 has no calibrated zero or
inherited PQQ bands**, so its sign is a representation warning, not a separately
validated binary decision. The repair effect is reference-independent. All five
other already-completed single-amide repairs shifted R positively by
2.815552–4.239884 kcal/mol; GGR is not exceptional in that comparison.

Printed endpoint terms account for the shift as follows:

| Bookkeeping contribution to repair effect on R | kcal/mol |
|---|---:|
| SCF total | +4.222280348 |
| Of which: printed CPCM dielectric | −3.576713161 |
| SCF minus printed CPCM | +7.798993509 |
| Dispersion correction | +0.018347750 |
| gCP correction | −0.000008785 |

The SCF total is the sum of the two rows immediately below it; do not add all
table rows together. Final energies close against SCF + dispersion + gCP within
5.3e-10 Hartree, consistent with printed precision. Dispersion/gCP do not explain
the shift; the printed dielectric term opposes it. **SCF minus CPCM is not a
vacuum calculation:** the density was polarized in CPCM. These terms cannot
uniquely separate direct electronic, density-response and solvent mechanisms.
No term was removed, and no correction or threshold was fitted.

## Experimental state and mechanistic clues

The primary Table I confirms Kd(Ca)=25±11 µM and Kd(La)=729±4 µM (SE): about
29-fold Ca preference, measured at 25°C, pH 6, 100 mM KCl/10 mM PIPES with
2.5 µM protein. Sugar was stripped during preparation; none appears in the assay
recipe. This establishes nominal depletion, not measured zero occupancy.
Smaller late lanthanides bind more tightly than Ca. Cavity deformation and
dehydration are proposed explanations, not separately measured contributions.
[Snyder et al., 1990, Methods and Table I](https://pmc.ncbi.nlm.nih.gov/articles/PMC2899690/).

Primary captures and verification notes are pinned in
`GGR_PRIMARY_SOURCE_INVENTORY.json`. The Tb-reference dialysis temperature is
inconsistent between Methods (23°C) and Figure 3 (22°C); the competition assay
explicitly specifies 25°C. The construct restores three missing C-terminal
codons, with no reported metal-site mutation.

Our source [1GLG](https://www.rcsb.org/structure/1GLG) contains Ca and
β-D-galactose, so the assay and crystal sugar states differ. No La-bound GGR
structure is present in the inspected benchmark evidence. We have not changed
the frozen source structure or weakened its experimental label.

A follow-up NMR study found local changes on Ca→Sr substitution and increased
local solvent exposure after metal removal; distant probes changed little.
Metal removal weakened galactose binding about twofold. These observations
support local response/hydration as hypotheses, but are not La-specific and
do not quantify a Ca/La sugar-state correction.
[Luck and Falke, 1991](https://pmc.ncbi.nlm.nih.gov/articles/PMC2899689/).

## Consequence and reproducibility

The original favorable sign may partly reflect cancellation involving an
incorrect donor representation. Restoring the amide is chemically justified
even though it moves the descriptor away from that sign. A frozen core also
does not account for metal-dependent protein deformation and changing hydration.
Neither statement identifies the unique cause of disagreement with affinity.

Keep both versions, the experimental Ca direction, and the existing default.
Do not rescue GGR by altering the assay label, water inventory, geometry or
reference. A future calculation must distinguish boundary sensitivity from
physical response with an explicitly agreed comparison; none was scheduled here.

`GGR_ARCHIVED_AUDIT.json` pins the four outputs, execution receipts, original and
repaired manifests, actual unrounded energies and coordinate differences.
`audit_ggr.py` verifies hashes, convergence, paired coordinates, unchanged
charge/protonation, component closure and score algebra using these real files.
It performs no electronic-structure calculation. From the repository root:

```bash
python diagnostics/benchmark_set_20260915/audit_ggr.py \
  --result diagnostics/baseline_benchmark_20260915/RESULT.json \
  --output workspaces/benchmark_set_20260915/ggr_audit_recheck.json
```

The output path must be new; the utility refuses to overwrite a prior record.
No production code, model, reference or default changed in this follow-up.
