# Accommodation adds discrimination beyond the consistent static pocket

The same strict native solvent solver and same prepared pockets give
**199 correct /1 wrong /8 inconclusive /17 unavailable** when held at their
origins, versus **207 /0 /1 /17** after the existing common-pool accommodation.
All208 available sources are shared. Seven inconclusives and one wrong call
become correct; no correct static call is lost. This demonstrates an incremental
benefit of geometry response on this consumed structural challenge.

The original released static method gives203/2/2/18. On207 sources shared with
the final candidate, that becomes206/0/1. The remaining inconclusive is C5AXV8
La sample3, which the released method called correctly; the improvement is not
uniform across every source. Its ambiguity is retained.

## Fair reference and physical comparison

Each method uses its own25-member canonical calibration under the same fixed
extrema/minimum-gap rule. The new origin-only reference was frozen before the
full225 transfer result was read. No transfer structure, crystal or PLM label
entered that fit. The strict adaptive reference remains unchanged.

All eight static-to-adaptive repairs also require a physical score change under
the *adaptive* bands: their static origins remain wrong or inconclusive under
those bands. The improvements therefore cannot be explained solely by using a
different calibration. Their accommodation shifts are:

| Source | Static → accommodated | Change in Ca-minus-La contrast, model kcal/mol |
|---|---|---:|
| A0A3F2YLY8 Ca0 | inconclusive → correct | +4.550410 |
| A0A3F2YLY8 Ca3 | inconclusive → correct | +6.978411 |
| A0A3F2YLY8 La0 | inconclusive → correct | +7.245735 |
| A0A3F2YLY8 La3 | inconclusive → correct | +9.491256 |
| A0ACD6B9F2 Ca4 | wrong → correct | +17.827801 |
| A0ACD6B9F2 La4 | inconclusive → correct | +21.659485 |
| MMOL_1770 Ca1 | inconclusive → correct | +16.092890 |
| Q88JH0 Ca0 | inconclusive → correct | +4.090336 |

These eight structures belong to four protein groups. They are not eight
independent biological discoveries, and their positive shifts do not establish
that every protein or model improvement must favor La. The finite candidate pool
contains the origin and both metal-specific proposals; it is not a thermal
ensemble or a converged binding free energy.

## Group results and coverage

Both arms retain23/25 complete La4 groups,22/25 Ca5 groups and21/25 balanced
groups; all available group decisions are correct. Static gives92correct,
2inconclusive and6unavailable three-source subsets; accommodation gives94correct
and6unavailable. The released method also has94correct subsets: the gain relative
to release is individual-source fidelity and one source's numerical coverage,
not newly correct protein-group medians. All original75 groups/100 subsets and
17 source preparation exclusions are preserved.

## Cost, artifacts and scope

This ablation launches **zero molecular calculations**. It reuses the strict225
matrix and the strict32 canonical origins. Costs of producing those energies and
the earlier searches remain in their actual reports; this analysis does not
claim them free or provide a new scanner timing.

The editable figure (SVG/PDF/PNG), exact plotted CSV and hash receipt are under
`workspaces/discrimination_transfer_figures_20260923/strict_ablation_v1/`.
The figure uses only a display offset from each frozen band midpoint, without
rescaling. Root visually inspected the generated PNG.

Full result:
`workspaces/strict_static_ablation_20260923/run_v1/COMPARISON.json`.
Input strict-transfer SHA256:
`6ed4c1f4acf4359010922ff3f088c5b36105caf52dc92c734e27dd7e35d0df46`.
Reference pins and exact commands are in [COMMANDS.md](COMMANDS.md).
Three real-artifact tests pass; no scientific executable is mocked.

**Recommendation:** retain accommodation in the candidate. Its practical
three-source pocket definition still needs the separately running buffered-pocket
test. Production remains unchanged.
