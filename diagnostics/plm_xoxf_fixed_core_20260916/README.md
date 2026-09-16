# Two existing PLM XoxF AF3 sites — fixed-core preparation

## Completed results — 16 September 2026

**Job1200300 completed/PASS: two pairs, four endpoints, zero failures/unsupported cases; 5m22s on64 CPUs.** Every endpoint has normal ORCA termination and converged SCF. Source/input/output/receipt hashes, core coordinates and the full preparation preflight pass after execution; contrasts and exact released band assignments were independently recomputed by the parent. The watcher delivered terminal status and exited.

| Protein | CN including PQQ | S (kcal/mol) | Fixed-core band |
| --- | ---: | ---: | --- |
| AG41 / Gemmatimonadota XoxF 38_54 | 7 | **17.21776065297137** | Indeterminate |
| Rokubacteriales XoxF 2928_2 | 8 | **29.971263921125484** | La-supported |

Rokubacteriales XoxF is **9.98714607833281 kcal/mol higher** than the AF3 ADH9 32301_3 from the same genome (19.984117842792674, indeterminate). This is a descriptive comparison of the existing fixed-core scores, not a measured affinity difference or a functional substrate assignment. No new ADH9 calculation was performed. Hydrogenophilalia ADH9 4380_6 remains model/preparation-sensitive (Protenix21.450083 versus AF3/repaired-H9.720464); its biological metal preference is unresolved.

AG41 is **not evidence of Ca dependence**: it lies between the supported bands. The two current XoxFs share the E/N/D/D/Arg core composition and charge; their 12.753503268154114 kcal/mol score separation therefore does not arise from different retained core composition/charge. Geometry and generated-H differences remain coupled, and no causal decomposition was attempted. Both inputs are predicted structures of phylogenetically curated PLM proteins, not experimentally verified positives. The original50-step H preparation is retained and not force-converged, although both actual core geometries pass; do not omit that limitation when interpreting transfer.

Raw energy pairs (Hartree): AG41 La−2536.085939966249, Ca−3182.136047549587; Rokubacteriales La−2536.103825453389, Ca−3182.133609035582. Reporting gauge A−646.0775458314704; factor627.509474. Exact fixed-core bands remain Ca-supported S≤14.857129202922806, La-supported S≥23.460061205609236. No calibration/default changes or further calculations launched.

Authoritative results: `results.json`, `results.tsv`, `completion.json`, `post_completion_review.json` in the PLM execution directory. All preparation and execution records remain unchanged.


Jacob approved this two-protein comparison with “Go for it”. The exact approval
is at `EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue/xoxf_fixed_core/authorization.json`.
Only the existing AF3 models listed below were considered. No folding, MSA
search, reference/control run, calibration refit or old-queue restart occurred.
This preparation agent submitted no jobs; the parent owns the endpoint executor.

## Selected structures

All six existing AF3 structures preserve their exact full protein sequences,
1-based residue numbering, one protein chainA, La chainB and PQQ chainC.
The existing input sequences have589aa (AG41) and562aa (Rokubacteriales).

| Protein/gene | Selected AF3 sample | Typed CN / protein CN | Protein–La iPTM | Explicit core | Asp–Arg contact |
|---|---:|---:|---:|---|---:|
| AG41, `PLM0_60_coex_jun17_scaffold_38_54` | 0 | 7 / 4 | 0.98 | E190/N274/D316/D318/R343 | 2.943638 Å |
| Rokubacteriales, `PLM2_30_coex_sep16_scaffold_2928_2` | 0 | 8 / 5 | 0.98 | E186/N252/D302/D304/R329 | 3.075643 Å |

Full target IDs are `PQQSEQ_242fa05e3ffc20087d42` and
`PQQSEQ_faa97386262eec4316fc`, respectively. Both use seed101/sample0.
Selection follows typed CN descending, protein–La iPTM descending, protein CN
descending, sample index ascending. All three samples of each protein tie on
the first three criteria, so sample0 wins. Admission remains confidence≥0.9,
CN≥7 at3.1 Å and at most2 direct N donors. Each has one direct N donor.
The catalytic cationic partner is checked only after selection; no alternative
model was selected to rescue preparation.

Existing alignment-derived `coordination_review/saved_residue_observations.tsv`
maps E/N/D/D+2 at columns3082/4141/4687/4708 to the four recorded full-protein
positions. Both catalytic Asp residues immediately follow Trp. Each chosen
structure has exactly one Arg/Lys side-chain contact within3.5 Å of that Asp:
Arg at D+27, with Asp OD2–Arg NH2 the closest pair. The next Arg/Lys lies
7.61 Å and9.14 Å away. Exact sequence windows, role evidence, all six metrics,
source hashes and candidate partner distances are in `selection_review.json`.

## Preparation and hydrogen checks

Both original frozen preparations pass core geometry; **neither required
hydrogen recovery**. The unchanged fixed-core wrapper and calibrated helpers
produced complete oxidized PQQ³⁻, the four conserved E/N/D/D residues and the
mapped cationic Arg. No waters, point charges, synthetic ligands or heavy-atom
relaxation were added. Both cores have80 atoms: C27H31N6O15 plus metal,
La charge−2 and Ca charge−3, closed-shell singlets. The two arms have identical
nonmetal coordinates and metal positions.

| Actual core check | AG41 | Rokubacteriales |
|---|---:|---:|
| Minimum H–H separation | 1.710199 Å | 1.665859 Å |
| Protein H–parent distance range | 1.172298–1.195718 Å | 1.181025–1.202092 Å |
| Minimum nonmetal atom separation | 1.010000 Å | 1.010000 Å |
| Source heavy-atom repairs | 0 | 0 |

The observer captured the original OpenMM System and minimizer states without
changing the call. Both original minimizers used all50 iterations; remaining
full-protein H force RMS values are40.2325 and42.2306 kJ/mol/nm. **These are not
force-converged hydrogen structures.** They retain the frozen calibration's
original preparation because the actual core geometry passes. No longer
minimization or alternative H preparation was silently applied. This differs
from the prior ADH9 repair, where severe overlapping H required an explicitly
approved numerical recovery. Observation traces and geometry checks are saved.

## Files and validation

All compute products reside under:

`/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/plm_xoxf_fixed_core_20260916/`

- `candidate_manifest.json`, `candidate_implementation_pins.json`,
  `selection_review.json`: frozen source selection, homology/coordinate evidence
  and the current approval.
- `adapter_inputs/`: byte-identical selected raw CIF/confidence copies with
  filenames expected by the unchanged preparation wrapper.
- `original_preparation/`: original normalized/protonated structures, heavy
  coordinate receipts, fixed-core XYZ/input pairs and original carve manifests.
- `hydrogen_audit/`: original captured Systems, full-precision positions,
  minimizer observations and actual core-geometry reports.
- **`prepared/prepared_pairs.json`**: both admitted cases for the existing
  manifested La/Ca endpoint executor.
- **`prepared/geometry_validation.json`**: PASS with per-case geometry and
  preparation mode. Each final carve also pins its individual validation.
- `prepared/implementation_pins.json`: source, original calibration helpers,
  approval, selector/preparation code and observed results.

Nine read-only software regressions pass in `test_preparation.py`: the two
actual cores pass; the prior overlapping ADH9 core fails; coordinate/identity
and charge invariants pass; synthetic overlaps/nonfinite values are rejected;
selection ties and CN priority are preserved. `software_validation.json`
records the completed checks. Tests do not generate additional predictions
or endpoint calculations.

The protocol remains `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`,
native ORCA r2SCAN-3c/CPCM(Water)/DefGrid3/NoAutostart. The existing release
defines Ca-supported S≤14.8571292029, La-supported S≥23.4600612056 kcal/mol,
with the intervening interval indeterminate. The parent will collect the
endpoint results. These are protocol-specific scores rather than measured
biological metal use, binding affinity or substrate specificity.

## Parent execution — 16 September 2026

Job1200300 is running on node-64-768g-16 with64 CPUs (four endpoints,16MPI ranks each). No new inputs or geometry changes. Parent preflight and independent executor review PASS. Manifest SHA a4dd2ad7e870be20ceec97b5931e7e89843fbf64d933869836032d6ea2657993; endpoint bundle `workspaces/plm_xoxf_fixed_core_20260916/execution/`. Same-session completion watcher3160464 is running. PLM execution/result root: `/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue/xoxf_fixed_core`. Final exported scope names XoxF; legacy internal schemas are retained from the frozen executor.
