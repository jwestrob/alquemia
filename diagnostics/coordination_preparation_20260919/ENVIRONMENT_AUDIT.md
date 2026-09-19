# What existing environmental work actually established

Read-only audit during the current ten-endpoint geometry pilot, before its DFT
outcomes. This is not a new model, altered threshold or permission to rerun a
failed construction. All cited cases have already been inspected.

| Construction | Actual finding | Limitation for another accuracy trial |
|---|---|---|
| Balanced fixed-charge PQQ embedding | About 63–68 kcal/mol boundary jumps; selective differences survive but generic MxaF gate fails | Correct charge closure did not make changing Asp's representation continuous |
| Matched APBS transfer | Partition jump remains about 10 kcal/mol with refinement | Direct, isolated reaction-field subtraction and target terms change together; grid refinement alone does not fix it |
| Global vacuum DFT + MBIS direct coupling + TABI reaction field | Partition jump about 18 kcal/mol; surface/rotation errors 1.4/1.7 kcal/mol | Both representation and declared numerical acceptance fail; no predictive trial completed |
| Exact-density permanent field diagnostic | Replacing MBIS improves the partition difference by 8.06 kcal/mol; self-consistent core response adds only 0.97 kcal/mol improvement | About 9.67 kcal/mol remains; core polarization by itself did not repair this partition |
| Uniform CHELPG | Better exterior potential/direct coupling on four real states; coupling errors −0.12/−0.53 kcal/mol | A useful charge representation is not itself a better biological score |
| Whole-protein native GFN2/ALPB | Whole-system diagonalization/startup expense and memory made the explored route unsuitable | No accurate cheap full-protein electronic model resulted |
| Intact POLAR + frozen monopole OBC-II | Both checkpoints reverse all three primary directions | Cheap solver and rigid-motion checks pass while prediction fails |
| Chemical-group POLAR + OBC-II | PQQ pair improves, but all six alpha/GGR directions remain wrong; solvent reverses them | Grouped vacuum has all seven directions but fails the existing grouping-sensitivity criterion; deleting solvent after inspection is not validation |
| Fixed complete second-shell native DFT | PQQ gap 30.57→28.18; alpha/GGR gaps 10.92→14.81 and 14.57→18.42 kcal/mol | Relative donor-context effects exist, but this adds no classification success and narrows the PQQ gap |
| Contextual physical H preparation | Strong alpha/GGR gain from water H; PQQ native gap 30.57→30.07 despite large energy reductions | Geometry can carry useful information; generic lower-energy preparation does not guarantee better separation |
| Old physical donor/metal fixed path | All native endpoint energies decrease; discrimination does not improve | Vacuum DFT plus masked whole-minus-core score and one tangent direction were tested; this is not a fresh idea |

Primary local evidence:

- `diagnostics/pqq_balanced_embedding_20260914/RESULT.md`
- `diagnostics/global_representation_20260915/{REPORT,ARCHIVE}.md`
- `diagnostics/global_electrostatic_20260916/REPORT.md`
- `diagnostics/density_embedding_20260916/{REPORT,CHARGE_FIT_RESULT}.md`
- `diagnostics/mace_global_benchmark_20260916/REPORT.md`
- `diagnostics/mace_charge_groups_20260918/REPORT.md`
- `diagnostics/second_shell_20260919/{REPORT,PQQ_H_REPORT}.md`
- `diagnostics/mace_site_response_20260918/COUPLED_PATH_REPORT.md`

The current geometry test asks a narrower unanswered question: do explicit,
complete donor neighbors guide a cheap iterative native MACE proposal whose
original-core CPCM DFT contrast is more useful? It avoids an environmental scalar
and fixed/moved quantum atom ownership. It does **not** supply full protein
dielectric response, conditional electronic polarization of the exterior,
alternative protonation/occupancy populations, or a binding free energy.

Any next environment test must distinguish these missing contributions and keep
its energy accounting explicit. The prior failures do not justify another
parameter sweep over radii, dielectric constants, charge schemes or thresholds.
The actual result of the current geometry trial should determine whether extending
this route has marginal utility. Other agents own water basins and chemical-state
competition; duplicating those would obscure rather than improve the evidence.
