# Fixed donor-displacement response figure

**Suggested caption.** Metal-dependent response to fixed, intact extra-Asp
carboxylate rotations in three consumed development contexts. Bars show
ΔR(q) = [E_Ca(q) − E_Ca(0)] − [E_La(q) − E_La(0)] for the two predeclared
χ₂ displacements, −0.2 and +0.2 rad, computed with native r2SCAN-3c/CPCM,
native MACE-OMOL, and the composite MACE + [native GFN2-ALPB(water) − native
GFN2-vacuum] model. Each method uses its own endpoint energy difference; a common
reference offset cancels. Positive ΔR indicates an electronic-contrast shift
toward La, through greater stabilization or smaller destabilization relative to
Ca. The 4MAE panel uses an expanded vertical scale; both PLM panels share the
same scale. All three contexts and both displacement signs are shown (six
differential responses, twelve displaced metal endpoints). Bars represent
individual fixed-geometry evaluations, with no fitted curves or estimated error
bars.

Native MACE reproduces the direction and approximate magnitude of the tested
native DFT differential response. The added solvent term increases the
differential error in both PLM examples. These fixed perturbations preserve
chemical states and all source coordinates outside the declared intact donor
motion; they establish neither biological metal preferences nor relaxed minima,
binding free energies or general force accuracy. In particular, this experiment
does **not** provide native DFT validation of the later adaptive-angle or
joint-metal optimized candidates. Both PLM proteins remain biologically
unlabeled. Full IDs are `PQQSEQ_83440678cbbd658047c9` and
`PQQSEQ_07ab500e3df76b30d71c`; each uses original AF3 sample0.

## Files

Products are under
`workspaces/nikasha_manuscript_20260922/fixed_displacement_response_v1/`:

- `fixed_displacement_response.svg`: editable vector graphic with text retained.
- `fixed_displacement_response.pdf`: vector PDF with embedded TrueType fonts.
- `fixed_displacement_response.png`: 220 dpi preview, visually checked.
- `figure_data.json` / `figure_data.csv`: complete unrounded plotted data.
- `PROVENANCE.json`: source, final analysis, rendering script and output hashes.

The renderer checked all six declared source rows, the twelve-endpoint
denominator, each Ca−La subtraction and each MACE+solvent component sum. It
verified the archived final analysis hash. No fit, new source selection,
molecular call or score modification was performed. Existing structural
comparison figures and the shared methods outline were untouched.

## Reproduce

Use a new output directory; the renderer preserves existing figures:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  diagnostics/nikasha_manuscript_20260922/render_fixed_displacement_response.py \
  --source diagnostics/accommodation_torsion_20260920/FINAL_COMPONENTS_v1.json \
  --report diagnostics/accommodation_torsion_20260920/FINAL_DFT_REPORT.md \
  --output workspaces/nikasha_manuscript_20260922/fixed_displacement_response_v2
```

The saved v1 figure was generated with Matplotlib 3.10.9. Exact editable-figure
paths and hashes are also recorded in [FIXED_DISPLACEMENT_ARTIFACTS.json](FIXED_DISPLACEMENT_ARTIFACTS.json).
