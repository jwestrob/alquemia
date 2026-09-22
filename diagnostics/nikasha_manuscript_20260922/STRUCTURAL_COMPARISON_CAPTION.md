# Editable figure: structural-repeat comparison

## Caption

**Figure X. La/Ca class discrimination across structural repeats of 25 previously
used PQQ reference proteins.** Each method retains its own frozen decision bands;
no threshold was fitted to these folds. **(A)** The declared challenge contains
225 additional models: four La-conditioned and five Ca-conditioned folds per
protein. On the identical 207 sources available to both methods, preserved DFT
classifies 198 correctly, 3 incorrectly and 6 as inconclusive; the released context
composite classifies 203 correctly, 2 incorrectly and 2 as inconclusive. The gray
segment retains the other 18 declared sources outside this common comparison.
Seventeen lack supported preparations, and one lacks the context solvent result;
the latter's otherwise available DFT result is also excluded from the matched
counts. **(B)** All four three-of-four subsets of the additional La-conditioned
folds are evaluated for each protein (100 subsets). Each score is the median under
the corresponding fixed protocol, and all three members must be available.
Both methods classify the same 94 available subsets correctly, with no incorrect
or inconclusive calls; six subsets remain unavailable. Folds from the same protein
and overlapping subsets are correlated, not independent biological observations.
These results assess structural robustness within consumed PQQ reference groups,
not broad independent affinity validation, calibrated probabilities or quantitative
binding free energies.

The preserved DFT method is the archived native r2SCAN-3c/CPCM core protocol.
The released context composite combines native MACE-OMOL vacuum context energy
with its GFN2 ALPB-minus-vacuum solvent contribution. This figure does not show the
subsequent proposal or adaptive-pool challengers and makes no hardware speed claim.

## Editable files and provenance

- SVG with editable text: [structural_comparison.svg](../../workspaces/nikasha_manuscript_20260922/structural_comparison_v2/structural_comparison.svg).
- Vector PDF with embedded TrueType fonts: [structural_comparison.pdf](../../workspaces/nikasha_manuscript_20260922/structural_comparison_v2/structural_comparison.pdf).
- PNG preview and machine-readable `figure_data.json` are beside these files.
- Generator: [render_structural_comparison.py](render_structural_comparison.py).
- Source: `workspaces/nikasha_recovery_20260922/proposal_comparison.json`.
- Source SHA256: `d2bd4a91d9653363f5f8d3998ad20009963e94f0303ac44a9b4ce92ef1d5c440`.

The generator reads the recorded outcomes, verifies the 225-source/25-group and
100-subset counts, replays the existing triple tallies, and records every excluded
source. No molecular calculation, fitting, probability estimate or new statistical
test is performed. SVG/PDF/PNG were rendered successfully with Matplotlib 3.10.5;
the final preview was visually checked for clipping and label overlap. Version 1
is retained as a layout draft; version 2 is the reviewed figure.

To regenerate into a new output directory:

```bash
/home/jwestrob/.pyenv/versions/3.12.0/bin/python diagnostics/nikasha_manuscript_20260922/render_structural_comparison.py \
  --source workspaces/nikasha_recovery_20260922/proposal_comparison.json \
  --output workspaces/nikasha_manuscript_20260922/structural_comparison_v3
```
