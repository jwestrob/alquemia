# Independent electronic-treatment diagnostic v1

Declared 2026-09-19 before new electronic outputs. Jacob explicitly requested
parallel investigation of electronic-treatment errors, chemical states and
protein environment; this agent owns electronic treatment. His blanket pilot
authorization remains active. Baseline/default and historical artifacts stay
unchanged. These are consumed development cases, not blind validation.

## Question and frozen first panel

Does an independent correlated wavefunction treatment change the relative
Ca/La interaction contrast, or the useful water-preparation effect?

Six endpoints: Ca/La for the complete existing repaired-amide GGR1GLG core
(52 atoms), original alpha1F6S core (40 atoms), and the same alpha core with
previously prepared contextual water H coordinates (40 atoms). Copy exact XYZ,
charge and multiplicity from their recorded native DFT tasks. No re-carving,
optimization, donor/proton/water-count changes or new experimental labels.
Only the existing water-prepared alpha endpoints have metal-dependent H
positions; method comparisons use each endpoint's identical stored geometry.
Alpha versus GGR is one condition-qualified biological direction comparison;
the two alpha preparations are a mechanistic contrast, not independent samples.

## Electronic model, fixed before outputs

ORCA6.1.1 DLPNO-CCSD(T1), TightPNO, TightSCF, RIJCOSX/DefGrid3,
CPCM(Water) with explicit PTES (`CPCMccm 2`), def2-TZVPPD on H/C/N/O/La,
cc-pwCVTZ on Ca, AutoAux correlation fitting and def2/J Coulomb fitting.
Keep native La46-electron def2 ECP; Ca remains all-electron. Explicit frozen
core: Ca10 and La46 total electrons (La46 are ECP electrons; no additional
explicit La electrons frozen). C/N/O1s stay frozen under ORCA's default.
Thus Ca3s/3p and La5s/5p remain correlated. Iterative triples use the program's
named-method settings. No added D3/D4/gCP, basis override of native baseline,
new aquo reference or inherited decision band.

This is a physically motivated **finite-basis correlated diagnostic**, not
CC/CBS truth: La basis/core-valence flexibility, residual basis/PNO errors,
single-reference character, approximate CC solvent response and inherited
core chemistry remain limitations. The explicit Ca core-valence basis avoids
using a valence-only basis for important semicore correlation. Def2's ECP keeps
the same scalar-relativistic metal-core representation as the native baseline;
this first comparison cannot diagnose errors common to that ECP.

Retain SCF, CCSD, triples, CPCM correction, T1 diagnostic, convergence and
basis/ECP/core counts from real output. T1 alone is not a universal test of
multireference character. CPCM's printed correlation term is already included
in CC energy; do not add it twice or treat it as the entire solvent effect.

## Frozen comparisons and execution

Report R=E_Ca-E_La and D=(R_alpha-R_GGR)*627.509474 kcal/mol for each alpha
preparation, separately for native baseline and CC. Positive D is the expected
relative direction. Also report (R_prepared-R_original)*627.509474 for both
methods. A method-specific common metal/aquo offset cancels within D; raw R
from different methods is not on the same calibrated scale. No fitted band,
affinity probability or post-result choice of favorable method/geometry.

Four concurrent16-rank tasks on64CPU/256GiB; maxcore3000MB/rank. The last two
tasks use the same16-rank policy. Six finite tasks, no project compute/time
budget; scheduler limits apply. CC throughput is not yet measured: anticipate
tens of minutes to hours, rather than the native endpoints' roughly2minutes.
Use existing manifest runner, immutable inputs and receipts. Failures remain
visible and do not become baseline substitutions.

After measuring small-core throughput, prepare the same-method PQQ1H4I/4MAE
four-endpoint guardrail. It is not silently included in this six-task manifest.
Any basis/convergence follow-on receives a separate version before execution.

## Prior work and primary sources

Legacy B97-3c comparisons retained the signs for XoxF, MxaF, tannase and
calexcitin, with one tier change. They used old preparations/references and
do not establish accuracy of today's native r2SCAN-3c protocol. See
`VALIDATION.md` section7 and `legacy/b97_3c_panel/` actual inputs/outputs.

- [ORCA6.1 coupled cluster](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/mdci.html): iterative triples and named PNO settings.
- [ORCA6.1 CC/CPCM](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html): PTES and energy accounting.
- [ORCA6.1 frozen core](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/frozencore.html): Ca10 default, ECP counting and semicore basis requirements.
- [ORCA6.1 basis sets](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basisset.html): cc-pwCVTZ Ca and def2/ECP availability.
- [Minenkov et al.2017](https://doi.org/10.1039/C7CP00836H): alkaline-earth subvalence correlation can matter; old noble-gas frozen cores are not adequate by default.
- [Preferential binding DFT study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8028316/): small-ligand CC benchmarks motivate an independent check but do not validate this assembled real-core protocol.
