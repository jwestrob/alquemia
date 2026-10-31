# Nikasha structural increment — frozen 2026-10-30, before new fits

Question: do predetermined topology and one empirical size-fit descriptor improve
La/Dy selectivity prediction beyond donor composition and the released graph model?
Repository HEAD: c3b250c1c3fbf6c8927d4def6175e0b3ebbc8bb3. No subsequent commits;
unrelated dirty working files are preserved. Historical donor result is reproduced
first in Slurm 31942, into a separate workspace, without overwriting old results.

Target: directly tabulated NIST SRD46 v8 ML=M.L, Ln(III), 1:1, 25 C, same ligand
species and ionic strength. y=-R*298.15*ln(10)*(logK_Dy-logK_La), kcal/mol;
R=0.00198720425864083. Negative means Dy preference. No partial-series labels.
Start with existing formula-checked direct pairs (113 including high-I records),
recover the 15 previously identified diamine tetraacetates with explicit relative
stereochemistry, then freeze a row-level curation/exclusion and grouping manifest.
All available matched conditions remain linked. One condition per ligand, nearest
0.1 M, is primary; other conditions are linked secondary records, not independent
chemical examples. No new ligand classes beyond this bounded recovery.
Ambiguous/censored constants, incompatible condition annotations, formula/identity
mismatches and unresolved exact stereochemistry are retained with exclusion reasons.
Chemical identity preserves stereochemistry/protonation metadata; overlap identity
removes stereochemistry, normalizes charges and tautomers. Neither is a feature.
NIST tables stay under the ignored workspace, per existing redistribution policy.

Models, no outcome-dependent feature expansion:
M0 training mean. M1 all 12 historical SMARTS counts plus historical macrocycle
indicator and charge of NIST's fully deprotonated L. M2 additionally: number of
potential 5-, 6-, and 7+-member chelate rings; potential donor excess over 8;
rotatable bonds on donor-to-donor shortest paths; ring bonds on those paths;
largest macrocycle size (>=9, else 0). Chelate ring size is donor-to-donor shortest
path bond count +2 (including the dummy metal); retain only edges of the donor
minimum spanning forest with paths <=5 bonds. One O per carboxyl group; amine,
pyridine, ether, alcohol/phenol and amide O donors are potential donors, not measured
coordination numbers. Sulfonates are spectator groups for the geometry rule.
M3 adds only s=E_organic(x_Dy*)-E_organic(x_La*).

Ridge alpha grid [0.01,0.1,1,10,100]; StandardScaler fitted only inside each train
fold, no feature selection. Outer leave-one-chemical-family-out, explicit table
frozen before fitting; broad families keep homologous aminopolycarboxylates together,
all diamine tetraacetates in that broader family. Inner leave-one-family-out,
ligand-weighted MSE for selection; lowest grid alpha wins ties. No random row CV.
Report an additional development-challenge fit holding the 15 tetraacetates and
all related aminopolycarboxylates out. Previously inspected residuals mean this
challenge is not pristine prospective validation. Outer folds identical by model;
M0–M2 full dataset and M0–M3 identical size-fit-success subset both reported.

Geometry protocol, independent of selectivity labels: RDKit ETKDGv3, four embedded
conformers, seed 20261030, one thread; MMFF94 on the deprotonated organic ligand
only. Dummy center has no atom/metal force-field parameters; optimize its position
and organic coordinates with analytic radial harmonic restraints. La-like CN9
Shannon radius 1.216 A, Dy-like 1.083 A (Shannon 1976,
https://doi.org/10.1107/S0567739476001551); empirical donor offsets 1.40 A for O,
1.55 A for N give targets O 2.616/2.483, N 2.766/2.633 A. Offsets are documented
radius-sum assumptions, not fitted metal–ligand parameters or crystallographic
bond-length validation. Up to eight selected donors; for over-dentate ligands use
at most eight deterministic subsets from combinations of eight potential donors,
with identical assignments and initial conformer/center for both radii. Disconnected
donor clusters, ambiguous charge/tautomer/relative stereochemistry, absent MMFF
parameters or failed minimization produce missing values, never zero substitution.
Select the minimum organic+restraint objective independently at each radius from
the same candidate pool, then difference ORGANIC energies only. Save every
candidate, objective, convergence status, selected donor assignment and coordinates.
Primary k=20 kcal/mol/A^2, L-BFGS-B maxiter=1000, gtol=1e-4 kcal/mol/A,
ftol=1e-10. Bounded sensitivity: same pool, k=10 and 40 plus first-two-conformer
subpool at k=20. No extra conformers/restraint tuning after results. Descriptor
convergence coverage and between-setting stability accompany predictive results.

Published reference: pin Chemprop 1.5.2 and released upstream commit, run a small
released application before evaluation; M1 and combined M4 best-model five-member
ensembles as available, with memberwise Dy-La differencing at matched T and I.
Audit every released training/validation/test input including M3/M4 protonation
rows for normalized exact and stereo-agnostic matches, other metals/conditions,
family exposure and close fingerprint analogues. Released overlaps are references
only. Clean comparison requires wholly external families with verified absence
from training/model selection. If none are supportable, explicitly unresolved;
a multi-family full published-model refit campaign is outside this small iteration.
Do not fine-tune checkpoints or infer superiority from contaminated predictions.

Metrics: OOF RMSE/MAE, per-family counts/errors, per-ligand paired differences,
family-cluster bootstrap 2000 resamples, seed 20261030, 95% percentile intervals.
Numerical targets unchanged: family-held-out RMSE <=0.6 and tetraacetate challenge
RMSE <=0.8 kcal/mol. Incremental support additionally requires a positive 95% lower
bound for the paired RMSE improvement over the preceding ablation, with the same
sign under the bounded sensitivity. Targets and incremental evidence are distinct.
Select one final decision from the supplied list. Do not expand on a failed gate.

Execution: numerical work through Slurm standard_expand, bounded CPU allocations,
no --mem/--time overrides. Independent descriptor tasks can use process workers
inside one allocation. Record scheduler and measured process CPU costs including
failures. No new electronic/MD/alchemical/force-field-fit/ML-potential work.
Reserved SpyCI-LAMBS and PLM outcomes remain unopened. Optional peptide illustration
only if exact terminal chemistry and assay conditions are verified, otherwise defer.
