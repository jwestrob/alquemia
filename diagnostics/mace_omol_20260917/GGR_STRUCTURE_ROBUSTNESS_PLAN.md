# Frozen masked-MACE test across the existing GGR crystal structures

Declared 2026-09-17 before preparing/scoring the two additional structures.
Jacob's active goal and blanket research authorization apply. This is a new
retrospective structural-robustness test; prior experiments stay immutable.

## Question and fixed inputs

Does the successful masked-MACE alpha/GGR ordering survive the two other
already consumed GGR crystal geometries, 2FW0 and 2FVY? The earlier ORCA study
showed source-associated shifts near 10 kcal/mol. Repeating more 1GLG calls
cannot answer this question. The new MACE results have not been inspected.

Use both archived chain-A protonated sources and repaired manifests under
workspaces/ggr_mechanism_20260915/stage_b_prepared_v1/{2fw0,2fvy}/.
Their exact source, selections, preparation, coordinates and evidence are pinned
before inference. Preserve all observed standard-protein heavy coordinates and
protonation identities. Apply the same ff19SB radial H-length normalization and
observed-terminal OXT policy as the qualified whole-chain model. Require valid
1.0–1.8 A peptide C–N connectivity, forcefield template closure, exact La/Ca
paired coordinates, recorded physical charges and singlet parity. No geometry
optimization, missing-residue construction or donor modification.

Both structures lack terminal residues A1 and S307/K308/K309 in their archived
sources. Keep those omissions explicit; model the observed polymer segment
with its recorded terminal chemistry. No internal gap may be bridged. Record
terminal distances from the selected site before scoring. This comparison is
not a controlled causal effect of sugar or a pure geometry-only perturbation.

Use the same declared protein-plus-selected-metal, zero-water representation
as GGR1GLG. Enumerate excluded sugars, citrate/malonate, sodium, glycerol and
source waters explicitly. Those species are not silently charged or added to
the input. 2FW0 is sugar-free/open; 2FVY glucose-bound/closed and has reported
Glu149 radiation damage; both remain qualified structural observations of the
same protein, not new independent affinity measurements. Any actual unsupported
polymer chemistry or incomplete nearby heavy atoms blocks preparation.

## Model, inventory and tests

Keep protocol mace_omol_intact_charge_feature_ablation_descriptor_v1, its pinned
100M OMOL checkpoint, float64, singlet embedding, zero total-charge feature,
exact 1024-edge/product execution, and fixed factorization reference V2.
Actual physical charges remain recorded; these are learned descriptor outputs,
not quantum electronic-state energies or binding free energies.

Four new forwards total: the bound La and Ca states for each of 2FW0 and 2FVY.
Use the two-call formula and already qualified fixed disconnected model terms.
Reuse actual 1GLG and alpha1F6S/6IP9 descriptor scores, with no recomputation.
Require existing component-accounting tolerance <=0.01 model kcal. Additional
rigid/repeat qualification is not needed for an unchanged execution method.

Report all three GGR scores in frozen order [1GLG,2FW0,2FVY], their range, and
both alpha scores minus each GGR score. Structural robustness passes only if
all six alpha-minus-GGR margins exceed the same pre-existing 0.02 model-kcal
threshold. No averaging away a failing structure, best-structure selection,
new affinity zero, PQQ bands, fitted coefficients or changed preparation after
seeing a score. Missing cases remain in the denominator and block an all-case
robustness claim. The prior successful development result remains as observed.

This adds structural coverage, not an independent biological pair: GGR remains
one direct Ca-favoring observation, alpha one qualified La-favoring observation.
They have different assay conditions. No broad affinity validation is claimed.

## Implementation and resources

Add a small source-row bridge to the existing whole-chain preparation function,
with explicit JSON source/evidence/expected-water inputs and source pins. It
must reproduce an existing real1GLG preparation exactly before new inference.
Use existing prepared-input scorer and executor, one A5000,16CPU,64474MiB.
Four forwards should fit a short allocation based on measured existing GGR
calls; record actual preparation, controller, model, failures and Slurm costs.
No project CPU/time stopping budget; scheduler QOS applies. No new DFT, solver,
training, gradient, dynamics, production rescore, default promotion or job
interference. Native H2001200809 and other users' jobs remain untouched.
