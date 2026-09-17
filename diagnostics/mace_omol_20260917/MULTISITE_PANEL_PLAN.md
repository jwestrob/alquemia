# Whole-chain multisite extension: parvalbumin and aequorin

Declared before new preparation or MACE inference, 2026-09-17. Jacob's blanket
authorization and active discriminator goal apply. The failed expanded GGR
test remains failed. This experiment adds family/site coverage; it does not
change that test, the model, or any production baseline.

## Inputs and physical model

Use all sites already mapped in the frozen
`diagnostics/nonpqq_direct_site_benchmark_20260915/panel_manifest.json`:
carp parvalbumin 4CPV ordered [CD,EF], and aequorin 1SL8 [EF1,EF3,EF4].
Use their archived pH7 protonated chain-A sources and raw deposited structures.
These structures/labels were consumed by earlier DFT work. New MACE predictions
are uninspected, but this is retrospective method development, not blind testing.

Retain every source chain-A calcium ion: two in 4CPV, three in 1SL8. Swap only
the named selected site between Ca and La; all background ions remain Ca2+.
Retain the union of previously selected site waters identically for every site
of each protein: 4CPV HOH166, 1SL8 HOH716/698/680. Other crystal waters remain
explicitly excluded. No water is added or selected from the new scores.
Record actual positions/protonation/atom inventories and physical charges.

4CPV has an actual covalently attached ACE0 N-terminal acetyl group. Preserve
its three deposited heavy atoms, archived hydrogens and ACE-C--ALA1-N bond;
validate the supported ff19SB template before inference. This is not a synthetic
QM cap. Do not strip ACE or replace the N-terminal chemistry. 1SL8 lacks ten
N-terminal residues; retain the observed segment with its archived terminal
protonation and report the omission. No internal gap may be bridged.

Use the existing ff19SB radial hydrogen normalization and terminal OXT policy;
preserve all observed heavy coordinates. Keep supported disulfides and reject
unknown cofactors/modified residues. Verify raw-source heavy-atom mapping before
filtering and check the complete retained ion/water inventory for every site.
New versioned multisite preparation policy; previous preparations unchanged.

## Fixed scoring and interpretation

Same masked OMOL checkpoint, float64, singlet, qualified 1024-edge/product
execution and factorization V2. Ten forwards: five selected-site La/Ca bound
pairs. Keep selected metal last in each source-mapped coordinate file for the
existing worker. Although all-Ca physical systems agree within each protein,
the atom order differs; do not reuse them without a qualified permutation check.
Compare all-Ca energies within each protein as an actual permutation/inventory
check, tolerance 0.01 model kcal. All endpoint accounting must pass 0.01.

Report ordered vectors, physical state and all unavailable cases. No canonical
PQQ bands or universal affinity zero apply. Parvalbumin is supporting cross-study
La-favoring evidence, not a same-assay gold label. Report both site scores minus
each of the three GGR scores: a supporting all-case direction check requires
all six differences >0.02 model kcal. No sub-kcal CD/EF affinity-order target.
These are two parvalbumin sites, one protein-level evidence group.

Aequorin's weak protein-level Ca-over-La comparison is not site resolved.
Report [EF1,EF3,EF4] without site labels, choosing a favorable site, or a binary
pass claim. These one-at-a-time substitutions with background Ca are not a
reconstruction of a cooperative/global metal titration. No average or fitted
threshold will be used to claim assay reproduction.

## Resources and failure handling

Use the existing prepared-input runner, OpenMM driver and separately pinned
MACE workers on one A5000,16CPU,64474MiB. Measured GGR bound forwards take about
13 seconds at4633atoms; these smaller proteins should require minutes for the
ten calls. Record actual total allocation and local preparation costs. No
DFT, solvent correction, gradients, relaxation, training or new environment.
Preparation failures remain explicit; do not remove troublesome ions/cofactors
or alter the rule after reading scores. Preserve frozen manifests and receipts.
