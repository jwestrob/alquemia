# Whole-protein simultaneous metal substitution

Declared2026-09-18 under Jacob's active autonomous MACE discriminator goal,
before new collective endpoint outputs. Run alongside unchanged coupled-path
DFT validation, without disturbing its jobs, records or criteria.

## Question and fixed cases

Does the current global masked MACE representation predict non-additive effects
when ALL modeled sites change from Ca to La, and does this improve domain-level
ordering over the mean of single-site substitutions with Ca background?
Use exactly the previously scored whole A0A7/HEW5/RTX author domains, parvalbumin
4CPV, and aequorin1SL8, with6/8/8/2/3 modeled metal sites respectively. Choose
the first site in each previously frozen site order solely to establish one
atom ordering. Retain every source atom, H/protonation state, cofactor exclusion,
explicit water, assembly and coordinate. No new structure, fold, carver or
geometry optimization. These are all consumed development cases.

A0A7/HEW5/RTX retain pH6 author Ca-conditioned models without experimental MPVP
scar. GGR references retain their original preparations including pH differences.
Protein-level La ITC versus Ca CD-folding comparisons remain cross-readout
supporting evidence, not same-assay affinity labels. Parvalbumin remains
supporting cross-study evidence. Aequorin remains a weak protein-level Ca
direction with unresolved individual sites. No individual label is invented.
The fully occupied model state is a computational counterfactual; experimental
stoichiometry, titration curves and cooperative free energies are not reproduced.

## Energy expression and invariants

Let n be the recorded metal-site count, C the all-Ca energy, L the all-La
energy, L_i each archived single-La/Ca-background energy, and A the same
archived isolated-atom Ca-minus-La reference in model eV. Define

    R_joint_per_metal = ((C-L)/n - A) * EV_TO_KCAL
    R_single_mean = ((n*C-sum(L_i))/n - A) * EV_TO_KCAL
    K = L + (n-1)*C - sum(L_i)
    R_joint_per_metal - R_single_mean = -K/n * EV_TO_KCAL.

Keep unrounded energies, n, all site values and K. This measures fixed-geometry
non-additivity of an energy-like descriptor, not biochemical cooperativity or
binding free energy. Same model/checkpoint, charge-feature mask, exact batched
native evaluator and float64; no added continuum, electrostatic, entropy or
mechanical term. All-Ca and every L_i come from matching actual archived receipts.
Total charge rises by n in all-La; verify every explicit metal index, elemental
inventory and even all-electron parity. No global neutralization or guessed
charges. The established one-selected-metal default guard remains unchanged;
collective state validation belongs only to this new opt-in protocol.

## Finite calculations and predeclared criteria

Exactly8 new scalar forwards: one all-La state for each of5 proteins; RTX
all-Ca replay; RTX all-La rigid rotation (the existing37degree rotation about
[1,2,3]); RTX all-La reversed atom ordering. Same source geometry in every
physical state; transformation indices must map exactly. Native readout closure,
replay, rotation and permutation tolerance0.01 model kcal. Direct algebra
closure1e-6 model kcal. Existing one-site GGR endpoint algebra must reproduce
the collective expression at n=1; no new GGR inference.

Compare all3 author-domain joint scores with all3 retained GGR scores:9 margins,
all must exceed0.02 for this declared supporting challenge to pass. Report the
old mean and joint score side by side, all failures and denominators;3 domain
groups versus1 GGR group, not9 independent observations. For parvalbumin report
3 supporting joint-score/GGR margins at the same0.02 rule, retaining both old
site scores. For aequorin report the joint descriptor and old ordered vector,
without a site label or new absolute threshold. No PQQ band or universal zero.
No alternative aggregation, site subset, microstate or cutoff after outputs.

Use existing run_pilot.sbatch, A5000/16CPU/64474MiB. The prior10 multisite
forwards used216GPU-allocation seconds; this pilot is on that minute-scale,
with actual cost recorded. No project CPU/time cap; finite inventory8, no DFT,
solver, trajectories or fitting/training. Preserve baseline/default and old
studies. Export actual results, runnable commands, tests and vault note. A
failure does not justify adding variants to this frozen experiment.
