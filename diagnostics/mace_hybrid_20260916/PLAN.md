# MACE / DFT hybrid discriminator: proposed experiment

Date: 2026-09-16. **Direction approved; planning recorded; no new calculations launched.**

Jacob: “i like it. i want to try it out. write down the plan, then real quick
i just want to make sure our baseline results are as good as I want them to be
before I go blasting using it on everything, since I've used it on several
XoxF-like proteins from the PLM hillslope and they've come back inconclusive
so far.” This authorizes documenting the hybrid direction and reviewing existing
baseline results now. Baseline review takes priority. The proposed run counts
and unresolved scientific choices below have not yet been discussed as an
execution manifest. Do not launch them from this note.

## Objective and architecture

Test whether a pretrained, responsive whole-protein energy model plus a local
DFT correction improves affordable La/Ca discrimination. Hybrid scoring is the
main track. Direct MACE scoring is a comparator from the same evaluations.
Structural response is a complementary, subsequent track.

Proposed engine: **MACE-POLAR-1**, initially the medium checkpoint, with explicit
long-range electrostatics, known total charge and spin. Pin the actual checkpoint,
code, dependencies, precision and electrostatic boundary settings before running.
Use an isolated environment. Do not replace existing environments or defaults.
No training or fine-tuning is proposed in the initial experiment.

For metal M, the vacuum additive/subtractive skeleton is:

```
E_hyb,vac(M) = E_MACE,vac(full, M)
             + E_DFT,vac(core, M) - E_MACE,vac(core, M)
R_hyb = E_hyb(Ca) - E_hyb(La)
```

MACE supplies internal scaffold energy and core–environment interactions; DFT
replaces its internal core approximation. Both core evaluations must use exactly
the same capped atoms, coordinates, charge and spin. The full structure must be
identical across numerical core partitions. Capping errors are not guaranteed
to cancel between different methods. This architecture still has a QM boundary.

An isolated DFT correction does **not** respond self-consistently to the protein
field. The coupling remains MACE-level in this initial architecture. Embedding
DFT in a field would be a separately specified model, not an implicit benefit of
the equation. Adding MACE energies to the old CPCM score is not this equation.

## Decisions to settle before execution

1. **Solvent:** MACE-POLAR has no automatic CPCM(Water) equivalent. Define a single
   coherent solution-phase expression and its reference. A whole-system solvent
   contribution must have an explicit density/charge source, cavity, dielectric,
   boundary convention and subtraction policy. No duplicate full solvation term.
   The archived PB/TABI model failed its gates and is not an accepted default.
   Vacuum engine tests below can establish limited capabilities, not affinity.
2. **Physical system:** freeze assembly, source/protonation, cofactor state,
   water inventory, global formal charge and spin. Audit compatibility of the
   existing full protein with each archived core; do not import cap nuclei into
   the physical full system. Fail missing/unsupported chemistry explicitly.
3. **Implementation:** checkpoint hash and element support, isolated software
   environment, nonperiodic electrostatics, precision, force units and mappings.
   Measure memory and inference time on the actual system; no speed claim yet.
4. **Acceptance:** freeze reproducibility/rigid-transform tolerances, partition
   tolerances and any force-comparison tolerances before those outputs are seen.
   Numerical tolerances and scientific effect sizes are different requirements.

## Proposed staged experiment and run inventory

These are proposed selections/counts, not an approved task manifest. Record any
agreed revision before execution. Reuse existing runners, receipts and caches.

### A. Reuse the consumed 1H4I development artifacts

Question: can the proposed engine evaluate the actual full protein affordably,
and does the subtractive correction behave sensibly when Asp303 crosses the
numerical core boundary?

- Reuse the four archived native **vacuum** r2SCAN-3c endpoints for qm33/qm36 ×
  La/Ca from `diagnostics/global_electrostatic_20260916/` and its workspace.
  Their protocol differs from the canonical PQQ fixed core; preserve that fact.
- Evaluate one matched full-protein La/Ca pair and both core pairs: **six MACE
  energy/force evaluations, zero new DFT endpoints** if the archive is compatible.
- Repeat the full pair, rigidly translate the full pair, and rigidly rotate the
  full pair using transformations frozen in the manifest: **six additional MACE
  evaluations**. Total proposed stage A: **12 MACE evaluations**.
- Report direct and hybrid raw contrasts, separate correction terms, partition
  change, global charge closure, repeat/transform errors, memory and actual cost.
  This is consumed development evidence with no new biological accuracy claim.

Do not mix the archived CPCM gradients with vacuum MACE forces and call their
difference a model error. Comparisons need matched physical energy definitions.
Do not substitute invalid failed EDA fragment energies for archived valid DFT.

### B. Test predictive behavior after the solvent expression is agreed

Proposed development cases: **1H4I, 4MAE, GGR 1GLG, alpha-lactalbumin 1F6S**.
These span the successful PQQ transfer and known generic limitations. All have
already been inspected. The PLM candidates remain unlabeled applications, not
new positive controls inferred from their XoxF-like annotation.

For each site, preserve the same physical preparation across direct MACE and
hybrid comparisons. Four full-protein pairs plus four matching core pairs imply
**16 MACE energy/force evaluations**. At most **eight new DFT endpoints** are
proposed, reduced wherever an exact compatible archived endpoint exists. Any
additional solvent/reference tasks must be enumerated before this stage becomes
runnable; these counts are not an all-inclusive solution-phase execution plan.

Compare the unchanged baseline, direct model and hybrid side by side. Report raw
contrasts and between-site differences until a compatible new reference or
designated calibration exists. Do not inherit baseline bands or universal zero.
Keep PQQ functional-class evidence distinct from condition-qualified affinity
directions. No threshold fitting on these development outcomes. Report every
invalid/unscorable case and the preparation denominator.

The question is whether the additional computation improves useful predictions,
not whether it reproduces baseline energies. Let direct MACE outperform the
hybrid if that is what the evidence shows. Retain cheap composition/geometry
comparators in the report without making them a separate large project.

### C. Structural response, separately agreed

Once the fixed-geometry energy model is understood, test whether its forces can
support limited metal-dependent movement and scaffold strain accounting. Freeze
physical displacement coordinates and a trust region before evaluation. Verify
forces/curvature against compatible existing derivatives and a specifically
enumerated set of new DFT checks if needed. Initial run count: **zero scheduled**;
no local relaxation or entropy correction is enabled by this plan.

No synthetic cap motions, independent fragment translations or arbitrary springs
as protein degrees of freedom. Minima reached with a cheap model are not
automatically minima or equilibrated ensembles of the hybrid model.

## Other angles retained

- **AMOEBA with published La-specific parameters:** an alternative responsive
  full-protein engine. Check compatible Ca parameters, POLPAIR support and PQQ
  coverage before a separately agreed comparison. No MD/FEP campaign implied.
- **Extended GFN-FF:** possible inexpensive geometry/scaffold comparator; not a
  validated selectivity potential. Do not confuse it with the failed whole-
  protein GFN2 electronic calculation.
- **CHELPG:** optional charge output for models needing it. It does not improve
  the existing DFT score by itself and is not a prerequisite for MACE scoring.
- Chemical-state ensembles and evolutionary/structural evidence combinations
  remain possible later directions, requiring explicit scope and evidence.

## Deliverables and operating rules

A pinned finite manifest, isolated reproducible environment, actual receipts,
paired component table, physical checks, measured costs and an accuracy-focused
decision. Candidate products go under `workspaces/mace_hybrid_20260916/` when
authorized; plans/reports remain here. No CPU/time stopping budget is imposed.
No new calculations have run for this planning note. Preserve baseline/defaults,
old artifacts, released thresholds and concurrent PLM work.

## Sources and relation to the earlier proposal

- [MACE-POLAR documentation](https://mace-docs.readthedocs.io/en/latest/guide/polar_mace.html).
- [MACE-POLAR preprint, including lanthanide isomer tests](https://arxiv.org/html/2602.19411v1).
  Those tests do not validate protein La/Ca exchange or whole-protein affordability.
- [Released configuration including Ca/La](https://github.com/ACEsuit/mace-foundations/blob/main/mace_polar_1/config-mace-polar-1.yaml).
- [La AMOEBA parameterization](https://doi.org/10.1021/acs.jpcb.2c07237).
- [Extended GFN-FF](https://doi.org/10.1021/acs.inorgchem.4c03215).
- Existing vault notes `2026-03-31_mace-potential-grant-strategy.md` and
  `2026-04-01_mace-training-data-strategy.md` describe a broader trained-potential,
  within-lanthanide program. This proposal first tests released pretrained weights
  for La/Ca; it does not activate those historical training/AIMD plans or treat
  their unverified claims as current results.
