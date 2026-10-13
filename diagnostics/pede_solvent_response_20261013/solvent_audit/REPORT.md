# PedE: solvent sensitivity is electrostatic, with large intrinsic compensation

The saved decomposition rules out a large surface-area or reference-offset
contribution to the previously identified PedE origin spread. The leading
question is the coupled electronic/continuum response of the fixed capped
context. No new parameter or accounting defect was demonstrated, and no
corrected score is justified by this decomposition alone.

Slurm3352 completed the declared36 saved endpoints: all five recovered two-Ca
sources and all four compatible historical La-conditioned sources, both target
metals, vacuum and water. All printed energy and component closures passed the
predeclared1e-9Eh tolerance. No molecular evaluation, population recomputation,
optimization or calibration ran. The original classification failures remain.

## What carries the origin contrast

Use the established sign `U_Ca − U_La`, where each
`U_M = E_GFN2,ALPB,M − E_GFN2,vac,M`. The already identified comparison is
source sample3 (correct Ca) to sample1 (wrong La after accommodation). It is
reported alongside every source, not used to filter the population.

| Saved term | Ca change | La change | Ca−La change |
|---|---:|---:|---:|
| Electrostatic Gelec | −210.973412 | −249.328633 | +38.355221 |
| Surface Gsasa | −0.069587 | +0.127283 | −0.196870 |
| Hydrogen-bond Ghb | −9.915747 | −11.618297 | +1.702550 |
| Reference Gshift | 0 | 0 | 0 |
| Intrinsic remainder | +164.229447 | +182.833719 | −18.604272 |

All values are kcal/mol. The component sum reproduces the saved solvent
double difference **+21.256628kcal/mol**, with residual−6.2741e-10kcal/mol.
The previous OMOL vacuum contribution to this same origin difference was
−0.150413kcal/mol. Thus the solvent correction carries the large origin spread;
it is not explained by unusually large accommodation in the wrong sources.
The separate calibration analysis already showed that correct sources can have
comparable raw accommodation shifts, while the calibrated bands also move.

For each endpoint the exact bookkeeping is

`U = Gelec + Gsasa + Ghb + Gshift + I`,

`I = E_ALPB − Gsolv − E_vac`.

The large opposing Gelec and I changes are important: the observable correction
is their compensated total, not printed Gsolv alone. I is the full intrinsic
energy change outside Gsolv, including any electronic-temperature and
density-dependent terms. It is not independently identified as pure
polarization, charge-transfer energy, or a removable error.

## What the mechanism does and does not establish

The implementation feeds the solvent potential back into the SCC electronic
calculation. Gelec uses the converged charges, and Ghb itself depends
quadratically on atomic charges. The split therefore does not separate a
changed dielectric cavity from a changed electronic distribution inside Gelec;
it also does not make Ghb a geometry-only nonpolar term.
[Pinned xTB6.7.1 implementation](https://raw.githubusercontent.com/grimme-lab/xtb/v6.7.1/src/solv/gbsa.f90).

The existing population analysis independently found substantial PQQ and
protein-fragment charge redistribution. Together, these observations support
source-sensitive electronic/continuum response as the leading mechanism to
test. They do **not** prove an erroneous electronic branch, incorrect oxidation
state or physical charge transfer. Strong redistribution and small orbital gaps
also occur in correctly classified controls. Agreement of two implementations
and SCC convergence qualify reproducibility of this model, not its accuracy.

The common192-atom composition, fourteen caps, Ca−3/La−2 charges,588 valence
electrons and medium-matched coordinates were retained. These invariants reject
an accidental change of molecular inventory. They do not establish that the
capped context is large enough or its interfragment electronic response is
physically appropriate. The independent cofactor audit owns redox/template
interpretation; crystal geometry is not assumed to be an oxidized minimum.

The raw command warns `Unknown reference state 'gsolv'`, but the actual setup
prints `gsolv [1 M gas/solution]` and the same Gshift0.001080759698Eh. The observed
reference contribution to this comparison is exactly zero. Preserve the syntax
warning without attributing the source dependence to it.

Printed electronic temperature is300K; the solvent thermodynamic setup is
298.15K. Separate electronic entropy is not printed as a scalar and remains
unavailable. The electronic smearing contribution belongs to the SCC energy;
it is not nuclear/configurational entropy.
[SCC source](https://raw.githubusercontent.com/grimme-lab/xtb/v6.7.1/src/scc_core.f90),
[documented smearing convention](https://xtb-docs.readthedocs.io/en/latest/sp.html#fermi-smearing).

## Decision

Proceed with the already prepared independent fixed-origin electronic reference,
comparing each metal's solvent response and the source double difference. That
will test whether a different electronic/continuum treatment reproduces this
sensitivity. Because it changes both electronic treatment and solvent model,
disagreement alone will not uniquely identify which component is responsible.

Do not replace U by Gsolv, drop the intrinsic remainder, remove solvent or apply
a fitted offset. The accounting is already coherent. A useful correction needs
evidence about the metal-dependent response, followed by transfer checks on the
retained reference population. No alternative-solver sweep is launched here.

## Artifacts and measured cost

- [Compact pinned evidence](EVIDENCE.json) retains all nine sources and both
  metals' decompositions. Exact values remain in the immutable result.
- Result: `workspaces/pede_solvent_response_20261013/solvent_audit_v1/execution/RESULT.json`,
  SHA256`b48d8dec8bd5895df435fa12dcf2f9d801bc105aa52be126249d158d344f052d`.
- Full table: the adjacent`COMPONENTS.tsv`; raw setup, warnings and nested SCC
  components remain in the JSON. No missing component was filled with zero.
- Analysis elapsed0.818110431s, processCPU0.703118s, peakRSS32608KiB.
  Scheduler:2s on1CPU,2allocatedCPU-s. Zero GPU and zero new molecular calls.
  This is saved-artifact analysis cost, not molecular inference latency.
- [Plan and command](PLAN.md), [implementation audit](AUDIT.md),
  [Slurm-only analysis](DECOMPOSE.py). The completed output directory must not
  be rerun or overwritten. Root owns the separate reference submission.

Production, thresholds, old absolute-parity failure and all wrong/inconclusive
results remain unchanged.
