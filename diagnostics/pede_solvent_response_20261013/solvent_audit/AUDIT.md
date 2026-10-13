# PedE solvent response: implementation audit and discriminating test

The current implementation forms the matched solvent correction correctly as
`E_ALPB − E_vac`, using separately converged endpoints at the same coordinates
and state. No demonstrated duplicate solvation term, silent charge change, or
recurrence of the native-ORCA La shell-order defect has been identified. The
important unresolved question is the **physical response of the charged capped
context**, not whether the completed scalar jobs returned success.

## Fixed implementation and state

- Actual primary: standalone xTB6.7.1; independent component check: tblite0.7.0.
  This is not the older native-ORCA GFN2 implementation. The used full parameter
  file is `/home/jwestrob/miniconda3/envs/lanm_bench/share/xtb/param_gfn2-xtb.txt`,
  SHA256 `f6f05c6c5264cb2f75c54f8ac6838044541c1e493b16f1a4d4ee4fa64de5f545`.
  Its La block declares `5d6s6p` and explicitly named KCNS/KCNP/KCND and
  POLYS/POLYP/POLYD fields. No reordered native-ORCA JSON is substituted into
  this route. Prior actual backend comparisons remain the executable evidence.
- Commands explicitly pass charge−3 for Ca,−2 for La, `--uhf 0`, electronic
  temperature300K, accuracy0.0001,500 SCC iterations and no restart. The shared
  finite192-atom state has588 valence electrons for either metal, complete PQQ,
  fourteen synthetic caps and no added solvent/protons. Saved preparation maps
  include every atom; physical motion coordinates alone are not the QM mapping.
- The source formal-charge ledger includes Asp220, Glu221, Asp358 at−1 each and
  Arg386 at+1; the remaining listed protein groups are neutral. PQQ is the fixed
  trianionic state and the target metal is Ca2+/La3+. This bookkeeping does not
  constrain the SCC charges of individual disconnected fragments.
- Actual water setup: internal GFN2-xTB/ALPB, dielectric80.2, P16 interaction
  kernel, GBOBC Born integrator, radius scale1.4744, zero offset, hydrogen-bond
  correction on, ion screening off,230 surface grid points per atom. The solvent
  thermodynamic setup prints298.15K, separately from300K electronic smearing.
  Neither temperature is a sampled protein/nuclear temperature.
- The command's `gsolv` token triggers an unknown-reference warning. The same
  output then explicitly prints the effective `gsolv [1 M gas/solution]` reference
  and Gshift0.001080759698Eh. It is a real syntax warning, not evidence of a
  source-dependent reference change. A shared constant cancels in within-metal
  source differences. Preserve the warning and effective setup together.

These facts refer to the actual nine-source PedE records. They do not validate
all parameter values for all metals, infer an oxidation state from fractional
population charges, or establish that PQQ's chosen electronic state is correct.
The separate cofactor branch owns the source/redox evidence.

## What the existing correction contains

The saved native summaries expose `Gelec`, `Gsasa`, `Ghb` and `Gshift`, whose sum
is `Gsolv`. That value is already included in the total energy, so the current
`E_OMOL,vac + E_ALPB − E_vac` expression does not add it twice.
[Official ALPB documentation](https://xtb-docs.readthedocs.io/en/latest/gbsa.html).

The solvent is self-consistent: its charge-dependent potential enters the SCC
Hamiltonian. In the pinned6.7.1 source, `getEnergyParts` computes electrostatic
energy from the Born matrix and current atomic populations; the hydrogen-bond
term itself is quadratic in those populations. Thus a changed Ghb can reflect
charge redistribution as well as accessibility. It is not a geometry-only
nonpolar term.
[Born implementation](https://raw.githubusercontent.com/grimme-lab/xtb/v6.7.1/src/solv/gbsa.f90).

The intrinsic SCC energy also contains density-dependent dispersion and
electronic-temperature entropy. The separate entropy contribution is not printed
as a scalar in the inspected endpoints. A Fermi level is not an entropy energy.
No missing entropy term is replaced by zero or interpreted as nuclear entropy.
[SCC source](https://raw.githubusercontent.com/grimme-lab/xtb/v6.7.1/src/scc_core.f90),
[documented smearing convention](https://xtb-docs.readthedocs.io/en/latest/sp.html#fermi-smearing).

Consequently, the exact bookkeeping identity is

`U = Gelec + Gsasa + Ghb + Gshift + I`,

where `I = E_ALPB − Gsolv − E_vac`. I is the entire remaining intrinsic-energy
change. Removing I or substituting printed Gsolv for U would alter the model;
the decomposition is not permission to do so. Nor does it separate cavity
geometry from charge response inside Gelec without additional information.

## Ranked hypotheses before the new decomposition

1. **Source-sensitive self-consistent electronic response is the leading model
   hypothesis.** The prior saved analysis localizes the21.106kcal origin contrast
   mostly to U and finds different PQQ/residue population redistribution. Small
   gaps make charge response plausible. These observations do not prove an
   erroneous branch, and similar features occur in correctly classified controls.
2. **The capped context may amplify a real geometry dependence.** Missing outer
   protein dielectric/mechanical context, exposed truncations and unconstrained
   interfragment charge transfer can alter electronic response. All sources use
   the same membership/caps, so an inconsistent atom inventory is not the
   explanation; the adequacy of the common representation remains unqualified.
3. **Surface/hydrogen-bond terms or numerical cavity sensitivity could contribute.**
   These are not yet separated from Gelec in the published raw spread. The saved
   decomposition is the minimal discriminator; no solver or grid sweep is needed
   to discover which terms carry the contrast.
4. **Reference, sign, unit or failed-SCC artifacts have lower support.** Exact
   algebra, raw charge/state checks and tight SCC admission passed previously;
   the known reference warning has a common printed effective convention. Retain
   failures and absolute-v1 parity evidence without treating them as causal
   explanations for this source-dependent contrast.

## Minimal runnable discriminator and decision

[DECOMPOSE.py](DECOMPOSE.py) and [PLAN.md](PLAN.md) read only36 saved outputs for
all five two-Ca and four compatible La-conditioned origins. The frozen configuration
is `workspaces/pede_solvent_response_20261013/solvent_audit_v1/CONFIG.json`.
Root runs it on one standard Slurm CPU. All terms, accounting residuals, effective
setup and missing entropy status are retained. No molecule or population is
recomputed. Actual results will be reported separately when the job completes.

If surface/reference terms are negligible, rule out that simple explanation and
compare the metal-dependent electronic/continuum response with the independently
prepared DFT reference. If a specific surface term dominates, investigate its
geometry/parameter definition before changing the electronic method. Either
result directs one targeted repair rather than an alternative-solver sweep.

No inexpensive corrected scalar is justified yet. In particular, a label-fitted
shift, deleting solvent, freezing favorable charges, or suppressing the intrinsic
response is not a repair established by these observations. This branch delivers
a runnable discrimination between mechanisms; the paired independent reference
is the qualification path for a physically supported correction. Production and
all calibrations remain untouched.
