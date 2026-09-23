# Standalone xTB: fixed PQQ pools and analytic solvent derivatives

Approved by Jacob's Sept23 discretionary overnight instruction; root assigned
this exact contained comparison after native self-continuation settled energies
but failed the declared derivative checks. Agent water_basins owns this separate
backend diagnostic. No production change, optimization, new MACE, DFT, reference
fit or larger benchmark is included.

## Frozen source and number of calls

Reuse all80 exact cells from native_pool_continuation_20260923: 1H4I,4MAE,
Q88JH5 and A0A3F2YLY8 Ca-conditioned sample1, each with the same five geometries
origin/proposal_Ca/proposal_La/adaptive_Ca/adaptive_La, both metals and media.
Reuse all32 physical displacement cells in native_solvent_force_20260923:
4MAE Glu172 chi3 and Q88JH5 Glu221 chi3, ±0.001/±0.0005 radians, both metals/media.
Exact atoms, coordinates, caps, charge/singlet, water and chemical states stay
fixed. Both numerical levels receive exactly the same112 cells: **224 calls**.
Eight pool origin cells per level also request analytic gradients; these are the
shared derivative centers, not additional evaluations. No retries/third accuracy.

## Installed backend and fixed physical accounting

Use the existing `/home/jwestrob/miniconda3/envs/lanm_bench/bin/xtb`, version6.7.1
(edcfbbe), unchanged. Verify binary/package/parameter hashes. A private parameter
folder contains an exact copy of installed `param_gfn2-xtb.txt`; XTBPATH names
only that folder, excluding user rc/parameter overrides. Versioned upstream
calculator source confirms that this parameter file is loaded before fallback
to compiled parameters. Isolated fresh task folders and `--norestart` prevent
wavefunction reuse. No installation or fallback engine is permitted.

Use `--gfn 2 --chrg <exact> --uhf 0 --etemp 300 --iterations 500`, eight threads,
and `--acc 0.2` or `--acc 0.02`. Explicit ALPB(water) `gsolv` reference state and
normal230-point surface grid are common to both levels; no added salt. The
installed default starting guess, radii and solvent parameters remain unchanged
and are recorded from output/version. Only accuracy differs between paired
runs. Do not pass the unsupported/unused `--xparam` switch or tblite.

Take each backend's printed total, including its own electronic, dispersion,
repulsion, electrostatic, surface/H-bond and solvent reference-state terms.
`E = archived native OMOL(vac) + standaloneGFN2(ALPB) − standaloneGFN2(vac)`.
Never mix a standalone term with native ORCA's other medium. This is a new
composite descriptor; matching method names do not imply matching Hamiltonians.
The `gsolv` convention adds no gas/solution concentration correction. No aquo or
binding free energy is fabricated. Solvent components, charge, SCC thresholds,
occupations and any reported electronic-entropy terms stay available for audit.

## Frozen decisions and numerical gates

Always report0.02, with0.2 as sensitivity, irrespective of classification/energy.
All112 cells must remain in each level's denominator. Failures make dependent
quantities unavailable; no successful native or loose-accuracy substitution.

Energy gates: |E0.02−E0.2|≤0.1kcal/mol per cell; same-geometry Ca−La solvent and
pooled composite contrast changes≤0.2kcal/mol. Pool means the same five candidates
for both metals, mathematical row minima and separate unchanged0.1kcal/mol
origin-retention rule. Missing required cells invalidate the whole case pool.
Report raw matrices even if numerical qualification fails.

For the four components, two metal solvent transfers and one Ca−La solvent
contrast at each derivative source, use centered differences h=0.001 andh/2.
Both refinement and analytic-vs-fine error must be≤max(0.2kcal/mol/radian,
0.05|fine derivative|). Project actual analytic Cartesian gradients through the
same source/cap Kinematics Jacobian atq0, converting Hartree/bohr once. Verify
order, coordinates and gradient-energy identity. Report both accuracies and the
existing moving-heavy RMS physical normalization. No numerical DFT gradients,
new step, smoothing or Hessian. An unavailable/mismatched gradient is explicit.

Assess utility alongside consistency: retain all four pooled contrasts/candidate
choices, original native primary and continued comparisons, known-class relative
ordering and old native/released-band transfer explicitly labeled incompatible
calibration transfer. No four-case calibration or claim of fresh blind accuracy.
Boundary candidates are finite search samples, not thermodynamic populations.

## Execution and documentation

One CPU-only64CPU/128GiB allocation, eight workers×eight threads, via existing
explicit-node GPU-partition sharing policy. No GPU. Record immutable224-task
manifest, attempts, raw outputs/gradients, receipts, wall/core time and failures.
Use pinned runner snapshot; report-only/parser fixes may consume existing output
without rerunning chemistry. Report actual consistency and useful ordering, cost,
commands and vault note. Any further experiment requires a distinct decision.

Primary verification: installed CLI/manuals and versioned source, plus
https://xtb-docs.readthedocs.io/en/latest/sp.html#accuracy and
https://xtb-docs.readthedocs.io/en/latest/gbsa.html .
