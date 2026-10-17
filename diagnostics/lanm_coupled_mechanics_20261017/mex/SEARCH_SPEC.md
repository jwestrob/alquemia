# Mex full physical response search

Module 2 implementation, declared before any new search. Root owns binding and
submission. All numerical preparation, qualification and optimization use Slurm.

## Prerequisites and unchanged model

Bind the exact completed MEX_CARTESIAN 34-composite qualification result, its
config/source, the common CARTESIAN_OPTIMIZER source, both exact old terminal
starts and this specification. Its 10 central-replay/directional checks must all
pass under the frozen numerical criteria before a new Mex search can run. The
parent plan additionally gates the main six-search comparison on the actual LBT
scout. No production change is made.

La and Dy each start from their admitted restricted-search terminal. All 2,049
physical atoms may move: backbone, sidechains, target, 167 waters and three Nd
spectators. Each own-metal search uses the same inherited electronic/MM/solvent
energy and analytic total physical gradient. No source-displacement envelope,
rigid donor groups or reference-geometry springs are introduced.

## Physical and state admission

The full native ff19SB/TIP3P XML supplies equilibrium covalent lengths, including
water O-H bonds. Require an exact join between all native force bonds and the
exported physical bond graph before any search. Native equilibrium lengths set
broad 0.65–1.35 bond admission intervals, enlarged only enough to retain any
pre-existing source value. This admits bond stretching and does not freeze
inherited lengths. Protein stereocenters and peptide cis/trans identity follow
the actual source. All physical pairs receive the common severe-contact check,
with bonded and 1–3 exceptions; inherited close contacts may not worsen beyond
the declared roundoff allowance.

Freeze each site's protein donor inventory from original source N/O atoms
within 3.5 angstrom of its ion. Freeze its retained-water inventory from original
water oxygens within 3.5 angstrom. Include the target and all three spectator
sites. This is a deterministic state bookkeeping rule, not a binding restraint
or an inferred coordination population. A site departure or loss of one of these
waters is labelled at the shared optimizer's 4.5 angstrom threshold. Such a
physically finite trial is evaluated and retained. A selected departed state
terminates the bound-state search and cannot enter its accommodation pool.
Coordination-number change alone does not fail a trial.

## Search and checkpoints

Use the common limited-memory Cartesian optimizer and its frozen step/Armijo
settings. Remove only instantaneous global translation/rotation from search
directions. Full raw physical forces determine convergence: maximum atom vector
<=0.10 kcal/mol/angstrom and atom-vector RMS <=0.03 kcal/mol/angstrom. Each of the
two searches has at most 80 accepted steps or 200 distinct composite evaluations,
including resumed work. Allocation interruptions preserve explicit Cartesian
checkpoint/history and the immutable component receipts; resume does not invoke
the old q7 admission path or reset the cumulative search budget.

Budget exhaustion and failed steps are explicit nonconvergence. No minimum,
affinity or biological preference is claimed from successful execution. Root
cross-scores both admitted terminal geometries in the common compatible source
pool after the searches; this driver never selects a favorable source or hides
a failed terminal.

## Configuration and command

A config contains `protocol=mex_coupled_cartesian_search_v1`, exact `scope` from
MEX_SEARCH.py, and pins `implementation`, `plan`, `specification`, `evaluator`,
`evaluator_config`, `optimizer`, and `qualification_result`. Optional `checkpoint`
is a pinned shared-optimizer checkpoint; null means a fresh start. Preflight may
bind qualification_result=null while the separate qualification is running;
`run` requires a complete admitted result. No metadata pin bypass is permitted.

```
PINNED_MACE_PYTHON MEX_SEARCH.py --stage preflight --config CONFIG.json --output FRESH_PREFLIGHT
PINNED_MACE_PYTHON MEX_SEARCH.py --stage run --metal La --config CONFIG.json --output FRESH_LA_RUN
PINNED_MACE_PYTHON MEX_SEARCH.py --stage run --metal Dy --config CONFIG.json --output FRESH_DY_RUN
```

Use the same pinned MACE environment as the evaluator, one allocated GPU and
scheduler-compliant CPU resources. Root arms actual terminal/health wakes and
records measured scheduler CPU/GPU cost. No new calculation is submitted by
writing this adapter or config.
