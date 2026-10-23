# Parallel tests of solvent response and reciprocal Dy accommodation

## Execution hold following strategy discussion

21 October 2026: Jacob requested a complete physical-model plan and will set a
new goal after reading it. **Both branches are paused and neither was submitted.**
The CPCM adapter is an unreviewed source draft with an inventory; CONFIG and READY
are unfinished. The Dy checkpoint was inventoried without creating an executable
bundle. Do not launch these branches from the earlier authorization below.
Read the [new selectivity plan](../lanm_selectivity_free_energy_20261021/PLAN.md)
for the proposed program and activation boundary. Existing calculation 16680
was left unchanged. The original finite plan follows for provenance.

## Original finite authorization

21 October 2026 UTC. Jacob explicitly authorized launching these two branches
with subagents, then returning to discuss the scientific mechanism. This is a
finite execution authorization, not automatic resumption of all earlier goals.
Root owns submissions and integration. No alternative MACE checkpoint is selected.

## Questions and interpretation

1. Does a separately evaluated solvent response support the deformation work
   assigned by the current fixed-force-field-charge ALPB contribution? Both
   existing MACE and PBE0 comparisons share that term, so their agreement does
   not validate it. At exact 7CCN La steps92→725, ALPB contributes −45.442067
   kcal/mol while the two vacuum models increase by about15–16kcal/mol.
2. Does completing Dy accommodation change the comparison when both metals can
   select from the same completed La and Dy geometries? The existing La725
   convergence does not establish an equivalent Dy outcome.

These tests separate solvent-model sensitivity from a sampling asymmetry. They
do not establish experimental affinity, physical state populations or a unique
cause of model failure. No model is fitted to a desired preference.

## Branch A: independent solvent response

Linnaeus implements; Arendt independently reviews; root executes. Freeze a
ten-cell inventory before any new native energy:

- Eight existing whole256 LBT cells:7CCO and7CCN, each before/after the previously
  referenced deformation, each with La/Dy. Reuse all eight admitted vacuum PBE0
  energies and existing fixed-FF ALPB records.
- Two additional existing La geometries: exact7CCN steps92 and725. Reuse their
  admitted vacuum PBE0 and fixed-FF ALPB records. This pair directly tests the
  late −45.442067kcal/mol solvent benefit.

New native calls are PBE0-D4/def2-TZVP CPCM(water) singlepoint energies on the
same physical atoms and coordinates. Preserve each parent's explicit basis/ECP,
charge, effective spin, proton inventory, grid, SCF controls and unseeded guess.
Removing EnGrad requests an energy-only endpoint; adding CPCM(water) changes
the solvent treatment. Do not import the current OMol-target recipe or its SCNL
settings. Pin the exact native ORCA6.1.1 water/cavity settings and verify them
from output, rather than assuming a literature dielectric value.

One one-CPU, zero-molecular-call preparation validates all ten source joins and
generated inputs. After its actual result passes, submit ten independent native
jobs with12MPI ranks each, one thread/rank, standard_expand, scheduler default
memory/time and allocation-derived MaxCore. No waiting collector allocation.
At most ten new solvent singlepoints; no new vacuum call, optimization, other
functional, numerical retry or automatic replacement for a failed endpoint.

Report solvent transfers, solvent changes along each deformation, total PBE0
solution work, and matched Dy-minus-La work where both metals are represented.
Retain all missing/failed cells. Differences from fixed-charge ALPB combine
electronic response, continuum/cavity and other solvent-model differences;
they cannot uniquely identify polarization error. The two late La cells give
La-path evidence only. Preserve the existing1kcal/mol reporting resolution,
which is not an uncertainty estimate or a new admission threshold. No solvent
forces or optimized solution minima are qualified by this energy-only panel.

## Branch B: reciprocal Dy accommodation

Harvey implements; Arendt independently reviews; root executes. Continue the
exact saved7CCN Dy own-metal full-coordinate search checkpoint, preserving its
history, evaluator, physical state, topology/admission guards and charges.
Retain whole MACE-OMOL100M plus the existing fixed-FF ALPB model. No reset to the
La geometry or replacement electronic model is allowed.

Use the same full-system force criteria as La: maximum0.10 and atom-vector
RMS0.03kcal/mol/Å. Bound cumulative work at1172accepted steps/2172evaluations,
accounting for historical work before assigning the new allowance. Stop at
convergence, a physical/numerical failure or the finite budget; an incomplete
relaxation remains incomplete. Do not extend automatically.

After a one-CPU zero-call preparation and actual review, use one H200 allocation
with12CPUs and a warm pinned model. Cross-score the resulting Dy endpoint and
the preserved La725 endpoint under both metals, reusing exact compatible cells.
Report source-compatible common-pool selections, component work, residuals and
all failures. This tests conditional mechanics under the current model; known
reference-force and solvent limitations prevent a physical-affinity claim.

## Execution and preservation

Both branches are independent of running16680, which retains its existing
completion/health watchers. Simultaneous maximum for the new molecular work is
132CPUs and one GPU, plus16680's24CPUs, within the observed standard_expand
500CPU/user and20job/user limits. Per-job limit is24CPUs; no memory/time override,
priority change or interference with another user's allocation is authorized.

All numerical science, including real-input checks, runs on Slurm. Source,
metadata, text/hash inspection and orchestration stay on login. Use fresh
versioned workspaces, preserve failures, record actual costs, arm tested
completion and alert-only native-health notifications, and avoid duplicate jobs.
No automatic production promotion or further campaign follows these tests.

Root will report the actual submitted jobs and return to the requested
bird's-eye discussion rather than making a new experiment depend on16680.
