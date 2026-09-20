# Hydrogen preparation defect: reproduced and numerically repaired

2026-09-20. Research diagnostic only. No default, historical preparation,
all250 population, fold selection or molecular score changed.

## Read-only audit:17 failures and four preselected successes

The completed source campaign's17 failures have two distinct causes:

- **16 off-site predicted ions:** C5AXV8 five Ca samples; Q88JH5 five Ca
  samples plus La sample1; Q9Z4J7 five Ca samples. Every selected ion is
  about25 Å from the closest PQQ atom, outside the declared4 Å test.
- **One failed H preparation:** MMOL1770 La sample4, with1,759 sub-0.45 Å
  H/H pairs throughout the protonated protein,10 already within its core.

None of the other16 failures or four fixed successful controls has a
sub-0.45 Å pair. The controls were MMOL1770 La0/La3 and A0A3F2YLY8 Ca0/La0,
chosen before inspecting their H metrics. All H in all21 files have exactly
one heavy-atom topology bond. This rules out absent H-parent topology bonds
as the explanation for the observed collapse. It is not an all233-success
global H-quality audit.

## Exact replays establish silent, nonconverged minimizer termination

Frozen plan: [HYDROGEN_DIAGNOSTIC_PLAN.md](HYDROGEN_DIAGNOSTIC_PLAN.md).
Slurm array1203793 ran the failed La4 and successful same-protein La3 replays
independently, with the pinned seed/pH/software/CPU1 policy. Both reproduce
the saved PDB coordinates exactly. A passive reporter records:

| Quantity | Failed La4, native | Successful La3, native | La4 repaired |
|---|---:|---:|---:|
| Accepted native iterations | 1 | 50 | 50 |
| H/H pairs below0.45 Å | 1,759 | 0 | 0 |
| Final movable-H force RMS, kJ/mol/nm | 3.04739e13 | 37.6422 | 13.5032 |
| Maximum H displacement from initialization, Å | 0.1105 | 1.8730 | 1.7896 |
| Fixed heavy-atom displacement, Å | 0 | 0 | 0 |

RMS uses only positive-mass movable H; it excludes all4,895 mass-zero atoms.
All runs have zero constraints. The failed minimizer exits after one accepted
iteration despite a huge residual gradient; it did not reach its50-iteration
limit or convergence. Its closest initialized H pair, Met22 HE1/HE3, is
**0.000176367 Å** apart, and initial generic potential energy is1.03462e17
kJ/mol. The successful control's closest initial pair is0.00991554 Å.

Installed `Modeller.addHydrogens()` places H near the same outward parent
direction, adds explicit H-parent bonds, then uses50 L-BFGS iterations to
spread them. OpenMM8.5.1's [reference minimizer source](https://github.com/openmm/openmm/blob/8.5.1/platforms/reference/src/SimTKReference/ReferenceMinimize.cpp)
calls L-BFGS without examining its return value. The observed early failure
is consistent with line-search/numerical difficulty from extremely close
initialized H, but the precise discarded return code was **not measured**.
Missing bonds, cap construction and heavy-atom motion are not the cause here.

## Same-potential repair passes the declared geometry checks

From the exact captured initialization and System XML,14 deterministic
gradient steps of at most0.02 Å per movable atom removed the severe initial
overlaps. Every step decreased the original potential; no backtracking was
needed. The unchanged native50 minimizer then produced the repaired column
above. No force term, seed, atom, bond, pH or protonation state changed.

Both output topologies retain the same9,582 atoms/9,700 bonds. All4,895 heavy
PDB coordinates remain exact. All H have one parent, no sub-0.45 Å pair
remains, and the frozen0.8–1.5 Å H-parent/2.5 Å displacement guards pass.
This repairs gross geometry;13.5032 RMS exceeds the requested1kJ/mol/nm,
so the50-step result is **not strictly force-converged**.

### The long H bonds belong to this generic preparation potential

PDB-coordinate bond lengths, Å (all protein H; no source waters/PQQ H):

| Parent type (H count) | Successful La3 min / median / max | Repaired La4 min / median / max |
|---|---|---|
| C-H (3,606) | 1.174 / 1.190 / 1.223 | 1.172 / 1.191 / 1.222 |
| N-H (984) | 1.164 / 1.185 / 1.209 | 1.175 / 1.186 / 1.208 |
| O-H (97) | 1.159 / 1.183 / 1.194 | 1.171 / 1.182 / 1.194 |

The captured nonbonded force has **zero exclusions**, including bonded pairs.
For an isolated H-parent pair, its actual terms are
`0.5*100000*(r-0.1)^2 + 100*(0.1/r)^4`, with r in nm and E in kJ/mol.
Their minimum is1.17705 Å, not the nominal harmonic target1.0 Å. This explains
the long distances in both repaired and successful controls. They are not
validated physical C-H/N-H/O-H lengths. The repair is compatible with the
existing approximate H-placement policy; it is not an independently accurate
hydrogen Hamiltonian and is not a molecular energy score.

## Measured cost and evidence

Failed replay plus repair:18.7664 s process wall; success replay:11.7733 s.
Slurm elapsed21/16 s, total actual CPU19.159/12.891 s. Although each task
requested1CPU/8GiB, cluster allocation assigned64/224 CPUs respectively:
**4,928 allocated core-seconds**,32.050 actual CPU-seconds. Slurm MaxRSS was
missing/zero, so no reliable peak-memory number is claimed. No GPU/DFT/MACE.

Primary artifacts under `workspaces/accommodation_controls_20260920/`:

- `hydrogen_preparation_audit_v1.json`: all21 coordinate/topology audits.
- `hydrogen_diagnostic_v1/manifest.json`, `submission.json`, `sacct.txt`.
- `hydrogen_diagnostic_v1/sample_{3,4}/receipt.json`: energy/force/iteration
  evidence, immutable System XML, unrounded coordinate arrays and PDB outputs.
- `hydrogen_diagnostic_v1/verification.json`: exact names/bonds/heavy atoms,
  initialization overlaps, monotonicity and displacement checks.
- `hydrogen_diagnostic_v1/bond_type_check.json`: per-parent bond lengths and
  captured zero-exclusion verification.

Read-only audit reproduction is available in
[audit_hydrogen_preparations.py](audit_hydrogen_preparations.py).
The actual numerical code is [hydrogen_probe.py](hydrogen_probe.py), with
[run_hydrogen_probe.sbatch](run_hydrogen_probe.sbatch). These scripts and plans
were copied and pinned before execution. The all21 audit was first executed
as a standalone read-only Python command; its reusable script preserves that
calculation without implying a second executed audit.

## Narrow convergence continuation: passed the declared force criterion

Separately frozen [HYDROGEN_CONTINUATION_PLAN.md](HYDROGEN_CONTINUATION_PLAN.md)
authorized one 500-iteration continuation of the repaired La4 state using
exactly the same System and CPU1 implementation. **Job1203809 completed**:

- 158 accepted iterations; mobile-H force RMS **0.806246 kJ/mol/nm**, below
  the declared 1 kJ/mol/nm threshold. Maximum single-H force is 27.4572;
  RMS convergence does not assert every component is below 1.
- Zero overlaps; all heavy coordinates unchanged; H-parent distances
  1.170993–1.207169 Å; maximum H displacement 1.771882 Å.
- Process wall 20.9412 s; Slurm elapsed 23 s; actual CPU 21.878 s.
  Scheduler allocated 224 CPUs despite the 1-CPU request: 5,152 allocated
  core-seconds. No reliable MaxRSS was reported.

The four numerical preparations together used **53.928 actual CPU-seconds**
and **10,080 allocated core-seconds**. Receipts and coordinates are in
`workspaces/accommodation_controls_20260920/hydrogen_continuation_v1/result/`.
No adaptation or scoring occurred.

**Conclusion:** an isolated numerical preparation failure was reproduced
and repaired under the same potential, preserving chemistry and source heavy
geometry. This offers a provisional coverage improvement, not evidence of
better La/Ca discrimination. A chemically realistic H-placement Hamiltonian
would be a separate development protocol; the current preparer remains intact.
