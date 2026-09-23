# Approved two-pass native continuation of four fixed PQQ pools

Jacob, 2026-09-23: “alright. keep on goin buddy, launch subagents to do what
experiments you like so they ping you on completion and keep you awake. i'm goin
to bed. email me when cool stuff happens.” Parent assigned this finite experiment
under that authorization. No production or reference promotion is included.

## Frozen question and population

Does the confirmed native restart remedy settle the existing solvent corrections
and change useful discrimination on four consumed sources? Use exactly 1H4I,
4MAE, Q88JH5 and A0A3F2YLY8 Ca-conditioned sample 1. Retain each existing
five-candidate pool: origin, proposal_Ca, proposal_La, adaptive_Ca, adaptive_La.
Twenty fixed geometries × Ca/La × vacuum/ALPB give 80 cells per stage.
All source GBW/xtbw pairs and state/geometry pins are in the completed
`diagnostics/native_xtb_restart_20260922/SEED_AVAILABILITY.json`.

Stage 1 starts from each exact archived matching GBW+xtbw pair. Stage 2 starts
from its own stage 1 pair. Both use native ORCA 6.1.1 GFN2, 300 K, MaxIter 500,
native mixer, unchanged vacuum or ALPB(water) settings, charge and multiplicity,
eight MPI ranks. AutoStart is activated by the matching runtime-basename pair;
actual XTBRESTART and native-mixer output are mandatory. No GBW-only fallback,
ordinary SCF, TightSCF assumption, geometry changes or new MACE/DFT calls.

Stage 2's eight origin cells for 4MAE/Q88JH5 × both metals/media request analytic
EnGrad under the identical energy Hamiltonian. Retain atom ordering, units and
actual native gradient-driver evidence; no numerical gradients. These outputs
are shared with the separate derivative-check experiment to avoid duplicate q0
calls. All other calls are single points.

Maximum 160 new logical GFN2 calls; no automatic third continuation. Failed or
unconfirmed stage 1 cells cannot seed stage 2 and remain missing in the fixed
80-cell denominator. Both stages and all attempts remain separate and visible.

## Reporting and qualification, fixed before outputs

Stage 2 is the reported endpoint regardless of which stage is lower or gives a
preferred classification. Preserve stage 1 as numerical sensitivity. Acceptance:

- Every stage 1→2 native cell energy change is at most 0.1 kcal/mol in magnitude.
- Every paired, same-geometry composite Ca−La change is at most 0.2 kcal/mol.
- Pooled Ca−La change is at most 0.2 kcal/mol (mathematical and operational).

The composite remains E = MACE_vac + GFN2_ALPB − GFN2_vac; R = E_Ca − E_La.
Reuse exact archived MACE values. Select only among the same five geometries
using the existing common-pool mathematical minima and separate operational
0.1 kcal/mol origin tolerance. Never choose a continuation stage by its energy.
Any missing cell makes its complete pool unavailable. Failed numerical checks
leave the diagnostic raw stage 2 matrix visible but the qualified score null.

Retain old, stage 1 and stage 2 components, row choices, geometric works and
contrasts. Report both existing released and adaptive bands only as explicit
old-band transfer; no fitting four cases, new absolute reference or affinity
probability. Labels are used only in reporting. Coverage stays four sources,
20 geometries, 80 cells/stage, with all failures shown. This is development,
not independent affinity validation or proof of a global electronic minimum.

## Execution

Existing ORCA executor, scoped continuation wrapper, finite immutable manifests.
CPU-only gpu-partition route: 64 CPUs, 128 GiB, eight concurrent eight-rank tasks,
explicit available host, no GPU. Always collect after execution errors. Record
pre/post seed hashes, actual receipts, allocation cost and analytic-gradient
availability. No separate project time/compute cap beyond this finite task set.
