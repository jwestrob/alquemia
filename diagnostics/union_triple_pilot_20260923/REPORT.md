# Three-fold context pilot retains calibration and the difficult-source repair

The context union constructed from three declared La-conditioned folds supports
25/25 canonical calls and all three consumed crystal controls. The separate
A0A3 Ca-conditioned stress source remains La-supported, with a 1.1024 model
kcal/mol margin, despite reducing its context from 193 to 169 atoms. This earns
the planned fold-transfer test; it does not establish its outcome.

## Actual result

| Population | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Original designated canonical calibration | 25 | 0 | 0 | 0 |
| Consumed crystal transfers, excluded from fit | 3 | 0 | 0 | 0 |
| Separate A0A3 Ca1 stress source, excluded from fit | 1 | 0 | 0 | 0 |

All three-member context selections were fixed before energies. Each protein's
canonical geometry uses its lexicographically first complete noncanonical
La-conditioned triple; skipped incomplete Q88 triples remain recorded. No Ca
source contributed to selecting the stress probe's context. Its La charge stays
−2, with 13 caps versus 15 in the ten-fold context. Composition and cavity change
together; this is not a uniquely identified hydrogen-bond effect.

The new 25-source reference independently recomputes Ca_max
−405463.71230032056 and La_min −405456.3869114709, a 7.325388849654701 model
kcal/mol gap. These happen to equal the ten-fold numerical candidate's edges:
unchanged canonical sources still determine both extrema. They were not copied
as the new calibration. Mathematical and operational pool policies give the same
calls here. Full unrounded components, old-band transfer and source/state pins are
in [the actual reference](../../workspaces/union_triple_pilot_20260923/REFERENCE_v1.json).

## What ran

Nineteen complete pools were reused only after exact graph/state/cap-map,
optimizer, checkpoint, component and receipt checks. Ten changed source contexts
received 20 native-MACE searches, using the existing four-mode selector,
ftol=1e−8 Hartree-equivalent, 200-iteration limit and unchanged ±0.8-radian/0.8-Å
physical domain. Both metals score the same origin/Ca-proposal/La-proposal pool.
All 20 searches and 20 cross-MACE evaluations completed. Eleven searches retain
boundary flags; these are finite proposals, not unconstrained minima.

Six origin pairs needed fresh energies/forces: 12 MACE and 24 native GFN2 calls.
Searches used 389 fresh MACE evaluations; candidate scoring added 20 cross-MACE
and 80 GFN2 calls. Total fresh work: **421 MACE and 104 GFN2**, no failures,
retries, DFT, new folds, protonation or water changes. Native GFN2 uses the frozen
ORCA native mixer, 300 K, MaxIter500 and no-autostart recipe. New scalar tasks use
one rank; exact reused eight-rank receipts retain their actual provenance.

The score remains E(Ca)−E(La), with E = native OMOL vacuum energy plus
GFN2 ALPB(water) minus GFN2 vacuum at the same coordinates/state. Unit conversion
is applied once per component. No force-field, entropy or probability term is
added. The large raw offset is a model/reference convention, not absolute affinity.

## Measured cost and verification

Four completed allocations used **14,368 allocated CPU-seconds and 296 allocated
GPU-seconds**. Jobs1211064/65/73/74 lasted141/86/155/67 seconds respectively,
all on node-224-2t-8gpu-1 with32CPUs; GPU stages requested one H200. Nested runner
receipts are not added again. Executor-only times were9.92s origin MACE,
111.60s searches,14.94s cross-MACE,52.73s origin GFN2 and44.54s candidate GFN2.
These exclude folding and most preparation; they are not a fresh scanner timing.
The separate preparation branch recorded433.57 local CPU-seconds including its
preserved technical attempts. Additional pilot-local audits and tests are not
included in allocation totals. See [costs](COSTS_v1.json) and [scheduler receipt](SACCT_FINAL.txt).

Ten real-artifact tests pass in14.324s, zero skips. They check exact reuse,
paired states, method settings, actual cell counts, component/pool algebra,
reference membership, boundary status and rejection of corrupted real manifests.
Actual molecular integration ran as listed above; parser tests are not substitutes
for it. See [tests](TESTS_final_v1.txt) and [commands](COMMANDS.md).

## Scope and next step

This is consumed-reference development evidence. Full100 three-of-four La-fold
transfer remains unexecuted:94 complete triples and six pre-existing unavailable
triples, represented by125 distinct source/context pairs. Prepare that exact finite
manifest and reuse inventory next, then review the counts before submission.
No arbitrary unknown-source API, production promotion, PLM rescore or broad
biological validation follows from these29 calls. Production and prior protocols
remain unchanged.
