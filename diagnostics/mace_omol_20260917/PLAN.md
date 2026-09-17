# MACE-OMOL candidate: numerical qualification and real descriptor benchmark

Declared before downloading or evaluating this checkpoint. Jacob's active goal
and discretionary pilot approval cover these contained stages. No production
default, historical output, evidence label or scientific reference changes.
Previous goal turn: progress (canonical and fixed-field candidates completed
and rejected); no MACE job remains live at this declaration.

## Why a different backend

Direct POLAR plus frozen solvent failed the complete canonical calibration, and
the DFT/CHELPG/POLAR-local mixture still failed its partition criterion. Those
failures do not establish the performance of another pretrained MACE model.
MACE-OMOL-0 uses local geometric features with explicit graph-level charge/spin
embeddings, rather than POLAR's iterated learned global multipoles. This tests a
different energy representation; it is not a retuned POLAR parameter.

Use only the maintainers' recommended 100M model:
`MACE-omol-0-extra-large-1024.model`, release `mace_omol_0`, asset digest
`9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`.
Do not switch to the 4M model after seeing outputs. Pin the downloaded file and
upstream release metadata. Reuse installed MACE 0.3.16 / torch 2.8.0, read only,
with the normal `mace_omol` API, head `omol`, native float64 and analytic forces.
No package upgrade, training, new pair kernel, long-range adapter or added
dispersion correction. Check actual checkpoint elements and embedding classes;
installed metadata says89 elements while the current release table says83, so
the checkpoint must establish Ca/La and all required atom support directly.

The dataset and release define spin as multiplicity; singlets use spin=1,
with the exact endpoint integer charge. Inspect the actual input batch and
embedding specification; do not assume a missing charge defaults correctly.
The learned target is OMol's wB97M-V/def2-TZVPD vacuum energy, not r2SCAN-3c/CPCM.
Its ASL license is retained in the provenance; this is academic research use,
not permission for commercial distribution or production promotion.

Primary sources: [release](https://github.com/ACEsuit/mace-foundations/releases/tag/mace_omol_0),
[training configuration](https://github.com/ACEsuit/mace-foundations/blob/main/mace_omol/mace-omol.sh),
[maintainer's multiplicity example](https://github.com/ACEsuit/mace/discussions/1120),
[OMol dataset definition](https://fair-chem.github.io/omol25/).
The installed code is also pinned. No claim of absent training overlap with PDB
or any calibration case is made; this is retrospective development evidence.

## Stage A: ten real endpoint evaluations

Use exact archived canonical 1H4I and4MAE cores from the completed canonical
inventory. Original atom inventory, caps, source H, zero waters, PQQ(3−),
coordinates and paired formal charges remain unchanged.

- 1H4I La/Ca: primary, exact repeat, fixed rigid rotation, fixed translation
  `(10,-7,3)` Angstrom:8 calls. Rotation is the existing deterministic runner
  matrix, frozen in each task. These transformed inputs are numerical checks.
- 4MAE La/Ca primary:2 calls. These are consumed engineering inputs; class
  directions do not decide whether Stage B launches.

Require every call to succeed with finite energy/analytic force, correct actual
charge/spin batch, unchanged model parameters, and complete source receipts.
For repeat/rigid energy errors and paired contrasts use0.01kcal/mol; for rotated
and translated forces use0.001eV/Angstrom. No numerical force substitution.
Missing density/atomic charges are explicit: OMOL need not produce them, and a
charge-conditioned input check must not masquerade as a predicted charge sum.
No GB/PB correction can be taken from an unrelated POLAR charge distribution.

## Stage B: conditional60 new calls, no fit to transfer cases

Only if Stage A numerical/state checks pass:

1. Score all25 exact canonical calibration cores plus consumed1H4I/4MAE/1KB0
   (56logical endpoints). Reuse the four exact Stage A primary endpoints;
   run52new. Medium/large POLAR caches cannot satisfy an OMOL task.
2. Score the exact four existing center representations from the completed
   mechanics preparation: GGR extended58 and connected111 atoms, alpha1F6S52
   and alpha6IP955 atoms, La/Ca each:8new calls. Their original source H and
   explicit waters are unchanged. They are existing qualified direct/supporting
   evidence, not PQQ functional-class labels. Both GGR representations are
   retained; neither is selected to make the result favorable.

Total across both stages:70 new OMOL energy/analytic-force evaluations,
zero DFT/charge extraction/solvent calls. No structural optimization, Hessian,
trajectory, or training. Existing exact mapping, preparation, state and receipt
validators must pass before submission. Other agents' jobs are untouched.

## Descriptor, interpretation and fixed criteria

`R_OMOL = (E_Ca - E_La) * EV_TO_KCAL`, converted once. Larger is La-like.
It is a vacuum electronic descriptor with no solution-affinity or free-energy
interpretation. Solvent, aquo reference, S, relaxation and entropy stay null.
It does not inherit the failed POLAR calibration or baseline thresholds.

Canonical calibration uses only its25 designated rows: U=maxCaR; L=minLaR.
Require all25 valid and L-U>0.02kcal/mol before releasing this candidate's own
Ca R<=U, La R>=L and open-gap inconclusive bands. No fitted weights, flexible
classifier, sign reversal or optimized threshold. Apply the frozen rule to
all three transfer cases, requiring each to reach its expected supported region
for the canonical operational gate. Failed calibration leaves transfer decisions
unavailable. Training separation is not independent predictive accuracy.

For the separate non-PQQ stratum, retain all four scores. Primary relative tests
are alpha1F6S minus GGRextended and alpha6IP9 minus GGRextended, each >0.02kcal/mol
for the previously frozen condition-qualified La versus Ca direction. Repeat
the same two contrasts with GGRconnected as a representation robustness test;
any disagreement prevents a robust-affinity claim. Both alpha structures belong
to one biological group; GGR's two cores are one structure/group. No threshold
is trained on these four cases and no PQQ band is applied to them. Passing these
consumed comparisons would motivate further testing, not prove broad affinity.

Preserve baseline results side by side. Group calibration homologues and both
crystal accession overlaps (1H4I/P16027,4MAE/I0JWN7). No unscorable case drops
from the denominator. Report every failure and the original H/capping limits.

## Resources and completion

One A5000,16CPUs and64474MiB per job, existing finite-manifest runner and wrapper.
The checkpoint is422242640bytes; actual inference memory/runtime is unmeasured.
Expect small-core rather than whole-protein cost; measure both startup and
evaluation, host/GPU peaks, complete allocations and failed attempts. There is
no project CPU/time budget. The scheduler's existing QOS limit remains required.
An OOM permits a technical memory recovery with the same model and inputs;
scientific parameters and planned cases cannot change to hide failure.

Deliver runnable preparation/dry-run/execute/collect/report operations, typed
caches, actual receipts and a compact paired report. End-to-end numerical
feasibility or this narrow trial alone does not complete the broader goal.
Retain useful components and continue from evidence; do not make a model win.
