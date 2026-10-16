# La/Dy resumed: target coupled donor/metal response

**Saved forces identify Asp85/Asp87 as useful structural-response leads, but do
not establish a better La/Dy classifier.** Moving infinitesimally toward the other
Hans protein conformation changes each conditional substitution contrast toward
reducing their disagreement. The entire EF3 source-directed displacement is
nevertheless uphill for both metals. A raw interpolation between crystals is not
an appropriate relaxation proposal.

## Actual result

All four accepted repaired Hans origin gradients assembled with their matching
classical terms. Both sources share 1,665 protein atoms and the same protein
charge map. Heavy-atom alignment RMSD is 1.07794 Å. Waters and neighboring-ion
states remain source-specific, and the exchanged metal is held fixed for these
protein-displacement projections.

Each row below points from its own source toward the other source after alignment.
Numbers are normalized directional **gradients**, in kcal/mol/Å; negative is
locally downhill. Differential means Dy minus La, not an affinity.

| Source and direction | Whole protein: La | Dy | Dy−La | EF3 alone: La | Dy |
|---|---:|---:|---:|---:|---:|
| Hans8DQ2 toward 8FNR | −4.821861 | −5.092169 | −0.270309 | +6.234828 | +5.466344 |
| Hans8FNR toward 8DQ2 | −6.711228 | −6.116692 | +0.594536 | +9.984232 | +11.675955 |

The opposite differential signs both locally narrow the existing source gap.
EF3 carries almost all of this differential projection: its contributions using
the whole-protein normalization are −0.270739/+0.595999, versus
+0.000430/−0.001463 for the remaining protein. Locally normalized residue slopes
are strongest at Asp85 (−3.612990/+6.685773) and Asp87 (−2.736180/+3.657369).
These are each residue's Cartesian source-difference directions, not measured
single-torsion responses. The completed Glu91 finite test remains closed.

Whole-protein downhill response largely reflects shared strain. Common heavy/H
contributions strongly cancel: +31.599/−36.556 for DQ2 and −21.107/+14.693 for FNR.
Differential force norms are predominantly heavy-atom, not hydrogen: protein
heavy/H norms are 39.853/0.526 and 74.921/0.992. The target metal itself has a
substantial differential force norm (13.798/12.894), but its motion is excluded
from this particular source-protein direction.

## What to do next

Prepare a small covalent-preserving **coupled donor/metal** candidate informed by
the actual Asp85/Asp87 source differences. Include the scaffold's mechanical cost;
do not move disconnected donor fragments or minimize the metal difference itself.
Inspect the mapped motion before declaring a finite paired endpoint matrix. The
existing whole-protein source vector is diagnostic only. No quantum batch is
queued from these slopes, and no relaxation energy is extrapolated across the
full crystal-to-crystal displacement.

The two source states still differ in hydration and neighboring ions. Near-zero
exterior differential loads do not show that scaffold mechanics are unnecessary:
most of this finite model's metal-dependent electronic response is local by
construction. The embedding uses formal spectator charges; the matched classical
ledger additionally includes ion-specific pure 12–6 terms. Frozen-f electronic
state limitations, missing bulk solvent and incomplete full hybrid derivative
qualification remain. These loads can guide a test; they cannot supply populations,
protein affinity, a LanM classifier, or proof of a unique cause of source reversal.

## Execution and preservation

Slurm **11788 completed** with exit 0:0, one standard CPU, 3 scheduler seconds
(**3 allocated CPU-seconds**), 3.363224 seconds measured application wall time,
2.992 seconds recorded CPU time and 108944K batch MaxRSS. No new quantum,
MACE, solvent or classical molecular evaluations; no GPU allocation. The previously
completed native and classical gradients were reused. Four physical assemblies
passed existing coordinate, energy, charge/state, source-ID, cap-Jacobian and
classical-exclusion checks. New per-residue additive sums and rotation checks
passed. These are artifact/algebra checks, not new molecular derivative tests.

Independent static review caught and corrected copied-XYZ path equality before
execution; comparisons now require verified matching content hashes. Root reviewed
the actual completed receipt. A completion wake was successfully queued. Delayed
notifications require acknowledgement only, never a repeat of this analysis.

Read [RESULT.json](RESULT.json), [per-residue loads](RESIDUES.tsv),
[execution receipt](EXECUTION.json), [accounting](ACCOUNTING.json) and the
[declared plan](PLAN.md). Full arrays and exact source pins remain in
`workspaces/lanm_source_response_20261016/profile_v1/execution/PROFILE.json`, SHA256
`aeeb796369328a31709895e8128a4e5c7e9851c5df92af6c12198e7cc940a059`.
The previous six-origin scientific result and every earlier failure are unchanged.

Optional saved-data reproduction, from the repository root, in a fresh directory:

```bash
sbatch --parsable \
  --output=workspaces/lanm_source_response_20261016/profile_v1/replay-%j.out \
  --error=workspaces/lanm_source_response_20261016/profile_v1/replay-%j.err \
  workspaces/lanm_source_response_20261016/profile_v1/RUN.sbatch \
  workspaces/lanm_source_response_20261016/profile_v1/PROFILE.py \
  workspaces/lanm_source_response_20261016/profile_v1/PLAN.md \
  workspaces/lanm_source_response_20261016/profile_v1/SOURCE_PROVENANCE.json \
  workspaces/lanm_source_response_20261016/profile_v1/replay_execution
```

No reproduction is needed to obtain this result. The next scientific experiment
requires a separately defined physical motion, not another unchanged force replay.
