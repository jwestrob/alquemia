# Response continuation: 8GY2 passes; GGR improvement is structure-sensitive

**The frozen radial hybrid retains PQQ fidelity and correctly transfers to
8GY2, but fixes only 2/6 alpha-versus-GGR structural comparisons.** DFT-only
fixes 0/6. The earlier 4/4 result included two partitions of 1GLG; the two other
GGR structures do not cross the desired ordering. This limits generality and
does not justify default promotion. The opt-in explicit-core scoring CLI works
end to end on real completed artifacts.

## Fixed method and actual execution

[Continuation plan](CONTINUATION_PLAN.md) was frozen before scoring. No refit:
use the exact all-25 PQQ DFT-only and DFT+radial models from evaluation_v2.
Reuse the same native MACE-POLAR-medium vacuum checkpoint/software/kernel and
physical donor mapping. 2FW0/2FVY use their existing 58-atom alpha-cap preparations,
matching the original GGR_extended/alpha chemistry policy. Coordinates, original
protonation and waters are unchanged. Four new TightSCF native DFT SPs match the
original direct gradient-center input policy; earlier NormalSCF SPs stay archived.

8GY2 uses its existing 71-atom canonical fixed-core v3 pair, Ca−2/La−1, dry
oxidized PQQ3− and unchanged canonical native r2SCAN-3c/CPCM input policy. Heme
and other chains are excluded under that existing core protocol. This is not a
heme-free whole-protein MACE preparation or a new assembly model.

All six MACE core force endpoints completed. All six DFT endpoints completed;
two earlier 8GY2 launches failed before SCF because the adapter removed the final
input newline. A separate two-task retry restored only that newline. Failures,
exact receipts and old inputs are retained; no successful endpoint was rerun.

## New PQQ association transfer

8GY2/O05542 is a Ca-associated membrane alcohol-dehydrogenase structural control,
not a measured La/Ca affinity direction. Its maximum sequence identity to the
previous panel is43.7068%, below the existing50% grouping cutoff; no independent
fold claim. It was prepared and evidence-qualified before this scoring. See
the [original curation](../benchmark_augmentation_20260918/REPORT.md).

| Quantity | 8GY2 result |
|---|---:|
| Native DFT E_Ca, Hartree | -2843.445775506783 |
| Native DFT E_La, Hartree | -2197.312201092884 |
| Baseline R, kcal/mol | -405454.93941420544 |
| Baseline reporting S, kcal/mol | -35.15846628858708 |
| Released baseline band | Ca-supported |
| Mean MACE radial Delta g, kcal/mol/A | 15.831441498979789 |
| Saved DFT-only logit / class | -1.5351023903027858 / Ca |
| Saved DFT+radial logit / class | -2.7623389479099507 / Ca |

Both fixed heads are correct at their unchanged zero-logit decision. The new
head does not inherit the baseline's energy bands or report an affinity
probability. A larger negative logit is not a calibrated confidence estimate.
This adds one PQQ structural-association success; the earlier grouped25/25 and
three consumed crystal successes remain unchanged, not newly independent tests.

## All direct orderings, with structural replicates grouped

Positive alpha-minus-GGR logit gives the experimental direction: La-associated
alpha above Ca-associated GGR. PQQ weights/scales are reused without fitting any
direct label. Absolute cross-target classes and probabilities are unavailable.

| Alpha structure | GGR structure/representation | DFT-only | DFT+radial |
|---|---|---:|---:|
| 1F6S | 1GLG extended | -0.412424 | **0.270990** |
| 6IP9 | 1GLG extended | -0.478508 | **0.105432** |
| 1F6S | 2FW0 alpha caps | -0.688068 | -0.548637 |
| 6IP9 | 2FW0 alpha caps | -0.754152 | -0.714195 |
| 1F6S | 2FVY alpha caps | -0.662700 | -0.346089 |
| 6IP9 | 2FVY alpha caps | -0.728784 | -0.511647 |
| 1F6S | 1GLG connected, partition check | -0.209415 | **0.514459** |
| 6IP9 | 1GLG connected, partition check | -0.275498 | **0.348901** |

Primary denominator: three GGR structures×two alpha structures, **2/6** versus
DFT-only0/6. The connected1GLG partition contributes two separate robustness
checks; do not combine it into a claim of four independent successes. There are
two biological groups and one biological comparison throughout.

Mean MACE radial loads are6.673648 for1GLG extended,−1.833071 for2FW0 and0.850778
for2FVY kcal/mol/A. The descriptor therefore varies strongly with GGR structure.
2FW0 is the sugar-free/open source;2FVY is glucose-bound/closed with reported
Glu149 radiation damage, with sugar omitted by its frozen core preparation.
These qualifications do not justify discarding their unfavorable outputs.
The present results do not isolate which structural/preparation difference
causes that sensitivity. No water-state mixing, geometry rescue or feature tuning.

## Runnable scorer and validation

[Scoring interface](SCORING_INTERFACE.md) exposes explicit prepared-core,
MACE-collection, DFT-manifest, frozen-model and baseline-release inputs. Three
real calls completed:2FW0,2FVY,8GY2. Baseline and new-head fields stay separate;
missing forces leave the head unavailable, never baseline-substituted success.
Generic/direct inputs export raw features without absolute affinity decisions.

Fifteen distinct real-fixture tests passed with zero skips:11 archive/mapping/
grouping/transfer tests in2.081s and4 scoring/receipt/failure/serialization tests
in1.221s. Actual native scientific integration was the six successful MACE and
six successful DFT endpoints, plus the preserved two failed input attempts.
No new Hessian, displacement, geometry optimization or stiffness model was used.

## Measured cost and artifacts

| Job | Outcome | Wall s | Allocated CPUs | Allocated core-s | GPU-s |
|---|---|---:|---:|---:|---:|
| 1202089 | 6 MACE gradients, H200 | 45 | 16 | 720 | 45 |
| 1202090 | 4 DFT complete;2 pre-SCF failures | 199 | 64 | 12736 | 0 |
| 1202092 | 2 DFT retries complete | 161 | 64 | 10304 | 0 |
| Total | all requested endpoints complete | | | **23760** | **45** |

1202089 moved from our pending A5000 request to an available H200 with the same
manifest/model,16CPU/200000MiB/oneGPU; no priority change or cancellation.
The retry requested32CPUs but Slurm allocated64, which is counted. MACE summed
model evaluation time7.839879s; maximum process RSS2436856KiB and CUDA allocation
516527616bytes. These small-core timings do not establish whole-protein runtime.
DFT batch peak RSS8639744KiB initially and5817084KiB on retry. Allocation figures
include failures/startup; no scheduler wait estimates are presented.

[Receipts](CONTINUATION_RECEIPTS.json) pin the actual outputs, task manifests,
execution events and accounting. Main workspace:
`workspaces/response_probe_20260919/continuation_v2/`:

- `8GY2_score_v1.json` (SHA256
  `f310cc1f6d4f581710ac721945a82a479abc3c9ae6c0094be9e383c8a1557253`).
- `GGR_2FW0_score_v1.json`, `GGR_2FVY_score_v1.json`, `direct_transfer_v1.json`.
- `mace/collection_job_1202089.json`; `dft/manifest.json` and
  `dft_8gy2_retry_v1/manifest.json`; `scoring_implementation_v1/`.

**Recommendation: retain the baseline default and the radial hybrid as an opt-in
research candidate, with its limited2/6 direct robustness clearly exposed.**
The tested tangential addition is rejected. This step produces useful partial
improvement and a runnable interface; it does not complete broad La/Ca validation.
