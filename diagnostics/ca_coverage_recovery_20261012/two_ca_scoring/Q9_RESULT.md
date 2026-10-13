# Q9/ExaA: five source exclusions recovered

All five new two-Ca Q9Z4J7/ExaA sources classify correctly as **Ca**, both at the
origin and after accommodation. This demonstrates usable source coverage for the
five original off-site-Ca exclusions. It does **not** demonstrate incremental
accommodation accuracy on Q9: origin already gives 5/5 correct.

The original 125 ledger remains unchanged. The separately annotated overlay now
contains **112 correct / 1 wrong / 2 inconclusive / 10 unavailable** under the
operational and mathematical policies; its origin counts are **110 / 1 / 4 / 10**.
It includes the earlier two MxaF and one MMOL recoveries. The remaining ten C5/Q88
results are not included in this Q9 checkpoint.

## Actual decisions and accommodation

All five samples belong to one consumed reference protein. `R = E_Ca - E_La` is
the raw composite contrast in kcal/mol, not an affinity. Each policy uses its
unchanged canonical bands; no new fit or triple policy was introduced.
Operational and mathematical selections agree. `Ca proposal` / `La proposal`
name the shared-pool geometry; work is the selected energy minus origin energy.

| Sample | Origin R | Accommodated R | Ca selected / work | La selected / work | All three calls |
|---|---:|---:|---|---|---|
| 0 | -405473.319861 | -405467.907129 | Ca proposal / -7.695491 | La proposal / -13.108223 | Ca, correct |
| 1 | -405478.932621 | -405463.458903 | La proposal / -4.978557 | La proposal / -20.452275 | Ca, correct |
| 2 | -405478.507987 | -405468.045479 | La proposal / -1.547800 | La proposal / -12.010308 | Ca, correct |
| 3 | -405479.872167 | -405462.658359 | Ca proposal / -8.842262 | La proposal / -26.056070 | Ca, correct |
| 4 | -405483.032252 | -405471.223117 | Ca proposal / -3.215799 | La proposal / -15.024935 | Ca, correct |

The machine record retains unrounded values and component/input pins. Ca selects
La's proposed geometry in samples 1 and 2; both metals still choose independently
from the same three geometries. No favorable-source filtering occurred.

## Preparation and numerical qualification

- All five spectator-aware preparations passed, retaining 192-atom local contexts
  and the exact frozen state anchor. The existing single-Ca preparation regression
  passed with unchanged coordinates. Both ions remain in each full source.
- Both seven-CPU lanes passed the two real qualification fixtures. Energy errors
  were exactly 0; maximum force errors were 3.9745984281580604e-14 (Ca) and
  9.547918011776346e-14 eV/Å (La), below the unchanged 0.01 kcal/mol /
  0.001 eV/Å gates. Both scientific preflights passed.
- All required six metal/candidate components per source were accepted. Every
  origin/Ca-proposal/La-proposal paired-coordinate error is exactly 0 Å.
- Corrected native xTB 6.7.1 GFN2 and independent tblite 0.7.0 passed the frozen
  matched-solvent-component gate (1e-6 Eh), population agreement (1e-5 e), and
  charge closure (1e-5 e). **Historical absolute-v1 parity remains failed**;
  component-v2 acceptance does not erase or pass that separate absolute test.

The MACE checkpoint/head/float64, four-mode force-selected SLSQP, shared pool,
0.1 kcal/mol origin-retention tolerance and corrected solvent settings are unchanged.
The canonical record remains SHA256
`929dedab7ba46c7088de7a4a695a8f19fc5af60066ecd3599645f687bd4e6896`.

## Measured work

| Stage / job | Allocation | Actual calls | Measured stage elapsed (s) | Allocated CPU-seconds |
|---|---|---|---:|---:|
| Native lane 0 / 2996 | 7 CPU | 97 MACE forwards | 1482.1968724459875 | 10375.378107121913 |
| Native lane 1 / 2997 | 7 CPU | 70 MACE forwards | 1074.963524626044 | 7524.744672382309 |
| Solvent adapter / 2999 | 12 CPU | 60 xTB + 60 tblite | 392.838965177536 | 4714.067582130432 |

Native counts include each lane’s two qualification forwards. Native model-call
CPU times were 8264.829313013 s and 5982.886364762998 s; no started calls lack terminal
receipts. Both lanes ran concurrently, so summed lane wall times are not latency.
Solvent backend CPU time was 4162.7366 s; summed individual call wall time was
4168.07093667984 s, with 112848 KiB maximum backend-process RSS. The inner solvent
executor through final collection measured 388.7263388633728 s /
4664.716066360474 allocated CPU-s; the table includes adapter preparation and overlay.

Source preparation 2969 separately measured 26.537397091975436 s wall and
81.25002799999999 s process CPU with 6 allocated CPUs and 5 source workers.
Its self/child peak RSS was 223436 KiB. It made no energy or optimization calls.
All four scoring/collector jobs completed successfully; the terminal watcher
recorded 1482 s, 1075 s, 0 s and 394 s for 2996–2999 respectively. This report does not
aggregate folding, scheduler waits, earlier failed/reused work or total campaign
cost, and does not extrapolate a GPU speedup or new-site latency.

## Interpretation and next handoff

The alternative full sources provide admissible, correctly classified local
pockets. The distant second Ca is omitted explicitly from the local Hamiltonian;
no global spectator response was computed. Folding settings were not established
identical to the historical sources, so the improvement cannot be assigned solely
to occupancy. Five correlated samples do not establish new biological transfer.

Keep the pending C5/Q88 batch and cumulative job 3049 on their frozen manifests.
Use their actual terminal collections for the full fifteen-source conclusion.
No production promotion, new optimization or retuning follows this Q9 result.

## Artifact record

- Full source: `workspaces/ca_coverage_recovery_20261012/two_ca_solvent_q9_v1/FINAL_COLLECTION.json`;
  SHA256 `ebf279eec368d4415ee7628a0e74ba7a59451bc7d54ed8167e1788b085b62dec`.
- Original 125 plus explicitly linked alternatives:
  `two_ca_solvent_q9_v1/FINAL_COLLECTION_ALTERNATIVE_SOURCE_OVERLAY.json`;
  SHA256 `fcad56728990524c4cbd0dd0d33687b82a99f7e24320f17903aa9d79ce199dc6`.
- [Machine-readable copied-field report](Q9_RESULT.json) pins actual collections,
  preparation, qualifiers, completion event, model and method settings.
- This documentation was assembled from saved Slurm outputs; no new scientific
  computation, reclassification or molecular call ran while writing it.
