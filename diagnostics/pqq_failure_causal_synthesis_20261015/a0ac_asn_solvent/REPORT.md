# Fixed Asn flip changes two solvent-inclusive pools

**All20 new endpoints completed and passed the declared scientific gates.**
Adding the fixed Asn290 flip to the full original three-geometry pool changes
the Ca-conditioned sample4 contrast by **+10.974999511265196kcal/mol** and the
La-conditioned sample4 contrast by **+0.6730578923597932kcal/mol**. The other
eight sources are unchanged. All ten sources and all60 original cells remain
in the comparison. No class bands, classification claims or new fit are applied.

The contrast is the selected local composite **Ca energy minus La energy**.
In both changed sources, Ca retains `adaptive_Ca` and La switches from
`adaptive_La` to `amide_flip`. Thus the positive change reflects lower selected
La energy. Mathematical and operational selections agree; the separately
reported origin policy is unchanged for every source.

## Complete source denominator

| Source conditioning | Sample | Common3→4 contrast change, kcal/mol | New Ca / La selections |
|---|---:|---:|---|
| Ca | 0 | 0 | adaptive_Ca / adaptive_La |
| Ca | 1 | 0 | adaptive_Ca / adaptive_La |
| Ca | 2 | 0 | adaptive_Ca / adaptive_La |
| Ca | 3 | 0 | adaptive_Ca / adaptive_La |
| Ca | 4 | +10.974999511265196 | adaptive_Ca / amide_flip |
| La | 0 | 0 | adaptive_Ca / adaptive_La |
| La | 1 | 0 | adaptive_Ca / adaptive_La |
| La | 2 | 0 | adaptive_Ca / adaptive_La |
| La | 3 | 0 | adaptive_La / adaptive_La |
| La | 4 | +0.6730578923597932 | adaptive_Ca / amide_flip |

The actual five-source Ca-conditioned range decreases
23.325984036317095→12.3509845250519kcal/mol; the La-conditioned range decreases
20.459434523945674→19.78637663158588. Coverage is5/5 in both groups. These are
correlated, consumed structural samples and raw local-energy contrasts.

## What changes the target result

For Ca-conditioned4, the component decomposition of the **change from the old
selected pool** is:

| Component | Ca−La change, kcal/mol |
|---|---:|
| MACE vacuum | +1.460981253638387 |
| Frozen-charge ALPB, total | +9.514018257640517 |
| ↳ Gelec | +9.280010684416428 |
| ↳ Ghb | +0.18693228335101925 |
| ↳ Gsasa | +0.047075289873041604 |
| ↳ Gshift | 0 |

The saved component-sum change is10.974999511278902kcal/mol; direct subtraction
of the stored pool contrasts gives10.974999511265196. Both are retained exactly
in [RESULT.json](RESULT.json). The final algebra check passes its unchanged
1e−8kcal/mol tolerance.

For La-conditioned4, the flip is **worse in vacuum than its old selected La
candidate**, giving a paired MACE contribution−8.883471120330942kcal/mol.
The solvent contribution+9.556529012658293 outweighs that penalty, producing the
smaller positive composite change. Its Gelec contribution is+9.301310644008552,
Ghb+0.21426276093739438, Gsasa+0.040955607712362146 and Gshift0kcal/mol.

These full-pool changes must not be confused with the much larger flip-from-origin
works. For Ca-conditioned4, the composite flip works are−1.4336755955376512 for
Ca and−35.09741659116525kcal/mol for La; much of the La stabilization was already
available in the old adaptive pool. For La-conditioned4 those works are
−1.4393188256742668 and−35.27323909668904. Exact per-metal components for all ten
sources are copied from the completed Slurm collection into RESULT.json.

## Geometry supports a specific, conditional explanation

The independent saved-geometry job4626 confirms that target Asn OD1 approaches
the metal,4.595274311725036→2.562017935679298Å, while ND2 retreats,
2.4460701134677234→4.718439745320599Å. La-conditioned4 shows the same qualitative
change. This supports sensitivity to the omitted terminal-amide orientation;
it does not isolate one atom's energetic contribution or establish a bad fold.

The target also has an included, nonbonded HD22–Gly292 backbone-H contact of
**1.3821908645320085Å**, versus4.943529710641982Å before the flip. The analogous
La-conditioned4 contact is1.3744771216817757Å. These contacts are represented
in the204-atom energy model and remain unrelaxed. The full source supplied
geometry/contact checks only. Dry oxidized PQQ3−, neutral Asn,14 caps and
Ca−4/La−3 singlets are unchanged. This single discrete180° candidate neither
establishes an equilibrium conformer nor tests every relaxed Asn basin.

See [the actual geometry report](../a0ac_asn_geometry_response/REPORT.md).
The old template starts remain rejected under their original local-domain
rules; the present declared discrete protocol does not retrospectively pass them.

## Scientific admission and preserved failures

The execution has20/20 available new endpoints and zero missing cells. All60
baseline cells were revalidated before new SCC work. Actual component-v2
acceptance, native water-diagonal replay, fixed-cavity checks and informative
vacuum-charge precision passed for every new and reused endpoint. Preflight
reproduced390 baseline pool comparisons with zero residual. Final collection
records420 algebra checks with maximum residual3.2443381314806174e−11kcal/mol.

**Absolute-v1 parity did not pass:** both saved vacuum and water absolute-point
flags remain false for every60 old and20 new endpoint. The accepted component
recipe does not relabel these historical failures. Actual raw receipts, partial
result records and all old source/solver failures are preserved. Scheduler
success is additional operational evidence, not the scientific acceptance rule.

The executed backend is published-parameter standalone GFN2 plus the qualified
internal native ALPB factory at vacuum charges, retaining Gelec/Ghb/Gsasa/Gshift.
The corrected full ORCA La JSON is lineage evidence; ORCA was not executed in
this stage. No default native ORCA parameter substitution occurred.

## Measured cost and evidence

Slurm4623 completed with exit0:0 on `node-64-768g-3`: **475 scheduler seconds,
7 allocated CPUs**, TotalCPU48:23.695 and batch MaxRSS1,449,900KiB. The stage
receipt measures474.2548175220145s and3319.7837226541014 allocated CPU-s through
collection. Actual execution records80 SCC starts/80 terminal backend receipts
and40 native component calls. Per-process resource receipts are pinned.

These costs cover this solvent stage only. They exclude the earlier20-call
MACE stage, geometry work and original60 baseline calculations. There were no
new MACE, DFT, GPU or optimization calls in4623, and no new scientific calculation
was performed to write this report.

Authoritative [FINAL_COLLECTION](../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_solvent_v1/execution/FINAL_COLLECTION.json)
SHA256 `212f8d0ffc82b070a3a7c78ef47e2441acaa76ce061825ffd3662c98f058202f`.
[RESULT.json](RESULT.json) pins configuration, execution, preflight, original
inventory, all source outcomes and geometry evidence. Frozen settings, thresholds,
production results and all historical records remain unchanged.
