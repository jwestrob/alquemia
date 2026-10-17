# Terminal force discrepancy persists with the expanded partition

17 October 2026 UTC. The saved four-case comparison is complete. Expanding the
existing electronic region reduces the excessive Phe1 CA force only modestly;
the classical contribution remains dominant. **The tested boundary change does
not repair the terminal force discrepancy against whole-peptide PBE0.**

## Same geometry comparison

All four cases use the exact admitted 7CCO/cartesian_01 starting coordinates:
baseline/expanded partitions × La/Dy. The reference is the before-La/before-Dy
pair from whole-peptide collector13424, with the same fixed-charge ALPB solvent.
The electronic region expands from residues3–14 to2–15. Phe1 CA stays outside
the electronic-real region and becomes an omitted cap anchor in the expanded
partition. This comparison therefore does not test direct electronic treatment
of Phe1 CA.

Force projections use the original frozen initial hybrid Phe1 CA direction.
Units are kcal/mol/Å. The MM column is its contribution to the total projection.

| Metal | Partition | Total hybrid | MM contribution | Whole PBE0 plus ALPB |
| --- | --- | ---: | ---: | ---: |
| La | Baseline | 51.063006 | 52.001900 | 3.341027 |
| La | Expanded | 47.742478 | 49.289919 | 3.341027 |
| Dy | Baseline | 51.064870 | 52.001979 | 3.345638 |
| Dy | Expanded | 47.764896 | 49.289997 | 3.345638 |

The enlarged electronic region changes the electronic/cap contributions, but
leaves a large classical terminal drive for both metals. The earlier exact
interaction inventory found no copying or double-counting defect; these force
comparisons still do not identify a unique faulty parameter.

## Whole system comparison

| Metal | Partition | RMS force difference from reference | Maximum atom difference |
| --- | --- | ---: | ---: |
| La | Baseline | 10.141916 | 48.447093 |
| La | Expanded | 8.656817 | 48.958287 |
| Dy | Baseline | 10.162397 | 48.453867 |
| Dy | Expanded | 8.680798 | 48.963662 |

These are atom-vector force differences in kcal/mol/Å. The maximum discrepancy
is always atom246, `A/17//C`: the Ala17 terminal carboxylate carbon, classified
as MM in both partitions. Global RMS improves while this maximum increases
slightly. Thus the residual problem includes both terminal regions.

Along the original normalized full-deformation chord, hybrid projections change
La19.928221→17.150118 and Dy20.194774→17.421113, versus reference3.638430 and
3.873822. This remains a starting-force comparison. No expanded after-state is
available, so its finite deformation work is untested.

## Evidence and retained replay failures

Job13463 compares the previously admitted, pinned POINT gradients from recovered
104-cell pool13013. All four cases pass exact physical-ID/coordinate/state and
common-ALPB joins, with full256 finite-array/internal component closure. The two
baseline cases additionally replay the original forces and projections exactly.
Result SHA256:
`e25076f56edffc9881427a6d7a0a13ac63906a873ee26d84b7575fa131a579c6`.
Full tables and actual reviews are in [the archived comparison](archived_comparison_result_v1/RESULT.json).

The separate raw bitwise reconstruction checks remain failed. Job13446 exposed
an incorrectly reassociated unit conversion in the new analysis code; its source
and missed static review are retained. After that one-expression correction,
13455 still failed raw component equality. Diagnostic13460 localized the
remaining differences to MACE/GFN cap-scatter entries, at most
3.552713678800501e-15 kcal/mol/Å. Stored component sums, MM and ALPB replay exactly.
The precise runtime cause is unproven. No tolerance was introduced to pass that
gate. The separately declared archived-result comparison uses the original
admitted forces as its explicit data authority.

All four jobs used one Slurm CPU and zero molecular calls. Scheduler elapsed
times were0,1,1,1seconds respectively:3 allocated CPU-seconds total, with actual
records in [SACCT.txt](archived_comparison_result_v1/SACCT.txt). A zero-second
scheduler entry does not mean no application work occurred.

## Decision

Keep the challenged LBT hybrid response closed. The existing expanded partition
does not provide the needed mechanical correction. Further model development
must explicitly address scaffold force balance and deformation cost, with
independent validation; this comparison supports no fitted parameter change.

Whole-reference SCF, vacuum-dianion, terminal-state, frozen-f and force-accuracy
limitations remain as stated in the [four-cell report](../lbt_whole_reference_continuation/FINAL_REPORT.md).
This is conditional partition sensitivity, with no new La/Dy preference claim
or production change. Mex13281's independent central-gradient/cost gate remains
pending at this checkpoint. The broader mechanics goal remains unmet.
