# Carve-free grouped-vacuum MACE transfer: 17/22 directions, qualification incomplete

**All55 evaluations completed.** The new whole-protein descriptor gets all
seven original directions right, including both alpha structures against all
three GGR structures. Supporting-family performance is10/15, unchanged from
the earlier masked-OMOL descriptor. RTX and parvalbumin CD remain failures.
The production baseline is unchanged. This is useful research progress toward
an empirical compatibility discriminator, without a validated replacement claim.

Protocol: `intact_POLAR_medium_typed31_group_vacuum_compatibility_v1`.
The source-based group policy and all comparisons were frozen in PLAN.md.

## Predictive results, separated by evidence

| Stratum | Passing / total | Interpretation |
|---|---:|---|
| XoxF4MAE over MxaF1H4I |1/1| Consumed functional-class development pair; margin57.141809 model kcal. |
| Alpha1F6S/6IP9 over GGR1GLG/2FW0/2FVY |6/6| Qualified alpha direction versus direct GGR affinity direction; not six independent biological observations. |
| Khoury all-site means over three GGR structures |6/9| A0A7 and HEW5 pass; RTX fails. LaITC versus CaCDfolding proxy, not a matched Ca/LaKd ratio. |
| Both parvalbumin sites over three GGR structures |4/6| Cross-study supporting evidence; CD fails against2FW0/2FVY. |

Nine biological groups are represented, eight with directional/supporting
labels. All structures/labels were previously inspected under other models;
none is claimed as blind.22 is a computational comparison denominator, not
22independent biological validations or a pooled affinity-classification accuracy.
Khoury means include every declared site:6A0A7,8HEW5,8RTX. Aequorin retains
EF1/EF3/EF4 as an ordered vector, with no sitewise experimental labels assigned.

| Remaining incorrect comparison | Margin, model kcal |
|---|---:|
| RTX_all_site_mean minus GGR_1GLG | -70.962794 |
| RTX_all_site_mean minus GGR_2FW0 | -129.662933 |
| RTX_all_site_mean minus GGR_2FVY | -120.507921 |
| PARV_4CPV_CD minus GGR_2FW0 | -17.732278 |
| PARV_4CPV_CD minus GGR_2FVY | -8.577266 |

R=E(Ca)-E(La), with exactly one eV conversion; larger R is more La-like.
Only differences between compatible sites cancel the shared metal-energy offset.
No aquo reference, inherited baseline bands, fitted zero or absolute class exists.
Aequorin raw R remains[-405478.821771,-405390.160027,-405364.905555] in EF1/EF3/EF4
order. Conditional single-site substitutions with other ions fixed Ca are not
experimental cooperative binding/folding free energies.

## What transfers, and what does not

The direct whole-protein model preserves the earlier vacuum component's six
alpha/GGR successes across all three GGR structures. Earlier bounded hybrid
response gains had been confined to1GLG. This result requires no new DFT or
mechanical correction. It remains consumed development evidence.

The uniform donor-contact grouping removes the primary dependence on a carved
QM region. PQQ separation grows from51.626321 in the earlier carve-defined
vacuum component to57.141809. The GGR and alpha scores reproduce the earlier
vacuum component; its declared solvent-corrected result remains1/7, unchanged.
No old failed total was reclassified as a successful vacuum model.

On the added supporting families, old masked-OMOL and the new descriptor both
pass10/15 with the same failures. Thus broader transfer adds evidence of scope,
but demonstrates no accuracy gain on those added families. No site, input state,
threshold or grouping radius was changed to repair an inconvenient result.

## Numerical consistency versus representation sensitivity

All61 numerical checks pass:55 endpoint charge closures,three rigid/permutation
energy checks and three force checks. Maximum rigid energy error7.818e-6 model
kcal and force error6.682e-8 eV/A. Reversed-order La energy reproduces exactly;
its force error is6.439e-15 eV/A. All-Ca source atoms, geometry, groups and charges
match across sites of each multisite protein, permitting one evaluated reference
per protein without replacing missing results.

All three GGR group-expansion shifts fail the unchanged2-model-kcal tolerance:

| GGR | Total shift | Electrostatic readout | Electron readout |
|---|---:|---:|---:|
| GGR_1GLG | -2.513855 | -1.866445 | -0.647410 |
| GGR_2FW0 | -3.116626 | -2.555801 | -0.560825 |
| GGR_2FVY | -3.275712 | -2.661396 | -0.614317 |

The learned short interaction readout is unchanged to numerical precision.
These components locate the representation effect algebraically, not a unique
biological cause. Every one of the22directional pass/fail outcomes is unchanged
under the declared alternate GGR grouping. Keep both facts: rankings survive
this test, while the frozen numerical score-sensitivity criterion fails.
Qualified comparison count remains0; the tolerance was not relaxed.

The model imposes formal group charges on learned density coefficients. This
is an unvalidated constraint, not exact charge partitioning or self-consistent
protein/solvent thermodynamics. Vacuum is an empirical descriptor choice; it is
not evidence that solvation can be omitted from an affinity calculation. No
relaxation, entropy correction or combined solvent gradient is supplied.

## Implementation and actual execution

New source-graph preparation covers all34 sites with complete residues,actual
amide C-N neighbors,disulfide joins,whole waters/PQQ and all metal sites at once.
Shared ligands merge groups. Heavy/H coordinates,source assembly,protonation and
water inventory remain unchanged. No caps, artificial bonds or force-field atom
charges enter the model. Same medium checkpoint,full learned readouts and
existing analytic/memory wrappers; baseline runner behavior is unchanged.

Job1201517 completed55calls in**1236GPU-allocation seconds (20min36s)**,
19776allocated core-seconds,1435.470reported CPU-seconds. Model inference sums
853.299751s. Peak GPU allocation10207527424bytes; peak reserved12043943936bytes.
No molecular failures,newDFT,solvent calls,training or trajectories.

Primary inference:alpha pair17.32–17.96s,GGR pair49.08–49.16s,PQQ pair112.99–
116.93s. Complete multisite primary sets cost36.04s(A0A7),60.84s(HEW5),89.70s(RTX),
21.96s(parvalbumin),54.60s(aequorin); these include one all-Ca reference plus
every single-La endpoint. Allocation totals include launch/validation overhead.

Two earlier preflights failed before inference:an unrelated gemmi import and
4e-16A NumPy norm rounding. Pure source-definition loading and1e-12A metadata
roundoff tolerance fixed them without changing contact/group/charge inventories
or any of55endpointXYZ hashes/settings. Both source versions/failures remain.
Two input preparations together47.951765wall/47.073242CPU seconds. Additional
local preflights,tests,reports and prior development costs are not claimed zero.
No matched production-overhead guarantee is established by this pilot.

Eight real-fixture tests pass96.968s,zero skips,including source/amide/multimetal
invariants,corrupted input rejection,exact donor-definition equivalence,actual
partial-output handling,full energy/mean/vector algebra and report replay.
Earlier six old-model regression tests passed45.119s. Scientific integration
used actual outputs from the55evaluations; parser tests fabricated no energies.

## Recommendation and next test

**Pursue the descriptor further; retain the production baseline.** It now
provides a reproducible,affordable whole-protein score with useful raw rankings,
but no absolute calibration and a retained grouping-sensitivity failure.
Before adding more physical corrections, test this unchanged model on the
existing25canonical PQQ preparations as a designated calibration attempt.
Those full-protein preparations already exist. Do not inherit baseline bands,
call the consumed crystal pair blind,or force a threshold if the classes overlap.
That follow-up has not yet been declared or launched by this report.

[Operations](COMMANDS.md) provide the exact report replay. Full raw endpoints,
source maps,execution receipts,partial reports,group components,costs and
standalone figures remain under workspaces/mace_group_transfer_20260918.
