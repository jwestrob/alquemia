# GGR mechanism investigation

**Complete: 38/38 endpoints, no failed attempts or retries; all 12 directional
energy/gradient checks pass.** This executes the
[approved plan](AGREEMENT.md). All cases are method-development observations;
the production baseline, PQQ inputs, reference and decision bands are unchanged.

## What the completed comparisons establish

GGR is sensitive to both local representation and source structure. Keeping
more peptide context moves its contrast by several kcal/mol; changing the
deposited source structure moves it by roughly ten. Extending a nearby amide
boundary alone does not resolve the alpha-lactalbumin–GGR ordering conflict.
These results do not identify a unique missing physical term or validate a new
affinity scorer.

### A: local representation

`R = E_Ca − E_La`; `S = (R − R_aquo) × 627.509474`, with R and R_aquo in Hartree
and S in kcal/mol. The unchanged
`R_aquo = −646.0775458314704 Hartree` is a reporting gauge for these new models.
Larger S is more La-like. The new models have **no calibrated zero or inherited
PQQ bands**. Differences in S equal converted differences in R; the common
reference cancels.

| GGR 1GLG representation | Atoms | S (kcal/mol) | Change from previous row |
|---|---:|---:|---:|
| Historical defective v2, archived | 50 | −1.183146 | — |
| Repaired formamide v3, archived | 52 | +3.057473 | +4.240618 |
| Extended amide, new | 58 | +4.269408 | +1.211935 |
| Connected native peptide, new | 111 | −3.074091 | −7.343499 |

The final row is not a rescued classification: its threshold is uncalibrated,
and the comparison changes electronic context and the continuum cavity together.
It is not a matched-system QM/MM partition test or proof of convergence to the
whole protein. No geometry was relaxed.

| Matched formamide → extended-amide change | ΔR (kcal/mol) |
|---|---:|
| GGR 1GLG | +1.211935 |
| Aequorin 1SL8 EF3 | +0.810092 |
| Alpha-lactalbumin 1F6S | +0.852523 |
| Alpha-lactalbumin 6IP9 | +1.445474 |

Under the common extended-amide policy, alpha minus GGR is **−14.918658** and
**−17.309166 kcal/mol** for 1F6S and 6IP9. The previously identified ordering
conflict therefore persists. These alpha structures are one protein observation;
aequorin EF3 has no independent site-specific affinity label here.

Evidence strata remain separate: GGR has a direct same-assay Ca-over-La affinity
direction (reported Kd Ca 25±11 µM, La 729±4 µM, 25 °C, pH 6); the alpha direction
is condition-qualified and does not supply a common-assay cross-protein affinity
scale. The models retain their frozen pH-7 preparation policy. Neither result is
an experimental physiological-occupancy or functional-use label. See the
[GGR evidence audit](../benchmark_set_20260915/GGR_FAILURE_AUDIT.md) and
[alpha audit](../benchmark_set_20260915/ALACTA_FAILURE_AUDIT.md).

Printed component bookkeeping for extended → connected GGR gives −7.567057
kcal/mol in SCF, +0.242704 in dispersion and −0.019146 in gCP, closing to
−7.343499. The printed CPCM contrast changes by +24.658529 while SCF minus printed
CPCM changes by −32.225587. Those large opposing terms forbid a simple
"the solvent caused it" attribution: SCF minus CPCM is evaluated on a polarized
density and is not a vacuum calculation. All components remain in the score.

### B: source-structure robustness

| Source | State | Formamide S | Extended-amide S | Extension ΔR |
|---|---|---:|---:|---:|
| 1GLG | Galactose-bound | +3.057473 | +4.269408 | +1.211935 |
| 2FW0 | Sugar-free/open | +13.326287 | +14.240438 | +0.914151 |
| 2FVY | Glucose-bound/closed | +12.146103 | +13.322872 | +1.176769 |

Units are kcal/mol, on the same explicitly uncalibrated reporting gauge. The
three structures remain **one GGR biological observation**. All recover the same
six donor residues, seven oxygen donors, ligand charge −3, and zero explicit
waters under the frozen selectors; La/Ca charges are 0/−1 singlets. All source
heavy-atom coordinates are preserved.

The range across the three sources is **10.268814 kcal/mol** for formamide and
**9.971030 kcal/mol** for extended amide, considerably larger than the
amide-extension effect. Sugar is outside these cores, and crystallization conditions differ;
2FW0 has sodium/citrate in the cleft and 2FVY has reported Glu149 radiation
damage. Thus this is a geometry-associated robustness result, not a causal
sugar effect or a comparison of metal-adapted La/Ca states. The experimental
GGR affinity assay used sugar-stripped protein, but that does not license
selecting whichever crystal score agrees with the assay.

## Preparation and provenance

| Protocol | Role/status |
|---|---|
| `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3` | Existing canonical baseline; unchanged inputs and released bands |
| `generic_peptide_amide_vertical_native_r2scan3c_v3` | Existing generic donor repair; archived comparisons reused, two new GGR sources scored |
| `generic_peptide_alpha_caps_native_r2scan3c_dev_v1` | New extended-amide development representation; no absolute calibration |
| `ggr_connected_segment_native_r2scan3c_dev_v1` | New connected 1GLG development representation; no absolute calibration |
| `ggr_local_sensitivity_tightscf_dev_v1` | New analytic-gradient/physical-displacement diagnostic; no mechanical score |

Existing source-graph, peptide repair, native ORCA runner, receipts, cache and
baseline collector were reused. New preparation IDs are
`generic_peptide_alpha_caps_native_r2scan3c_dev_v1` and
`ggr_connected_segment_native_r2scan3c_dev_v1`. The gradient recipe is separately
versioned as `ggr_local_sensitivity_tightscf_dev_v1`.

The connected segment retains Asp138–Gln142, with native Lys137 CA/C/O and
Phe143 N/CA terminal atoms. Overlapping source atoms merge before capping;
charged Lys137 side-chain chemistry is not imported or neutralized. Independent
inspection found no new coordinating donor (nearest added typed oxygen 5.605 Å)
and no obvious cap overlap. Source bonds, cuts, caps, charge/electron ledgers and
exact paired mappings accompany every preparation. The unaffected 25-case PQQ
panel's 50 endpoint input/XYZ pairs retain their pinned hashes.

See [preparation A](STAGE_A_PREPARATION.md), [independent audit](STAGE_A_INDEPENDENT_AUDIT.md),
[preparation B](STAGE_B_PREPARATION.md), [unrounded comparisons](COMPARISON_AB.json)
and [component records](COMPONENTS_AB.json). Immutable endpoint outputs and receipts
are linked from those records under `workspaces/ggr_mechanism_20260915/`.

## C: physical derivatives

Four analytic-gradient centers and sixteen signed physical displacements used
the approved TightSCF recipe. All actual `.engrad` files pass element/order,
coordinate, unit, ECP, charge/multiplicity, energy, source-map and execution
checks. ORCA reports analytic SCF/CPCM/D4/gCP components and the La ECP derivative;
no numerical-gradient fallback occurred. Both normal/TightSCF bridges pass the
frozen 0.05 kcal/mol criterion: ΔR is +0.000119177 for extended and +0.000117475
for connected, far below the representation shifts.

For each direction, `[E(+h)−E(−h)]/2` was checked against `h × dE/dq`, separately
for La, Ca and R. All **12/12** comparisons pass the frozen
`max(0.02 kcal/mol, 5% × |predicted odd change|)` tolerance. Maximum absolute
residual is **0.000586557 kcal/mol**. No block triggered C3; **zero** half-step
calculations were needed or executed.

| Representation / physical coordinate | dR/dq | Actual odd ΔR (kcal/mol) | Residual (kcal/mol) | Even ΔR (kcal/mol) |
|---|---:|---:|---:|---:|
| Extended / metal, ±0.02 Å | +7.781294 kcal/mol/Å | +0.155185 | −0.000440 | −0.063085 |
| Connected / metal, ±0.02 Å | +3.473536 kcal/mol/Å | +0.068914 | −0.000557 | −0.065077 |
| Extended / peptide, ±1° | +6.691518 kcal/mol/radian | +0.116421 | −0.000368 | +0.016815 |
| Connected / peptide, ±1° | +5.938247 kcal/mol/radian | +0.103427 | −0.000215 | +0.016556 |

The metal moves toward the frozen Gln140 carbonyl oxygen; the peptide rotation
uses the frozen CA140→CA141 axis. Both preserve the selected donor/water/charge
state. Maximum heavy-atom movements are 0.02 and 0.030148 Å, within the frozen
0.05 Å limit. The two motions leave cap anchors fixed, so this validates these
physical derivatives, not an end-to-end moving-cap force calculation.

The peptide derivative of R is similar across representations; the metal
derivative changes appreciably. The tested motions produce sub-kcal/mol local
effects. Neither measures the protein's restoring forces or supplies a relaxed
affinity. Even terms are curvature-sensitive, not tests of anharmonicity or
extra score corrections. In particular, La's peptide even terms are negative
(−0.010722 and −0.007299 kcal/mol); a stable mechanical basin has not been
established, and a second derivative along this curved coordinate is not a
Cartesian Hessian eigenvalue.

Full endpoint projections, even/odd terms, residuals per coordinate unit, source
gradients and center bridges remain in the sensitivity collection. No relaxation
or entropy score is enabled: **`response_model_not_validated`**, numerical
corrections unavailable. No Hessian, force-field stiffness or uncertainty
covariance was invented. See [executed integration](GRADIENT_INTEGRATION.md).

## Execution and measured cost

| Stage/job | Completed endpoints | Batch elapsed | Allocated core-seconds | Batch peak RSS |
|---|---:|---:|---:|---:|
| A / 1199770 | 10/10 | 516 s | 33,024 | 14,231,164 KiB |
| B / 1199802 | 8/8 | 292 s | 18,688 | 8,130,340 KiB |
| C / 1199805 | 20/20 | 2,224 s | 142,336 | 20,967,644 KiB |

Total allocated quantum-job cost is **194,048 core-seconds = 53.902222
core-hours**, with zero GPU time. The jobs overlap; their elapsed times must
not be added and called end-to-end turnaround. No failures/retries occurred.
Summed endpoint receipts contain 11,241.978805 wall-seconds and 179,871.660880
assigned rank-seconds; allocation accounting also includes idle ranks and
wrapper overhead.
The conditional half-step branch was not triggered. **45 distinct real-artifact
software/integration tests passed**, with four repeated after all gradient
centers were available; scientific endpoint and derivative-check counts are
reported separately in [TESTS.md](TESTS.md).

These are actual Slurm accounting values, not scheduling estimates. Each job
used four 16-rank workers on 64 CPUs; GPU time is zero. Endpoint-pair sums differ
from concurrent batch elapsed time. GGR's extended pair used 248.632615 summed
endpoint-seconds / 3,978.121840 assigned rank-seconds; connected used 948.715536 /
15,179.448576. The archived repaired pair used 264.288805 / 4,228.620880 on a
64-CPU node with the same 16-rank endpoint layout. The connected pair's observed
endpoint-cost ratio is about 3.59×; identical processor hardware and contention
were not independently established, so this is not a controlled hardware
performance comparison. These results do not establish an equally cheap
replacement. This is research accounting, with no spending or time stop.

Stage A preparation consumed 3.904089 wall-seconds across five models. Stage B
source protonation/preparation consumed 11.641096 wall-seconds and 11.247494
process CPU-seconds; peak process RSS was 213,520 KiB. Source acquisition,
software tests and Stage C packaging were not independently timed; those costs
are unavailable, not zero.

## Recommendation and next scientific decision

**Retain the baseline as production default.** The larger fragment supplies a
useful sensitivity finding, not sufficient evidence for promotion or broad
La/Ca affinity discrimination. This is not a rejection of further development.

The most informative next proposed test is the connected representation on the
two already selected alternative GGR sources, 2FW0 and 2FVY: **four endpoints**,
using their existing pinned protonated geometry, the same topology rule, native
normal-SCF r2SCAN-3c/CPCM policy, and unchanged donors, charge and waters. This
would complete the source-by-representation comparison and determine whether
the −7.34 kcal/mol connected shift is peculiar to 1GLG. No structure or sign
would be selected as the winner. Those four calculations are outside the
current agreement and have **not** been prepared or submitted; they require a
new scientific decision. No curvature model, geometry optimization, hydration
change or global-model retry is implicitly included.

## Delivered artifacts and replay

- [Derived 68-record benchmark ledger](../../workspaces/ggr_mechanism_20260915/report_v1/benchmark/benchmark_manifest.json): all historical fields preserved; nine new development pairs appended to three existing protein groups; twenty sensitivity endpoints retained separately on GGR, with null corrections and full execution accounting.
- [All Stage C energies and physical checks](../../workspaces/ggr_mechanism_20260915/report_v1/collection_c.json), [actual execution cost](../../workspaces/ggr_mechanism_20260915/report_v1/execution_cost.json), and [figure/data manifest](../../workspaces/ggr_mechanism_20260915/report_v1/plots/plot_manifest.json).
- Standalone PDF figures: [representation ladder](../../workspaces/ggr_mechanism_20260915/report_v1/plots/ggr_representation_ladder.pdf), [source comparison](../../workspaces/ggr_mechanism_20260915/report_v1/plots/ggr_source_representation.pdf), [alpha minus GGR](../../workspaces/ggr_mechanism_20260915/report_v1/plots/alpha_minus_ggr.pdf), [physical derivative checks](../../workspaces/ggr_mechanism_20260915/report_v1/plots/physical_gradient_checks.pdf). Matching PNGs and exact-value TSVs are alongside them.
- [Compact final status](RESULT.json), [tests](TESTS.md), [benchmark export details](BENCHMARK_EXPORT.md), and [exact replay command](COMMANDS.md).

The full replay script ran successfully into `report_v1`, without new quantum
calculations. Quantum energies and artifact hashes reproduce exactly. Recomputing
mapped/projected gradient arithmetic on the login versus compute node differs
only at floating-point last bits (maximum absolute numeric-field difference
1.07×10⁻¹⁴); all decisions are identical. Both independently pinned collections
are retained. No values were rounded to force matching.

Numerical credibility is established for the executed endpoint arithmetic and
these two physical derivative directions. Predictive improvement from a larger
fragment is **not established**. Its observed cost is higher, and no production
change is justified by the current evidence. The prior environmental/global
routes remain archived; no new environmental scalar was computed here.
