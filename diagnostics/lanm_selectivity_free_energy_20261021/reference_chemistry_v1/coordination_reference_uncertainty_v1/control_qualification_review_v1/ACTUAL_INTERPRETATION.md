# Eighteen local-correlation controls: admitted calculations, failed qualification

The saved review admits all 18 raw calculations and reproduces all 35 original
comparisons. **Three individual CP comparisons fail the unchanged 0.25 kcal/mol
gate. The full-cluster scout remains unreleased.** All paired comparisons, all
18 RHF comparisons, and all six same-parent AuxC checks pass. No new molecular
calculation, changed threshold, or replacement result is included here.

## Where the discrepancies enter

Every number below is local minus canonical, in kcal/mol. CP is the complex
energy minus the two matching ghost-fragment energies. Components are extracted
from the actual, hash-verified local and canonical outputs, using all 18 cells
from each method. The printed CCSD/triples sum closes to the reported total
within 0.000000628 kcal/mol.

| CP comparison | CCSD contribution | Triples contribution | Total error | Original gate |
| --- | ---: | ---: | ---: | --- |
| Water La | +0.003926 | +0.014468 | +0.018393 | Pass |
| Water Dy | +0.025232 | +0.017540 | +0.042773 | Pass |
| Acetamide La | +0.075925 | +0.136483 | +0.212409 | Pass |
| Acetamide Dy | +0.117603 | +0.236178 | +0.353780 | **Fail** |
| Bidentate acetate La | +0.098360 | +0.183924 | +0.282284 | **Fail** |
| Bidentate acetate Dy | +0.187585 | +0.256112 | +0.443698 | **Fail** |

Triples contributes 58–67% of each failed comparison; CCSD contributes the
remainder. Every RHF difference is exactly zero at the saved printed precision.
The positive errors mean weaker local CP binding relative to the corresponding
canonical result. The error is larger for Dy in each donor family. The La
acetate failure exceeds the original threshold by only 0.032284 kcal/mol.

| Paired Dy-minus-La CP error | CCSD | Triples | Total |
| --- | ---: | ---: | ---: |
| Water | +0.021306 | +0.003072 | +0.024380 |
| Acetamide | +0.041678 | +0.099694 | +0.141372 |
| Bidentate acetate | +0.089225 | +0.072188 | +0.161414 |

The paired donor-minus-water errors are +0.116992 for acetamide and +0.137034
for acetate. Their passing cancellation does not override any individual
failure or establish the accuracy of a W8/W9 cluster contrast.

## Source, localization, and triples interpretation

The saved admission verifies the executed geometry, ECP, orbital basis, ghost
roles, effective state, NoFrozenCore electron counts, AuxC, declared local
cutoffs, and original/new producer identities. Positive localization completion,
LMP2 convergence, final CC residuals at or below 1e-8, and iterative-triples
convergence all pass. This provides no evidence of a wrong source, omitted
semicore electrons, accidental RHF change, unfinished localization, or simply
unfinished CC iterations. Matching RHF energies alone would not establish
identical orbitals; the other input/admission checks remain necessary.

The actual method used iterative T1. Its native triples CP discrepancies versus
canonical are smaller than the printed local semicanonical T0 intermediates:

| Failed case | Native local T0 component error | Actual iterative T1 component error |
| --- | ---: | ---: |
| Acetamide Dy | +0.560818 | +0.236178 |
| Acetate La | +0.573138 | +0.183924 |
| Acetate Dy | +0.746254 | +0.256112 |

These are saved component comparisons, not newly admitted T0 calculations.
ORCA's printed “Conventional (T0)” belongs to the local algorithm; it is not the
archived canonical CCSD(T) calculation. The official manual describes iterative
T1 as an improvement to the approximation of canonical perturbative (T), and
documents separate triples-natural-orbital truncation scales. Thus the labels do
not imply that T1 intentionally targets a different higher-order correlation
energy. [ORCA 6.1.1 MDCI manual](https://orca-manual.mpi-muelheim.mpg.de/contents/modelchemistries/mdci.html).

The remaining differences are consistent with combined local-domain/PNO,
triples-space, and correlation-fitting approximations. The accounting does
**not** identify a unique cause or demonstrate that tightening one cutoff will
fix it. All-source agreement and converged iterations do not prove a unique
localization solution or eliminate fitting error. Neither the CCSD nor triples
offset can be ignored. Native SCF printed residual flags remain as recorded;
the method is not thereby physically qualified. The Dy-acetamide T1 diagnostic
is 0.021198447 locally and 0.023312701 canonically, so that warning is not a new
artifact unique to the local calculation and does not establish multireference
character by itself.

## One proposed next reference discriminator — unexecuted

If root continues this branch, use **three new single points** for the worst
failed comparison: Dy bidentate-acetate complex, ion with ghost ligand, and
ligand with ghost ion, all at the already declared tighter local tier
(`TCutPairs=1e-6`, `TCutPNO=1e-8`, `TCutMKN=1e-4`, `TCutDO=5e-3`). Preserve
iterative T1, geometry, ECP, orbital/AuxC basis, NoFrozenCore, effective state,
SCF/CC convergence tolerances, and four-rank default-grant resource rules.

Declare this as an approximation-sensitivity scout before execution, recording
all actual triples-specific cutoffs/scales as well as the four principal
cutoffs. The generic tighter tier must not be called a proof of converged
triples. Compare both CP components against the saved Tight and canonical
values; do not add a fitted correction. A reduction mainly in the triples
component would support triples-space/domain truncation sensitivity; a CCSD
reduction would identify sensitivity in the singles/doubles approximation.
Little movement would leave fitting and other local approximations unresolved.
No outcome alone qualifies the other failed controls or releases the clusters.

The completed Tight three-role cycle used 2,753.56 / 83.16 / 1,330.21 seconds
on four ranks: **4.630 CPU-h** by runtime times allocation. This is a measured
baseline, not a tighter-method cost bound. No such new inputs, jobs, or release
are created by this note. The original 18-control qualification remains failed.

## Evidence and actual cost

- Collector RESULT: `9d17e658f751af6ddbfeeabe443afc3b17a0001757374fa2536269ea477ad2a6`.
- Independent saved raw/35-comparison REVIEW: `0d056bbcfb44514cc9a7685b0246520f536f89ab6729cdc0422d1a7e7ab6149e`.
- CCSD/triples COMPONENTS: `e1b12c5bb8457676c931b5e71a510e587ac9ec476406d0327fc805f730a02cac`.
- Native T0/T1 TRIPLES_DETAIL: `eaca7d218678662d4ca494da545dc42de8562765fbd84e7c80af315b74f3d958`.

The structured receipts reside in the matching workspace
`control_qualification_review_v1/actual_v1/`. The original 18 local calls cost
64,168.333 runtime-times-allocated-rank CPU-seconds (17.825 CPU-h), distinct from
Slurm accounting. Collector: 20.664 s wall / 18.937 CPU-s. Saved review:
10.430 s wall / 10.108 CPU-s. Component accounting: 1.623 s / 1.596 CPU-s;
triples detail: 2.662 s / 2.648 CPU-s. All saved analyses used one local thread
and zero native calls. Historical native-parser failures and their saved-only
recoveries remain preserved.
