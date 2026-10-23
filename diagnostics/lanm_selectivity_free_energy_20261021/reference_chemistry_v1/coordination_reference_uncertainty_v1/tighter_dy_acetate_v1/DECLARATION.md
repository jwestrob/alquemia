# Three-role tighter Dy-acetate attribution scout

Root authorized this separate preparation after the completed 18-cell TightPNO
qualification failed three original individual CP comparisons. Those results,
the 0.25 kcal/mol diagnostic limit, and the unreleased W8/W9 cluster stage remain
unchanged. This scout cannot qualify the full 18-control panel or release a cluster.

## Question and finite inventory

Does the previously specified tighter local tier reduce the CCSD and/or triples
contribution to the largest observed CP discrepancy, Dy bidentate acetate?
Exactly three new single-point energies are declared: the existing complex,
ion with ghost ligand, and ligand with ghost ion. No geometry, state, ECP,
orbital/AuxC basis, canonical target, frozen-core convention, iterative T1,
SCF/CC/iterative-triples convergence threshold, D4 convention, or parser
recovery convention changes. No forces, optimization, new donor, or cluster.

The existing Tight settings were TCutPairs=1e-5, TCutPNO=1e-7,
TCutMKN=1e-3, TCutDO=5e-3. Use the existing tighter tier:
TCutPairs=1e-6, TCutPNO=1e-8, TCutMKN=1e-4, TCutDO=5e-3.
TightPNO and DLPNO-CCSD(T1) keywords remain. Preserve automatic triples-specific
behavior; record its actual TNO, strong/weak domain, iterative-selection and
scale settings explicitly. Do not assume that generic tighter cutoffs prove
converged triples or silently introduce a new triples setting.

## Resources and execution

Root alone releases and submits the exact three inputs after source/preparation
review. Eight MPI ranks, one thread each, standard_expand, default memory/time;
no memory/time overrides. MaxCore=floor(0.75*actual granted MiB/8), with positive
actual Slurm memory required. The 32,000 MiB/3,000 MiB-per-rank template is only
a resource example; execution records the actual grant. The original three
Tight jobs used four ranks and 4.630 CPU-h total. That is a measured baseline,
not a runtime/memory bound for this new tier. No retries or automatic extension.

Each role has one immutable execution directory and at most one native call.
Missing/failed roles remain in the three-slot denominator, and no partial CP
value is admitted. Root retains all release and subsequent scientific decisions.

## Saved comparison and interpretation

CP = E(complex)-E(ion with ghost ligand)-E(ligand with ghost ion). Report the
new total, CCSD and triples CP components, their changes from the preserved
Tight values, and their differences from the unchanged canonical values.
Retain RHF equality at the original 1e-6 Eh gate and identical same-parent
AuxC assignments. Report the original |CP error| <=0.25 kcal/mol diagnostic;
do not reinterpret this three-role diagnostic as the full qualification gate.

Movement concentrated in triples supports sensitivity of that approximation;
movement in CCSD supports sensitivity in the singles/doubles approximation.
Little movement leaves correlation fitting and other local approximations
unresolved. Neither outcome identifies a unique mechanism or guarantees a
useful cluster result. Preserve native SCF/core/localization/T1 warnings and
the lack of physical-state, basis-limit and reference-accuracy qualification.
