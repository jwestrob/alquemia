# A bounded correlated check of the consumed W8/W9 contrast

## Decision

Use a **26-single-point DLPNO-CCSD(T1) panel**, conditional on numerical and
approximation controls. Do not launch a full-cluster canonical calculation by
default. This tests the raw fixed-geometry contrast
`C = E(DyW9) - E(LaW9) - E(DyW8) + E(LaW8)` in the existing orbital basis and
frozen-f ECP representation. It does not replace the PBE0 targets or qualify
forces, collective B, sampling, or physical accuracy.

| Stage | New energies | Purpose |
|---|---:|---|
| Exact archived water/acetate-bidentate/acetamide, La/Dy, complex and two CP fragments | 18 | Compare local correlation with the already completed canonical results |
| Four consumed full W8/W9 clusters, TightPNO | 4 | Raw C method sensitivity |
| The same four clusters, fixed tighter local cutoffs | 4 | Cluster-local approximation sensitivity |

The first cluster is W9 Dy (434 orbitals, 98 explicitly correlated electrons),
after all small controls pass. Root reviews its actual runtime, memory, scratch
and native admission before the other seven cluster calculations. The scout is
one of the 26, not an additional cell. No automatic fallback, retry, threshold
selection or full canonical extension follows a failure.

## What is already available

The archive contains six admitted central canonical-water CP energies, twelve
canonical donor CP energies, and 24 displaced canonical-water energies. None is
a W8/W9 cluster coupled-cluster result. The canonical water geometry differs
from the literal cluster dimers (metal–oxygen distance about 2.288 versus 2.500 Å);
it is an approximation control, not a substitute for a consumed endpoint.

Saved canonical/PBE0-D4 CP donor-versus-water paired contrasts differ by
−1.401006 kcal/mol for acetamide and −0.254895 for acetate. The latter barely
misses the old 0.25 working resolution. These results motivate a matched
cluster-method check but provide neither a universal correction nor an error
bound on C. T1 reaches 0.023312701 in the saved donor set; retain this warning.

The original raw PBE0-D4 C is −0.534206667 kcal/mol. Its original D4 contribution
is +0.106638333. The common-parent-basis PBE0 C/P/B is
−0.209483747/−9.245663037/+9.036179290; the corresponding model B is +2.797567940.
Those calculations test basis-partition sensitivity, not correlation-method
uncertainty. New CC energies must never be combined with old PBE0 fragments and
called a CC CP or B decomposition.

## Why local correlation, and its limits

ORCA6.1.1 documents `DLPNO-CCSD(T1)` as iterative triples; plain
`DLPNO-CCSD(T)` uses the semicanonical variant. Local correlation requires an
AuxC fit basis. TightPNO uses iterative local-MP2 preparation. Its cutoffs and
the tighter tier are explicit in SPEC; they are sensitivity settings, not a
canonical-limit certificate. Conventional RHF integral handling remains
available. NoFrozenCore includes every electron remaining outside the ECP;
native core/valence partitions still require inspection. [Official MDCI manual](https://orca-manual.mpi-muelheim.mpg.de/contents/modelchemistries/mdci.html).

Use light-atom def2-TZVP/C and the existing explicit metal AutoAux expansions,
changing only the auxiliary-slot label from AuxJ to AuxC. This preserves their
exponents and contractions on real and ghost centers. AutoAux is documented for
correlation fitting, and AuxC is distinct from AuxJ. The old J-fitting admission
does not itself establish the new correlation-fit accuracy; the canonical
controls test the combined RI/local approximation. [Official basis manual](https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/basisset.html).

The installed ORCA driver, MDCI serial/MPI and localization modules are pinned
in EVIDENCE.json. Their presence and official options establish implementation
feasibility; the first actual local result must establish native execution and
parser compatibility. No installed molecular executable was invoked here.

## Cost evidence and decision rule

The 12 saved canonical donor jobs cost 57,948 allocated CPU-seconds (16.10 CPU-h).
The 201-orbital, 40-electron La acetamide scout took 2,303 scheduler seconds on
four ranks; MaxRSS was 16,383,484 KiB, close to the 16,000 MiB grant. This is not
evidence that larger canonical or local jobs fit that grant.

A stage-wise engineering extrapolation of the measured donor non-triples MDCI
and triples timings, using `o²v⁴` and `o³v⁴`, gives about 53–56 four-rank hours
per W8 full complex and 104–111 per W9. Four full canonical cells total about
1,295 CPU-h; a 12-role canonical CP cycle about 2,307 CPU-h. These are rough
cost projections, not measured runtimes or memory bounds, and are why the
canonical shell is not the recommended next test.

There is **no measured DLPNO runtime in this archive**. A compact hydration
cluster may retain large domains, so a claimed 100-fold speedup would be
unsupported. Start with the small controls under the existing four-rank default
grant; use the largest small control's actual timings and memory to decide
whether the W9 Dy scout is affordable. Its actual cost then bounds the decision
to release the remaining seven, not their runtime mathematically. Preserve
failed/time-limited attempts and do not substitute a cheaper Hamiltonian.

## Scope of any result

Four same-method full energies are sufficient for raw C because matched free
water and own-ion terms cancel. A CP C at that method needs 12 roles. A
common-parent-basis C/P/B needs 80 roles; no exact extra cancellation is known.
An own-basis B would instead need 40 energies and answer a different basis
partition question. None of those larger panels is included in these 26 calls.

Small-control agreement and target threshold stability support a bounded method
sensitivity statement. They do not bound residual domain/RI errors, basis
incompleteness, omitted explicit f/SOC/multiplet physics or single-reference
failure. Even a stable contrast cannot establish a unique charge-transfer
mechanism. Report electronic PBE0 and PBE0+D4 separately against CC without D4.

Exact archive pins, all 18 control identities and measured timings are in
EVIDENCE.json. Numerical decisions and finite scheduling order are in SPEC.md.
