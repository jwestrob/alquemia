# Strict native energy stopping on the exact five failures and Ca partners

Root explicitly authorized20calls after reviewing the completed uniform32
continuation: “five failing La cells plus their exact Ca partner at the same
geometry/media =10cells. For each, run TolE=1e-10Eh native mixer/300K/MaxIter500/
rank1 from(a) its exact cold GBW+xtbw seed and(b) its exact uniform pass2 seed:
20calls total, no third pass ... compare starts with frozen0.1kcal/component
tolerance ... no score/reference repair from only these cells.”

## Fixed question, cells and method

Does a much tighter energy stopping condition remove the transient plateaus
seen in the existing native-SCF traces? Selection is all five cells failing the
previous predeclared0.1kcal pass1→2 gate, plus each exact Ca counterpart:

| Source | Shared candidate geometry | Medium | Metals |
|---|---|---|---|
| c5b120-pqq-la_model |adaptive_La|ALPB(water)|Ca,La|
| bbl57595.1-pqq-la_model |origin|vacuum|Ca,La|
| p15279-pqq-la_model |adaptive_Ca|vacuum|Ca,La|
| p38539-pqq-la_model |adaptive_Ca|vacuum|Ca,La|
| q4w6g0-pqq-la_model |origin|vacuum|Ca,La|

Every one of these10cells gets two independent strict calculations, seeded
from its actual cold and uniform-pass2 GBW+xtbw pairs. No cross-metal,
cross-geometry or alternative-state seed. Source coordinates/state are copied
exactly; Ca/La nuclear correspondence obeys the existing1e−12Å mapping policy.
No water, protonation, fragment, charge, multiplicity or parameter change.

Use the same ORCA executable and nativeGFN2 mixer, electronic300K, MaxIter500,
one rank, matched-basename AutoStart. Add **only `%scf TolE 1e-10`** to the
existing continued recipe. Require actual printed TolE1e−10Eh, native mixer
and XTBRESTART; a requested but ignored setting is a failed capability check.
No ordinary-SCF backend, ConvergenceTight bundle, different smearing,
additional pass, automatic retry or geometry search. The older six-context
nativeTight test used1e−8Eh and did not contain these exact failing cells.

## Frozen comparison and interpretation

Retain both strict energies, full iteration traces, charge vectors, parameter
and electronic-state audits, and compare to each exact cold/pass2 result.
Acceptance is |E(strict from pass2)−E(strict from cold)|≤0.1kcal/mol for every
one of10cells. Report all10 and any unsupported/nonconverged attempts. A failed
start leaves its two-start comparison unavailable. Do not choose the lower
or more favorable result. Report signed Ca−La component contrasts for each
shared geometry/medium separately; these are not complete solvent corrections
or a full classifier score.

The native mixer implements its own stopping logic and overrides other SCF
settings. Printed density limits do not automatically define its operative
criterion. Retain density residuals and actual energy-stop behavior; even
repeatable strict energies do not qualify forces or prove a unique electronic
ground state. No labels, reference, classifier pool or default are changed.

## Execution

Twenty immutable native scalar tasks,20 concurrent one-rank workers/20CPUs,
48GiB host memory, CPU-only existing explicit GPU-partition host, noGPU.
Exact saved starts are available; no cold recomputation. Recent native calls
were seconds each; tighter stopping may require more iterations and can fail
at the unchanged500-iteration bound. Measure actual allocation/runtime and
retain failures. Zero new MACE/DFT/optimization and no project budget cap.
Stop after these20attempts. A broader numerical-policy decision is separate.
