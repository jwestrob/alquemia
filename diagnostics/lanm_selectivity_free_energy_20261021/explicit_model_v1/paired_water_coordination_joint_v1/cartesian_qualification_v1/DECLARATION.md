# Conditional joint-water Cartesian qualification

Source preparation only; root owns release and submissions. No current17676
trajectory/result is opened. Bind only its immutable terminal candidate passing
all ten energy and fifteen rigid-force development controls and the exact
independent selected-ten review produced by the queued17683 reviewer. This
adapter relies on that review and does not repeat its raw-ten/score machinery.

Four consumed clusters: aquo8_central_La/Dy (25 atoms each),
aquo9_central_La/Dy (28 atoms each). Copy the exact selected XYZ, key and parameter
bytes. No geometry generation, parameter fitting, donor or validation geometry.

Use the same pinned CPU testgrad and unchanged gradient parser admitted by17612:
`testgrad system.xyz Y Y 0.00001`. The input specifies total two-sided separation
1e−5 Å, each side5e−6 Å. Require every analytical/numerical atomic component and
maximum discrepancy≤1e−3 kcal/mol/Å. Require total energy replay against the
selected cell≤1e−4 kcal/mol and analytic force replay≤1e−3 kcal/mol/Å. Preserve
all four independent gates; no paired-cancellation substitute or relaxed tolerance.

The native source testgrad.f calls gradient once, then energy twice per active
Cartesian coordinate. All atoms/components are required:6×(25+25+28+28)=636
internal displaced-energy evaluations, plus four analytic gradient/energy
baselines, in at most four testgrad program invocations. This is a requested
upper bound; a failed program may complete fewer internal calls, which are not
invented. Count completed636 only from complete admitted native outputs.

Exact identical W8 inputs with already admitted17612 Cartesian evidence reuse
that pinned output; no duplicate native program. W9 has no such old Cartesian
receipt. Four single-thread workers in one four-CPU standard_expand task,
default memory/time, no GPU. First native/admission failure stops further starts;
retain completed siblings and full four-slot denominator. No retry/fallback.

A numerical pass checks these exact cluster forces, not physical accuracy,
coordination transfer, bulk sampling, or fitness of a combined donor model.
Candidate/review/root release remain null in staged preparation. Root separately
binds actual inputs and releases execution only after the stated terminal gate.
