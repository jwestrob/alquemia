# Uniform SLSQP precision: canonical fidelity and two numerical failures

Root explicitly authorized this bounded expansion under Jacob's overnight
instructions on September23. Same single numerical change as the completed
four-source pilot: optimizer ftol1e−8 Hartree-equivalent, rather than
3.6749322179347735e−11. Native objective, analytic forces, matched four-mode
selector, q0 start,200 iterations, ±0.8rad and final0.8Å heavy displacement,
1e−7Å geometry tolerance and all chemistry guards remain unchanged. There is
no extra start, tolerance sweep, retry, relaxation of geometry or scorer change.

Fixed34 sources: the exact25 canonical calibration sources and3 consumed
crystals from the original UNION reference; all four noncanonical sources in
common8; Q92WY9 Ca-conditioned sample2 and Q60AR6 La-conditioned sample0. The
last two were selected because their original searches exhausted200 iterations,
not because a new outcome is favorable. Every source is consumed development.
The canonical set alone calibrates; the other9 cannot enter class extrema.

Reuse the completed exact four precision pools (Q89, Q9,1H4I,A0A3 Ca1).
30fresh sources mean60 bounded searches, at most60 cross-MACE and240 native
GFN2 endpoints. Reuse all existing q0 native energies/forces and both solvent
endpoints. Zero q0, DFT, standalone xTB, numerical DFT gradient or new structure
calls. Fixed common pool: origin and both new precision proposals, evaluated
under both metals. No failed original proposal is recycled as a successful one.

Same numerical gates as the original pilot, frozen before outputs: each matched
successful native proposal energy differs by at most0.001kcal/mol; each matched
native GFN2 medium cell by at most0.1kcal/mol; common-pool operational R by at
most0.2kcal/mol. Report all actual atom displacements, choices, mathematical
versus operational minima, decisions, boundary flags, calls and allocation even
when a numerical gate passes. These gates do not prove biological accuracy.
For old failures, successful-old comparisons remain null, and any energy or
geometry comparison against the archived final unsuccessful iterate is labeled
as such, never counted as numerical equivalence to a valid old result.

Report old-reference transfer first, then freeze a separate canonical25-only
reference using the existing extrema/minimum-gap rule independently for both
pool variants. Preserve original reference and full225 results. No promotion,
posthoc thresholds, favorable-source selection or automatic full225 rescore.
All34 denominators and all failed calculations remain explicit.

Use the same private engine as the qualified pilot and immutable source
manifests. Warm1H200/32CPU/200000MiB execution; native ORCA eight8-rank workers
on64CPU/128GiB. Native GFΝ2 MaxIter500 and all other electronic settings stay
identical. Count loading, preparation, failures and allocated overhead. The
four-source pilot required75GPU-s plus95s on64CPUs; molecular search cost scaled
roughly with cases is a development estimate only, not a promised runtime.
The concurrently queued MMOL1770 repaired-H source is a separate approved
experiment, with its original optimizer policy; no cells overlap these34.
