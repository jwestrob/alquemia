# Uniform SLSQP stopping precision: approved bounded numerical test

Jacob's overnight discretionary authorization covers this contained test; root
specified the four sources and one numerical change on September 23. Freeze
this plan before new model evaluations. Existing full225 union/adaptive work
retains its original settings and is unaffected.

Question: can a stopping tolerance appropriate to molecular score precision
remove a pathological boundary-settling tail without changing useful proposals?
The completed union-context Q89GY2 La search used 134 iterations/1304 evaluations,
with its last 30 angle changes below 7e-7 rad and identical printed native MACE
energies. Its final heavy displacement passes the existing 0.8+1e-7 A rule.
This suggests excess numerical settling, not proof of a chemical basin problem.

## Exact sources and unchanged physical problem

In order: canonical Q89GY2, canonical Q9Z4J7, crystal 1H4I, and
A0A3F2YLY8 Ca-conditioned seed-1 sample1. Use the exact completed UNION contexts,
source coordinates, charge/singlets, common four-angle selector, source mappings,
archived paired q0 energies/forces and electronic settings. Q89 comes from the
canonical continuation; the other three from the original common8 union pilot.
All old proposals, failures, traces, solvent cells and references stay immutable.

Eight searches (one q0 start per metal/source), up to eight cross-metal MACE
single points and 32 native GFN2 cells (two metals x two new candidate poses x
two media x four sources). Reuse all q0 cells and own-metal proposal energies.
No DFT, new chemical states, extra starts, reweighting, tuning or retries.
Native GFN2 uses the existing paired ALPB-minus-vacuum policy and MaxIter500;
no standalone backend substitution or continuation is added.

## One numerical version

Load a private instance of the existing completion optimizer. Change only
`optimizer_ftol` from 3.6749322e-11 to exactly 1e-8 Hartree-equivalent units
(6.2751e-6 kcal/mol). The shifted native MACE objective, analytic derivative,
constraints/Jacobians, angle bounds +/-0.8 rad, final maximum heavy displacement
0.8 A (+ unchanged 1e-7 A numerical tolerance), maxiter200 and all final geometry
and chemistry guards remain unchanged. SLSQP's ftol controls several internal
checks; no independent physical-force convergence or unconstrained minimum is
claimed. Intermediate infeasible trials remain recorded; no clipping/penalty
or fabricated energy at another coordinate is introduced.

## Predeclared comparison and acceptance

Retain all eight endpoint results even if any fail. Require the original final
geometry/pair invariants, and compare actual new-vs-old native proposal energy
at every endpoint: absolute difference <=0.001 kcal/mol. Report full coordinate
and angle changes, boundary flags, iterations, evaluations and timing.

For every matched named new proposal cell, report the change in each native
GFN2 vacuum/ALPB energy relative to its same archived q0 (q0 cancels in the
new-old difference). Each difference must be <=0.1 kcal/mol. Report solvent
transfer and MACE components separately; no compensating cancellation conceals
a failed component gate. Require pooled operational R difference <=0.2 kcal/mol
for all four pairs. Show mathematical and operational selections, all candidate
contrasts and decisions under the existing frozen union/adaptive reference.
Old bands are an explicit numerical-transfer check, not calibration inheritance
or promotion of this numerical version. No threshold is fitted here.

Tail removal requires Q89 La to finish with fewer evaluations and lower actual
search wall time than its archived run. Report all other search costs; timing
is hardware/runtime-dependent and excludes no failed attempts. An accuracy or
geometry failure is a negative result, not grounds for another tolerance.

Use existing warm 1H200/32CPU/200000MiB allocation for proposals/cross scoring
and 64CPU/128GiB, eight8-rank native ORCA CPU workers for finite solvent cells.
Record allocated and molecular cost separately. No production/default change.
