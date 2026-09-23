# SLSQP precision pilot passes and removes the long boundary tail

**The single stopping-tolerance change preserves all four pooled decisions and
candidate choices while removing Q89GY2 La's pathological optimizer tail.** All
frozen geometry, native-energy, solvent-component and score gates pass. This
supports the numerical policy on these four consumed sources; it does not
increase biological accuracy or change the existing full225 protocol/default.

## What changed, and why

A private instance of the existing optimizer uses `optimizer_ftol=1e-8`
Hartree-equivalent units instead of `3.6749322179347735e-11`. Nothing else in its
settings changed. The objective remains shifted native MACE energy with its
analytic gradient. The four common selected angles, q0 start, 200-iteration
limit, +/-0.8 rad box, coupled final 0.8 A displacement bound, 1e-7 A final
numerical tolerance and chemical/final-geometry guards are identical. Shared
modules/settings and pinned earlier inputs were untouched.

The installed SciPy `_slsqp_py.py` documents that ftol governs multiple stopping
checks: constrained optimality, constraint violation, step length and objective
change. It is not an energy-only convergence switch. At the original Q89 La
endpoint, a displacement of 0.80000000018 A already met the independent physical
acceptance rule, yet the saved tail contains tiny angle changes at identical
printed MACE energies. The original in-flight audit pins 30 tail evaluations
with angle ranges below 7e-7 rad and zero saved-energy spread. This is consistent
with excess numerical settling at the active constraint; the test does not
identify a unique internal SLSQP failure mechanism.

The new Q89 La endpoint is still boundary-limited (0.80000000420 A), within the
unchanged physical tolerance. It is not an unconstrained minimum. Four of the
eight old and new endpoints have the same boundary flags; no boundary is hidden.

## Actual search changes

| Endpoint | Old evaluations | New evaluations | Old wall, s | New wall, s |
|---|---:|---:|---:|---:|
| Q89GY2 Ca | 16 | 15 | 4.179 | 7.695 |
| Q89GY2 La | 1304 | 16 | 364.093 | 4.217 |
| Q9Z4J7 Ca | 17 | 16 | 3.887 | 3.960 |
| Q9Z4J7 La | 15 | 14 | 3.124 | 2.850 |
| 1H4I Ca | 19 | 17 | 9.756 | 6.169 |
| 1H4I La | 19 | 18 | 6.552 | 6.075 |
| A0A3F2YLY8 Ca1, Ca endpoint | 76 | 16 | 21.914 | 4.202 |
| A0A3F2YLY8 Ca1, La endpoint | 13 | 11 | 3.608 | 2.810 |

Eight search wall times sum to **417.114 -> 37.976 s**. Unique actual proposal
MACE calls fall **1471 -> 115**; each search additionally reuses its exact q0.
Evaluation counts include that q0 lookup. The new pilot uses eight further
cross-metal MACE calls. Wall times include actual runtime variability and the
first warm-worker evaluation; this is not a controlled universal speed factor.
Q89 Ca is slower despite one fewer evaluation, illustrating that distinction.

Maximum new-versus-old own-endpoint native energy difference is
**1.76312e-6 kcal/mol**, below the frozen 0.001 gate. Maximum actual Cartesian
atom displacement is **0.000585934 A**, at 1H4I Ca; all eight mapped coordinate
and angle differences are retained. Q89 La changes by only 3.57454e-5 A with
native energy difference 6.04035e-8 kcal/mol. All final mapping, fixed-atom,
charge, cap, bond and overlap guards pass.

## Full scoring comparison

All 32 fresh native GFN2 cells and eight cross-MACE cells completed. Each score
uses the common q0/adaptive_Ca/adaptive_La pool and the same native
MACE + GFN2(ALPB - vacuum) expression. All q0 values were reused. No standalone
xTB, additional continuation, DFT, chemistry change or reference refit occurred.

| Source | New minus old pooled R, kcal/mol | Old/new decision |
|---|---:|---|
| Canonical Q89GY2 | -0.0124666080 | La-supported / La-supported |
| Canonical Q9Z4J7 | -0.0002018480 | Ca-supported / Ca-supported |
| Crystal 1H4I | +0.0012942858 | Ca-supported / Ca-supported |
| A0A3F2YLY8 Ca-conditioned sample1 | -0.0005092325 | La-supported / La-supported |

All pooled differences are below the frozen **0.2 kcal/mol** gate. Every named
new-versus-old GFN2 vacuum/ALPB cell difference is below **0.010202 kcal/mol**,
passing the independent **0.1 kcal/mol** component gate. Subtracting the shared
q0 makes these exactly the corresponding work differences. The largest solvent
transfer difference is 0.0125291 kcal/mol; no component cancellation conceals a
failed cell. The complete MACE/vacuum/ALPB decomposition is preserved.

Mathematical and operational selections remain unchanged. Operational Ca/La
choices are respectively: Q89 adaptive_Ca/adaptive_La; Q9 adaptive_La/adaptive_La;
1H4I adaptive_La/adaptive_La; A0A3Ca1 adaptive_Ca/adaptive_La. Original frozen
union/adaptive bands are used only as an explicit numerical-transfer check.
This new numerical version has not acquired a newly fitted calibration.

## Execution, testing and disposition

Proposal/cross job **1210561** used 75 s x 32 CPUs, with one requested H200 and
200000 MiB host memory. Native solver job **1210572** used 95 s x 64 CPUs,
128 GiB host memory and zero GPU. **Total: 8480 allocated core-seconds and
75 requested GPU-seconds.** There were 115 new proposal MACE calls, eight
cross-MACE calls, 32 native GFN2 calls, eight reused native q0 endpoints and
16 reused q0 solvent cells. No failed molecular attempts or retries occurred.

The solvent job was routed while still pending, from the explicit 128-CPU
GPU-partition host (only 52 CPUs free) to the established H200 CPU-sharing host.
The same job ID, tasks, rank count, memory, zero-GPU request and priority stayed
unchanged; pre/post scheduler records are preserved. No duplicate job or active
allocation was changed. Actual executor times and scheduler RSS are in COSTS.
The proposal worker reports peak CUDA allocation 5,580,981,248 bytes (5.20 GiB)
and peak host RSS 1,862,476 KiB; its cross-scoring worker is slightly smaller.
These process readings complement the scheduler RSS records. Local setup/
reporting and historical reused calculations are additional costs.

Four real-artifact tests pass in 5.268 s, zero skips. They verify exact input/
selector reuse, private module isolation with only one changed setting, the
actual pathological source and final comparison algebra. The scientific result
test was explicitly unrun before execution. No fabricated successful fixture.

**Recommendation:** retain this versioned stopping policy as a qualified small-
pilot numerical improvement for a deliberate future integration. It removes
wasted search work without material score changes here. Do not rewrite prior
results or change the running full225 protocol. No additional tolerance,
restarts, calibration or production change is made in this branch.

Protocol: `union_four_angular_native_OMOL_SLSQP_ftol1e8_v1`.
Full results: `workspaces/slsqp_precision_20260923/COMPARISON_v1.json` and
`SEARCH_AUDIT_v1.json`; all exact endpoint receipts and coordinates remain pinned.
[Plan](PLAN.md), [commands](COMMANDS.md), [compact result](RESULT.json).
