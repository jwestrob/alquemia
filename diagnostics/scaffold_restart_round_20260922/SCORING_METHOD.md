# Collective proposal comparison

The new opt-in scorer protocol is
`nikasha_collective_scaffold_geometry_native_OMOL_GFN2_ALPB_v1`.
It adds admitted source-connected Cartesian proposals to the existing shared
pool; it does not change the released scorer or copy a historical reference
onto a newly validated numerical scale.

For each of the same eight sources there are three required donor targets:
original, Ca-adaptive and La-adaptive. Each target has its own protein-only
ff19SB proposal. All three physically admitted proposals are cross-scored by
both metals, with exact existing charge, multiplicity, cofactor/proton inventory,
context membership, cap mapping, native float64 OMOL and native GFN2 MaxIter500.
The native-restart diagnostic remains separate; its electronic seeds are not
mixed into these unchanged scoring endpoints.

E_M(q) = E_OMOL,vac,M(q) + E_GFN2,ALPB,M(q) − E_GFN2,vac,M(q)
R = min_q E_Ca(q) − min_q E_La(q)

The existing operational0.1 model-kcal selection tolerance and mathematical
row minima are both retained. There is no new tolerance chosen for this experiment. Larger R
is more La-like on this particular scale. No aquo-subtracted score, affinity,
probability or equilibrium population is inferred.

The adapter verifies the actual mechanics receipt and its frozen implementation,
source parent, complete atom order, fixed donor target, physical geometry checks,
nonincreasing measured parent work, and actual initial/final force evaluations.
It independently reconstructs source coordinates and offset caps. Final Ca/La
geometries differ only in metal identity. Collective proposals carry their full
parent coordinates, never fictitious old angular coordinates. The ff19SB energy
is a proposal objective and contributes zero terms to the composite expression;
that is an explicit model choice, not a missing correction filled with zero.

Primary comparisons require all three targets. Failed required energy cells leave
the expanded pool unavailable; the archived adaptive result remains separately
named. Exact coordinate duplicates reuse prior energies. No favorable partial
pool, refit, threshold padding or changed chemical state is introduced.

Report, in order:

1. Each target's paired contrast and Ca/La component work versus the same archived
   target before scaffold response. This includes the original-target control.
2. Target response minus original-target response, a descriptive comparison,
   not an exact force-field bias correction.
3. Expanded-pool selection, raw contrast, released/adaptive band transfer and
   selected-source coverage. There is no new calibrated decision model.
4. The two predeclared structural pairs: A0A3 Ca-samples1/3 and A0AC Ca/La-sample4.
5. Actual mechanics, MACE and GFN2 work, failed attempts and allocation costs.

The parent model omits direct PQQ/metal force-field coupling and solvent. Donor
constraints and overlap checks cannot replace those forces. Boundary-limited
proposals are not stationary minima. Improved parent energy alone is not useful
classifier evidence. These repeatedly inspected sources are development cases.
