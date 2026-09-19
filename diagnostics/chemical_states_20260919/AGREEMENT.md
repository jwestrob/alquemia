# Balanced water/carboxylate proton-location test

Jacob authorized the competing chemical-state track in parallel. Root approved
this concrete test: all three source water/carboxylate-O contacts within 2.60 Å,
two metals and four new points per contact, 24 native analytic-gradient endpoints.
No further permission is required within this scope. Existing baseline remains
unchanged. No new biological labels or broad scanner rescore.

Question: does the existing neutral-water/deprotonated-carboxylate assignment omit
a competing proton location with a metal-dependent electronic stabilization?
The exact reaction is `COO− + H2O <=> COOH + OH−` inside the same prepared complex.
Atoms, total charge, water-derived oxygen inventory and metal are conserved.
Thus the proton/water reservoir coefficients are exactly zero. This experiment
does not estimate a pKa, pH occupancy, entropy or full chemical-state population.
Native CPCM(Water) settings remain unchanged; no externally assigned pH is used
to weight states. The source protein protonation protocol remains frozen.

## Freeze before endpoint energies

Use the full-water contextual 1F6S and 6IP9 native endpoints already included in
the completed occupancy table. Both are bovine alpha-lactalbumin, one protein
group. Hans-LanM's present three cores have no waters and are not altered.
PQQ zero-water inputs and default decisions are unaffected.

- Select every variable-water O / retained ASP or GLU carboxyl O pair with
  source O–O distance <=2.60 Å. No selection by energy or desired discrimination.
- Three contacts result: 1F6S A211→Asp87 OD2; 6IP9 A310→Asp82 OD1 and Asp88 OD1.
- Select the donor H nearest that acceptor in each actual endpoint; ties resolve
  by source atom index. Equivalent H labels can therefore differ across metals.
- The proposed acid H lies 0.98 Å from the acceptor toward donor-water O. Accept
  only source C–O–O_water angles 90–130 degrees, acid H at least 0.70 Å from
  another H and at least 0.75 Å from nonparent heavy nuclei. These are geometric
  input gates, not fitted energetic acceptance thresholds.
- A straight Cartesian interpolation can cross the water oxygen when the
  starting H points away. Instead interpolate radial distance linearly about
  donor O while following the shortest spherical arc to the final direction.
  Native path fractions are 0.25, 0.50, 0.75, 1.00; fraction zero reuses its
  actual archived energy/analytic gradient. Only this one physical H moves.
- Retain all 24 expected new endpoints in the denominator. Unsupported geometry
  is explicitly recorded and is not replaced by another contact or path.
- Use native r2SCAN-3c / CPCM(Water) / DefGrid3 / TightSCF / EnGrad, singlet,
  original endpoint charges, 16 ranks per endpoint, four concurrent endpoints.
  No MACE model, numerical DFT gradient, Hessian or optimization is required.

Report each path's native energy change, analytic projected slope and components,
and `ΔE_transfer,Ca − ΔE_transfer,La`. Positive differential transfer stabilization
means this proton arrangement is relatively more favorable for La. Report all
points; no favorable-state picking or calibrated classifier call. Gradients are
derivatives, not forces. Retain their original Cartesian units and path mapping.

This fixed-scaffold path can identify a candidate proton arrangement or a high
cost along the prescribed route. It cannot rule out a separately relaxed acid/OH
basin. A low-energy candidate requires a separately defined physical relaxation
and actual selectivity test, not automatic occupancy assignment. Bound-state
quantum/thermal/nonpolar corrections remain missing.

Compute estimate from the prior full-context native pair and wider-basin checks:
six batches of four 70/76-atom analytic-gradient endpoints, approximately 15–40
minutes on 64 allocated CPU cores; actual receipts supersede this estimate.
No arbitrary project compute/time cap; the finite manifest and scheduler policy
bound execution. Failed tasks remain visible and are not silently replaced.
