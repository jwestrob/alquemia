# Anchored strain and PME mesh phase — source-only finding

The existing metal-anchored strain changes absolute fractional coordinates.
This can add a translation derivative to the reported cell virial when the
finite PME discretization is not exactly invariant to continuous translation.
It is a concrete explanation to test, not an established cause from source alone.
No arithmetic, molecular program, preparation, build or new solver was run here.

## Exact path and chain rule

`displaced_preparation_v1/PREPARE.py:93–103` defines, for each diagonal strain,
`r(q) = r + ε S(r−a)` and `L(q) = (I+εS)L`, with `ε=q/L0`, `L0=1Å` and
`a=(2.5405,1.8185,−1.86)Å`. For yy strain:

```
y(q) = (1+ε)y − ε a_y
Ly(q) = (1+ε)Ly
y(q)/Ly(q) = y/Ly − ε a_y/[(1+ε)Ly]
```

`pmestuf.f:251–293` computes the spline/grid position from absolute reciprocal
coordinates, `nfft*(w−anint(w)+1/2)`, with no compensating metal-anchor shift.
`grid_mpole` (726 onward) and `grid_uind` (861 onward) use those spline weights.
Continuous translations generally change the discrete spline-assigned structure
factor; integer-cell invariance does not imply continuous translation invariance.

The origin-affine path `r0(q)=(I+εS)r` instead keeps fractional coordinates fixed.
The two paths differ only by the uniform translation `−ε S a`. With positive
Cartesian energy gradients `g_i=∂E/∂r_i` and `G=Σ_i g_i`, the exact chain rule at
q=0 is:

```
D_anchor = D_origin − (S a)·G/L0
```

If the native virial supplies `D_origin = (V:S)/L0`, the anchored comparator is
`(V:S − (S a)·G)/L0`. The correction has a minus sign. For all three existing
diagonal directions it is `−a_axis*Σ_i g_i,axis/L0`, not just the yy correction.
The inherited collector `projection` at lines131–144 uses only `V:S/L0`.
Native reciprocal virial accumulation in `empole1.f:3781–3855,3909–3980` and
`epolar1.f:8617 onward` contains lattice, multipole and frame contributions;
there is no supplied metal-anchor parameter. This supports checking the missing
path term, while not proving the complete reciprocal virial correct.

## Shortest next check

First use the already saved central full-gradient tokens to evaluate the anchor
correction for all three diagonal directions, with the original two steps and
0.001 kcal/mol/Å gate unchanged. Preserve the original failed result. Propagate
printing uncertainty from both virial and `Σg` (triangle sum weighted by `|S a|/L0`),
and identify any ordinary arithmetic allowance separately. A corrected comparator
would qualify this discrete derivative identity, not physical translational
invariance, active-response accuracy or dynamics. Total gradients cannot assign
separate mpole/polar anchor corrections without corresponding saved term gradients.

If the saved check is inconclusive, one bounded diagnostic is four new yy states
at the same ±1e−4 and ±5e−5Å steps using `y0(q)=(1+ε)y`, identical strained boxes,
cutoff12, alpha0.4,64³ mesh,order5 and inactive-response input. Reuse the existing
central state and four anchored outputs; use the admitted term-observer binary
without a rebuild. This requires four gradient programs, no new factorization.
At each step compare anchored versus origin-affine term slopes and the total
chain-rule prediction. Removal of the nearly constant electrostatic offset only
on the fixed-fractional path would support the mesh-phase/path explanation;
persistence there would direct attention to the underlying virial derivative.
No states or execution are authorized by this note.
