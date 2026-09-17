# Installed-source rotation audit

The pinned graph-longrange 0.4.4 implementation constructs each dipole from
three positive-axis displacements and a compensating central charge:

```
q_a = p_a / h at r + h e_a, for a = x,y,z
q_center = q - sum_a p_a/h at r
```

The checkpoint retains h = 0.1 Angstrom for field features and 0.02 Angstrom
for electrostatic energy. Fields also use forward differences at displaced
receiver positions. The Gaussian radial kernel and retained self terms are
unchanged by the memory implementation.

This matches the primary [tagged backend source](https://github.com/WillBaldwin0/graph_electrostatics/blob/v0.4.4/graph_longrange/realspace_electrostatics.py),
classes `RealSpaceFiniteDiffereneEnergy` and
`RealSpaceFiniteDifferenceElectrostaticFeatures`. Actual installed-source hashes
are pinned by the software manifest/backend inventory reused in this audit.
MACE 0.3.16 explicitly permutes Cartesian coordinates to [y,z,x]; its
`compute_total_charge_dipole_permuted` maps coefficient columns [3,1,2] back
to Cartesian [x,y,z]. The audit uses that same convention for density and
feature rotations, and maps rotated Cartesian forces back with F_rot @ R.

## Why this can depend on orientation

For a smooth potential f, the one-sided displaced representation gives

```
sum_a p_a [f(r + h e_a) - f(r)] / h
= p . grad(f) + (h/2) sum_a p_a d_a^2 f + O(h^2).
```

The leading dipole term has the expected rotational behavior. The finite-h
remainder introduces higher moments tied to the chosen axes. Rebuilding the
same dipole along fixed laboratory axes after rotating a molecule generally
produces a different finite charge cloud. A pure coordinate rotation of the
whole auxiliary cloud, including its axes, should preserve the module result.

This is a mathematical explanation of a possible source of error, not a
numerical attribution for the complete neural model. The paired core and
frozen-density tests establish how much of the observed failure it explains.
The physical geometry is unchanged by the coordinate transformation; no
protonation, charge, water or donor inventory is changed.

## Scope of a remedy

A future analytic Gaussian monopole/dipole kernel could evaluate derivatives
of the same radial kernel without axis-dependent finite displacement. That
would change this approximation and needs separate validation of feature
normalization, self terms, epsilon regularization and learned model behavior.
Shrinking h blindly risks subtractive cancellation and changes the model's
field inputs. Neither remedy is implemented or scored in this investigation.
