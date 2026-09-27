# Native selected-direction force qualification — 27 September 2026

Declared before any new endpoint evaluation. Reuse the actual completed Ca/A
and La/A scout energies, analytic core gradients, point-charge gradients and
execution receipts. This tests the declared embedded electronic component,
not a full protein Hamiltonian or a relaxation coordinate search.

Exactly 20 new native r2SCAN-3c/DefGrid3/TightSCF/EnGrad endpoints:
for each Ca and La, two modes at ±0.001 and ±0.0005 Å (eight endpoints),
plus one exact repeat and one joint rigid transformation (two endpoints).
All states, charges, field charges, atom inventories and unrelated coordinates
remain unchanged. Root owns submissions, behind PQQ 1219217–1219220 (refresh
actual queue before launch). Requested explicit layout: 20 concurrent workers,
17 MPI ranks each, on the 344-CPU full-memory allocation; the existing runtime
renderer derives MaxCore from actual allocated memory and reserves 25%.

MM mode: source Thr159 HG1 rotates around its actual CB→OG1 axis. Arc coordinate
h has units Å; angle is h divided by the hydrogen's perpendicular radius to the
axis. The analytic derivative is the point-charge gradient dotted into the unit
rotation tangent at A. No hydroxyl normalization or inventory change.

Boundary mode: the actual omitted Glu177 CA moves along a unit vector obtained
by projecting the least-aligned Cartesian axis perpendicular to CB→CA. Rebuild
its one mapped H cap using the recorded 1.09 Å source-based map, retaining the
constant vector offset between exact map and archived serialized cap at A.
Only a genuinely represented MM CA would also move; this source removes its
MM charge. The analytic derivative combines the cap gradient through the
omitted-source Jacobian and any actual MM CA gradient. This is a physical
coordinate derivative of the mapped component, not an admitted independent
protein conformational degree of freedom. No synthetic H moves independently.

Rigid test: rotate every core and field coordinate by +0.37 rad about Cartesian
z, then translate by (0.173, 0.117, 0.231) Å. Transform both core and field
analytic gradients with the same rotation. Repeat uses byte-identical geometry
and point-charge files. Native scientific input stays identical.

Frozen tolerances inherited from ../PLAN.md: repeat energy 1e-7 Eh, component
gradients 1e-6 Eh/bohr; rigid energy 1e-5 Eh, transformed component gradients
1e-4 Eh/bohr. Each centered finite-difference derivative and step-to-step change
must agree within 0.05 kcal/mol/Å + 0.5% of absolute analytic projection.
Report signed derivatives and residuals, never only pass/fail. Unavailable
required cells block qualification. No full numerical gradient, optimization,
new biological label, or model comparison is implied.
