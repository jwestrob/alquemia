# Native derivative qualification — declared 27 September before new energies

Question: are the native analytic derivatives of the **embedded electronic
component** consistent with its finite energy changes and physical cap mapping?
No full QM/MM energy, relaxation, new classifier or ML comparison is claimed.
Use both completed 1H4I environment-A endpoints from scout1219207, with all
source preparation and the original Hamiltonian unchanged. The successful
Ca/La response contrast alone did not qualify derivatives.

## Finite matrix: twenty new analytic-gradient endpoints

For each Ca/La endpoint:

1. Four signed hydroxyl displacements: ±0.001 and ±0.0005Å arc length along
   the actual Thr159HG1 rotation around CB→OG1. Rotate the H exactly; normalize
   the angular step by its actual perpendicular distance to the axis. This
   preserves the inherited bond length and angle. Compare the central energy
   derivative with native point-charge gradient dotted into the unit tangent.
2. Four signed boundary displacements: ±0.001 and ±0.0005Å movement of the
   actual omitted Glu177CA. The direction is the Cartesian axis least aligned
   with CB→CA, projected perpendicular to that bond and normalized. Recompute
   its link position through the source map and preserve the original constant
   XYZ-serialization offset. Project the core link gradient using the omitted
   atom Jacobian, plus any actual represented point-charge contribution.
   The source boundary map removes CA's charge; verify this instead of assuming
   it. This is a derivative with respect to a physical source coordinate, not
   an independent artificial-cap motion or a proposed relaxation coordinate.
3. One exact repeat of A, and one simultaneous rigid transform of every core
   atom and external charge: rotate around z by0.37radian and translate by
   (0.173,0.117,0.231)Å. Compare energy and transformed core/external gradients.

Total2metals×(4+4+1+1)=20. Reuse the two original centers and their actual analytic
gradients; do not recompute them under a hidden new reference. All calculations
remain native ORCA6.1.1 r2SCAN-3c/DefGrid3/TightSCF/EnGrad/NoAutostart, physical
singlets and original charges/ECP. No numerical full DFT gradients or Hessians.

## Acceptance: unchanged from initial PLAN.md

- Repetition energy residual≤1e−7Hartree; gradient maximum residual≤1e−6
  Hartree/bohr, retaining both core and external derivatives.
- Joint transform energy residual≤1e−5Hartree; maximum rotated gradient
  residual≤1e−4Hartree/bohr.
- Each directional finite difference and step-refinement change must agree
  within0.05kcal/mol/Å +0.5% of the analytic projected derivative magnitude.
- Record all residuals and unavailable cells, not only pass counts. Passing
  these checks validates derivatives of this electronic component, not the
  missing classical interactions or its biological predictive value.

No constant-potential shift is approximated with an arbitrary cloud of charges.
The native finite-point-charge interface fixes zero at infinity; a separate
scalar-potential backend gauge test remains unavailable until that interface
exists. The algebraic Q×delta_phi test is separately labeled.

## Resources, execution and wake

Target one exclusive344CPU node,344actual MPI task slots and `--mem=0`.
Twenty independent17-rank workers use340slots. Per-rank MaxCore derives from
full scheduler node RAM with25%operational headroom. Original16-rank scout
versus17-rank repeat also tests the relevant numerical parallel-layout change.
This is resource use, not a measured scaling advantage. Count all attempts.

At declaration, other-session PQQ chain1219217→1219218→1219219→1219220 is pending.
Keep this research behind that chain, recheck before submission, and retain
the startup guard against newly arrived PQQ work. No test-partition science.
Record the finite manifest and actual dependency/submission receipts. Attach
an afterany collector and the demonstrated `codex queue` completion watcher
to both owned worker/collector IDs before ending the turn. No email needed.

If native checks fail, diagnose the specific residual before any expanded-region
or LanM energy calculation. Preparing an expanded-region map is not executing
that dependent molecular experiment. Exact MACEPOL-EF weights/engine remain
blocked as documented in upstream/AUDIT.md.
