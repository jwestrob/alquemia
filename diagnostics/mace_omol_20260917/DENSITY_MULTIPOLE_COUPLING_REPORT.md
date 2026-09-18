# Actual density coupling to native permanent multipoles: numerical checks pass

Eight saved-density ORCA potential queries completed as job1201063; three
native initialization/rotation-only exports also passed. No DFT, MACE,
induced-response or force calculation was added. These are component checks
on consumed GGR/alpha development structures, not biological predictions.

At every real exterior atom, evaluate the native vacuum endpoint potential
at its center and at the declared 36 offsets. Spatial potential Hessians use
0.01 and 0.005bohr steps; these are not nuclear Hessians. Actual quantum electric
fields come from the separately qualified finer field experiment. The direct
coupling is `sum(q*phi - mu*E + Q:Hessian(phi))`. Convert native e-A moments to
atomic units once. Native `kmpole` already divides quadrupoles by three; retain
the full symmetric tensor contraction without another factor.

All eight endpoint and four paired refinement checks pass. Maximum endpoint
quadrupole change is **0.000206909396kcal/mol**, maximum paired change
**0.000074223742kcal/mol**, versus declared 0.02/0.01 tolerances. Center potentials
reproduce the earlier utility within1.5231e-13au; rigid contraction residual is
4.441e-15kcal/mol. Native local/global moment reconstruction, charge and axis
checks pass. Actual rounded native quadrupole traces reach9.334284e-7eA²;
preserve and report their small trace contraction rather than silently forcing
traceless tensors. The native point-kernel comparison accounts for this known
convention difference; no exact quantum-density/native-kernel identity is claimed.

The unchanged exterior residue charges in this component experiment are not
charge-closed QM/MM boundary charges. Do not treat the component sums as a
complete correction. The following boundary preparation supplies that missing
piece without rerunning these density observations.

## Actual cost and tests

Job:846wall seconds on8allocated CPUs, **6768allocated core-seconds**,
**2951.569reported CPU-seconds**,1551296KiB peak batch RSS,zeroGPU. All8utilities
succeeded first try. Exact per-utility wall/CPU/RSS receipts and paired component
values are in [DENSITY_MULTIPOLE_COUPLING_RESULT.json](DENSITY_MULTIPOLE_COUPLING_RESULT.json).
Preparation23.501669s wall/19.442915s processCPU; native build0.499663s
wall/0.463864s CPU, and the three native exports are separately receipted.
Costs are additional to earlier DFT/charge/field work, not total production cost.

Five real-artifact/algebra tests pass10.237s, no skips. The analytic-monopole
finite-difference test originally demanded factor-four improvement even when
already roundoff limited. Its assertion now includes the independently derived
floating-point cancellation bound. This changed no scientific gate or output
and required no scientific rerun. No fabricated quantum output was used.

Protocol `vacuum_density_native_AMOEBA_permanent_multipole_coupling_diagnostic_v1`.
Manifest0747bd0b0c58996dbaf659ea29bdc8b13d718101a96cdeed8517995a09fa6d80.
Products: `workspaces/mace_omol_20260917/density_multipole_coupling_v1`.
Baseline unchanged; no new reference, calibrated decision, force or relaxation.
