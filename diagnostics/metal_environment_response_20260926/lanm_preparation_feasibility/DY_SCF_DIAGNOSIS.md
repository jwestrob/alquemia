# Dy molecular SCF fails before a qualified electronic state exists

Read-only diagnosis, 27 September 2026. No molecular calculation or restart was performed here. Source: `workspaces/metal_environment_response_20260926/lanm_ef3_hans_scout_v1/`; exact inspected output hashes are in DY_SCF_DIAGNOSIS_EVIDENCE.json. Interrupted Dy outputs are not terminal converged energies; original collection remains unchanged.

## Finding and next decision

Do not extend the same large Dy runs or treat their final GBW files as qualified seeds. Both independently initialized A/B calculations diverge during the first few molecular DIIS steps. A two-start isolated, real 50-atom Hans EF3 diagnostic with the same native r2SCAN-3c Dy ECP28/basis and physical sextet is a useful smaller test. It can distinguish failure that already exists locally from difficulty introduced by the much larger embedded problem. It changes region/environment together, so success alone cannot uniquely assign the large failure to embedding rather than region size.

Use PModel and HCore only as diagnostic independent molecular starts; neither is a promised remedy. Preserve fixed chemistry and inspect spin/local occupation and convergence before any subsequent use. No exact native-r2SCAN 50-atom duplication was found in the inspected consumed-transfer records: the previous small-core transfer used MACE/GFN2, not DFT. This is a targeted search, not an exhaustive claim about every historical file.

## Direct evidence

- Actual large Dy input is charge −1, multiplicity 6, native r2SCAN-3c/DefGrid3/TightSCF/EnGrad, DoEQ false, external point charges, NoAutostart. Output reports 833 explicit electrons and 2,719 basis functions. No saved molecular guess was read.
- PModel successfully creates its neutral-Dy atomic fitting density. The internal calculation uses ROHF/def2-SVP with ECP28 and HCore, atomic multiplicity 5 / 4f10 6s2; it converges in 51 cycles (−881.638158107720 Eh). This auxiliary neutral atom is **not** the molecular Dy(III) sextet and does not validate that state. It rules out a failed atomic fitting subprocess here.
- Dy_A first molecular energy is −6552.580269036807 Eh; at cycle5 it has risen to −5661.684477334602 Eh. Cycle4 MaxDP is 2.42e3, far from a small residual. Dy_B shows almost the same immediate trajectory (cycle5 −5660.822038691345 Eh).
- AutoTRAH begins after 50 ordinary iterations. Its Dy_A first macro energy is −4192.145207601669 Eh; microstep predicted changes reach approximately −1e26. These are internal optimizer diagnostics, **not physical computed energies**. Large negative orbital gaps and orbital-gradient errors persist. No final molecular orbital-energy, spin-population, or S² result exists to diagnose a particular 4f occupation.
- The smallest overlap eigenvalue is 1.551e−5 versus the 1e−7 cutoff; zero eigenvectors are removed. The observed report does not identify linear-dependence removal as the trigger.
- La is not trivially convergent either: it also invokes AutoTRAH, and prints SCF convergence after 362/245 cycles for A/B. Its early density errors are far smaller, and its TRAH microproblems converge rather than exploding. Actual final La energies are −5730.072520781592 and −5730.072802298282 Eh. These alone give no La/Dy contrast.
- The saved Dy GBWs are about116MB each, associated with unconverged/interrupted iterates. They are evidence to retain, not successful restart states. PMIX warnings also appear, but both molecular computations advance for many iterations; those warnings do not explain the electronic pathology by themselves.

## What is and is not diagnosed

The immediate failure is molecular SCF/orbital optimization, not absent Dy ECP coverage, a pre-SCF syntax failure, or ordinary almost-converged density. A poor molecular initialization is plausible. Unsupported functional/explicit-4f behavior, orbital occupations, region/environment effects, and optimizer behavior remain possible causes. The present incomplete outputs cannot distinguish them. The promolecular fitting-density integral is not used to infer the final molecular electron count; the explicit output count is the actual state check.

The [ORCA initial-guess manual](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html) documents HCore as molecular one-electron diagonalization and warns that it often gives excessively compact orbitals. PModel instead constructs a model potential; changing the guess does not change the target Hamiltonian. A successful small-core guess would still need a documented mapping or fresh full-region initialization before an embedded calculation; merely possessing a different-size GBW is not a validated transfer.

Legacy `benchmarks/hans_lanm_dy_qmmm_correction_v1/qmmm_calibration/PATOM_ENGINE_CAPABILITY_AMENDMENT_PRECOMMIT_2026-08-05.md` records a real ORCA6.1.1 PAtom abort: stored Dy atomic orbitals unavailable. `BASIS_ECP_AUDIT.md` permits HCore/PModel in its separate PBE0/ECP28 calibration. Do not repeat PAtom or import that protocol's different Hamiltonian as if it were this test. No spin, ECP, iteration ceiling, or scientific acceptance threshold was changed by this diagnosis.
