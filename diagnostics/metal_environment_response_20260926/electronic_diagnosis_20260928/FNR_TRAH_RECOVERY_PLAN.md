# Same-Hamiltonian Hans8FNR numerical recovery

Declared before recovery. Both repaired source-A PModel/default-SCF endpoints in job1220334 diverged together: energy rose by hundreds of Hartree over the first five iterations; maximum density changes exceeded1000 and DIIS error remained~1.2 rather than entering its .2 startup domain. These are not molecular energies and no biological score has been read. Stop this owned pair rather than increasing iterations or letting unchanged fixed-point divergence consume hours. Original manifests/output/failure/actual costs stay preserved; snapshot residuals are in FNR_DEFAULT_DIVERGENCE.json.

Read-only checks found no obvious core atom overlap (shortest core contacts are actual .9572angstrom O-H bonds), no exceptional cap/field coincidence (closest~1.286angstrom like other sources), and smallest overlap eigenvalue2.093e−5 above unchanged1e−7 cutoff. This does not establish a unique electronic cause. Do not move atoms, delete waters, alter boundaries/charges, or loosen thresholds to fix numerical convergence.

Execute one separately versioned two-endpoint origin recovery La_A/Dy_A with explicit TRAH from the same PModel initial construction. All repaired Hans8FNR v2 coordinates, point charges, charge−1/816explicit electrons, physical states, effective restricted valence, PBE0-D4 basis/ECP/AuxJ, grids and VeryTightSCF convergence remain unchanged. No old unconverged wavefunction is reused. No iteration-limit increase or blind cascade of solvers. The other source jobs remain untouched.

TRAH is documented for restricted Kohn-Sham and RIJCOSX; it replaces the failing fixed-point/DIIS convergence route with orbital trust-region optimization. It is not a different functional, a nuclear Hessian, guaranteed global electronic state, or validated affinity. Require actual TRAH iterations, converged states and complete gradients. The already running repaired Hans8DQ2 TRAH residual improved substantially but is not yet accepted; count both attempts independently.

Two112MPI endpoints on the released224CPU node,mem0, same existing runner and completion/health watches. Preserve original primary outcomes as cancelled/unavailable. Report this numerical recovery separately in the paired table, never overwrite the primary manifest or present fallback as original success. No displaced geometry, additional metal or model threshold change.

Reference: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html .
