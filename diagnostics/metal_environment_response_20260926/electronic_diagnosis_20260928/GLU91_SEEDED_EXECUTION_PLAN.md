# Repaired Glu91 finite response: same-metal seeded execution

Declared 28 September2026 after accepted repaired origins, before displaced energies.
Question: does ±1degree intact Glu91 carboxylate rotation produce the differential energy response predicted by the repaired physical gradients, and can common-pool selection change the conditional metal contrast?

The existing glu91_motion_v1 real source and four-cell manifest define coordinates, membership, caps, environment, proton/water/spectator states. Reuse exact admitted repaired origins. Evaluate La/Dy at both displacements; no optimization or new geometry. Repaired normalized differential gradient2.7169243kcal/mol/angstrom supports the previously selected direction. Individual angular gradients1.3512089La/11.3711441Dy kcal/mol/radian remain reported, not only their difference.

Create a new MORead manifest with initial orbitals from each metal's accepted repaired origin. Validate exact electronic recipe, ordered physical/cap maps, state and checkpoint provenance using existing orbital-seed tooling. Keep unchanged default convergence thresholds and four-endpoint executor. Do not reuse origin energies as displaced results or silently recover failed cells. Source Dy origin used TRAH, target ordinary same-metal MORead; this is explicit numerical initialization, not a Hamiltonian change.

Four concurrent86MPI tasks on idle344CPU high-memory node, exclusive mem0, runtime MaxCore from allocation; normal priority. Existing completed original-H matrix took5314s, but nearby converged seeds may shorten this; no guaranteed completion before shutdown. Start before10PDT cutoff; preserve any partial outcomes. No repeat origins or third conformer. Completion and health watchers plus dependent collector must be armed.

Compare centered energy derivative at ±1degree with actual origin gradient, individual works and Dy-minus-La response. Use existing numerical derivative tolerance where applicable and report raw residuals; this finite-step motion is not a stationary minimum or qualified full hybrid force surface. Classical contributions must be matched before finite-energy conclusions. Common pool comprises origin, minus, plus with identical composition; no candidate populations or affinity. FNR source-reversal test remains separate and required.
