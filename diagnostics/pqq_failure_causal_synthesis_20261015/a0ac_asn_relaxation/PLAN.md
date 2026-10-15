# Finite Asn290 alternative-basin relaxation

## Question and immutable scope

Does bounded physical relaxation around the already admitted terminal-amide flip change the fixed-flip response or relieve its represented short contacts? This is a developmental basin test, not a class-targeted adjustment. Preserve all ten A0ACD6B9F2 sources in the existing inventory order, all20 completed original-orientation searches, and all80 actual original/common3 plus fixed-flip endpoint cells. No labels enter construction, mode selection or optimization.

Exactly20 new searches: ten flipped origins × Ca/La. Each source retains its exact original four angular mode IDs/order for both metals; La-conditioned samples0/2 therefore retain χ2 without adding χ1. The frozen `slsqp_precision.engine()[0].completion.optimize` uses the original SLSQP binary, maxiter200, Hartree-equivalent objective/gradient and ftol1e−8. Bounds remain ±0.8rad and maximum0.8Å physical-heavy displacement, now explicitly relative to the admitted flipped start. This is a separate domain, not a claim that the flip lies within the old source domain. No additional starts, adaptive mode selection, clipping, springs, parameter fitting, geometry rescue or automatic retries.

## Preparation and gates before energy calls

Create a derived physical map by grafting only OD1/ND2/HD21/HD22 into `context.positions_A` and setting `context.core_positions_A` to the exact admitted204-atom flip. Preserve all source identities, H parents, bonds, mode axes/order/membership, physical exterior, caps and cap-link/Jacobian definitions. Keep the old immutable map/source and explicit graft provenance. The existing Kinematics must reproduce the admitted start, fixed atoms/caps and paired Ca/La coordinates; analytic Jacobian and covalent checks run on Slurm before optimization. Reject unsupported moving boundary members or cap anchors.

Use all20 cached actual4617 MACE energies/forces for q0 and project those forces with the derived Jacobian. No q0 model repeats. Replay all80 prior endpoint input/state/model/force receipts, the20 prior search settings and four-mode identities, and actual runtime qualification. Full-source identities/H remain in the same geometry-only union as preparation4614; no new atoms, protonation or exterior energy. Retain existing relative severe-overlap policy1.0Å heavy/.55Å H, source/cap/bond checks and final convergence/nonincreasing-native-energy gates. Apply the source-union check to proposed physical coordinates without a new radius or steric threshold; preserve inherited contacts and report final named H contacts. No independent H or synthetic-cap coordinates.

## Execution order and common pools

Root submits all numerical work. First run preparation/preflight on Slurm. Then execute only the first declared source (Ca-conditioned sample0), both metal searches and required opposite-metal cross cells, under one unchanged seven-thread warm CPU worker. Root inspects actual execution, map/gradient fidelity, contact admission and costs before releasing the remaining18 searches. Each stage is immutable and exclusive-created; interrupted calls/searches are retained and cannot silently restart.

Reuse each successful own-metal final result and cross-score its exact physical geometry under the other metal: at most40 new final endpoint cells, with model trial calls recorded separately. Retain common3 (`origin,adaptive_Ca,adaptive_La`), common4 (plus`amide_flip`) and common6 (plus`relaxed_flip_Ca,relaxed_flip_La`) side by side. Use the existing chooser, original-first tie order and strict−0.1kcal/mol operational work rule; mathematical minima remain separate. No partial-pool minimum or favorable-source subset. All ten rows and every failed proposal/cell stay visible.

Later native frozen-charge ALPB is a separate stage, not submitted here: only actual complete/credible common6 candidates and a meaningful sign-blind common6-versus-common4 paired change warrant root review. Use the existing1kcal/mol absolute mathematical-change criterion for that review flag; never use biological correctness or the already-observed common3→4 change to trigger it. No calibration/default change.

## Cost and interpretation

No fixed optimizer-forward count is promised: maxiter200 does not bound line-search evaluations. Record every model request/forward, cached call, optimizer termination, boundary flag and actual worker/scheduler resource receipt. Existing4617 measured20 forwards/329.435s worker time; old same-source CPU searches used12–15 function evaluations in roughly154–199s. These are planning anchors only. No MACE minimum implies a composite stationary point, equilibrium rotamer, folding correction or validated class improvement.
