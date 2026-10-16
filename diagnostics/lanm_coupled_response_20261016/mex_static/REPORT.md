# Mex whole-system model assembled; response qualification is next

The repaired Mex8FNS EF3 source now has an executable local-electronic,
protein-mechanical and solvent energy with gradients for both La and Dy.
Job12685 completed in13s on12CPUs and oneH200; application11.566s.
Exactly two new local MACE energy/force calls and two global fixed-charge ALPB
calls ran. Four accepted local GFN2 and two classical origin records were reused.
No new SCC, classical evaluation, DFT or optimization occurred.

The finite system contains2049 physical atoms,209 local electronic atoms
including four caps,1840 external field sites and three retained Nd spectators.
Its energy is MACE_vac + GFN2_emb - GFN2_vac + MM_ledger + ALPB(qFF,R).
Read ../MEX_STATIC_PLAN.md for the interaction ledger, charge/state distinctions
and boundary rules. All geometry, composition, source and raw reuse gates pass.
The initial12673 failure and saved-only12680 diagnosis are retained separately.

This is capability, not a preference prediction. Unlike-metal raw MACE totals
contain element offsets and must not be read as binding selectivity. Total-force
directional checks and the matching archived PBE0 response comparison are next.
No stationary structure, stable electronic state or whole-protein affinity is
established. The model is conditional on fixed EF3 substitution and three Nd
spectators; the whole-protein solvent uses fixed force-field charges, not
self-consistent quantum populations or a fully polarizable environment.

Actual receipt and all raw component pins are in RESULT.json. Component timings
and allocation accounting are separate in the workspace and ACCOUNTING.txt.
No production, PQQ calibration or prior failure was altered.
