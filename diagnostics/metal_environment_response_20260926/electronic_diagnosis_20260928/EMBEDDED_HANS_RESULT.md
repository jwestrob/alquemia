# Embedded frozen-f Hans donor response: four complete endpoints

The frozen-f model now completes both metals in the 195-atom Hans8DQ2 electronic region and actual finite protein field. This overcomes the execution failure of the earlier explicit-f Dy treatment, but does not establish a correct LanM preference model.

The prescribed +2 degree Asp85 motion changes the electronic energy by **−0.04353110 kcal/mol for La** and **+0.06425619 kcal/mol for Dy**. Dy-minus-La response is **+0.10778728 kcal/mol**. This small one-source response is not an affinity, a classifier improvement, or a test that changing the exterior alone produces that effect: the donor itself moves inside a fixed field.

| Endpoint | SCF cycles | Torsion energy derivative, kcal/mol/radian |
|---|---:|---:|
| La A |50|−5.06783955|
| La B |49|+2.56093048|
| Dy A |61|−4.76638776|
| Dy B |60|+5.99101580|

All four terminated normally, with native analytic electronic and point-charge gradients and checked effective electron/ECP states. The force projections change sign along the prescribed motion; no stationary geometry or numerical full-hybrid derivative check is inferred from that observation. Original source C-bound hydrogen direction defects remain present in this experiment. The repaired-H/TRAH scout is separate and cannot be pooled with these endpoints.

## Exact scope and cost

PBE0-D4, def2-TZVP CHNOS, matched lcecp1-TZ La ECP46/Dy ECP55, explicit AutoAuxJ, RIJCOSX, DefGrid3, VeryTightSCF, PModel, finite embedding, DoEQ false. Both states have charge −1 and 806 explicit electrons, restricted valence multiplicity 1. Physical DyIII f9 sextet is metadata distinct from this frozen-f approximation. No spin-orbit validation, continuum solvent, or classical total energy is included here.

Worker1220312 completed in5314s on344allocatedCPUs; collector1220313 used2s on1CPU: **1,828,018 allocated CPU-seconds**, zeroGPU. This is reference-development cost, not affordable scanner throughput or measured CPU utilization. No failed restart was required. A PMIX diagnostic appeared during Dy B postprocessing, but native normal termination, successful execution receipt and actual gradient collection all completed; preserve the diagnostic rather than interpreting it alone as a failed molecule.

Authoritative collection: `workspaces/metal_environment_response_20260926/lady_frozen_embedded_hans_v2/FINAL_COLLECTION.json`; compact hashes, energies, projections and scheduler records are in `EMBEDDED_HANS_RESULT.json`. Existing inputs and outputs remain immutable.

## Decision

Do not repeat this matrix or promote it. Use its actual gradients with the already implemented physical cap/environment mapping and matching classical terms to diagnose shared scaffold motions. Finish the independently running repaired Dy scout before choosing further repaired electronic cells. The compact Hans/Mex source reversal remains unresolved; no calibrated within-series discriminator is claimed.
