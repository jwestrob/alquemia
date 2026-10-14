# First ExaA reference: matched vacuum donor work

**Use the existing native r2SCAN-3c vacuum treatment first.** It directly tests whether a second electronic model supports the large MACE-vacuum energy decrease along the exact saved donor motion. CPCM is a useful subsequent solvent sensitivity, but comparing MACE vacuum with DFT/CPCM changes the medium as well as the electronic method. Such a disagreement cannot be assigned to MACE's local electronic surface alone.

This is an independent reference, not a reproduction of MACE's training Hamiltonian or a ground-state oracle. Root owns the finite source/geometry manifest and execution; this review submitted nothing.

## Exact executable starting point

Reuse the successful PedE scalar recipe without adding a guessed basis, functional, solvent correction or empirical restraint:

```text
! r2SCAN-3c NoAutostart DefGrid3 TightSCF
* xyzfile -3 1 core.xyz
```

Use charge −2 for La; the shown −3 input is Ca. Preserve the actual saved 192-atom coordinates, atom ordering, dry oxidized-PQQ3− hypothesis, fragment/cap membership and both physical singlets. Native r2SCAN-3c supplies def2-mTZVPP, def2-mTZVPP/J, D4 and gCP; the actual La treatment uses Def2-ECP replacing 46 electrons. Calculate and check the ExaA electron inventory from its own source during Slurm preparation; do not copy PedE's 838/828 explicit-electron numbers just because both contain 192 atoms. No GFN2 parameter JSON belongs in these DFT inputs.

Pinned ORCA6.1.1 executable:
`/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca`,
SHA256 `38b5f057452fef275c0a1b98d270ad03d0411dab70453820d1c678bd732f6c83`.
Reuse the PedE pinned runner/runtime-library policy, 12 actual MPI slots and memory derived from the shared allocation. The old observed MaxCore3100 reflected a 48,000MiB allocation; it is not a new fixed memory request.

For each declared motion, report each metal's `E(q)−E(q0)` and their Ca-minus-La difference. A common displaced geometry evaluated under both metals isolates a metal-dependent response. Comparing each metal's separately selected endpoint instead measures the actual two-path accommodation contribution; name those different quantities explicitly. The current ambiguous sample3 selects adaptive_Ca for Ca and adaptive_La for La, whereas crystal A selects adaptive_La for both. A four-endpoint common-La-proposal test alone therefore does not reproduce sample3's whole operational pool score. No new DFT optimization or classifier calibration is needed.

Start with the finite target/control energy endpoints and inspect the first actual SCF results before expanding. If a coherent vacuum response is unavailable because of an electronic-state or convergence failure, record that result. Do not silently substitute CPCM and call it a matched vacuum qualification. A separately named CPCM counterpart uses the same inputs plus `CPCM(Water)`; its work difference then diagnoses environmental modification of the response. The complete CPCM endpoint already includes its dielectric term, which must not be added again.

## Existing evidence and practical limits

- Historical donor validation used **r2SCAN-3c/CPCM(Water)**, explicitly frozen in `diagnostics/accommodation_torsion_20260920/PLAN.md:21`. Its agreement with MACE on 4MAE and two compressed PLM paths is useful precedent, but not a matched vacuum qualification. Its 16 new endpoints cost 589,696 allocated CPU-seconds; those different fragments/NormalSCF inputs are not fresh ExaA predictions.
- The recent eight-cell PedE vacuum/CPCM reference completed at a molecular cost of **1,074,804 allocated CPU-seconds**, using 12 MPI ranks. Per-endpoint receipts show roughly one hour for the easier sample3 water states and over twelve hours for sample1 Ca vacuum. Job3372's **50,971 seconds includes both water and vacuum**, not vacuum alone. No fast new-site runtime should be inferred from a successful inexpensive xTB replay.
- PedE Ca vacuum outputs have positive occupied frontier energies. Several final water states are nearly degenerate, and sample3 La vacuum has an apparent non-Aufbau final ordering. These are warnings about the conditional electronic reference, not proofs of electron detachment or a unique preparation defect. Normal termination does not establish electronic-state uniqueness, strict density convergence or physical redox assignment.
- Preserve actual ending-solver evidence. The repaired collector is `diagnostics/pede_solvent_response_20261013/reference_preparation/recovery/COLLECT_V2.py:35`: it uses the final TRAH macro residual against the actual printed gradient threshold and does not require an optional/stale `Last Orbital Gradient` summary. Its original failed parsing remains preserved. SOSCF's native energy/gradient stopping conditions must likewise not be mislabeled as all-criteria convergence.
- Native `MaxIter125` was not an overall 125-cycle work cap: the unchanged long PedE attempt executed 50 ordinary DIIS iterations, TRAH through macro98, and printed 843 SCF cycles. Do not describe the same default as a guaranteed finite cost bound or claim that a rank/time/iteration override occurred. Monitor actual progress; repeated inner-solver warnings alone are not accepted endpoint energies.

## Reuse and gradients

`EVIDENCE.json` records a bounded metadata search: none of the **18 exact saved ExaA XYZ hashes** from two-Ca samples0/3 and crystal A, both metals and all three candidate geometries, appears in eight inspected native-DFT manifests. This is not an exhaustive archive claim. Those manifests cover the known environment, donor-torsion, coordination, fold and current PedE references.

All eight PedE final GBW files remain present, but none is an exact ExaA cache. They should not be used as same-state ExaA restarts merely because the atom count matches. PedE was energy-only and contains no `.engrad` artifact. Historical analytic r2SCAN-3c vacuum/CPCM support exists in `scripts/mace_omol_vacuum.py`, for other fixtures; its parser's element-number table omits sulfur, so it is not a drop-in ExaA gradient collector. Native analytic D4/gCP/ECP gradient checks are implemented there, but a new exact ExaA gradient requires an actual new EnGrad evaluation and compatible full-element/state checks.

The smallest decisive first test is therefore **fixed-geometry scalar work**, using actual saved MACE energies and native DFT totals. It can identify an excessive MACE finite-motion response without new optimization, a full force campaign or any fitted correction. If it agrees, that supports the tested local response only; scaffold constraints, source populations, redox and solvent remain separate questions. If it disagrees, inspect electronic-state admissibility before attributing the entire discrepancy to MACE.

Read-only review completed14October2026. No molecular imports, calculations, installations, submissions, scientific recomputation or threshold changes were performed.
