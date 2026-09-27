# Remaining consumed LanM EF3 sources: preparation only

Hans8FNR and Mex8FNS are ready for a separately authorized electronic-response comparison after the first Hans8DQ2 scout is interpreted. No molecular endpoints, Slurm jobs, optimization or additional source search occurred here. Frozen Hans8DQ2 inputs and implementation remain untouched.

The existing pre-scoring `lanm_preparation_feasibility/RESULT_v2.json` determines every selected atom, water, peptide boundary and perturbation. `PREPARE_TRANSFER.py` generalizes that exact source export and the reviewed residue-local heavy-neighbor charge policy; `CHECK_TRANSFER.py` verifies the resulting real fixtures. Native ff19SB/TIP3P topology export makes no Context or energy call, adds no atoms, and attaches charges to archived normalized coordinates.

| Source | QM atoms | Field rows | Source/QM waters | QM charge | Field charge | Total charge | Actual spectator identities |
|---|---:|---:|---:|---:|---:|---:|---|
| Hans8FNR |198|2008|179/3|−1|+9|+8|Dy EF1,EF2,EF4|
| Mex8FNS |209|1840|167/7|0|+2|+2|Nd EF1,EF2,EF4|

Spectators retain source positions and formal+3 monopoles in both exchanged-metal states. They are not changed to La, deleted, assigned invented ion force-field parameters, or treated as responding electrons. Target LaIII is singlet; target DyIII is the physical sextet hypothesis. All-electron counts La/Dy are862/871 for Hans,878/887 for Mex. ECP and spin-orbit/reference-state qualification remain the runner's responsibility.

Hans uses the same selected loop83–94, flanking peptide units and supporting65–66 amide as 8DQ2; waters314/345/369 are retained in QM. Mex uses homologous loop84–95, C/O83 and N/H96 boundaries, supporting66–67 amide, and waters320/331/358/373/376/394/434. Full-loop formal charges are−4 and−3 respectively; adding target+3 gives the declared QM charges. Both regions use four fixed-length C–H1.09Å/N–H1.01Å source caps with mapped analytic Jacobians.

The boundary rule removes selected QM charges plus the four boundary CA field charges, retains those real CA atoms in the physical mapping, and restores each partial residue's original charge among its two actual exterior heavy neighbors. Hans increments reproduce the prior local rule exactly (within floating arithmetic). Mex increments per recipient are−.01120e for Leu83 N/CB,−.07260e for residue96 C/CB,−.01120e for Leu66 N/CB, and−.09230e for residue67 C/CB. Native source residue identities and all charges are retained in the exact export. This is local charge conservation, not dipole conservation or a native ORCA charge-shift equivalence.

A preserves each normalized source; B applies the same +2° CA→CB rotation to Asp85 (Hans) or homologous Asp86 (Mex), moving CG/OD1/OD2/HB2/HB3 only. The A/B external field is byte-identical within each source. Each source has its own composition and state; absolute electronic totals must not be minimized or directly contrasted across source regions.

Seven real-artifact checks pass per source. Maximum bond-length changes are1.6e−15/5.6e−15Å; maximal displacement.08499/.08360Å. Closest moved nonbonded contacts excluding1–2/1–3 neighbors are2.22765/2.31059Å, both longer than their source contact. Analytic cap Jacobian residuals are9.3e−11/3.9e−10. Source coordinates, state parity, all physical-atom bookkeeping, spectator identities and fixed fields pass. These are preparation checks, not electronic force or biological qualification.

## Inputs and commands

- `workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Hans8FNR/INPUTS.json`, SHA256 `28ac59e2dc444502cef21b2c11852ee25a59f89edb5d010e4cc7e64af9ad6be1`.
- `workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Mex8FNS/INPUTS.json`, SHA256 `235708eb3a6460cb592fda434484baa5f44bd8698f1ff9f928df20d7277f11a2`.

Reproduction uses a fresh output directory from repository root:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/metal_environment_response_20260926/lanm_ef3_preparation/PREPARE_TRANSFER.py --repository "$PWD" --source-id Hans_8FNR --output workspaces/metal_environment_response_20260926/lanm_ef3_transfer_reproduction/Hans8FNR
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/metal_environment_response_20260926/lanm_ef3_preparation/PREPARE_TRANSFER.py --repository "$PWD" --source-id Mex_8FNS --output workspaces/metal_environment_response_20260926/lanm_ef3_transfer_reproduction/Mex8FNS
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/metal_environment_response_20260926/lanm_ef3_preparation/CHECK_TRANSFER.py --repository "$PWD" --inputs workspaces/metal_environment_response_20260926/lanm_ef3_transfer_reproduction/Hans8FNR/INPUTS.json --output workspaces/metal_environment_response_20260926/lanm_ef3_transfer_reproduction/Hans8FNR_CHECK.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/metal_environment_response_20260926/lanm_ef3_preparation/CHECK_TRANSFER.py --repository "$PWD" --inputs workspaces/metal_environment_response_20260926/lanm_ef3_transfer_reproduction/Mex8FNS/INPUTS.json --output workspaces/metal_environment_response_20260926/lanm_ef3_transfer_reproduction/Mex8FNS_CHECK.json
```

Root owns any reference execution and continuation decision. Complete classical cross interactions, responsive spectators, full solvation/free energies, dimerization and state populations remain unavailable. No inference about which metal these proteins prefer follows from prepared inputs.
