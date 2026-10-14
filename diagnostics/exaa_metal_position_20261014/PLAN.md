# ExaA: test the contribution of frozen metal placement

Authorized under Jacob's autonomous discriminator-improvement scope. Root owns
execution and integration. Production, calibration, previous results and other
sessions remain unchanged. This plan precedes preparation and new evaluations.

## Question

Does the source-dependent position of the fixed ion relative to PQQ explain a
substantial part of the different donor-motion responses in ExaA samples 0 and 3?
The completed eight-cell native DFT comparison reproduces MACE's larger target
response and smaller control response. Adding the previously selected packing
fragments had little differential effect. Neither result establishes faulty
folding or correct experimental chemical-state representation.

The donor-only search holds both PQQ and the metal fixed. The target's markedly
uneven Glu contacts and different PQQ contacts motivate testing that constraint.
An earlier joint-metal search on different cases added no class benefit; we will
test this specific mechanism before implementing another optimizer.

## Frozen construction and denominator

Use only the consumed two-Ca ExaA sample 0 control and sample 3 target from
`exaa_donor_response_20261014/inventory_v1/INVENTORY.json`, with its pinned
192-atom origins and saved `adaptive_La` donor configurations. Retain the current
dry oxidized-PQQ3− hypothesis, source/cap mapping, proton inventory and remote
spectator exclusion. Charges remain Ca −3 and La −2, physical singlets.

1. Identify the same 24 named PQQ heavy atoms through the actual source maps.
   Reject incomplete, duplicated or inconsistent identity mappings.
2. Align donor-source PQQ to host-source PQQ by the unweighted rigid
   least-squares proper rotation (Kabsch, no reflection). Use all 24 named atoms;
   do not choose a subset or weight atoms based on response.
3. Transform the donor source's original ion point into the host frame. Change
   only the host metal coordinate, applying that identical point to the host's
   saved origin and `adaptive_La` arrangements. All nonmetal coordinates remain
   exactly those of each original host arrangement. No geometry averaging or
   donor transplant occurs. Repeat the reverse transfer.
4. Require the existing source/cap, bonded/overlap and physical-displacement
   checks; the transferred metal must stay within the existing 0.8 Å limit from
   its host origin. The unchanged donor arrangement must also retain its
   existing admission. Do not relax a gate to admit a transplant. Alignment
   residuals and actual displacements are reported, not fit to a desired score.
5. Score both metals for both arrangements in each host: **eight new MACE
   vacuum cells**. Reuse the exact eight compatible original cells. A failed
   required cell leaves that host's differential unavailable. Preserve the
   declared eight-cell denominator and all failures.

These are modeled perturbations of real prepared structures, not observed
structures, equilibrium samples or experimental specificity labels. Metal
substitution is the only element change. Paired-coordinate equality, charge/
spin conventions and unchanged host membership must be checked before scoring.

## Energy and interpretation, frozen before evaluation

For each metal and placement, W_M = E_M(adaptive_La) − E_M(origin).
D = W_Ca − W_La is a La-directed change in the Ca-minus-La contrast when positive.
Report each host's delta_D = D_transplanted − D_original, and the target-minus-
control contrast before/after transplantation. Do not compare raw energies across
host compositions or Hamiltonians, attach class bands, or change calibration.

The narrow reciprocal-placement hypothesis gains support if the target's
delta_D is below −1 kcal/mol and the control's delta_D is above +1 kcal/mol.
The 1 kcal/mol scale is the existing interpretive resolution, not a validated
error bar or a biological accuracy criterion. Report magnitudes and individual
works regardless of this descriptive outcome. Near-zero effects close the
narrow hypothesis; asymmetric effects are mixed evidence. Failed geometry or
matrix cells leave the causal comparison unresolved. No desired class sign is
required, and either transferred position may still be physically inappropriate.

## Execution and next gate

All numerical preparation, physical checks and model calls run on Slurm. Reuse
the proven seven-thread persistent CPU MACE worker, exact model/runtime/precision
and archived qualification from packing job 4091. No extra qualification forwards
are planned when those identity pins match; incompatible pins block reuse. The
previous larger-context worker used 243.43 s for 14 calls, including 1.90 s model
load; eight calls are expected to take minutes, not measured until execution.
Use current shared scheduler limits and actual allocation. No GPU, GFN2, DFT,
new environment, optimization or state enumeration is part of this first stage.

Root reviews actual preparation and independent implementation checks before
executing the eight cells. Preserve exact commands, model/source/software hashes,
wall time, allocation and execution receipts; arm completion wakes. Successful
matrix collection is not evidence of a better classifier.

If reciprocal placement is implicated, the next gate is the qualified corrected
frozen-charge solvent response on the same admitted cells, before any bounded
joint ion/donor search. That continuation needs a separately recorded finite
manifest and decision based on these results. If this test rejects the narrow
hypothesis, do not enlarge motions or tune thresholds to rescue it.

## Ownership

- `prepare/`: Curie, exact mapped perturbations and existing physical gates.
- `scoring/`: Heisenberg, pinned warm MACE evaluations and saved-cell comparison.
- `review/`: Gauss, independent mapping/gate/energy review.
- Root: this plan, integration, Slurm submissions, result interpretation and logs.

Candidate products belong under `workspaces/exaa_metal_position_20261014/`.
