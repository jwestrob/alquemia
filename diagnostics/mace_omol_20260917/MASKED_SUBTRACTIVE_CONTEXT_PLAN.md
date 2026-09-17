# Static DFT core plus masked-MACE protein context

Declared after the completed masked response screen, before calculating this
candidate's partition contrast or new whole-system energies. Active-goal autonomy
applies. This is a new static descriptor; no optimization or failed curvature
correction is enabled. Earlier POLAR/GB hybrid failures remain distinct results.

## Energy definition and purpose

Keep DFT's local chemistry and let the masked MACE descriptor supply only the
change from a capped core to its source-matched intact protein:

H_M = E_DFT,CPCM(core,M) + [T_mask(full,M) - T_mask(core,M)].
R_H = H_Ca - H_La.

Convert native Hartree DFT energies and model eV values once to their declared
kcal scales before addition. No new full solvation, permanent-charge term,
relaxation, entropy, fitted multiplier or calibrated threshold is added.
This is a fixed-weight hybrid descriptor, not a rigorous QM/MM Hamiltonian or
binding free energy. The residual CPCM-core/vacuum-learned-context approximation
and artificial core caps remain explicit. The model is neither a PB replacement
nor the old POLAR electrostatic backend.

The charge-feature mask removes the known pathological global charge alias from
the learned model. DFT replaces its local core contribution, including the poor
local response exposed by the latest screen. That motivates a static test; it
does not prove that either component's geometry gradient is a validated hybrid
force. No gradient or mechanical correction is calculated here.

## Reuse and exact state accounting

Use the same pinned mechanics preparedV2 record as the response screen. Reuse
all eight actual masked core-center results from `masked_response_v1` and their
native r2SCAN-3c/CPCM DFT references. The four representations are GGR_extended,
GGR_connected, ALPHA_1F6S and ALPHA_6IP9. They share three physical structures;
GGR representations are one biological observation, and the alpha structures
remain one qualified cross-study direction comparison. All are consumed cases.

The archived mechanics `physical_cases/.../center` full coordinates preserve the
original source H positions used by DFT, unlike the H-normalized direct-score
benchmark. Keep those exact coordinates, original charge/spin, assembly, waters
and source graph. Existing source-H bond defects remain a declared limitation;
do not silently normalize H and combine the result with old DFT cores. Verify
shared source atoms, cap mappings, source physical inventory and paired states.
Do not reuse differently H-normalized whole-chain outputs.

## Stage1: cached partition test, zero new model calls

The same GGR full-system term appears in both partitions and cancels. Therefore:

R_H,connected - R_H,extended
 = (R_DFT,connected - R_DFT,extended)
   - (R_mask_core,connected - R_mask_core,extended).

Compute this directly from the four actual GGR core endpoints, with source hashes
and each component retained. Acceptance: absolute partition shift <=2kcal/mol,
consistent with the earlier hybrid development target and small compared with
the baseline canonical separation. Do not change it after seeing this contrast.
This test does not establish predictive accuracy.

If this fails, retain the result and do not run Stage2 for this candidate. Missing
or incompatible source/receipt state is unsupported, not a zero correction.

## Stage2: six conditional whole-chain energy calls

Only if Stage1 and source mapping pass, evaluate Ca and La for each of the three
actual whole structures using the pinned float64 masked model, existing qualified
energy-only edge/product adapter and existing runner. One A5000,16CPUs,64474MiB.
Measured comparable whole forwards cost about13seconds for GGR and less for
alpha; expect a few minutes including startup, and record actual costs.
No new DFT, solvent solver, gradients or training. Six is the complete initial
new inference inventory; preparation and any failure attempts remain recorded.

Report baseline core R, learned core R, learned full R, context contribution and
R_H side by side. Compare alpha-minus-GGR directions for both alpha structures
and both GGR partitions: all four must be positive to pass this consumed
relative-order development gate. Preserve all denominators and components;
no PQQ bands, universal zero or broad affinity claim. A pass would motivate a
new uniform H preparation with fresh matched DFT endpoints and additional real
structures, not production promotion. A failure rejects this fixed candidate;
no per-site weights, cap choices or label changes to make it win.
