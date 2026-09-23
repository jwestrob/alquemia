# Approved nine-source six-angle warm-start pilot

Frozen before new molecular evaluations, 2026-09-23. Parent authorization:
“Proceed with this exact9-source six-mode pilot under Jacob's discretionary
overnight authorization.” This implements the previously inventoried proposal;
it does not change the original four-angle protocol or production defaults.

Exact inputs/selection: `workspaces/c5ax_response_20260923/SIX_MODE_FEASIBILITY_v1.json`
SHA256 `d580f4bf6ede0b5fd312602475832b4e29c6b481c7e5475dcbac26702e6b18af`.
All five C5AX La sources, then A0A3 Ca1, A0AC Ca4, A8 La1 and A8 La3,
in that inventory's exact order. All are consumed development cases.

Question: does adding the next two independent paired-origin force-selected
angles yield useful differential accommodation beyond the existing four-angle
result? Continue the identical alternating differential/individual normalized
force rule and Gram–Schmidt threshold1e-5; keep the original four-ID prefix.
No forced donor identity or label-based mode choice.

Each Ca/La search starts from its own archived four-angle proposal, added angles
zero. Reuse both q0 and warm-start MACE energies/full Cartesian forces, projecting
all six gradient components. Preserve source chemistry, context, metal/PQQ/waters,
maps/caps and all other coordinates. SLSQP/native-vacuum-MACE objective and exact
gradient; Hartree-equivalent scale, ftol1e-8,200iterations, one start. Every angle
remains within ±0.8rad relative to originalq0; final maxheavy displacement0.8Å
plus existing1e-7Å numeric tolerance. Intermediate infeasibility recorded; source
bond/cap/overlap guards retained. No clipping, artificial energies or retries.
Final success needs optimizer success, physical checks and native nonincrease
against originalq0 under the existing1e-7eV allowance; warm-start added work is
reported separately, without a new outcome gate.

The common pool retains originalq0, oldCa, oldLa, newCa and newLa, scored under
both metals. Exact108old strict scalar cells and all oldMACE cells reused.
At most18new single-start searches,18crossMACE and72fresh nativeGFN2 cells;
search MACE evaluations are data dependent. No newq0,DFT,entropy,water inventory,
protonation, fold, solvent-force search or calibration. Use freshNoAutostart
TolE1e-10/300K/MaxIter500/nativeGFN2/rank1 with matchedvacuum/ALPBwater from the
qualified strict32 recipe. No seeds or favorably selected continuation stage.

All required five-geometry cells must succeed for an extended pool result;
missing candidates/cells remain unavailable, with the old result separately
visible. Report mathematical minima and existing0.1kcal-origin operational
selection separately. Transfer frozen strict32 bands (REFERENCES.json
SHA171bd31d466ff97ef6073ce23286af8519ef23690f7d6ad23bb444cab30ecd5c), with no new
nine-source calibration. Record each metal's added native/solvent/composite work,
selected candidates, score changes, classes, geometry extents, residuals and costs.
Row minima cannot increase by discarding old candidates because none are discarded;
that does not guarantee beneficial change in their difference.

Finite resources: oneH200/32CPU/200000MiB warm worker; CPU-only32one-rank workers,
32CPU/64GiB for the72cells. Existing18four-angle searches took65.55summedseconds/
260objective evaluations;108strict scalar receipts5001.52summedrank1seconds.
This is a minute-scale batched development test by prior measurements, not an
assured wall time. Every actual attempt and allocation is recorded.
