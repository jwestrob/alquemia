# Matched fixed-FF solvent comparison: both source signs agree

Slurm12715 completed all four declared LBT3 cells. Both source-specific balanced
exchanges favor **Mex–La + LBT–Dy** within the declared conditional model.
The common two-source pool selects7CCO for both metals, with no exact ties. It
does not conceal a source sign reversal, but the large magnitude difference
between sources remains.

R = (E_Mex,Dy − E_Mex,La) − (E_LBT,Dy − E_LBT,La), in kcal/mol.
Positive R favors Mex–La(+2)+LBT–Dy(−2) over Mex–Dy(+2)+LBT–La(−2).
Values below are rounded copies of the saved comparison; no new numerical
analysis or molecular evaluation was performed for this report.

| LBT source/menu | MACE vacuum | GFN field difference | MM ledger | Fixed-FF ALPB | R |
|---|---:|---:|---:|---:|---:|
|7CCO, La-deposited model1|−1.145156|−0.041319|+0.066701|+3.652327|**+2.532554**|
|7CCN, Lu-deposited model1|+70.904103|−0.041319|+0.066701|−2.828247|**+68.101238**|
|Common two-source menu; both select7CCO|−1.145156|−0.041319|+0.066701|+3.652327|**+2.532554**|

The7CCO sign depends on the fixed-FF solvent contribution overcoming its
oppositely signed vacuum exchange. The7CCN result is dominated by vacuum MACE;
its solvent term opposes that direction. Thus agreement of the final signs is
useful consumed-source evidence, while quantitative source robustness remains
unestablished. No favorable source was removed or chosen as the only report.

## Complete denominator and common-menu admission

|Source|Metal|Status|New fixed-FF ALPB energy, kcal/mol|
|---|---|---|---:|
|7CCO|La|complete|−321.73800250543655|
|7CCO|Dy|complete|−332.3736897620542|
|7CCN|La|complete|−283.88818198116354|
|7CCN|Dy|complete|−288.04329450919903|

Preflight verified identical ordered physical membership, terminal state, native
FF charges and XML topology, covalent graph, whole-Q partition and electronic
recipe across the two LBT sources. All four cells were available to selection.
The two-source minimum is a finite common menu, **not an equilibrium population**.
The original source-specific results remain authoritative alongside that menu.

The four native MACE12500 rows and two Mex12685 static assemblies were reused.
Only four256-atom fixed-FF native ALPB calls were new. LBT uses its existing255
native ff19SB charges plus formal+3, with whole256 as Q and partition-defined
zero field/MM terms. No peptide FF vacuum energy was added to whole MACE.
The native solvent functional/build matches Mex; only wrapper atom-count
admission changed. No old self-consistent LBT solvent offset was used.

## Limits

This is conditional fixed-EF3 Mex with three Nd spectators versus a whole LBT3
monomeric peptide tag. Mex retains its source waters; LBT is the declared dry
peptide. Local capped209-Q Mex and whole256-Q LBT use the same energy construction
with different partitions; partition invariance is not established. LBT’s free
NH3+/COO−, peptide−5/holo−2 hypothesis and the primary SI terminal-charge ambiguity
remain. Physical MACE La1/Dy6 states and the frozen-f GFN field approximation
retain their existing limitations. This is no protein Kd, matched experimental
affinity, independent blind test or classifier qualification.

## Execution and preserved failure

|Job|Outcome|Scheduler seconds|CPUs|New molecular calls|
|---|---|---:|---:|---:|
|12711 V1|failed during source-map admission|3|1|0|
|12715 V2|complete,4/4 cells|5|1|4ALPB|

Total incremental scheduler cost is8allocated CPU-seconds, zeroGPU time. V2
application time was5.029727724962868seconds. Parent computations are reused and
excluded from this incremental cost. V1 remains preserved: it confused a native
map key (`A/1//H`) with the deposited provenance name (`A/1//H1`). V2 corrects
that join and additionally requires exact deposited-record identity, retaining
the original coordinate tolerance and all scientific settings.

## Pinned evidence

- [Actual result](../../../workspaces/lanm_coupled_response_20261016/lbt_matched_model_v2/execution/RESULT.json): `2644f6d8bf3dff9b0d9e4e7358c904b0642d60402ab41c9c03d7d3c0f6fe32ec`.
- [Preflight](../../../workspaces/lanm_coupled_response_20261016/lbt_matched_model_v2/execution/PREFLIGHT.json): `02ffd7465217b3e313ceff6c9218f1deedebc0c37df45f9b2a2b8ec1b9daa023`.
- [Plan](../LBT_MATCHED_MODEL_PLAN.md): `c1ffb2b13e354f1283e70b60fc7038d68180c9f8e7bcab8903c1ad5493a2c1c8`.
- [V2 source](../LBT_MATCHED_MODEL_V2.py): `563b5e1a32fd3bf94544a42c244c54270c7d5bb497d60747db52d4415d4ff22e`.
- [V2 config](../LBT_MATCHED_MODEL_V2_CONFIG.json): `3f330329eac5cbc5e3935333a83f7b5f62dd4081d1a3e5b8c949be5007d70a34`.
- [Recovery note](../LBT_MATCHED_MODEL_RECOVERY.md): `af6c7cba99c3845f958f68b21dee4bc07236579894ac39fe8256f9157c07d66c`.

`RESULT.json` beside this report copies the saved comparison and parent pins;
`ACCOUNTING.json` preserves the observed terminal scheduler records.
