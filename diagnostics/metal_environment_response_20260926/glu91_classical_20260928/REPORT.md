# Glu91 classical response is common to both metals

The four prescribed repaired Hans8DQ2 Glu91 ±1° configurations completed. Both La and Dy classical energies increase by **0.033739411802 kcal/mol from A(−1°) to B(+1°)**. The classical La-response-minus-Dy-response difference is exactly **0** in the saved calculation. At each geometry, the classical Dy-minus-La contrast remains **+0.058069241205 kcal/mol**. Root can add the matched classical arrays/energies to the electronic endpoints while retaining the electronic-only result separately.

The changing component is donor–MM LJ. Target-metal position and its outer contacts are unchanged; donor cross LJ uses the same native donor parameters for both metal substitutions. MM–MM interactions and retained classical bonded terms do not change. Purely QM Glu91 bonded terms remain excluded by the existing ledger, for the electronic calculation to represent. This cancellation is a property of this finite model and motion, not a biological conclusion or evidence that the exterior has no influence.

## Exact reuse and output

The original repaired `scaffold_origin_DQ2_v1/Hans8DQ2/LEDGER.json` and all pinned parameters, exceptions and retained terms are reused without alteration. Coordinates come directly from each `glu91_motion_v1/physical_source_{A,B}.json`; only the declared CD/HG2/HG3/OE1/OE2 move. Core XYZ, physical caps, source identities, charge/state metadata, field and actual La/La/Na spectators pass direct comparisons. No coordinate search or new parameterization occurred.

Full results:
`workspaces/metal_environment_response_20260926/glu91_classical_20260928/RESULT.json`
plus `{La,Dy}_{A,B}.json` with1891 physical source IDs, exact coordinates, all component energies, full component and total Cartesian gradients (negative native forces; kcal/mol/angstrom), input/state and ledger pins. Physical target ID is `A/203//LA` for both metals. Caps are not independent physical particles.

Classical expression remains retained native bonded + MM LJ + shifted-charge MM Coulomb + QM–MM LJ. No QM–MM classical Coulomb, C4, cap LJ or implicit solvent. The full-model qualifier remains false; missing bulk solvation, unqualified cross parameters and previous numerical limitations are unchanged. No origins or electronic endpoints were repeated.

## Tests and actual cost

Four actual-input/result checks pass: source geometry, cap mapping, unchanged charge/proton/occupancy state, matched paired coordinates, exact common-component energies/gradients, finite1891×3 arrays and their component sum. These are implementation tests, not a new quantum or full-energy derivative test.

Job1220361 COMPLETED: **3CPU-s** (3s×1CPU), shared allocation, mem=0. Sixteen component queries, zeroQM/GPU/optimization. Script wall2.146840s, process peakRSS166844KiB. Scheduler MaxRSS0 is not taken as a memory measurement.

Recheck without molecular evaluation:
```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/metal_environment_response_20260926/glu91_classical_20260928/TESTS.py
```
The batch file records completed work and must not overwrite this output. Root owns electronic collection and combination. No active job remains in this subtask.
