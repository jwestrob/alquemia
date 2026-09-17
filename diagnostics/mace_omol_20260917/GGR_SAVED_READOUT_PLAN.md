# GGR saved-readout diagnosis

Declared after the robustness results were inspected: 1GLG 23.98742146,
2FW0 44.88412736, 2FVY 45.97550041 model kcal. The predeclared robustness
test fails. This diagnosis cannot rescue or replace that result.

Use the actual bound La/Ca node and embedding arrays for all three structures,
the fixed disconnected-metal terms and their exact physical atom mappings.
Decompose the existing score by source residue and distance shells
metal, (0,6], (6,12], (12,18], >18 Angstrom. Require closure within the existing
0.01 model-kcal numerical tolerance. Enumerate nearby O/N/S atoms within 3.5 A
and compare all three GGR structures in their original fixed order. Record
full charge, terminal and water differences rather than attributing everything
to a single changed bond. These are learned bookkeeping terms, not observable
atomic energies or proof of a causal physical mechanism.

No new model calls, altered inputs, biological labels, calibration, fitting,
optimization or DFT. This can identify where the learned readout changes and
inform a separately declared response experiment; it cannot demonstrate that
an untested relaxation or environmental correction will improve discrimination.
