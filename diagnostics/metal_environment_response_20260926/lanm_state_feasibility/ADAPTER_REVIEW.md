# Explicit-state adapter review

Reviewed the new research-only `metal_environment_electronic_state.py`, preserving the old production parser.

Concrete fix: ORCA's contributor credits mention “stability analysis”; the initial broad evidence regex matched those credits and saved only spin-population headings, not actual values. The revised parser now records actual S2 values, spin-pure comparison and deviation, per-atom Mulliken/Loewdin charges and spin populations, integrated alpha/beta counts, and explicit missing/unqualified state. Atomic indices and elements must match the real geometry. No absent S2 or population becomes zero. Stability headings/verdicts are anchored to exclude credits; absent actual analysis remains `not_run_or_not_printed`, and even a recognized printed block is not automatically a stability pass.

Executed RHF/UHF is now checked against the requested physical state. Both ECP gradient component gates remain active. Nonfinite gradients fail. Input multiplicity must be an actual integer. This adapter still never declares electronic-state or full-hybrid qualification, and does not extract f-orbital occupancy or claim spin localization from total multiplicity.

Validation:3tests passed,1explicit Dy execution skip. Six actual Ca/La endpoints retain bit-identical energies and gradients; six actual LanM XYZ inputs retain correct state counts. Real archived ORCA6.1.1 Er output (25-atom aquo fixture) exercises open-shell format: S2=3.752366, atom0 Mulliken spin3.017912 and Loewdin spin3.016722; an explicitly corrupted atom index is rejected. This old Er calculation is a parser-format fixture only, not a new supported metal or evidence of new Dy reference convergence. Exact fixture hashes are in PINS.json. No molecular calculation ran.

Command: `python -m unittest discover -s tests -p test_metal_environment_electronic_state.py -v`.
