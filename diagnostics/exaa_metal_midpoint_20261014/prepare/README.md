# One fixed midpoint, preparation only

The frozen manifest prepares eight cells: samples0/3 × original saved
adaptive_Ca/adaptive_La arrangements × Ca/La. It reuses the actual admitted full
transplant metal point and computes `original + 0.5*(full_transplant-original)`
once per host. No alignment is recomputed. Both donor arrangements and both
metals receive the identical midpoint.

The original source maps, state, cap mapping and full donor q vectors are pinned.
The unchanged frozen `adaptive_metal_proposals.final_geometry` validates each
saved and midpoint geometry against the original host baseline, including the
0.8 Å displacement limit and source/cap/covalent/overlap guards. No threshold is
relaxed. Only the metal XYZ line changes; every nonmetal line is preserved
exactly. Midpoint serialization uses17 significant digits and requires exact
floating-point roundtrip. Kinematic coordinate replay retains the existing
1e−8 Å tolerance.

`PREPARATION.json` contains both rows and all eight cells, even on failure. A
failed required geometry invalidates all four cells of its host; the root
scorer requires both rows admitted. Each row/cell pins its geometry admission,
original XYZ/MACE/solvent result and donor-candidate identity. Row metadata also
pins the original and full-transplant origin XYZ files and reports all three
metal points. New scoring labels are owned by the scorer, while the preparation
matrix retains `adaptive_Ca` / `adaptive_La` candidate keys.

[READY.json](READY.json) records the exact one-CPU Slurm command and manifest
hash. Metadata freezing and syntax inspection are complete; scientific
preparation is unexecuted until root submits it. No model/solvent/optimizer calls
belong to this script. Root owns submission, wake monitoring and interpretation.
