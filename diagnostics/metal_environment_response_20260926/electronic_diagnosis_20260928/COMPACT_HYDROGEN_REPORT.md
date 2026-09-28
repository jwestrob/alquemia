# Real hydrogen strain also affects the compact LanM cores

**Yes: the compact comparisons retain the severe real-source hydrogen-angle defects. They also retain the older overlong H bonds.** This is relevant to interpreting source-sensitive energies and metal loads. It does not invalidate their exact paired coordinates, demonstrate which metal is preferred, or uniquely explain the SCF failures.

Read-only audit of the actual compact manifests/XYZ, amide-v3 source maps, hashed original protonated PDBs and ff19SB XML hydrogen-parent bonds. No molecular calculation, geometry repair or new label was introduced. All selected source coordinates exactly match their archived QM maps. Protein real hydrogens and synthetic boundary caps are identified separately; retained crystal-water hydrogens are not included in the protein-H table below.

| Actual compact source | Real-source angle | Degrees | Source heavy neighbor retained? |
|---|---|---:|---|
| Hans8DQ2 | CA83–CB83–HB2 |166.711397|no; replaced by collinear cap|
| Hans8DQ2 | CG83–CB83–HB3 |161.333748|yes|
| Hans8DQ2 | CG94–CB94–HB3 |148.640588|yes|
| Hans8FNR | CA91–CB91–HB3 |164.927825|no; replaced by collinear cap|
| Hans8FNR | CG91–CB91–HB2 |153.328931|yes|
| Mex8FNS | CB95–CG95–HG3 |130.337161|yes|

The corresponding compact cap–CB–real-H angles are166.711396913° at Hans8DQ2 CB83 and164.927825431° at Hans8FNR CB91. Across all audited cap substitutions, the angle agrees with the original missing-heavy direction within3.81e−9degrees. Thus these are inherited real-H orientation problems, not synthetic caps being optimized or placed in a new direction. Retained-heavy-neighbor examples independently establish that the issue is not restricted to boundary angles.

| Source | Audited real protein H count | Parent–H distance range, Å |
|---|---:|---:|
| Hans8DQ2 |17|1.178017–1.203213|
| Hans8FNR |17|1.170975–1.191392|
| Mex8FNS |11|1.174053–1.201502|

These are the exact older compact sources used by `lady_compact_exchange_v1`, not the subsequently normalized whole-chain/new embedding inputs. Radial length normalization in those newer inputs does not change hydrogen angles and therefore cannot repair the orientation defect. Neither version should be described as an angularly relaxed hydrogen network.

## Interpretation for the running comparisons

Each source's La/Dy pair has identical positions and inventory; no unequal H preparation was introduced between the metals. But equal strain in the coordinates does **not** guarantee an equal electronic penalty: metal-dependent polarization/coordination can couple to it. Likewise, the different hydrogen distortions between source structures can contribute to apparent source dependence. We have identified a concrete confound, not quantified its energy contribution or proved it is the dominant cause.

Large metal or ligand gradients can now reflect both actual coordination mismatch and preparation strain. The predeclared force diagnostic should remain unmodified and report the original results. A separately versioned, uniformly applied H-angle preparation comparison would answer a different and useful question; it must preserve heavy geometry/occupancy and compare both metals on the same corrected geometry. No repair or such calculation is performed here, and no previously inspected outcome is treated as a fresh test.

Reproducer: `audit_compact_hydrogens.py` in this directory. Exact source/manifest/topology/XYZ pins, all tested real-H bond lengths, heavy-neighbor angles and cap correspondences are in `COMPACT_HYDROGEN_AUDIT.json`. The script exclusive-creates that file; use an explicit new output path in a copied/versioned analysis if rerunning, rather than overwriting evidence.
