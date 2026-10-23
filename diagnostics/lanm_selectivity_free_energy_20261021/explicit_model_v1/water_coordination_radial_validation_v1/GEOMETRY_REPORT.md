# Fresh W9 radial geometries frozen

Four geometry-only cells are prepared: `aquo9_minus_La`, `aquo9_minus_Dy`,
`aquo9_plus_La`, and `aquo9_plus_Dy`. Each contains28 atoms: one metal and nine
intact waters. Every water was translated rigidly by−0.10 or+0.10Å along its own
original metal-to-oxygen unit vector. The metal remains fixed. La and Dy share
identical numerical coordinates and water graphs for each displacement.

The original typed central XYZ files match the exact inputs admitted for17670
byte for byte. Both ordinary and typed XYZ outputs use the same coordinates,
rounded once to12 decimal places. The largest radial-step discrepancy is
6.08×10⁻¹³Å; the largest translation-component discrepancy is4.57×10⁻¹³Å.
All27 intrawater pair distances are preserved within3.34×10⁻¹⁶Å. All graph,
state and arithmetic checks pass under the declared2×10⁻¹²Å numerical tolerance.
Charge remains+3; physical multiplicities are La1 and Dy6. Metal/water types and
bonds remain unchanged. No model key, candidate parameters or reference method
is bound by this geometry preparation.

## Contact changes

Distances inÅ; La and Dy values are identical. These are geometric observations,
not selection filters or an energy-based acceptance test.

| Metric | Central | −0.10Å | +0.10Å |
| --- | ---: | ---: | ---: |
| Metal–O shell radius | 2.500000 | 2.400000 | 2.600000 |
| Closest metal–H | 3.177364 | 3.080335 | 3.274572 |
| Closest interwater O–O | 2.505029 | 2.404828 | 2.605231 |
| Closest atoms on different waters | 2.308322 | 2.217946 | 2.399479 |

No distinct fragments become coincident. The9-water connectivity, internal water
geometry and orientations are preserved; water–water contacts change as intended
when the shell moves radially.

## Frozen evidence and execution boundary

Manifest: workspace `preparation_v1/MANIFEST.json`, SHA`a34e9d05`.
Source: `GEOMETRY.py`, SHA`6ec4f936`. Configuration: `GEOMETRY_CONFIG.json`,
SHA`2d1d5d58`. The manifest contains all ordinary/typed XYZ hashes, original source
pins, per-water displacement vectors, bonds/types/bodies and contact metadata.
The root declaration is unchanged.

Preparation passed22 checks in0.089s wall/0.050s CPU on one thread. Zero native
calls; no fit trajectory, fitted candidate, reference energy/force or validation
outcome was read. Molecular release remains conditional on terminal17676 passing
all development controls, independent selected-ten admission and a separate
root source/config release. Until then this panel remains prepared and unbound.
