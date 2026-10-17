# Local matched-core force diagnostic

Status: complete.

Fixed Cartesian target vectors in eV/Å; signed difference is MACE minus PBE0.

| State | MACE target (x, y, z) | PBE0 target (x, y, z) | Difference (x, y, z) |
|---|---|---|---|
| 7CCO_La | -0.349341, -1.023436, 1.095643 | -0.412464, -0.962867, 1.191963 | 0.063123, -0.060569, -0.096320 |
| 7CCN_La | -1.402865, 2.966764, -0.682817 | -1.431322, 2.818149, -0.802792 | 0.028458, 0.148615, 0.119975 |
| 7CCO_Dy | -0.704868, -0.187859, 1.808232 | -0.753351, -0.165119, 1.892354 | 0.048482, -0.022740, -0.084121 |
| 7CCN_Dy | -1.442759, 1.850268, -0.678615 | -1.579696, 1.780410, -0.782653 | 0.136937, 0.069858, 0.104039 |

RMS atom-vector differences (eV/Å), with no acceptance threshold:

| State | All atoms | Target | Other real atoms | Synthetic caps |
|---|---:|---:|---:|---:|
| 7CCO_La | 0.153731 | 0.130118 | 0.164303 | 0.073341 |
| 7CCN_La | 0.182303 | 0.193107 | 0.194809 | 0.075533 |
| 7CCO_Dy | 0.159797 | 0.099720 | 0.171673 | 0.072577 |
| 7CCN_Dy | 0.162747 | 0.185623 | 0.173151 | 0.073339 |

All 48 atom vectors, maps, native gradients, group means and SCF evidence are retained in RESULT.json.

GFN force arrays are unavailable from these exact energy-only core receipts. No new calls or substitute geometries are used.

PBE0 force consistency and electronic stability remain unqualified. Its frozen-f effective singlet differs from physical Dy sextet metadata used by MACE; spin–orbit accuracy is unqualified. These raw cluster forces, including synthetic caps, are descriptive and do not establish whole-peptide physical forces, affinity or transferable preference accuracy.
