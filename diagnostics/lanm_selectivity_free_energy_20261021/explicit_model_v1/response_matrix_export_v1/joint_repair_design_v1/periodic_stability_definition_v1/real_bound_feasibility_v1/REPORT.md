# Fixed aqueous real-space bound: shortcut insufficient

The two-box diagnostic completed under the exact declared0.4/A Ewald split and
9A real cutoff, using baseline H=0 and unchanged geometries. Both fixed block
partitions give negative lower bounds in both La and Dy boxes:

| Partition | Lower bound | Upper bound | Block attaining lower bound |
|---|---:|---:|---|
| Single atom | −1.73567019 | 3.35051677 | Cl atom11956 |
| Intact water / singleton ion | −1.13055495 | 3.16774658 | Water atoms6578–6580 |

These are ordinary-double exploratory bounds for the real+self part. They
cannot certify the full periodic operator, and they do not show that either
operator is indefinite. Do not use them as failed physical-model evidence.

The limiting chloride is35.7997A from the metal; the water's three sites are
20.8601–21.7773A away. All lie beyond the9A real cutoff. With the fixed intrinsic
alpha normalization, a correction only to the metal's diagonal stiffness leaves
these nonmetal bound rows unchanged. Retuning that correction therefore cannot
make this particular global certificate positive. Close the atom/water block
Gershgorin shortcut; do not scan metal coefficients or Ewald settings to rescue it.
Another stability method requires a separately defined implementation.

Each box contained1,815,015 real atom pairs and262,878 cross-molecular blocks.
Actual execution cost93.844wall seconds/93.768childCPU seconds on one authorized
login thread, zero scheduler allocation and zero native/reference/fit calls.
The independent saved review passed50 provenance/arithmetic/geometry checks
without rerunning the bound. Source and synthetic controls are retained.

Artifacts in the matching workspace: actual_v1/RESULT.json(7033c722),
PROCESS_TERMINAL.json, ROOT_RELEASE.json(b55ae176), and
independent_actual_review_v1/REVIEW.json(82c54077). The current component-fit
campaign and production implementation are unchanged.
