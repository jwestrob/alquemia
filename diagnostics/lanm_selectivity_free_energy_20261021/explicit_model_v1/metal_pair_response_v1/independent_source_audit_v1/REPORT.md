# Prototype source audit

**No source blocker remains in source_v3.** The null-value parser issue found
in the first draft is fixed: both type integers start at zero and both deltas
start at IEEE quiet NaN before each list-directed read. Missing/null values
cannot reuse a previous record. Earlier source versions remain preserved.

The actual patch matches the independently reviewed derivative:

- `kpolar` resolves the existing POLPAIR table before applying the two role
  multipliers. Only unordered type pairs 500/36 and 501/36 are admitted;
  duplicates, nonfinite values, deltas ≤−1 and invalid resulting coefficients
  fail. Both matrix directions receive identical coefficients and selection.
- `damptholed` selects D; `dampthole` selects M, including the existing mutual
  preconditioner calls. No global alternate-Thole flag is enabled.
- `epolar1a` uses D factors and matching radial derivatives before direct
  contractions and permanent-frame torque accumulation. It independently
  recomputes M before the cross-induced `uscale` contractions, including the
  direct-only case. It retains the exact pair coefficient rather than the
  stock alternate-direct `min(thole_i,thole_k)` fallback.
- Absent/all-zero declarations leave `npairdm=0`; unmatched pairs return
  `selected=false` without changing the damping output array. The original
  factor and derivative branches then execute. Changed water dipoles through
  the coupled solve are not mistaken for edits to water–water kernels.
- The nonzero-feature guard runs at `induce` entry, after mechanic setup and
  before PB/GK/TCG dispatch, field construction, mutual solves or preconditioning.
  It rejects the declared unsupported modes and requires positive finite site
  polarizabilities. Energy/analysis callers reach this guard before their pair
  loops. Absent/all-zero declarations bypass the new mode guard.
- Only the declared five source files change. The original energy formulas,
  native solve, exclusion scales, force/virial accumulation and torque
  redistribution remain intact. The modified double-loop derivative path is
  protected from the unmodified list/periodic/solvent alternatives.

This is a source audit. Root owns the isolated compile. No compiler, linked
kernel probe, molecular calculation or numerical derivative test was run by
this reviewer. Executable behavior, induction stability, energy-route equality
and the declared finite force tests remain unqualified until actual execution.
The dual-field energy is stationary; positive coefficients do not by themselves
prove a positive-definite response matrix or physical accuracy.

The corresponding workspace `REVIEW.json` pins the actual source_v3 manifest,
patch, reviewed five-file pairs, unchanged callers, implementation map and
independent derivation. It records the bounded saved pin/patch checks.
