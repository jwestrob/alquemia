# Independent static review: fixed LBT scaffold transition

Reviewed 2026-10-19T07:42:10Z by `/root/mex_health_review`.

**PASS for the declared bounded Slurm diagnostic. No outstanding static blocker.**
This review inspected source, archived schemas, text and file hashes only. It
performed no scientific imports, numerical execution, model calls or submission.
Actual numerical closure and the complete four-cell result remain untested.

## Frozen bundle

| File | SHA256 |
|---|---|
| ANALYZE.py | `79d311d8ea1fe24dd7582efd568b9ab107877e04677189a6f3578d1c3b5f2a13` |
| SPEC.md | `5bbfeef09afc0763c5e7884e16b81580ef7bee61d0f007161b9ff1e44c75e7f0` |
| CONFIG.json | `8406085218c8f2b9136341ca05fa9fe28c8f7d840302ca7c1b1610b046e686bd` |
| RUN.sbatch | `81c78eb865ab370024e9db23b8348e935096cec82fec9c6213e2916db44dfdcd` |
| ROOT_SCOPE_v2.json | `e7493639861f02a5309dda574f7780147d81f2fb8ff209b780eeb06e79ba6a40` |
| READY.json at review | `862c259455b64da0879e073a78ab3e211a3b91b0447e7aa49dcf5d8f5ef71ee5` |

All path/SHA records embedded in CONFIG passed byte-hash verification. The pinned
old AUDIT.py (`a7154cc7…`) and PREPARE_V2.py (`0d02c79f…`) interfaces match their
uses here; neither helper's main program is invoked.

## Findings

- **Exact four-cell joins:** physical parents, positions, ordered IDs, charge,
  physical/effective state distinctions, native task signatures and individual
  input/XYZ/receipt/output/gradient pins are checked. Before-Dy correctly uses
  point023 at the shared La-before coordinates; the own-Dy point024 is not reused.
  Inline reference task metadata is handled as an object, not as a file pin.
- **Force signs and availability:** saved gradients are negated to forces; the
  reference composite is `-(PBE0 gradient + fixed-FF ALPB gradient)`. POINT
  component closure and archived comparator forces are checked exactly. The
  actual before-Dy energy-only force status is recognized. After-Dy solvent uses
  its completed row/raw receipt and full-reference admission, rather than the
  earlier preparation's unavailable placeholder. Before-Dy whole-MACE forces and
  after-Dy full-hybrid forces/work remain unavailable.
- **Complete terminal inventory:** the unchanged class counts are 34 bonds,
  61 angles, 90 periodic torsions and 2 CMAP terms: 187 unique interactions
  touching any of the 33 Phe1/Ala17 atoms. Every selected native/retained term and
  ledger identity is checked. Shared terms enter the union once. The corrected
  aggregate supersedes the preserved draft transcription of 197.
- **Splitting and closure:** CMAP subsets preserve all 16 maps, original map
  indices, atoms and periodic flags. Other subsets use the pinned parameter-copy
  helper. Each batch retains every original interaction exactly once, including
  both group3 containers. Full256 group/class/batch energy-force closures and
  term-force sums at every terminal atom use the original 1e-8 kcal/mol energy
  and 1e-5 kcal/mol/Angstrom Cartesian-force thresholds. The helper supports the
  33-atom terminal arrays. After-Dy requires admitted after-La control before
  querying; common groups0–2 use the pinned comparison-only inactive-sigma proof.
- **Finite execution:** original groups/total require 5 queries, classes/total
  5, and seven term batches 194: 204 queries and 9 contexts per cell, capped at
  816 queries and 36 contexts overall. All four input bindings precede queries.
  A failed cell remains unavailable; remaining declared cells may be attempted
  once, with no retry and no complete transition interpretation unless all four
  pass. The launcher requests one CPU on standard/banfield, disables requeue,
  fixes numerical threads to one and leaves scheduler memory/time defaults.
  Source requires Slurm, the pinned interpreter and OpenMM8.5.1 Reference.

## Interpretation limits retained

The result separates endpoint force errors, their changes, and same-model
energy differences along the exact saved deformation. Both frozen directions,
full256 arrays, named terminal vectors and component contributions are retained.
The descriptive interaction ranking adds no queries. Interaction energies are
whole-interaction contributions; overlapping terminal memberships are explicit.

This cannot uniquely infer stiffness/equilibrium parameters, quantum bond
energies, physical accuracy, preference or a parameter repair. New after-Dy MM
does not fill the missing full hybrid. Native SCF/state evidence and conditional
termini, vacuum-dianion, frozen-f, nonstationary-geometry and unqualified
whole-reference-force limitations remain. Prior failed raw reconstructions and
the closed Mex derivative candidate are unchanged.

Root may proceed with the frozen bounded run; interpretation requires the
actual Slurm result, all original closure gates, availability and measured cost.
