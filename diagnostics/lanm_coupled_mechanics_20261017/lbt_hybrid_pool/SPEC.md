# Fixed complete LBT pool under both hybrid partitions

This finite panel is declared before its hybrid pool energies are examined.
Root owns all execution. No geometry, optimizer, electronic method, parameter,
state, solvent charge or MM ledger changes in this stage.

## Gates and complete denominator

Require all of the following actual receipts before any new molecular call:

1. The unchanged local hybrid evaluator's nine-cell 7CCO/baseline/La scout
   passes both directions and both finite-difference steps.
2. Its eight original cells (two sources × two partitions × two metals) are
   complete and replay from raw component evidence.
3. The remaining three 7CCO cases pass all 24 displaced checks under
   `QUALIFY_EXPANDED_SPEC.md`, using the same frozen gates.
4. The Module2 LBT pool is complete: all four broad Cartesian terminal slots
   are admitted, all original 22 candidates / 44 states replay, and every
   required paired state in the expanded Module2 pool is available.

Use **every** candidate in that actual Module2 pool, in its recorded order.
Its exact-coordinate, source-local deduplication is authoritative. Let
`N = N_7CCO + N_7CCN`, where `22 <= N <= 26`; each source has its inherited
11 candidates and at most two additional broad terminal geometries. Require
exact original source membership and retain all aliases. Do not select only
the old energy minima, rank candidates by score, add search trial geometries,
or combine the two source menus.

The required hybrid denominator is **4N cells**: both metals and both regions
at every candidate. All eight original cells are reused. The new-composite
ceiling is `4N - 8 <= 96`, reduced by any other exact compatible hybrid cache
hits. Freeze the actual candidate pins, all required slots and missing-call
count in `PREFLIGHT.json` before any model initialization. Use fixed order:
source 7CCO then 7CCN; region baseline then expanded; metal La then Dy; then
the parent candidate order within that source.

## Unchanged model and exact reuse

Use the unchanged `lbt_hybrid/EVALUATOR.py` SHA256
`a482559ca4febd240a40ec7d6042f4fa24ec4358e4f5ba05770bcfaacea4b784`.
The Hamiltonian remains local MACE vacuum plus embedded-minus-isolated GFN2,
the full MM ledger, and whole-system fixed FF ALPB. Baseline Q172/PC84 and
expanded Q210/PC46 retain their prepared maps, caps and charge shifts. All 256
physical coordinates come directly from the admitted parent pool. Physical
La1/Dy6, effective GFN2 UHF0 and all charge/state limitations are unchanged.

Reuse the existing evaluator's full-cell and global-solvent cache keys. Replay
the admitted scout, eight origins and 24 qualification cells. A hybrid cell
can be reused only for exact source/region/metal/coordinates and model identity.
Do not reuse whole-peptide MACE as local MACE or transfer local components
between different partitions.

Every Module2 paired cell already has global fixed FF ALPB evidence. Import
that native energy and full native gradient into the existing solvent cache
only after raw artifacts, binary identity, return code, actual XYZ/symbols,
exact unshifted charge vector, charge/rounding gates and saved component energy
are verified. Require the same vector in both regions. This panel makes
**zero new ALPB calls**. An absent or incompatible receipt makes preparation
unavailable; it does not authorize a replacement calculation. Retain this
distinctly from a later electronic/MM backend failure.

Maximum new work is 96 MACE forwards, 192 standalone GFN2 calls, 96 MM
configurations / 384 MM group queries, and zero ALPB calls. Actual ceilings
are recomputed from the fixed `4N - compatible_cache_hits` denominator. One
warm model is reused; component execution/threading remains the frozen
evaluator's 12-CPU / one-GPU policy. No new caching engine, worker scheduler,
retry, continuation or partial-run recovery is implemented.

## Analysis declared before evaluation

For candidate c, source s, metal m and region r, report stable native-component
work `W_r,m(c) = E_r,m(c) - E_r,m(original_s)`. Select each metal's mathematical
minimum work from the **complete** source-local menu; preserve exact ties in
parent order. Report individual La/Dy works, their common mean and Dy-minus-La
accommodation, component contributions, selected IDs and all per-geometry
works. Also report the inherited 22-candidate submenus explicitly, with their
fixed membership. They are diagnostic submenus, not fallback success.

Report same-metal source strains separately for original and selected
geometries: `E_r,m(7CCN) - E_r,m(7CCO)`, followed by their Dy-minus-La
difference and absolute magnitude. These comparisons stay within one region.
Keep individual-source outcomes and the parent Module2 comparison visible.

At every identical physical geometry, compare expanded-versus-baseline
**accommodation work**, component by component, and compare each hybrid total
work with the parent whole-peptide model's work at that same geometry. Also
retain the hybrid works at the parent model's selected geometries. Never use
`E_expanded - E_baseline` as an absolute energy comparison or mix their minima
in one pool. No threshold fitting, biological labels, affinity, populations or
convergence claim follows from this fixed finite panel. The qualification is
directional at 7CCO, not a proof of electronic accuracy or all-gradient accuracy
at every candidate.

## Failures and output

All required cells are declared before evaluation. A component/backend failure
stops new calls; remaining slots stay `unavailable_not_attempted`. No failed
exact cell is retried. A source/region menu needs every paired cell in its
declared membership. Overall completion requires both regions and both sources;
an unavailable expanded result must never appear as successful baseline
fallback. Retain successful partial cells, missing-cell reasons, raw request
pins, exact parent costs and new attempted costs separately.

Run a real saved-cell missing-data check on the same menu routine: remove one
actual original La cell from a test-only copy and require an unavailable menu.
This makes no molecule or alternate production result.

CLI: `--stage preflight|run --config CONFIG.json --output FRESH_DIRECTORY`.
Both stages require Slurm. Preflight makes zero new molecular calls and does
not initialize MACE. Config contains `protocol`, `implementation`,
`specification`, `plan`, `qualification_implementation`, `qualification_result`,
and `module2_pool_result` pins. The qualification and Module2 receipts bind
their actual configurations and frozen helpers. Root binds only completed
parents, reviews preflight, and launches the declared panel.
