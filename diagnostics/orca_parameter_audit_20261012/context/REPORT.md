# Context versus accommodation: exact finite factorial prepared

12 October 2026. **The missing matched control is now geometrically executable.**
Both boundary proteins support transferring their actual saved donor motions
between original and expanded contexts without inventing coordinates for added
atoms. No MACE, GFN2, DFT, optimization or scoring call ran in this branch.

The independent corrected audit already completed 2,796/2,796 original points.
Its original accommodated reference separation fails; this preparation does not
repeat or reverse that result. Its separate expanded-context positive gap changes
both composition and the searched donor motions. The finite experiment here asks
which of those changes contributes, without fitting a replacement threshold.

## Exact sources and preparation

| Protein | Pinned parent fold | Original atoms/caps | Expanded atoms/caps | Shared entries including metal |
|---|---|---:|---:|---:|
| q9z4j7-pqq-la_model | La-conditioned, seed-1/sample-0 | 145 / 9 | 192 / 14 | 142 |
| c5axv8-pqq-la_model | La-conditioned, seed-1/sample-1 | 172 / 10 | 192 / 14 | 171 |

These are the two already-consumed canonical boundary cases, not new biological
validation. The exact parent protonated PDB hash is identical between each
protein's contexts. Both retain Ca charge −3, La charge −2 and multiplicity 1.
All shared origin coordinates, including hydrogens and the retained PQQ atoms,
are exactly equal. Every original physical atom exists in the expanded context;
the only disappearing entries are synthetic caps: three for q9 and one for c5.
Their retained/omitted source bonds remain explicitly mapped.

Additional q9 physical atoms belong to A139/A140, A211/A212, A282, A351/A352 and
A442/A443; c5 additions belong to A130/A131, A316/A317 and A409/A410. The exact
atom and formal-charge records are in PREPARATION.json. This changes fragment
chemistry, caps and solvent cavity together; a matched-coordinate context effect
must not be called exclusively an electrostatic effect.

Original candidate pools were recovered from
`workspaces/slsqp_precision_20260923/pool_v1/manifest.json` (q9) and
`workspaces/slsqp_precision_expansion_20260923/pool_v1/manifest.json` (c5).
Expanded pools are in
`workspaces/motion_envelope_20260923/pilot34_pool_v1/manifest.json`.
The preparation artifact pins exact original/expanded XYZ files, proposal
receipts, protonated parents, cap maps, physical kinematic maps, selector results
and implementation files. Original cap construction comes through
`consistent_context_20260922/prepared_v1`; expanded maps come through
`motion_envelope_20260923/run_v1/prepared`.

## An additional confound is now explicit

The force selector changes the fourth active torsion when context expands:

| Protein | Three retained modes | Original fourth | Expanded fourth |
|---|---|---|---|
| q9 | A213 chi1/chi2/chi3 | A350 chi1 | A350 chi2 |
| c5 | A319 chi1; A275 chi1/chi2 | A317 chi2 | A198 chi1 |

Thus the original comparison varied the active search subspace as well as the
energy context. A new search is unnecessary to test the saved hypotheses:
cross-score both contexts' actual proposals under both energy contexts first.

## Prepared factorial and checks

For each source, both metals see the same five physical arrangements:
origin, original-context Ca proposal, original-context La proposal,
expanded-context Ca proposal, expanded-context La proposal. Each arrangement is
represented in both contexts. This supplies 20 context/geometry rows, 40 metal
energies, and 80 logical GFN2 vacuum/water cells before reuse.

`PREPARE_FACTORIAL.py` replays each saved full torsion vector by physical mode ID
through the target context's existing `mace_site_kinematics.Kinematics`. It does
not graft disconnected Cartesian fragments or optimize new angles. Replaying
the donor's own map first verifies the archived candidate. Then the target map
retains its own cap construction and full source bond graph.

The actual preparation completed in about 9.7 seconds using one numerical thread,
within Jacob's newly authorized login allowance and this agent's two-thread
allocation. All 16 nonorigin replay rows passed:

- Shared target/source candidate coordinates match exactly (maximum error zero).
- Every context-exclusive atom and cap remains fixed.
- The worst retained physical bond error is 3.78×10⁻¹⁵ Å.
- Existing analytic-Jacobian, fixed-atom and 0.8 Å displacement gates pass.
- No new severe overlaps occur under the existing H-containing 0.55 Å / other
  1.0 Å admission thresholds. These coarse checks are not steric-energy validation.
- Ca/La endpoint coordinates are identical, with only the metal symbol and
  corresponding total charge changed.

These are real-fixture mapping tests, not successful molecular integrations.
No unavailable energy or correction is set to zero.

Actual generated coordinates and all checks are under:
`workspaces/orca_parameter_audit_20261012/context_factorial_preparation_v1/`.
`FACTORIAL_PLAN.json` pins its `PREPARATION.json`, contains endpoint paths for
every row, and provides an explicit fresh-output preparation command.

## Finite next work and interpretation

Twelve context/geometry rows already have their own-context archive calculations.
After binding and validating those exact corrected raw receipts, only eight
cross-context geometry rows need new scoring: **16 MACE evaluations and 32 GFN2
single points**, plus up to 32 independent tblite checks under the integrator's
qualified per-input policy. No optimization, DFT, fold generation, full-panel
rescore or calibration campaign is included.

Keep `R(c,g) = E_Ca(c,g) − E_La(c,g)` and retain all MACE, water, vacuum and
solvent-transfer components. Report context change at matched geometry,
`R(expanded,g) − R(original,g)`, and the interaction with motion:

`[R(expanded,g) − R(expanded,origin)] − [R(original,g) − R(original,origin)]`.

Report all original and expanded origin/adaptive_Ca/adaptive_La ablations under
their named source protocols. Additionally report mathematical five-pool minima
and the explicitly retained operational rule (origin unless improvement is
strictly greater than 0.1 kcal/mol; recorded tie order). A missing required cell
blocks its pooled score. Candidate counts are not populations. No new class
decision is supported by these two inspected references alone.

The frozen geometric comparator and its source poses stay intact. Existing
corrected-original failure and corrected origin-only outcomes remain separate.
Any later reaccommodation experiment must freeze a common search policy first;
it is not silently appended to this finite cross-score matrix.

## Execution gates at preparation time

The integrator owns backend qualification and new endpoint execution. The
selected vault snapshot contains reports and selected receipts, but not all
corrected per-candidate raw endpoint archives; exact reuse still needs binding
to verified cluster equivalents or retrieval of the unchanged saved artifacts.
Do not populate a cell from summary prose or automatically repeat 48 old GFN2
cells because their raw files were not copied. MACE model/head/dtype and corrected
GFN2 state/numerical pins must match before reuse or execution.

The cited `independent_review/strict208_corrected_edge_audit.json` is absent from
the 129-file synced selection. The mapping results above instead come from the
actual local manifests, source maps and saved geometries. Historical reports'
Mac paths remain provenance, not cluster executable paths.

Before launching, read current `~/.claude/CLAUDE.md` and synced cluster-operations
instructions in that same turn. This branch used the explicit new login
allowance; it does not authorize background jobs outside the current shared
thread allowance or prescribe stale Slurm resource layouts.

## Subsequent finite execution — 12 October2026

After the missing corrected per-candidate archives were confirmed, root authorized
the separately recorded80-input Linux runtime-transfer/factorial matrix. It is
now complete:80 inputs in each of two independent backends, all40 actual
same-geometry solvent components qualified under the explicitly versioned v2
gate. The original absolute-energy parity failure and its causal dispersion
diagnostic remain preserved. See [the solvent report](../context_GFN2/REPORT.md).
The prepared geometry set above was unchanged; no unavailable old result was
filled from prose. Root separately owns MACE evaluation and combined selection.
