# Electronic-model sensitivity on the same geometry pools

**The large LBT source gap changes substantially with the electronic model.**
Saved-only job **12823** admitted all **36 LBT cells, 10 Mex cells and 15 menus**.
The previous MACE-based comparison replayed exactly. No molecular calculation,
geometry generation or optimization was performed.

The substituted energies are GFN2 vacuum plus fixed-FF-charge ALPB for the
whole LBT peptide, and GFN2 embedded plus the existing MM ledger and fixed-FF
ALPB for Mex. The original MACE term is replaced; the isolated-GFN2 field
correction is cancelled consistently. Geometries, composition, mechanics and
solvent models are unchanged.

## Actual result

Here, R = (Dy−La)Mex − (Dy−La)LBT. Positive R is the conditional Mex-La/LBT-Dy
exchange direction. Values are kcal/mol; these are not binding affinities.

| Source / menu | Original MACE-based R | GFN at the same MACE selections | GFN reselection in the same pool |
|---|---:|---:|---:|
| 7CCO origin | 1.042026 | 0.851761 | 2.856602 |
| 7CCN origin | 66.610710 | 2.384443 | 4.389283 |
| 7CCO candidate / union | 1.019378 | 4.043404 | 2.636567 |
| 7CCN candidate / union | 29.281571 | 3.434845 | 3.809642 |
| Global candidate / union | 1.019378 | 4.043404 | 3.809642 |

The origin source gap falls from **65.568684 to 1.532681** at fixed selections.
The candidate/union gap falls from **28.262193 to 0.608559** at fixed selections;
GFN reselection gives **1.173075**. Thus, most of this gap change is already
present before changing the selected geometries.

All 15 GFN reselected exchanges are positive and outside the predeclared
±1 kcal/mol interpretation interval. The fixed-selection 7CCO/global origin
value, 0.851761, remains inside it. Reselection changes both models' selected
geometries: for the global union, GFN chooses LBT `geometry_17` for both metals
and Mex `Dy_terminal` for both. All original outcomes, component balances,
selections and source gaps remain in [RESULT.json](RESULT.json); the full
[15-menu table](../../../workspaces/lanm_coupled_response_20261016/electronic_substitution_v2/execution/REPORT.md)
retains the prior and control menus as well.

## Scientific decision and limits

Prioritize a small, predefined independent electronic reference for the same
two-source La/Dy contrast before expanding the optimizer or promoting GFN.
The natural four-cell question is both original LBT sources × both metals,
with a matched chemical domain and explicitly declared electronic states;
feasibility and settings must be frozen separately before execution. The
decision should be whether an independent treatment supports the source
contrast within the declared resolution, regardless of the exchange sign.

This result establishes model sensitivity, not that MACE is wrong or GFN is
accurate. GFN uses an effective singlet, f-in-core treatment; MACE's physical
La multiplicity 1 / Dy multiplicity 6 metadata does not make those electronic
treatments equivalent. Spin-orbit and physical Dy-state accuracy remain open.
These are finite, previously generated pools. LBT has two differently refined
NMR sources; Mex still has one source and its fixed occupancy/context. The
experimental assays are not matched, and no affinity, population or robust
protein classifier is established. New matched-objective response geometries
are outside this archived 46-cell comparison.

## Execution and preserved failure

Job 12823 completed in **8 scheduler seconds on 1 CPU**, with 6.262 seconds
of measured CPU time; application time was 7.714 seconds. **Zero new molecular
calls.** [ACCOUNTING.json](ACCOUNTING.json) pins the complete result and receipts.
The first attempt, job 12770, failed its energy-only receipt schema check;
its unchanged output and 10-second cost are retained in
[ATTEMPT_V1.md](ATTEMPT_V1.md). V2 recovered the actual archived charge evidence
at the unchanged 1e-5 e limit, with no molecular reruns.
