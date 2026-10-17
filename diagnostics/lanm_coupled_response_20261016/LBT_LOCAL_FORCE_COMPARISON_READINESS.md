# Local target-force comparison: saved data unavailable

**None of the four exact matched-core MACE cells from job 12471 retains a
force vector.** The declared inspection covered both original sources and
both metals, with target index 0 and fixed Cartesian components. No case,
direction or tolerance was selected from DFT outcomes.

The actual pinned MACE adapter used `torch.no_grad()` and
`compute_force=False`, and explicitly required `v['forces'] is None`.
Each separate 48-atom cell receipt exactly matches its combined-result row;
all contain completed energies and native-state metadata, with no saved
forces. This is intentional energy-only execution, not lost force output.

| Declared cell | Target index | Saved MACE target force |
|---|---:|---|
| 7CCO / La | 0 | Unavailable |
| 7CCO / Dy | 0 | Unavailable |
| 7CCN / La | 0 | Unavailable |
| 7CCN / Dy | 0 | Unavailable |

Existing whole-peptide forces cannot substitute for these charge −1,
48-atom capped-core states. No force comparator was implemented, no unfinished
DFT output was read, and no recalculation or job was requested. The inspection
stops here under the declared missing-data rule. Native DFT force consistency
remains unqualified; no validation conclusion follows.

[LBT_LOCAL_FORCE_COMPARISON_READINESS.json](LBT_LOCAL_FORCE_COMPARISON_READINESS.json)
records the exact result, config, adapter, four cell and XYZ pins plus their
saved field inventories. Zero new molecular calls or numerical analysis.
