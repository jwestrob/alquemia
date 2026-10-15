# Actual rebased-start preparation: ten sources admitted

**Recovery4821 admits all ten sources and their twenty flipped-start search
tasks. No new energy or optimizer call occurred.** The five original failures
were caused by comparing computed cap coordinates with raw XYZ coordinates using
bitwise equality. Actual saved-map diagnosis proves that rebasing did not move
the caps. Scientific inputs, cap definitions and tolerances are unchanged.

## Full source denominator

Ca/La below denotes source conditioning, not an assigned biological preference.
Each source has both exchanged-metal tasks.

| Source | Original4814 | Shared old/new cap residual in4819 (Å) | Recovery4821 |
|---|---|---:|---|
| Ca0 | prepared | 0 | prepared |
| Ca1 | rejected: Rebased caps changed | 4.440892098500626e−16 | prepared |
| Ca2 | prepared | 0 | prepared |
| Ca3 | rejected: Rebased caps changed | 2.220446049250313e−16 | prepared |
| Ca4 | rejected: Rebased caps changed | 8.881784197001252e−16 | prepared |
| La0 | rejected: Rebased caps changed | 8.881784197001252e−16 | prepared |
| La1 | prepared | 0 | prepared |
| La2 | prepared | 0 | prepared |
| La3 | prepared | 0 | prepared |
| La4 | rejected: Rebased caps changed | 1.1102230246251565e−16 | prepared |

Original4814 remains immutable and partial5/10. Diagnostic4819 compared all14
caps in all ten sources: raw XYZ caps, raw map caps, links and physical anchors
are exactly unchanged; old/new computed q0 caps are also bitwise identical.
The rejected predicate reproduces the original five failures. The nonzero
residuals above were already present in the old map's cap algebra.

## Narrow recovery and actual gates

The independently reviewed v2 compares raw-before/raw-after caps exactly, old/new
raw map caps exactly, and old/new computed caps exactly. Both original and rebased
full-coordinate q0 replay checks retain the existing1e−12Å tolerance. There is no
coordinate rounding, cap repair, new tolerance, mode change or optimizer change.
[RECOVERY_V2.md](RECOVERY_V2.md) records the frozen technical rationale.

Actual4821 passes all three exact cap checks for all ten sources. It also passes
paired q0 replay, cached flip-force projection through the rebased Jacobian,
covalent/source/contact admission, fixed atoms/cap anchors and retained physical
hydrogens. Recorded Jacobian residuals are about1e−9Å/unit, within the existing
1e−7 gate; bond residuals are zero or2.220446049250313e−16Å. The cap residuals
match4819 exactly. The full denominator retains twenty old searches and eighty
prior endpoint cells as provenance/reuse inputs.

These are actual preparation gates, not evidence of a relaxed basin or successful
new optimizer. The first declared Ca-conditioned0 pair is now available for the
separate reviewed warm scout. Remaining searches still require the staged root
decision; no search was launched by this preparation or reviewer.

## Immutable evidence and cost

- Preserved original: `workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_v1/preparation/PREPARATION.json`.
- Actual diagnostic: `workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_caps_v1/execution/RESULT.json`,
  SHA256 `b59de647f1a0764836d04a1f5e17187978d5bda6b6e5470c781b1bcd33510603`.
- Recovery configuration: `workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_v2/CONFIG.json`,
  SHA256 `c355ba19e7b0a808ed3ba4e49e48bb39db4da5a0a1cf490c9ec89969c055b831`.
- Actual recovery: `workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_v2/preparation/PREPARATION.json`,
  SHA256 `6405de39b9d3ffa88aea7f99f895a685ddcd8ec5b1b37529766b72c365d9d57d`.

4819 used one Slurm CPU and0.546393s process wall.4821 used one Slurm CPU,
21.310315s process wall and19.692892 CPU-s. All three stages report zero new
energy calls. The recovery result and source-specific admissions were read here
without recomputing numerical quantities; root owns submissions and subsequent
scientific execution.
