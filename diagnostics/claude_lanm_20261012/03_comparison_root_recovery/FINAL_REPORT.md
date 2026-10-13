# LanM six-cell comparison: source reversal persists

**All six electronic endpoints are now available, and the two Hans structures
still give opposite conditional preference directions relative to Mex.** The
missing calculation did not resolve the discrepancy. The current fixed-state,
fixed-geometry model is not a robust LanM discriminator.

| Hans source versus Mex8FNS | Embedded electronic | Plus matched classical terms |
|---|---:|---:|
| La-conditioned 8DQ2 | +11.698489 | +11.689857 |
| Dy-conditioned 8FNR | −18.276098 | −18.269452 |

Units are kcal/mol. The frozen contrast is
`D = [E_Dy − E_La]_Hans − [E_Dy − E_La]_Mex`.
Positive means greater **conditional relative La preference in Hans**; negative
means the opposite ordering. Neither value is an absolute affinity or calibrated
biological classification. Each metal substitution has fixed stoichiometry;
the common elemental reference cancels in this protein–protein exchange.
Source states remain different, so this is not a comparison of their populations.

## What the result establishes

The fuller electronic region with the actual finite protein field still gives
source-dependent ordering. The separately evaluated matched classical terms
shift these two contrasts by only −0.008632 and +0.006646 kcal/mol and preserve
the reversal. Completing the numerical solve closed the missing-data question.
It does not establish which physical approximation causes the disagreement.

This is conditional EF3 substitution in a monomer, with fixed environment and
coordinates. It has not tested coupled scaffold relaxation, folding or dimer
response. The sources also retain different water and spectator states:

| Source | QM waters | Total source waters | Fixed EF1/EF2/EF4 spectators |
|---|---:|---:|---|
| Hans8DQ2 | 2 | 74 | La / La / Na |
| Hans8FNR | 3 | 179 | Dy / Dy / Dy |
| Mex8FNS | 7 | 167 | Nd / Nd / Nd |

Spectators are represented as formal monopoles. Thus geometry, hydration,
neighboring occupancy/identity, finite-field representation and electronic
method remain confounded. Do not label this uniquely a geometry problem or
select the favorable Hans source. Bulk solvent and state-population accounting
are absent. The full hybrid force qualification remains incomplete.

## Actual execution and checks

La 2665 and Dy 2666 both converged with native electronic and point-charge
gradients. Their final energies are −5804.400869806574 and
−5809.917822123823 Hartree, respectively. Their absolute unlike-metal totals
are not directly interpreted as affinity. La converged after 52 cycles; Dy
converged after 70, with final orbital residual 1.6467e−7 below 2e−6 and final
energy change −3.7289e−11 Eh below the 1e−9 tolerance. Final energy matches the
Dy gradient header; 198 QM / 2008 PC counts and the state/seed hashes agree.

Both targets have charge −1 and 816 explicit electrons. La uses ECP46. Dy uses
ECP55 with effective singlet valence, while its physical 4f9/sextet hypothesis
is recorded separately. This is a spin-free frozen-f approximation, not an
explicit open-shell Dy sextet or spin–orbit calculation. Small-gap warnings
remain; SCF convergence does not establish electronic ground-state stability.
The final Dy RIJCOSX grid correction is +0.000040984 Eh. Physical-state and
complete-force qualification flags remain false.

The original eight archived collector tests passed in Slurm 2664. Independent
read-only review checked raw Dy convergence and 52 comparison artifact pins.
All six collection energies/states agree, and the saved contrast algebra/sign
matches the pinned comparator. Four electronic endpoints and all six classical
cells were reused; only the two missing FNR electronic endpoints were recovered.
No scientific imports, structural recomputation or molecular evaluations were
performed on login for this report.

Job 2669 wrote `comparison_complete`, 6/6 accepted cells, no unavailable cells,
and zero new molecular calls. At review Slurm still labels that job RUNNING on
`node-48-256g-18`, while its worker reports `node-64-768g-18`. Scientific output
completion and scheduler termination are recorded separately. The terminal
watcher remains armed; no cancellation, rerun or collector duplication occurred.

## Cost, preservation and recommendation

| Completed continuation | Wall time | CPUs | Allocated CPU-seconds |
|---|---:|---:|---:|
| La 2665 | 9h04m00s | 12 | 391,680 |
| Dy 2666 | 12h35m50s | 12 | 544,200 |

These two continuation jobs total **935,880 allocated CPU-seconds**, zero GPU.
The comparator records 3.471861 seconds of actual work; its scheduler allocation
is not yet terminal. These costs exclude the earlier failed/interrupted runs
and four reused origins; those receipts remain preserved. This is reference
work, not demonstrated screening throughput. No runtime cap was loosened and
no convergence threshold changed to accept the endpoints.

**Close this static-source comparison as a negative robustness result.** Do not
launch another unchanged SCF retry or promote the model. Any next scientific
branch should distinguish source-state differences from coupled structural
response before claiming within-series selectivity; this completed comparison
does not itself justify a large relaxation or library campaign. No new molecular
job or production change follows from this event. La/Ca recovery is separate
and remains 119 correct / 3 wrong / 3 inconclusive across 125 inspected models.

[FINAL_RESULT.json](FINAL_RESULT.json) is the exact completed comparison;
[FINAL_PROVENANCE.json](FINAL_PROVENANCE.json) pins both continuations, raw inputs,
state inventories, accounting and the original output. The actual workspace is
`workspaces/claude_lanm_20261012/comparison_root_recovery_v1/`.
To inspect this delivery from the repository root:

```bash
cat diagnostics/claude_lanm_20261012/03_comparison_root_recovery/FINAL_REPORT.md
```
