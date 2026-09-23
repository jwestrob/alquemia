# Qualified restart restores the three missing summaries

The separate two-start numerical diagnostic completes the missing La-at-Ca-
proposal vacuum cell without changing its geometry, state, Hamiltonian or
convergence tolerances. Both actual native restarts agree within
1.44×10⁻⁷ kcal/mol. The origin-seeded value was designated before execution;
the other start is an agreement check, not a source of a favorable replacement.

The read-only pooled sensitivity retains every original record and the frozen
envelope reference. It adds one explicitly identified recovered cell and
recomputes the affected shared geometry pool and strict three-member medians.

| Result | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Original primary, 104 distinct source/context pairs |99|0|4|1|
| Separate recovery sensitivity, same 104 pairs |100|0|4|0|
| Original primary, 100 declared triples |91|0|0|9|
| Separate recovery sensitivity, same 100 triples |94|0|0|6|

The six remaining exclusions are the original unsupported triples. These100
overlapping triples belong to25 consumed reference protein groups; they are not
100 independent biological observations. Released and strict-tenfold methods
already classify the same94 available triples correctly. This restores practical
reference fidelity, not an accuracy gain over those group summaries.

The restored source is A0ACD6B9F2 La-conditioned sample0. Its recovered R is
−405445.69509985024 model kcal/mol, La-supported. Both metal rows select the
Ca-generated proposal after evaluating the complete common pool. Thus the missing
cross-evaluation could not legitimately have been ignored. Triples004/005/006
change from unavailable to correct; all other source and triple values are
unchanged. No threshold is fitted or transferred from a different model.

Limitations remain: four individual envelope pairs are inconclusive versus one
under the strict-tenfold comparator; the separate A0A3Ca3 pilot probe also loses
its earlier confident call. Across the91 primary complete triples, structural
ranges are not uniformly smaller: relative to same-context origins,39 decrease
and52 increase, with median range7.180→8.941 kcal/mol. The two PLM compression
examples show large reductions, but that observation is not a universal property.

## Status, artifacts and reproduction

This is a separately named post-failure sensitivity, not an overwrite of the
primary experiment or an implemented automatic retry policy. Current generic
execution v2 still uses fresh strict scalar inputs. Production and DFT remain
unchanged. The recovery required two scalar calls,21s allocation wall,
42allocated CPU-seconds, zero GPU; the original failed attempt is retained in
the primary campaign cost. See [the actual recovery report](REPORT.md).

Primary comparison SHA256:
`d58fa76be5e08971a8c837795c9a4363beff3665589e51088c36cc53d7248648`.
Recovery sensitivity:
`workspaces/motion_envelope_transfer_20260923/RECOVERY_SENSITIVITY_v1.json`,
SHA256`5874ace7232b768f86c7ec865bb40d3d3735b52ec1df867c47716e6c138cb59c`.
It pins the original comparison, qualified recovery collection, unchanged
reference, code, reconstructed matrix and all104 source/100 triple summaries.

Three actual-artifact tests passed in2.250s, zero skips: only the designated
source/three triples change while original bytes remain unchanged; a corrupted
copy selecting the other seed is rejected; a copy missing the agreement start
is rejected. No fabricated scientific values or new molecular evaluations.

This command reproduces only the algebra into a new output file:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/envelope_recovery_sensitivity.py \
  --primary workspaces/motion_envelope_transfer_20260923/COMPARISON_v1.json \
  --recovery workspaces/native_failed_cell_recovery_20260923/run_v1/COLLECTION.json \
  --output workspaces/motion_envelope_transfer_20260923/RECOVERY_SENSITIVITY_replay_v1.json
```
