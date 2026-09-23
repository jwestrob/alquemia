# Solvent-guided finite probes: no additional class repair

**Close this tested probe rule.** The eight-source pilot completed, but its
expanded geometry pool changes no calls on the frozen adaptive bands:
**7 correct / 1 wrong before and after**, with no missing scores. The remaining
error, A0A3F2YLY8 Ca-conditioned sample1, is unchanged. This is an inspected
development comparison, not independent validation or a calibration of the
expanded search. Production and all existing references remain unchanged.

## What the experiment establishes

Archived candidate rankings justified testing whether solvent preference could
suggest missing geometries. Across 204 complete noncanonical pools, native MACE
and composite selections differed for 93 Ca and 29 La rows; 90 Ca and 24 La
composite regrets exceeded 0.1 kcal/mol. All 225 sources, including 21 incomplete
pools, remain in that replay. The existing pipeline already selects by composite
energy; adding solvent selection itself was not a new method.

The frozen rule made two short, geometry-limited probes per metal around its
archived composite winner. Both metals then competed on the same expanded pool.
All four source-connected angular coordinates, chemical states, heavy-atom
displacement limits and numerical selection settings were preserved. Exact rule
and prior authorization: [PROBE_PLAN.md](PROBE_PLAN.md), [PLAN.md](PLAN.md).

Of 32 declared probe slots, 27 were geometrically admissible; five unsupported
directions remain recorded. These are partial search coverage, not failed
molecular calculations. All **54 native MACE and 108 native GFN2 calculations
succeeded**, with no retry, new DFT, optimization, fold, CPCM calculation or
threshold fit. All eight source scores are available.

## Actual discrimination and structural spread

R is E(Ca)−E(La), in model kcal/mol. The table transfers the existing adaptive
bands without changing them; the changed search has no independently calibrated
decision. Labels entered reporting only.

| Source | Expected class | ΔR versus adaptive pool | Adaptive-band call, before → after |
|---|---|---:|---|
| 1H4I | Ca | −0.160404 | Ca → Ca |
| 4MAE | La | −0.054573 | La → La |
| Q9Z4J7 | Ca | −0.244967 | Ca → Ca |
| Q88JH5 | Ca | 0 | Ca → Ca |
| A0A3F2YLY8 Ca sample1 | La | 0 | Ca → Ca |
| A0A3F2YLY8 Ca sample3 | La | +0.158465 | La → La |
| A0ACD6B9F2 Ca sample4 | La | 0 | La → La |
| A0ACD6B9F2 La sample4 | La | 0 | La → La |

Six of 16 endpoint selections adopt a new probe, each accepting higher native
MACE energy in return for a lower solvent contribution. The largest absolute
score change is 0.244967 kcal/mol. The two predeclared A0A3F2YLY8 source scores
spread slightly farther apart, 8.201646→8.360111 kcal/mol. The selected A0ACD6B9F2
pair stays 15.035742 kcal/mol apart. These are selected pairs, not estimates of
whole-protein conformational distributions or independent biological repeats.

The older released bands give 5 correct / 1 wrong / 2 inconclusive before and
5 correct / 0 wrong / 3 inconclusive after: Q9Z4J7 becomes inconclusive on that
different scale. That transfer is retained separately and does not establish
an accuracy gain for the new protocol. All new-protocol validated decisions
remain explicitly null.

## Numerical qualification

A parallel source-start experiment exposed a roughly 4.77 kcal/mol native GFN2
contrast discrepancy at nearly identical Q88JH5 coordinates. Its numerical audit
is owned separately. Consequently, these solvent energies and structural effects
remain descriptive; this pilot does not establish their physical accuracy.
That limitation does not erase the observed absence of class repair.

Across all 54 new endpoint cells, the largest solvent change from its actual
generator center is 1.882467 kcal/mol. It occurs at Q88JH5 Ca on a 0.10-radian
probe with 0.272857 Å maximum context-heavy displacement, accompanied by
2.030286 kcal/mol native work. This is a finite displacement, not a repeat at
nearly identical coordinates. It neither reproduces nor resolves the parallel
SCF concern. No extra probes or retries were introduced after seeing results.

## Measured execution cost and reproducibility

| Job | Actual allocation | Elapsed | Allocated core-seconds |
|---|---|---:|---:|
| 1210100 | 32 CPU, one H200, 200000 MiB host | 35 s | 1120 |
| 1210101 | 64 CPU, 128 GiB host; eight 8-rank GFN2 workers | 267 s | 17088 |
| Total | 54 MACE + 108 GFN2 calls | — | **18208** |

GPU allocation was **35 seconds**. Inner executors account for
17126.579533 core-seconds and 21.348035 GPU-worker seconds. Local preparation and
reporting were not metered and are additional. This is no matched-hardware speed
comparison. Both Slurm jobs completed with exit0 and all scientific cells valid.

The immutable final selection is
`workspaces/solvent_guided_20260922/pool_v1/after_solvent_0_1210101.json`.
`comparison_v1.json` retains every source, selection, component, frozen-band
transfer and selected-pair spread. `probe_component_audit.json` retains the
all-cell displacement/component check. `COSTS.json` and scheduler/worker receipts
retain actual allocation and attempts. [ARTIFACTS.json](ARTIFACTS.json) pins these
files; [COMMANDS.md](COMMANDS.md) provides a report-only replay and focused tests.
All nine focused real-artifact tests pass, with zero skips: archived ranking and
signs, missing-result handling, physical subspace/backoff, label-independent
selection, and actual completed/MACE-only report paths. These parser/geometry
tests are separate from the 162 successful scientific calls above.

Recommendation: retain the current scorer and close this finite line-probe rule.
It found solvent-favored coordinates but earned no extra classification or
structural-robustness benefit. No wider rescore or production integration follows.
