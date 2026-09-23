# Saved native restart seeds: all 80 cells available

Read-only follow-up, 2026-09-23 UTC. All four proposed five-candidate pools
have complete matching GBW and xtbw files for both metals in both media.
**80/80 compatible seed pairs; zero missing. No new molecular calls or inputs.**

| Consumed source | Available / required |
|---|---:|
| 1H4I | 20/20 |
| 4MAE | 20/20 |
| Q88JH5 | 20/20 |
| A0A3F2YLY8, Ca-conditioned sample 1 | 20/20 |

Each set comprises `origin`, `proposal_Ca`, `proposal_La`, `adaptive_Ca`
and `adaptive_La`, crossed with Ca/La and vacuum/ALPB. Selection comes from
the existing common-panel INPUTS; no candidate or energy selection was added.

## Compatibility checked

The actual source task, successful receipt, output, input and XYZ pins agree.
Coordinates match their existing pool candidates within 1e-12 Å; Ca/La
coordinates and corresponding vacuum/ALPB coordinates are identical. The
actual native GFN2 path, ORCA 6.1.1, charge/spin, parameter export, 300 K
recipe, eight MPI ranks and recorded endpoint energy also agree. Both saved
files are nonempty and come from the same completed task/runtime basename.

The source iteration cap is **125 for 48 cells and 500 for 32**. This numerical
difference is retained explicitly; the physical method and remaining primary
settings agree. Source artifacts span the earlier static, proposal, common
pool and angular pool protocols. Their compatible endpoint recipes are
verified individually, rather than assuming protocol names establish identity.

GBW/xtbw hashes were captured in this inventory. Where older receipts did not
hash these files, this is a current pin of the matching saved pair, not a claim
of a historically recorded seed hash or previously demonstrated consumption.

## Next scope remains proposed

The archive supports all 80 proposed self-continuations without generating
replacement seeds. This inventory does **not** execute or validate that
rescore. The separate restart experiment demonstrated activation for its
Q88JH5 La test; actual XTBRESTART/native-mixer evidence would still be required
for every new continuation. No new scores, classifications or calibration
are reported here.

Exact per-cell paths, hashes, states and source caps:
[SEED_AVAILABILITY.json](SEED_AVAILABILITY.json). The JSON pins the read-only
inventory implementation in
`workspaces/native_xtb_restart_20260922/seed_availability_v1/inventory.py`.
