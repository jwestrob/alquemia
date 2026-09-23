# Strict stopping preserves the accuracy gain and removes the Q4 regression

**Use the fresh strict scalar setting for the next opt-in candidate qualification.**
It retains both repaired released errors, removes the Q4 numerical regression,
and adds no classification regressions across the full supported panel. This
supports the numerical change. Root's now-completed, separately frozen static
ablation also isolates accommodation:199 correct /1 wrong /8 inconclusive on
the same208 sources versus207 /0 /1 after movement. All eight repairs require
physical score changes even under the adaptive bands; see the separate
[strict-static report](../strict_static_ablation_20260923/REPORT.md). These eight
structures represent four protein groups. Production stays unchanged.

The final225 ledger is207 correct /0 wrong /1 inconclusive /17 unavailable.
On the same207 sources available to the released static method, strict scores
206 correct /0 wrong /1 inconclusive versus203 /2 /2. On the same208 sources,
prior precision scores206 /0 /2 and strict207 /0 /1. All17 original preparation
exclusions remain. All94 available La triples remain correct; no gain is claimed
where aggregate accuracy was already complete.

The sole changed decision is Q4W6G0 Ca-conditioned sample2, inconclusive→Ca,
with ΔR−7.987209876 model kcal/mol. **This is a previously inspected strict32
reuse, not a newly discovered repair in the200 fresh pools.** All200 fresh pools
retain their prior decisions. A0A3F2YLY8 Ca1 and A0ACD6B9F2 Ca4 retain their
correction from the released method's wrong Ca calls to La-supported.
C5AXV8 La3 remains inconclusive, R−405457.848334666. No threshold adjustment,
source replacement, additional search or automatic rescue was used.

For Q4, the La-at-adaptiveCa vacuum component decreases11.155581210 kcal/mol,
while matched ALPB decreases only0.001161289. The solvent transfer contribution
therefore increases11.154419920, and the selected La candidate changes from
adaptiveCa to adaptiveLa. Geometry and MACE are unchanged. This establishes a
numerical-electronic solution effect on selection, without establishing a unique
physical electronic ground state. All225 classification decisions are unchanged if transferred
through the old precision bands; mathematical and operational policies also
produce identical decisions. Thus the repaired call is not a calibration trick.

Strict32 previously qualified fresh versus original-cold-seeded starts on384
cells. This transfer used one fresh start per new cell and does not repeat that
agreement test on all2496 cells. Scalar repeatability is not analytic-gradient
qualification. These are consumed PQQ protein/fold controls, not new independent
La/Ca affinity measurements or physiological-occupancy evidence.

## Actual throughput and overhead

All2400 new scalar calls and96 exact reused cells pass collection; zero molecular
failures or retries. New execution uses no MACE, search, DFT or preparation
chemistry. Four32CPU/64GiB jobs completed in1316,1316,1314 and1315s, respectively.
Total168352 allocated core-seconds,18960.007 accounted CPU-seconds, zero GPUs.
Maximum batch MaxRSS is6451544K. No undocumented CPU-time ceiling was applied.

The observed span from first molecular launch to last completion is1182.501771s:
**2.029595 new scalar calls/s across the four simultaneous batches.** Summed
executor wall is4779.779032s, representing152952.929011 allocated core-seconds.
Full allocation adds15399.070989 core-seconds for validation, setup, collection
and scheduler/accounting boundaries. Do not equate elapsed ORCA time with CPU
consumption. The source/manifest audit adds150.387s local wall; implementation,
tests and reporting have additional unmetered local cost.

For context only, the earlier single32-worker fresh qualification completed384
calls in226.770720s of executor time (1.693340 calls/s). The source mixture and
host load differ. These measurements cannot establish128-versus32 scaling or a
hardware-controlled speedup. Existing32-worker routine policy is unchanged.
Historical MACE, strict32 qualification and the separately owned four-pool
comparator costs are excluded from new transfer totals and remain recorded.

Six final real-artifact tests pass in3.507s, zero skips, including full prior-ledger
preservation, exact matrices/bands, four disjoint shards and explicit
missing-cell behavior using a corrupted copy of a real result. Actual scientific
execution is separately evidenced by the2400 successful ORCA receipts.

Full report: [REPORT.md](REPORT.md). Full source/matrix/aggregate data and exact
cost/timing receipts are pinned in [ARTIFACTS.json](ARTIFACTS.json). The next
scientific decision should retain accommodation and use the separate
motion-envelope qualification for the practical three-source interface, without
expanding or tuning C5AXV8
merely to force its remaining abstention into a class.
