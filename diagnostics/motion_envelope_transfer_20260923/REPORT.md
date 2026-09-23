# Motion-envelope full100 transfer: fidelity retained where complete, coverage reduced

**All 91 complete three-source summaries are correct.** The new envelope recovers
both earlier A8 threefold abstentions and matches released/strict-tenfold decisions
on the same 91 triples. One real SCF failure makes three additional triples
unavailable: the primary result is **91 correct / 0 wrong / 0 inconclusive / 9
unavailable**, versus 94 correct / 6 unavailable for strict tenfold scoring.
This is a useful practical three-source result with a numerical coverage limitation;
it does not establish equivalent coverage or justify default promotion.

## Frozen comparisons

| Method, each with its own canonical-only reference | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Released static | 94 | 0 | 0 | 6 |
| Tenfold strict adaptive | 94 | 0 | 0 | 6 |
| Original 3.5 Å threefold adaptive | 92 | 0 | 2 | 6 |
| 4.3 Å envelope origin only | 92 | 0 | 2 | 6 |
| 4.3 Å envelope adaptive | 91 | 0 | 0 | 9 |

On exactly the common 91 complete triples, envelope adaptive and strict tenfold
are both 91 correct; envelope origins and old threefold each give 89 correct and
two inconclusive. Search repairs origin abstentions in A0A3 triples 000/002. Both old
A8 abstentions (008/010) are correct. The old 3.5 method used the historical looser
scalar policy, so its difference is not a pure context-radius effect. The strict
tenfold comparator uses the same current scalar profile. Mathematical and
operational envelope variants have the same classification counts.

The reference was frozen from the original 25 canonical geometries before these
transfer energies. No new thresholds, member selection, label fitting or production
change occurred. All 100 original memberships remain. These are correlated subsets
of 25 consumed proteins, not 100 independent biological validations or affinity tests.

## Structural variation on the matched triples

On the same 91 complete triples, envelope accommodation narrows the range in
39 and widens it in 52 relative to its own origins (median range 7.180465 →
8.940608 model kcal/mol). Against released static scoring, 42 ranges shrink and
49 grow (7.633593 → 8.940608). Against strict tenfold, 49 shrink and 42 grow
(7.950702 → 8.940608). The median paired change versus strict tenfold is −0.107064;
it differs from the change between the two marginal medians. There is no uniform
spread improvement. See [all matched ranges](SPREAD.md) and [unrounded values](SPREAD.json).

## Individual sources and coverage

There are 104 distinct source/context pairs from 98 physical source structures and
300 declared triple-member occurrences. These denominators are reported separately.
Unique pairs give 99 correct / 4 inconclusive / 1 unavailable, versus 103 correct /
1 inconclusive for strict tenfold evaluated on those same 104 pairs. The
inconclusive pairs are C5AX La-sample3 (already inconclusive), Q60AR6 La-sample3 under
two different envelope selections, and Q88JH5 La-sample3. Their complete triple
medians are correct without removing or replacing members. Three newly inconclusive
pairs are a limitation even though aggregation resolves their calls.

The six old preparation exclusions remain: triples 037–039 (MMOL1770) and
076/077/079 (Q88JH5). The new numerical failure invalidates A0AC triples 004–006.
The separate Ca-conditioned A0A3 sample3 pilot regression remains inconclusive;
a favorable La-conditioned transfer does not erase it.

## Actual numerical failure

La evaluated at the Ca-proposed geometry of A0ACD6B9F2 La-sample0 vacuum fails
native SCF after 500 cycles. Its last energy changes remain around 10^-4–10^-3
Hartree, well above TolE1e-10. The native process returned 0 but explicitly reported
error termination; its receipt correctly records no normal termination or SCF
convergence. No energy is accepted and the required whole pool remains unavailable.
The cell took 201.459371s. It used one rank, maxcore 2000MB; the outer allocation had
64 GiB/32 CPU and batch peakRSS approximately 6.37 GiB, which is not a per-cell memory
measurement. There is no OOM evidence. See FAILURE.json for exact pins and accounting.
The separately owned [two-seed restart diagnostic](../native_failed_cell_recovery_20260923/REPORT.md)
is excluded from this primary result.

## Execution, costs and integrity

Three exact completed pilot pools were reused; 101 new pools were executed in 51/50
source shards. All 202 searches yielded admissible proposals. There were 3,408 new
MACE calls (202 origins, 3,004 search evaluations, 202 cross evaluations) and 1,212
strict native GFN2 attempts (1,211 accepted, one failed). No DFT, new folds, states,
waters or protonation changes occurred. Both metals score the same candidate pool.
97 final proposals touch a bound, with no unconstrained-minimum claim. The unchanged
SLSQP engine records 641 evaluated trials outside final admissibility; none is an
admitted final proposal. The motion-envelope argument applies to admitted original
anchor displacement against omitted fixed atoms, not every optimizer trial.

All 21 outer allocations, including empty/failed technical stages, cost 82,389 allocated
CPU-seconds and 1,571 requested GPU-seconds. GPU time is requested allocation time,
not measured utilization; this cluster omits GPU from AllocTRES, so pinned wrapper
requests supply that count. Local preparation and read-only collection/recovery time
are additional and not included as metered allocation cost. No nested runner costs
are added twice. Search executors took 503.780/466.984s for 102/100 searches, with
worker times 497.090/459.265s and model loads 1.174/1.185s. This is prepared-batch timing,
not folding-inclusive or isolated three-source latency.

Two metadata defects were recovered without repeating molecules. First, exact-byte
runner/renderer copies at different paths caused receipt rejection; recovery retains
all actual paths and checks the original receipts. Second, a missing scalar cell join
field was reconstructed from exact case/metal/candidate identity plus physical state.
Original inputs, receipts and rejected collections remain intact. Separate JOIN files
record the recovery. Ten failed allocation statuses comprise nine technical/reporting
failures plus the one allocation containing the real SCF failure.

Nine focused real-artifact tests pass: five scope/recovery checks and four final
matrix/median/denominator/old-context-join checks. Raw scores, components, works,
spreads and every member occurrence remain in the full comparison. The original 3.5 Å
context is joined separately for each triple where larger envelopes merge older
selections. No favorable old context is chosen.

[Compact result](RESULT.json) · [triples](TRIPLES.csv) · [pairs](PAIRS.csv) ·
[commands](COMMANDS.md) · [failure](FAILURE.json).
Full comparison: `workspaces/motion_envelope_transfer_20260923/COMPARISON_v1.json`.
