# Alternative starts: no new basin; native solvent sensitivity found

**The eight searches return the existing adaptive basins. This policy adds no
demonstrated discriminatory utility.** It does reveal a concrete numerical concern:
Q88JH5's native GFN2 solvent transfer shifts 4.805 model kcal/mol when the optimized
coordinates change by less than one millionth of an ångström. Exact-input repeats
reproduce both outcomes. That difference is not evidence of better geometry.

The production/default scorer and every older output remain unchanged. This is
a consumed-reference development experiment, with no new calibration or independent
biological validation. The full eight-source denominator is retained.

## What was tested

Two reference-derived starts were specified for each of eight target sources,
using real Ca-associated 1H4I and La-associated 4MAE dihedrals. Self-protein targets
use the frozen canonical Q9Z4J7/A0A3F2YLY8 replacements. Only identical role/residue
chemistry is mapped, through each target's existing four selected angular modes.
Absent chemistry stays unchanged and explicit: the Ca template supplies no extra
Asp, so the admitted 4MAE start changes only its two mapped Asn angles. Carboxylate
O-name equivalence is used only for the archived deprotonated identical groups.
No atoms, protons, waters, states or source coordinates are transplanted directly.

Only **4/16 starts** fit the unchanged ±0.8 radian / 0.8 Å domain and chemical guards:
1H4I ← La, 4MAE ← Ca, Q9Z4J7 ← La, Q88JH5 ← La. Both starts for all four difficult A0A3F2YLY8/
A0ACD6B9F2 folds are rejected: required heavy motion is approximately 1.14–2.42 Å,
with some angles exceeding the limit by more than threefold. Their branch status
is `no_admitted_distinct_start`; their existing baseline/adaptive results remain
available separately. No bound was enlarged, template clipped or replacement
chosen after seeing a score. Homology-independent transfer is not claimed.

The eight admitted metal/start searches use the unchanged native OMOL checkpoint,
source-based physical map and scaled four-angle SLSQP policy. All eight converge
successfully. Their own-metal energies differ from the original adaptive search
by at most 2.04e−10 eV; largest final angle difference is 2.18e−5 radian. Thus these
starts return the same local basin to optimizer accuracy. The common scorer then
cross-evaluates all eight new geometries for both metals, retaining every old
candidate. All 8 cross-MACE and 32 GFN2 calls succeed normally.

## Paired outcomes

Changes below are relative to the already completed adaptive common pool, using
`R = min E_Ca − min E_La` and unchanged
`E_M = E_OMOL + E_GFN2,ALPB − E_GFN2,vacuum`.

| Source | ΔR, model kcal/mol | Transfer under unchanged adaptive bands |
|---|---:|---|
| 1H4I | −0.000013791 | Ca remains Ca |
| 4MAE | +0.002537201 | La remains La |
| Q9Z4J7 | +0.000003770 | Ca→inconclusive at exact calibration edge |
| Q88JH5 | +4.769598660 | Ca remains Ca; numerically unqualified change |
| Four difficult folds | unavailable | No admitted start; no claimed improvement |

Mathematical and operational pool results agree. Literal old-band transfer is
3 correct / 1 inconclusive / 4 unavailable. Q9's 3.77e−6 shift crosses its exact fitted
Ca boundary; the record is preserved without numerical padding. It is **not a
substantive predictive failure**, nor a reason to widen the bands. Q88's much
larger shift is almost wholly from the solvent subtraction. Neither result
supports adding this search stage to production.

## Q88 component and reproducibility diagnosis

The old La-adaptive and new La-template-search geometries differ by only
6.4623e−7 Å maximum atom displacement. Atom order, charge −2, singlet state,
438 electrons, 300 K, native GFN2 parameters and actual input controls agree.
The actual runtime inputs are byte-identical for each medium: NoAutostart,
UseXTBMixer true, MaxIter 500 and 8 MPI ranks. Parameter exports also agree exactly.

| La endpoint at nearly identical geometry | Old, Ha | New, Ha | Change, kcal/mol |
|---|---:|---:|---:|
| Vacuum GFN2 | −252.116197616031 | −252.108540025437 | +4.805210646 |
| ALPB-water GFN2 | −252.570227444230 | −252.570227099694 | +0.000216200 |

Native OMOL is identical at printed energy precision. A higher vacuum energy is
subtracted in the composite and therefore lowers the new La score. The Ca vacuum
change at this geometry contributes a further −0.0354 kcal endpoint shift, giving
the pooled +4.7696 change in R. This component accounting does not identify a unique
physical or solver cause; it rules out presenting the shift as an established
new native-MACE basin. [Actual audit](Q88_AUDIT_v1.json).

Eight further approved calls used old/new exact XYZ × vacuum/ALPB × two fresh
identical-primary repeats on one host. All eight report normal convergence.
Each pair of identical-input repeats gives **exactly the same printed energy**.
The old XYZ on this host reproduces the prior other-host vacuum value within
6.26e−7 kcal/mol; the new XYZ reproduces exactly. The matched La transfer difference
remains −4.804993821 kcal/mol. The observed sensitivity therefore survives the
same-host reproducibility check; scheduling/host variation does not explain it
at this scale. No favorable repeat replaces an older result.

The old vacuum output's printed MAX/RMS density changes exceed their displayed
thresholds; both ALPB outputs do as well, despite normal native-SCF termination.
The new vacuum passes those printed diagnostics. These are recorded observations,
not an assertion that generic ORCA density criteria govern the native mixer.
[ORCA 6.1 documents a special native-xTB SCF path that overrides other SCF settings](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/semiempirical.html#scf-with-native-xtb).
Consequently, adding generic Convergence Tight and seeing a smaller printed tolerance
would not establish tighter native convergence. The proposed Tight manifest was
preserved **unsubmitted**; it was replaced before numerical execution by the
exact-input repeats. Native tightening remains unavailable, not a successful
qualification. No ordinary SCF retry was run.

## Cost, reproducibility and next step

| Job | Work | Allocated CPU-s | Requested GPU-s |
|---|---|---:|---:|
| 1210098 | 8 searches, 287 native MACE evaluations | 3,200 | 100 |
| 1210103 | 8 cross-MACE singlepoints | 768 | 24 |
| 1210104 | 32 primary GFN2 singlepoints | 5,952 | 0 |
| 1210115 | 8 exact-input GFN2 repeats | 1,792 | 0 |
| **Total** | **295 MACE + 40 GFN2; 0 DFT** | **11,712** | **124** |

Costs include reserved idle capacity and in-job collection. Local preparation,
tests and reporting are additional unmetered work. These are development costs,
not measured source-inclusive production latency. Peak measured CUDA allocation
in the search worker was 6,056,950,784 bytes. Twelve real-artifact tests pass with
zero skips; parser/algebra/replay tests are separate from the actual molecular
runs above. Explicitly corrupted real copies test missing role/convergence fields.

[Commands](COMMANDS.md) validate and replay the saved results. [RESULT.json](RESULT.json)
pins preparations, actual collections, comparison, costs and tests. Earlier
preparation drafts remain unexecuted; the final source implementation is pinned
in prepared_v3 and molecular executors are snapshotted in each actual run.

**Recommendation:** close this bounded alternative-start policy. A wider real
rotamer domain requires a separately defined physical model and source-wide
collision checks; it must not be obtained by clipping these templates. Before
attributing additional composite-search improvements to geometry, investigate
native GFN2 electronic-solution continuity. An eight-call native .xtbw cross-start/
self-start test is [proposed separately](RESTART_PROPOSED.md), with actual saved
compatible files; it has not been launched. It may identify a numerical remedy,
but no convergence fix or new classifier has been claimed here.
