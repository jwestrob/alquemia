# PedE: solvent response, with substantial ligand charge redistribution

**The strongest clue is the solvent contribution at the starting geometries.**
Between correct two-Ca sample 3 and wrong sample 1, the origin contrast changes
by +21.106215 kcal/mol: OMOL contributes −0.150413 and the matched GFN2
ALPB-minus-vacuum component contributes +21.256628. The wrong sources also show
stronger PQQ/nearby-residue charge redistribution in ALPB. This identifies a
specific electronic/solvent response to investigate; it does not establish that
the folds are bad, that a discrete electronic branch switch occurred, or that
this response is physically wrong.

## Actual scope and accounting

Slurm 3341 decomposed every saved metal/candidate cell for all five recovered
two-Ca PedE sources and four exact-method La-conditioned controls (54 cells).
Slurm 3345 read all 36 origin points in both saved backends, including actual
charge files, execution receipts, SCC tables and complete preparation mappings.
There were **zero new molecular evaluations, force calls, optimizations or fits**.
These are consumed, correlated structures of one protein, not independent blind
observations. Historical La sample 1 is absent from the compatible energy panel;
the separate geometry audit retains its original offsite-metal exclusion.

The scored energy is `E = V + U`, with `V = E_OMOL,vac` and
`U = E_GFN2,ALPB-water − E_GFN2,vac`. The solver is the corrected-method
standalone xTB 6.7.1, checked against tblite 0.7.0, not the old native-ORCA
parameter assignment. The matched subtraction includes electronic redistribution
between vacuum and continuum endpoints; it is not an isolated, fixed-charge
protein electrostatic energy. `R = E_Ca − E_La`; larger R is more La-like on this
protocol's calibrated scale. These quantities are not binding free energies.

## The origin effect is larger than the distinctive accommodation effect

Both metals become more stabilized by U in sample 1 than in sample 3:
Ca changes by −56.729299 and La by −77.985927 kcal/mol. Their difference gives
the +21.256628 kcal/mol contribution to R. Thus this is not simply a bad calcium
endpoint or a constant shift applied only to La.

| Two-Ca sample | Final call | Origin U_Ca | Origin U_La | Accommodation change in R: OMOL / solvent / total |
|---|---|---:|---:|---:|
| 0 | correct Ca | −371.573404 | −255.598041 | +6.835944 / +2.833149 / +9.669093 |
| 1 | wrong La | −424.382214 | −326.721608 | +6.451441 / +2.967590 / +9.419031 |
| 2 | wrong La | −403.982195 | −300.070382 | +9.482000 / +0.721345 / +10.203345 |
| 3 | correct Ca | −367.652915 | −248.735681 | +7.637878 / +2.945342 / +10.583220 |
| 4 | correct Ca | −381.366121 | −268.515115 | +5.325786 / +0.646427 / +5.972213 |

All energies in the table are kcal/mol. Samples 1/2 begin inconclusive; 0/3/4
begin Ca-correct. Correct sample 3 actually has the largest accommodation shift
among the five. Against the canonical PedE origin, samples 1/2 have OMOL changes
−2.513729/−0.721848 but solvent changes +11.641803/+5.390595 kcal/mol. Their
starting contrast is already displaced toward La before donor motion.

Calibration matters. The independent [frozen-band analysis](../CALIBRATION_RESULT.json)
finds that the Ca edge itself moves +8.189098 kcal/mol between static and
accommodated protocols. The +9.419/+10.203 raw changes therefore worsen the
expected-class margins by only 1.230/2.014 kcal/mol; correct sample 3 worsens by
2.394. The canonical abstention gap contracts from 10.297679 to 3.472171 kcal/mol.
The failures cannot be attributed solely to an unusually large relaxation work.

All four compatible La-conditioned PedE sources remain correct under both
policies. Their origin contrast range is 13.960892 kcal/mol, with 10.911753 OMOL
and 3.049139 solvent contributions using the same two range-defining sources.
Canonical La sample 2 changes by +8.761966 kcal/mol on accommodation. The full
per-control values and every cross-scored cell are retained in the generated
energy report, not pooled across differently conditioned sources.

## The saved populations support a charge-redistribution hypothesis

The table contains actual xTB population charges, in elementary-charge units.
PQQ includes its three fixed hydrogens. Arg386 is the retained local residue
fragment; synthetic cap charges are separate. These are population-analysis
values, not formal oxidation states or measured residue charges.

| Two-Ca sample | PQQ Ca: vacuum → ALPB | PQQ La: vacuum → ALPB | Arg386 Ca / La in ALPB |
|---|---:|---:|---:|
| 0 | −2.5972 → −3.0607 | −2.2162 → −2.9291 | +0.8903 / +0.9863 |
| 1 | −2.6398 → −3.7524 | −2.2568 → −3.6392 | +1.4065 / +1.5320 |
| 2 | −2.4932 → −3.4183 | −2.1308 → −3.3274 | +1.2932 / +1.4271 |
| 3 | −2.5147 → −2.9235 | −2.1258 → −2.7813 | +0.8822 / +0.9225 |
| 4 | −2.5889 → −3.1083 | −2.1996 → −3.0306 | +1.0216 / +1.0525 |

Metal populations vary much less: for example, ALPB Ca is +0.60284249 in wrong
sample 1 and +0.61263377 in correct sample 3; La is +0.72426776 and +0.75078717.
The larger differences involve the cofactor and ligand/context groups. All
mapped residue groups, both backends and all individual atom charges remain in
the saved result and corrected TSV exports.

This is a useful correlation, not proof of a specific donor-to-PQQ electron
transfer mechanism. Correct canonical La sample 2 also has strongly negative
ALPB PQQ populations (−3.4468 Ca/−3.3161 La), with a different distribution among
nearby residues. Population sums alone cannot diagnose the correct electronic
solution, determine whether charge delocalization is excessive, or justify a
new rejection threshold.

## Convergence passed; electronic-surface accuracy remains unqualified

Every required component was previously admitted. Across all 54 cells, the
largest saved component discrepancy is 1.2900415e−7 Eh against the frozen
1e−6 Eh gate; the largest population discrepancy is 8.3335735e−8 e against
1e−5 e. Raw population arrays exactly reproduce the stored admitted arrays.
All 192 atoms map to explicit preparation records, with the metal handled
separately; this avoids confusing 189 physical motion coordinates with QM order.
Paired coordinates, chemical-state signature, −3/−2 charges, singlet states,
model and runtime/settings compatibility pass unchanged.

Actual xTB SCC iterations, listed as vacuum/ALPB:

| Two-Ca sample | Ca iterations | La iterations |
|---|---:|---:|
| 0 | 45 / 41 | 68 / 46 |
| 1 | 67 / 38 | 86 / 44 |
| 2 | 53 / 39 | 87 / 44 |
| 3 | 49 / 39 | 81 / 46 |
| 4 | 45 / 42 | 77 / 43 |

Both backends report convergence throughout. For example, sample 1 ALPB Ca's
last printed SCC energy change is +1.76215e−12 Eh with residual 6.57e−10;
ALPB La is −5.68434e−14 Eh with residual 5.74e−9. The actual SCC tables, including
tblite counts and final residuals, are retained rather than inferring iterations
from mixed orbital-population output.

Wrong sample 1 has small ALPB HOMO–LUMO gaps (0.004696 Ca/0.012594 La eV) and
fractional occupations at the existing 300 K electronic temperature. Those are
not wrong-only features: correct canonical La sample 2 has gaps
0.022685/0.007722 eV, and another correct historical control has 0.005845 eV for
La. These observations make state/response accuracy worth checking but do not
prove a discontinuous branch switch. Tight convergence and cross-implementation
agreement do not validate the underlying Hamiltonian.

Water stdout retains the warning `Unknown reference state 'gsolv'`. Its subsequent
setup explicitly prints `gsolv [1 M gas/solution]` and Gshift 0.001080759698 Eh
(sample-1 La example). The warning is real; the effective printed convention is
also known. A common coordinate-independent reference offset cancels in the
within-metal source differences. There is no evidence that this warning explains
the observed 21.26 kcal/mol source spread. The old absolute-v1 parity failure
also remains failed; component-v2 acceptance does not erase it.

## Proposal selection and interpretation

For all five sources, La selects its own proposal. Ca selects its own proposal
for 0/1/4 and the La proposal for 2/3. The latter two Ca choices differ from their
vacuum minima because the actual solvent-containing energy selects the candidate.
The complete six-cell matrices retain both cross evaluations. Negative endpoint
works are not evidence of composite stationary minima: proposals use vacuum
forces, and no composite gradient or unrestricted minimum was validated.

All saved searches pass their unchanged numerical and geometry admission rules.
The two boundary-touching Ca proposals belong to samples 2 and 3 and are not
selected. Wrong sample 1 is interior throughout; selected La proposals are
interior for every sample. No wrong-only boundary failure is established.
The separate [geometry report](../geometry/REPORT.md) finds no gross defect or
confidence rule that isolates the wrong sources.

**Recommendation:** test the metal-dependent solvent/electronic response at these
fixed origins before increasing accommodation freedom or discarding folds.
The current evidence localizes sensitivity to a specific component and plausible
charge redistribution, but cannot yet separate a legitimate structural response
from a model error. Do not retune thresholds, alter protonation, or remove wrong
sources based on these inspected outcomes. No further calculation is launched by
this report.

## Artifacts, export repair and measured analysis cost

[EVIDENCE.json](EVIDENCE.json) pins the actual output/receipt files and contains
copied per-source values. Complete workspace outputs are under
`workspaces/pede_failure_diagnosis_20261013/`:

- `energetics_v1/execution/`: `RESULT.json`, `CELLS.tsv`, `WORKS.tsv`, generated
  `REPORT.md` and `EXECUTION_3341.json`.
- `energetics_charges_v1/execution/`: authoritative `RESULT.json`,
  `ATOM_CHARGES.tsv`, `ORIGIN_SOLVENT.tsv`, `EXECUTION_3345.json`.
- `energetics_charges_export_v2/execution/ORIGIN_CHARGES.tsv`: corrected 72-row
  group-charge table; `EXECUTION_3348.json` records export-only work.

The original charge TSV had a column-name collision: numeric group `metal`
overwrote element identity during flattening. Its JSON was unaffected. Job 3348
copied the immutable JSON into a new table with `charge_e:` prefixes and explicit
Ca/La identity checks. The old TSV remains preserved and **superseded**; it must
not be used as an identity-bearing export. No charge analysis or chemistry was
repeated for this repair.

| Job / work | Process wall seconds | Process CPU seconds | Peak RSS KiB |
|---|---:|---:|---:|
| 3341 energy decomposition | 0.643151 | 0.368029 | 63512 |
| 3345 saved charge/SCF analysis | 1.352954 | 0.846636 | 34460 |
| 3348 export-only repair | 0.022937 | 0.085623 | 22908 |

These are measured saved-artifact analysis costs, not new-site scoring latency or
the cost of the previously executed molecular calculations. Every job used one
actual Slurm CPU. No production protocol, calibration or biological label changed.
