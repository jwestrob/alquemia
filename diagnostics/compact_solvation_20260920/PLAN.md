# Compact MACE plus matched GFN2 solvation transfer

## Approval and question

2026-09-20 Jacob approved the discussed experiment: “go for itl”, followed after
an interruption by “ok. go for it. sorry for the interrupt.” Standing discretionary
classifier-improvement and parallel work permission remains in force.
Question: can inexpensive solvent screening preserve the complete-context MACE
alpha/GGR gain while restoring robust PQQ class separation? The completed native
CPCM DFT charged-context comparison motivates the test: three Ca-family context
shifts are+1.111/+2.242/+3.278kcal/mol, versus native vacuum MACE+49–51model-kcal.
This suggests solvent response but does not uniquely identify the discrepancy.

## One frozen physical expression

At exactly the archived geometry, composition, total charge and singlet state:

E_composite,M = E_nativeOMOL,vac,M + (E_nativeGFN2,ALPBwater,M − E_nativeGFN2,vac,M).

R_composite = R_nativeOMOL + (delta_solv,Ca − delta_solv,La).

Use actual final native GFN2 single-point totals consistently in both environments;
record printed components and default electronic smearing. This is a low-level
solvation transfer for a vacuum learned endpoint, not self-consistent QM/PB,
a protein-shaped dielectric, an addition to CPCM DFT, or binding free energy.
ALPB treats the whole supplied cluster boundary. Atomic radii, element/solvent
parameters are installed ORCA6.1.1 native defaults and never fitted to labels.
No changed waters, protonation, caps, core membership, coordinates, label or source.
GFN2 all-electron-equivalent parameters have their own valence convention; verify
actual Ca/La parameter exports/electron counts, not an assumed DFT ECP count.

## Fixed inventory and run order

An explicit hash-pinned inventory contains25canonicalPQQ,3consumed crystal
transfers(1H4I,4MAE,1KB0),2alpha structures and3GGR structures:33cases. Preserve
both core and complete-context representations where their matching real native
MACE receipts exist. The33×2representations×2metals×2environments define264distinct
low-level single points; missing/unsupported source data remain explicit.

1. Syntax/cost/convergence pilot: context1H4I andq9z4j7-pqq-la_model, bothCa/La,
   vacuum/ALPBwater =8nativeGFN2calls. Charged versus neutral context selection
   follows the question, not a new score. All MACE energies are reused.
2. Numerical check: repeat those same8calls with ordinary ORCA SCF and TightSCF,
   instead of the native xTB-specific mixer. Same GFN2 Hamiltonian/parameters,
   nuclear state and300K electronic smearing. This is a numerical comparison,
   not a favorable-score backend choice. Default mixer remains primary.
3. If primary computations are supported/converged and actual cost permits useful
   throughput, execute the remaining256primarycalls, reusing pilot8. No prediction
   or desired sign is a condition for executing the full panel. Failure of the
   numerical check is reported alongside preliminary discrimination; it cannot
   be hidden by substituting a selected backend or baseline output.

Nominal initial resources:64CPU/128GiB,8concurrent8-rank native ORCA endpoints,
no GPU, no new DFT, no optimization/training/folding. Initial expected low-level
cost is far below native DFT but must be measured before expansion. No project
runtime/compute-budget stopping rule. Scheduler policy and finite task manifests
remain mandatory. Native SCF iteration safeguards remain algorithmic convergence
criteria, not a wall-time budget. All failed attempts count in reported cost.

## Numerical and predictive checks, fixed before outputs

- Native versus ordinary TightSCF endpoint solvation transfer difference <=0.10
  kcal/mol; Ca−La correction difference <=0.20kcal/mol. These target numerical
  uncertainty below roughly10% of the current2.315model-kcal PQQ gap and below
  the4.56kcal/mol weakest DFT alpha/GGR context margin, not a classification sign.
- Actual normal termination/convergence, finite totals, exactly matching atoms,
  charges/spins, one native parameter set and parsed solvent settings. Charge
  sums must close within5e−4e (printed charge precision); absolute atomic charges
  exceeding4e are flagged as outside this diagnostic's plausible ionic range.
- Solvation identity algebra uses an actual endpoint on both sides; no fabricated
  electronic output. Retain full endpoint terms so no correction is double counted.
- MACE vacuum, GFN2 vacuum, GFN2 ALPB, transfer and composite contrasts remain
  separate. No compatible aquo reference exists; S isnull. All cases are consumed
  development/calibration; no new blind accuracy claim.
- Report25case calibration gap, complete pair ordering, class spread, frozen
  class-extrema calibration on only those25, three transfer calls, and all6
  alpha/GGR structural comparisons. Six comparisons represent one biological
  direction, not six independent observations. No label fitting of the correction.
- Success requires preserving PQQ fidelity and alpha/GGR ordering with improved
  separation/robustness and measured useful throughput. An unfavorable outcome
  stays a result. Neither new raw score nor calibration is automatically promoted.

## Software basis and provenance

Verified ORCA6.1 manual, native GFN2/ALPB sections:
https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/semiempirical.html#native-gfn-xtb-and-gfn2-xtb
Native xTB defaults to its own mixer (overrides ordinary SCF settings),300K
smearing, and implicitALPBwater via keyword. The alternate numerical check
explicitly disables that mixer before requestingTightSCF. Each task starts in
an empty directory; no xtbw/gbw restart. Parameter exports use installed defaults.
Existing9141-atom nativeGFN2 attempts never produced converged energies and are
not reusable endpoint data; this compact test is scientifically distinct.

Root owns preparation/execution/collection and this plan. The second-shell agent
owns pinned source inventory and comparison. Use existing ORCA manifest runners,
preserve other jobs, and store compute products only under workspaces/. Baseline,
water-preparation release, thresholds and inbox/watchers stay unchanged.
