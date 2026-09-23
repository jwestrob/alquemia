# Approved restart and collective-pocket experiment

## Agreement

After the completed parallel pilots, root recommended a targeted native solver
restart test alongside collective pocket motion. Jacob replied:
**“i authorize the next round of experiments you suggest.”**

This authorizes the two contained experiments below, including their implementation,
routine numerical fixes, actual execution and comparison. It supersedes the
proposed-only status of these two designs, not their archived records. Standalone
xTB, redox, new DFT, wider rescoring and production promotion remain outside this
round. Production, historical outputs and unrelated working-tree changes stay intact.

## A. Native electronic initialization

Question: does the Q88JH5 energy discontinuity depend on the electronic initial guess?
Use the exact two archived La geometries and their actual saved native charge/
multipole solutions from the preceding eight-repeat diagnostic.

Execute the eight calls defined in
`diagnostics/structure_informed_starts_20260922/RESTART_PROPOSED.md`: old/new
destination geometry × vacuum/ALPB × same/opposite geometry seed. Native GFN2,
300K, charge−2, singlet, MaxIter500 and8MPI ranks remain unchanged; remove
NoAutostart consistently and copy only the pinned native.xtbw seed to the actual
runtime basename. Keep self controls, seed-before/after hashes, full iteration
traces, charges, energies and positive native-restart evidence. No guessed
tightening, ordinary-SCF substitution or favorable-solution selection.

One64CPU/128GiB allocation with eight8-rank concurrent calls; no GPU requested.
Submit directly to the supported gpu CPU-sharing route rather than inherit the
standard partition's exclusive-node setting. Count all actual attempts. A restart
that cannot be confirmed remains unavailable. This diagnoses initialization; it
does not automatically recalibrate or repair the scorer.

Khoury owns `native_xtb_restart.py` and its versioned diagnostic/workspace.

## B. Collective response to matched donor targets

Question: does allowing neighboring protein atoms to respond improve the existing
metal-accommodation result, beyond ordinary protein force-field relaxation?
Use the exact common eight consumed sources and archived adaptive pools pinned
in `diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json`.

The agreed design is
`diagnostics/collective_scaffold_20260922/PROPOSED.md`: three starts/targets per
source (original, Ca-adaptive, La-adaptive), with donor coordinates fixed to that
target. Move complete protein residues having a heavy atom within8Å of the
canonical-core heavy atoms, core-contributing residues and their actual bonded
peptide neighbors. Choose this set from original coordinates before energies.
All other protein atoms, metal, PQQ, waters and chemical states remain fixed.
Preserve actual source bond lengths, H inventory, chirality, peptide cis/trans,
cap mapping and0.8Å maximum heavy displacement. Reuse exact previously supported
parent parameterizations; unsupported parents remain explicit, without fresh
terminal/protonation repairs or invented metal/cofactor parameters.

Protein-only ff19SB supplies proposal forces. It is not added to the final score.
Use sparse constrained projected L-BFGS with Newton/LSMR bond retraction, an
Armijo line search, at most200 accepted steps and maximum0.05Å atom step. The
mechanics owner freezes all numerical details in EXECUTION_PLAN before forces.
Required acceptance: bond residual≤1e−6Å (retraction target1e−8), fixed-coordinate
error≤1e−12Å, heavy displacement≤0.8+1e−7Å, unchanged source chirality/cis-trans,
and no new severe overlaps (1.0Å heavy/0.55Å involving H), including fixed
cofactor/metal contacts. Initial feasibility precedes any force call.

Each last accepted physically feasible endpoint may be scored if parent energy
is no higher than its own starting energy+1e−6kcal/mol. Preserve iteration/bound/
line-search status; these finite proposals are not called stationary minima unless
projected maximum atom gradient≤0.1kcal/mol/Å and Cartesian RMS≤0.03kcal/mol/Å.
Zero-progress cases remain explicit and exact-coordinate duplicates are reused.
All three declared targets are required for a source's primary comparison.

Maximum24constrained searches,48new compact MACE evaluations and96new native
GFN2 singlepoints before exact reuse. Each metal scores the same admitted pool.
Keep the existing source-specific context membership, native float64 OMOL and
native GFN2 ALPB-minus-vacuum expression. Do not mix restart results into this
branch's unchanged numerical recipe or fit a new threshold on these eight cases.

Report matched-target responses and original-target controls separately, then
expanded-pool selection, raw contrasts, explicitly transferred bands, source-pair
spread, failures and incremental cost. Unknown PLM proteins are not targets.
Lower protein FF energy alone is not classifier benefit. Omitted protein–PQQ/metal
forces and absent matched scaffold-energy subtraction remain explicit limitations.

Water_basins owns exact parent preparation/source maps; second_shell owns mechanics;
root owns physical candidate admission, existing MACE/GFN2 runners and comparison.
Only one owner runs each molecular task. Record platform/precision and actual FF
calls in the first declared search; no separate performance campaign is needed.

## Delivery and boundaries

Candidate products belong under workspaces; compact plans, receipts and summaries
under diagnostics. Use finite manifests and current allocations, with no arbitrary
project-total compute/time budget. Capture actual allocated CPU/GPU cost and
source/collection overhead separately. Every agent writes a vault note; root
coordinates the shared index for scoped commits. No push or default change.

An electronic numerical failure and a negative scaffold result are distinct
outcomes. The restart test can identify a numerical remedy worth qualifying;
the scaffold test must demonstrate improved discrimination or structural robustness
to justify further classifier work. Neither is presumed to win.
