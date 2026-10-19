# Fixed 7CCN whole-peptide reference, version1

## Question and frozen selection

Does whole-system MACE plus the same fixed-charge ALPB describe the energy change and available full forces along an existing, independently selected 7CCN deformation? The geometry pair is exactly the admitted start and 80-step terminal of the whole-MACE 7CCN-La Cartesian search12968. Both La and Dy are scored at those same coordinates. Selection predates these new reference outcomes and does not depend on their sign.

INVENTORY.json freezes four cells, in order before_La, after_La, before_Dy, after_Dy. Common pool12974 supplies geometry_20 and cartesian_03. The La start POINT is50a597c8 and terminalcf59ef7a; their Dy cross-scores are6fc810c6 anded6a2ac. The exact full hashes and all children are in the inventory. No own-Dy terminal substitution, alignment, geometry repair, source rebasing or new optimization.

This is Lu-conditioned7CCN, the same17-residue peptide as consumed7CCO, with a distinct fixed deformation. It is not a blind independent-protein validation. Search12968 was physically admitted and bound but stopped at80 accepted steps, nonconverged (stored maxforce9.739697567 and RMS1.987798073 kcal/mol/A). Preserve the existing source/terminal topology and state records. The fixed physical inventory is255peptide atoms plus target255, charge−2, FIDTNNDGWIEGDELLA, no waters, spectators, caps or external field. Free NH3+/COO− and five deprotonated acidic side chains remain conditional.

## Reference and unchanged native gates

Reuse only low-level functions from exact old whole REFERENCE.py c22c9896: whole state(), asset_checks(), admit_alpb(), JSON canonicalization and the original native parser/runtime/runner assets. Do not invoke its 7CCO inputs(), physical_inputs(), prepare() or scientific_binding(). This adapter has its own exact7CCN inventory and signatures.

Each native input is unchanged PBE0 D4 def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF EnGrad, Guess PModel, vacuum charge−2/effective multiplicity1. Same pinned ORCA6.1.1, La/Dy basis/ECP46/55 and auxiliary files. Physical multiplicities remain La1/Dy6; effective frozen-f singlets are separate. No CPCM, SMD, external charges, changed grid/SCF/MaxIter/temperature or parameter sweep. Preflight derives all4 electron counts and checks parity, atom order, native asset identity and exact input coordinates; no copied7CCO state receipt substitutes for that check.

Admission preserves native normal termination/convergence, original analytic ECP/D4 evidence, charge/electron/core/RHF checks, finite256 gradients, original1e−8Eh energy and1e−6Å coordinate parser tolerances, exact task/manifest/runner/runtime receipt binding and original SCF convergence semantics. Printed SCF residual flags remain reported even when false. No independent whole-system derivative qualification, force-accuracy threshold, physical spin qualification, electron-binding or global ground-state claim.

## Reuse, units and controls

Reuse all4 original whole-MACE vacuum energies and all4 exact native fixed-FF ALPB energy/gradient receipts. Both La endpoints have archived MACE and composite forces; both Dy endpoints are energy-only and their MACE force comparisons remain unavailable. There are no new MACE, ALPB, GFN or MM evaluations.

Replay native model state, charge, exact coordinates, IDs, target255, recipe/cache identity, model weights/build/library pins, all4 native scalar energies and the two available complete gradients. Preserve original producer arithmetic: whole-MACE gradient is `-forces_eV_A * EV + raw_ALPB_gradient * (HA / BOHR)`; original reference solvent admission is `(raw_ALPB_gradient * HA) / BOHR`. These are separate saved representations, not bitwise equated after reassociation. Constants areEV23.06054783061903,HA627.509474,BOHR0.529177210903. Both originate from the same admitted raw solvent receipt. Exact four-cell solvent and source pairing remains required.

Preflight freezes one normalized full256×3 before→after coordinate chord, without alignment or changed coordinates. Its source is the pinned prior search, not a new reference force. Every available reference retains full256 force arrays, per-atom identities, maximum/atom-vector RMS force, this fixed chord projection, and available La MACE-minus-PBE0 force errors. Dy reference force arrays remain available when admitted; missing Dy model forces are not synthesized. Small projections cannot conceal the saved full-atom errors. Projections are descriptive and not integrated as work or interpreted as a Hessian.

## Finite stages and resource limits

All stages enforce SLURM_JOB_ID. No scientific imports, numerical replay or tests on login. Root owns scheduler submission.

* `preflight`: one saved-data allocation, zero molecular calls; create and validate all4 exact input decks and states, archived components, paired coordinates and chord. One CPU is sufficient; no GPU. The bound result is mandatory before any run.
* `run --task before_La`: at most one native invocation/cell, first scout. Config allowlist must be `["before_La"]`; future scout/decision inputs must remain null. Twelve actual MPI tasks, one CPU per task, one node; original helper reads actual AllocTRES and derives memory without changing requested defaults. No GPU, no requeue, scheduler/default3-day limit and native SCF limits unchanged.
* `run --task after_La|before_Dy|after_Dy`: each is a separate singleton allocation after raw scout admission and root's explicit measured-cost decision are pinned. Config allowlist is exactly those three, forbidding a repeat before_La. Each run has one native invocation/cell ceiling, same12MPI/native settings, no automatic retries or adaptive continuation.
* `collect`: saved-only, zero new molecular calls. Re-admit each actual singleton raw result and its preflight/state/config/receipt identity. Missing or failed slots stay explicit; failures never yield a survivor work comparison. All4 are required for the work comparison. The collector may retain admitted force rows from available endpoints.

Four new DFT cells maximum for this declared panel. Root binds the actual first scout before deciding whether to execute the remaining three. Existing consumed7CCO cells each took23058–24919s on12CPUs (total1154004 allocated molecular CPU-seconds); these actual costs motivate a separately measured scout, not a runtime estimate or a claim of cheap execution.

## Required future decision schema

The config initially has `preflight_result`, `scout_result`, `continuation_decision` null and `endpoint_results` empty. Future actual pins must replace these; no placeholder success is accepted.

A continuation decision must have protocol `lbt_whole_7ccn_reference_v1`, status `approved_remaining_three`, exact `scientific_binding` copied from this panel's actual result, exact `scout_result`, `preflight_result`, `scout_native_receipt`, `approved_tasks` = `["after_La","before_Dy","after_Dy"]`, boolean `measured_cost_reviewed:true`, nonempty `cost_review_reason`, and `scout_scheduler_cost_receipt` pin. That factual cost JSON must contain actual `job_id`, `state:"COMPLETED"`, `allocated_CPUs:12`, finite positive `elapsed_seconds`; root retains its underlying scheduler evidence. This decision is not created by the adapter. New molecular stage outputs must be fresh directories; existing endpoint pins disallow execution/retry.

## Comparison and limitations

Compute native differences before unit conversion: each metal's `E_after−E_before` separately for wholePBE0 vacuum, MACE vacuum and unchanged fixedFFALPB. Report each composite sum and the Dy-minus-La difference of deformation work. Report signed MACE-minus-PBE0 work errors and whether their magnitude is within the existing1kcal/mol reporting resolution, per metal and for the differential. This is a discrepancy report, not a native admission, accuracy classification, target preference or fitting criterion. No absolute cross-Hamiltonian total comparison.

A vacuum dianion in a finite basis and frozen-f Dy remain conditional references. Post-added ALPB is not self-consistent electronic solvation. Native convergence is not proof of electron binding or a stable physical state. Nonstationary endpoint force differences and a single deformation do not establish all-basin stiffness or transfer accuracy. Preserve all native residual flags, unavailable Dy model forces, failed attempts, partial outputs and actual per-allocation costs.
