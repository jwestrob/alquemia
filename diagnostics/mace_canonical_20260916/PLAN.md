# Direct canonical PQQ MACE scorer: calibration and retrospective transfer

Declared after the coupled response trial, before any new canonical MACE outputs.
Jacob's active goal and discretionary pilot authorization cover this contained
candidate. Baseline defaults and all historical experiments remain unchanged.

## Question, inputs and scientific scope

Can the already working local MACE+GB calculation provide an affordable
operational PQQ class discriminator across the complete canonical panel?
The earlier medium local calculation orders 4MAE above1H4I correctly; the full
protein contribution reverses it. That narrow local result justifies this test,
but is already consumed development evidence, not blind validation. The coupled
mechanics study validated small Ca response energies but cannot provide paired
La/Ca corrections; it is not added to this scorer.

Use the25 exact released canonical fixed-core calibration preparations and
labels in diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json.
Use exact archived1H4I/4MAE crystal cores and the exact external1KB0 fixed core
as three retrospective transfer cases. They have already been inspected under
other methods;1H4I/4MAE have also been inspected with medium local MACE. Never
call them fresh blind tests. Do not open new external candidates for this pilot.
The panel is PQQ functional-class evidence, not direct La/Ca affinity data.
Preserve protein/sequence groups from the source records; no claim that every
structure or homologue is an independent biological observation. No statistical
confidence interval based on treating the25 structures as independent samples.

Keep original frozen coordinates, atoms, caps, PQQ state, zero-water inventory,
formal charge/multiplicity and source versions byte-exact. No hydrogen repair,
geometry search, cofactor substitution or revised selection after scoring.
Preparation validates the actual archived manifests/outputs, not reconstructed
energies or inputs from prose. Unsupported/missing artifacts remain explicit.

## Candidate, numerical settings and finite inventory

Primary candidate: pinned MACE-POLAR-1-medium with the validated analytic
multipole/memory adapter, float64, vacuum plus frozen-monopole OBC-II reaction
energy. Large checkpoint is a predeclared sensitivity calculation with the same
recipe, reported separately; do not select the better checkpoint after scoring.
Use the existing solvent model exactly (protein dielectric1, solvent78.5,
zero salt/SA, matched radii including identical Ca/La radius1.8Angstrom), frozen
endpoint MACE monopoles. This is direct MACE+GB, not an additive correction to
DFT/CPCM. No second Coulomb energy and no double CPCM solvation. The charge/radius
approximation retains its existing limitations and no self-consistent solvent
feedback is claimed. No response or entropy term is included.

28 preparations x2 endpoints x2 checkpoints =112MACE+112GB logical endpoints.
Reuse the four exact medium1H4I/4MAE archived-core MACE+GB outputs only if input,
charge/spin, checkpoint, backend and solvent settings match their receipts.
Expected new work:52medium+56large MACE calls and the same108GB calls. If a reuse
compatibility check fails, stop preparation and document that specific mismatch
before changing the declared new-call inventory. No new DFT calls. Reuse the
existing finite-manifest runner, isolated environments and allocation policy:
oneA5000,16CPU and64474MiB perGPU job. Expected small-core inference is seconds
per endpoint; model startup adds overhead, so budget scale is minutes per job.
There is no project CPU/time stopping budget; record actual complete costs.

Preserve existing numerical gates, charge closure1e-5e, exact paired coordinates
and scientific cache keys. Those solvers have already passed real rigid/grid
checks; no new unvalidated backend is introduced. Failures remain visible and
cannot be substituted from DFT or the other checkpoint. All source method and
implementation hashes, receipts, component energies and memory/timing are saved.

## Calibration, decisions and frozen success criteria

Primary raw descriptor R_MACE = (E_Ca,vac+GB_Ca)−(E_La,vac+GB_La), in kcal/mol.
Convert MACE eV once; GB is already in kcal/mol. Larger values are La-like.
No sign reversal, old aquo reference, baseline threshold or universal zero.
Retain vacuum and short components for accounting, not alternate winners.

Calibrate each checkpoint on the25 designated calibration rows only:
U = max R over Ca-class rows; L = min R over La-class rows.
If L−U<=0.02kcal/mol (twice the existing0.01 numerical tolerance), no supported
bands are released and the calibration-separation criterion fails. Otherwise
Ca-supported means R<=U, La-supported means R>=L, and the open interval is
inconclusive. This simple rule is fixed now; no optimized threshold or fitting
of weights/features. All25 calibration outputs, invalid counts and denominators
remain reportable, but calibration separation is not an independent accuracy test.

Apply the resulting bands unchanged to all three named transfer cases. The
operational criterion is all three in their expected supported region, with
inconclusive and failed cases counted separately and neither counted as success.
Report the primary medium result regardless of the large sensitivity outcome.
Compare against the saved baseline classifications/raw ordering on the same
exact inputs. Report cost, checkpoint sensitivity and protocol-specific bands.
A successful result is a narrow, inexpensive canonical PQQ candidate; it does
not fix the known alpha/GGR direct-affinity failures, establish incremental
information beyond motif/core composition, or complete broad affinity validation.
Motif and core charge already perfectly confound with the calibration labels.
No production promotion or rescore follows automatically.
