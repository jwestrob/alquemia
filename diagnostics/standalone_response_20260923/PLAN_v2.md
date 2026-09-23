# Six-source bounded solvent-aware accommodation pilot (approved)

2026-09-23. Parent first approved the four-source read-only proposal, then explicitly
expanded the new execution scope before any molecular call: retain Q9Z4J7,
C5AXV8,Q88JH5,A0A3F2YLY8 Ca-sample1 and add A0ACD6B9F2 Ca-sample4/La-sample4.
The original four-source PROPOSED_INPUTS.json remains untouched. Source selection
uses inspected development failures, not a blind evaluation. Jacob's blanket
continuation authorization is recorded in overnight_discrimination CURRENT.md.
Parent: "Approve this contained ... experiment ... Same uniform search/state/settings"
and then "six sources/twelve searches ... MaxsearchGFN960+24cross=984 plus<=16baselinecells,
maxMACE480 ... Confirmprecisecount/reuses/readiness beforelaunch."

## Question and fixed sources

Can the numerically qualified standalone solvent gradient guide useful accommodation
that vacuum-only proposals miss? Static standalone was worse than the native baseline
on the full structural challenge. Reliable derivatives therefore justify this research
question, not replacing the production backend. Test both strongest static standalone
and native-proposal standalone three-candidate pools, and retain native static/adaptive
values for context. Two Ca and four La structural samples represent five protein groups;
the two A0AC samples are the same protein. All source/physical states remain fixed.

For each source and metal choose the minimum own *standalone composite energy* among
origin/adaptive_Ca/adaptive_La, in that tie order. This uniform start rule is frozen now.
Fifty-six existing standalone cells cover the original four complete pools plus two
added origins. Compute the16 missing standalone cells for the added adaptive geometries
before selecting their starts; never substitute native ORCA solvent energies.

## Same physical model; new bounded search

E(q)=MACE_native(vac,q)+GFN2_standalone(ALPB,q)-GFN2_standalone(vac,q).
Project exact Cartesian gradients through the existing source/cap Jacobian:
`gq = J^T[-F_MACE*eV_to_kcal + grad_ALPB - grad_vac]`.
Use only the existing common four angular modes for each source. All other modes,
source atoms, charge, spin, protonation, PQQ, waters and source context are unchanged.
Bounds remain +/-0.8radian and final heavy displacement<=0.8A relative to the sourceq0,
not relative to a shifted start. Intermediate SLSQP constraint violations are recorded;
broken bonds/new severe overlaps or unsupported geometry fail the endpoint.

Existing xTB6.7.1 exact parameter file,accuracy0.02,GFN2,500SCC max,electronic300K,
ALPBwater/gsolv,230surface points,no ionic screening,no restart; native OMOL checkpoint
unchanged. Both low-level media come from standalone xTB. Reuse verified actual MACE
start forces/energies, but request both standalone analytic gradients at every start.
Check start gradient-evaluation energies reproduce old standalone cells within0.1kcal
per medium. No numerical DFT gradients or approximation to the composite derivative.

Twelve independent SLSQP searches,maximum20iterations and40distinct evaluated points
per endpoint, including the start. Shift the objective by own starting E and divide
energy/gradient by627.509474, so ftol=0.001/627.509474 Hartree-equivalent. This is50times
the observed standalone accuracy-pair energy change and100times smaller than the
existing0.1kcal origin-selection scale; no stationary-minimum claim follows.

## Admission and failure rules, before output

On successful termination or the declared iteration/evaluation limit, offer the lowest
own energy among actually evaluated, complete, physically admissible points; tie goes
to earliest evaluation. The start participates, allowing a legitimate zero-motion
result. No label or Ca-minus-La contrast enters this selection. This is a finite best-
feasible proposal algorithm, not convergence evidence or a thermal distribution.

A molecular/native/SCF/gradient/geometry failure or other optimizer failure makes that
endpoint unavailable, preserving the actual attempts. No retry, fallback to an earlier
successful trajectory point after molecular failure, alternate start or radius change.
Search-limit termination is reported separately from convergence and scientific failure.

Offer BOTH final proposals to BOTH metal rows alongside all original three candidates.
Only exactly identical coordinates may deduplicate; retain aliases. Own final energy/
gradient evaluations supply their own cells; opposite-metal cells require new MACE and
both standalone media. Any required failed cell makes the extended pool unavailable;
old pools remain separately reportable. Apply the existing mathematical minimum and
operational0.1kcal source-origin retention rules to the same common pool.

Primary utility: raw known-class ordering/gap and changes versus static and prior3pool,
all six rows and five biological groups, no training. Static bands are transferred
checks only for changed candidate methods. No six-case calibration or inherited
compatible reference claim. No Hessian/entropy, new modes, DFT, protein folding or
proton/water changes; no automatic broader rescore or promotion.

## Finite execution and cost

At most16 missing baseline cells +12*40*2search cells +24final cross cells =1000new
standalone calls. At most12*39new search MACE +12cross =480newMACE, assuming all12
start-force reuses validate. Count all failures and unused line-search evaluations.
The40-point cap defines the scientific proposal algorithm, not a project time budget.

One32CPU/200000MiB/oneH200 allocation. Two concurrent endpoint searches, each with
vacuum/ALPB eight-thread calls in parallel; one serialized existing warmGPU worker.
This occupies at most32CPU solver threads, without MPI or oversubscription. Existing
standalone calls took0.85–7.02s under eight-worker qualification; modest hundreds of
calls should take minutes, but report actual costs, not a matched-speed claim.
Manifest, actual reuses and preflight are frozen before submission. No job is authorized
outside this explicit scope by this document. Retain final receipts/results and vault.
