# Active MACE discriminator checkpoint — 2026-09-17

**Goal active; no broadly validated MACE replacement.** Jacob authorized autonomous
analyses and available resources. Home/project AGENTS now remove the old
per-analysis permission gate. Preserve baseline/default, immutable studies,
other agents' edits/jobs. No push/promotion. No project CPU/time/token budget.
Recover GOAL.md and inspect live jobs before resuming.

## ACTIVE: charge sampling complete; isolated ddX build ready

The sampling experiment completed 40 native NoIter high-level property replays
and 40 actual saved-density potential queries. No SCF optimization or scientific
retry. All energy, default-charge, density, fit and projection gates pass.
All eight finer-to-finest sensitivity checks pass, but five of eight fitting-
extent checks fail the frozen 0.5 kcal screen: overall sensitivity gate FAIL.
GGR 2FW0 minus 2FVY source-self differences (kcal): default 18.12902095,
finer 16.38528087, finest 16.15398128, extent 14.96728109. Sampling contributes
up to 3.16174 kcal sensitivity; the large discrepancy persists. No setting was
selected per case, no full hybrid rescore, calibration or baseline change.

Reports: CHELPG_SAMPLING_STABILITY_REPORT/RESULT and QUALIFICATION_REPORT in
mace_omol_20260917. Full output: chelpg_sampling_report_v1/result.json.
Population: chelpg_sampling_panel_v1/population_job_1201177.json.
Density queries: chelpg_sampling_panel_identity_v1/collection_job_1201211.json.
Reporter snapshot: chelpg_sampling_reporting_source_v1. Read-only collection
chelpg_sampling_panel_recollection_v1 reproduces all rows without new calls.

Jobs 1201173/74/77/1211 used 1030 summed job-wall seconds, 65920 allocated
core-seconds, 53718.534 reported CPU seconds, zero GPU. The largest job reports
CPU only to whole seconds. Local report: 18.373 s. All 11 actual-fixture tests
pass in 9.732 s, no final skip. Maximum density identity error 1.03001e-13 au;
energy error 3.23739e-8 Hartree. The explicit NoIter contract retains raw generic
SCF failure flags: 1201173 Slurm FAILED, but native property calls terminated
normally. The shared scorer/runner/default parser remain unchanged.

### Next question and available backend

Does a resolved PCM protein boundary show the same source-self discrepancy as
GK on identical physical cavities and fixed default charges? Declare the exact
physical setup, convergence checks and solve inventory before the next pilot.
No ddX scientific calculation has run yet. This question tests the solvent
approximation; agreement with GK or a desired biological direction is not a gate.

DDX_CAPABILITY_NOTE records the source audit. Exact-density coupling requires
both surface potential phi and a density integral psi; exact phi plus fitted-
charge psi would mix representations. Upstream Psi4 uses Becke quadrature.
An ORCA density bridge remains unvalidated. Projected-charge source-only tests
can examine the GK approximation first, with their limitations explicit.

The isolated ddX 0.9.0 build/import passed in job 1201203:
`ddx_software_v2/attempt_0001/receipt.json`. Interpreter:
`ddx_software_v2/attempt_0001/venv/bin/python` (under the common workspace root).
Python 3.11.15, NumPy 1.26.4, SciPy 1.17.1, pybind11 3.0.1,
setuptools 82.0.1, wheel 0.45.1, packaging 25.0, CMake 3.31.6.
No existing environment changed. Build: 89 wall seconds, 5696 core-seconds,
76.595 CPU seconds. Initial job 1201202 failed the CMake executable hash check
before installing or compiling (1 wall second, 64 core-seconds, 0.474 CPU s).
The new build uses a pinned CMake wheel; both attempts remain intact.
Plans: DDX_BUILD_PLAN/RECOVERY. Source archive in ddx_capability_source_v1,
SHA256 1f8faf41d61483f063d29934ea969747ea5e060442eed188adbf629e2a945561.

No study/build job remains live. Inspect the queue before continuing and
preserve other agents' work, PLM 1200796 and pending H200 1200809.
No project compute/time budget or per-pilot permission gate.

### Preceding findings remain frozen (commit 2bbc8a1)

The five-site full hybrid passed 39 numerical checks but failed all six supporting
parvalbumin-minus-GGR comparisons. GGR localization passed 1072 checks.
Native source-only/empty job 1201172 matched analytical source-self energies
within 6.324e-13 kcal. Symmetric attribution of the 18.129021 kcal shift:
charges +18.452793, distances -0.372800, Born radii +0.049028. The common 30 A
Born bound is verified in native/OpenMM code; that alone does not establish a
bug. No broadly validated hybrid predictor or absolute reference exists yet.
The goal remains active. Production baseline is unchanged.

## Latest: frozen hybrid fails GGR structural transfer; multisite preparation works

Job1201162 completed27static/24field/30response calls. All18numerical/identity/
rigid/radius checks pass, but both2FW0 comparisons fail. Alpha1F6S minus2FW0
=-14.764159019kcal; alpha6IP9=-15.638938500. Both2FVY comparisons pass+8.927577642/
+8.052798162. All six crystal comparisons retain1GLG,2FW0,2FVY; GGRrange24.891221454.
No old result overwritten, threshold change, case selection or baseline change.
These remain2consumed groups; no absolute reference or broad validation.

2FW0-minus2FVY=23.691736661kcal: intrinsicQM+0.242251473,MACEshort+0.634064156,
directdensity+13.421787726,GK+12.600209216,induction-3.206575910. Environment accounts
for22.815421032. This identifies dominant terms, not a unique cause. Trialdensity
responds toff19SB, evaluated withAMOEBA/GK; no currentfunctionalselfconsistency.
Projected CHELPG remainsGKproxy. Do not rescue by selectingfavorablecrystal.
Report/actualpins/costs: TRIAL_DENSITY_GK_EXPANSION_REPORT.md / RESULT.json.
Full products trial_gk_expansion_native_v1/collection_job_1201162.json and
trial_gk_expansion_comparison_v1/result.json. New panel interface/scripts:
mace_density_panel.py and mace_density_panel_report.py; existing nativeexecutor/
collector reused. Scientificprotocol unchanged5f7d164; runnerprotocol
'declared_source_graph_responsive_density_GK_POLAR_panel_v1'.

Fresh4DFT1201154,8MACE1201158,4CHELPG+4queries1201161 allcomplete. Source/reference
8initializations pass. Earlier preflights1201149/1201152 failedbeforeDFT, preserved;
qualifiedparserrecollection solvedconditionalORCAwarning withoutDFTrerun.
Recorded allocationtotal45208core-s/17650.217382CPU-s/112GPU-s,1261s jobwall sum
(includesdevelopmentcontrols/failures; notproductionlatency). Native409wall,
26176core-s,2812.966CPU; kernels100.439151s. Newpanel7distincttests allpass across
preflight6pass102.894s and integration1pass15.474s; oldhybrid6pass33.062s. Exact
legacyreplay scores/decisions; fourrigid dipole contractions differ<=3.553e-15,
test allows1e-12onlyonthese intermediates, remainingfields exact. No scientific
rescue. Expansionnojoblive. Email20260918T040505Z acceptedrelay withthisfinding.

Independent next-family preparation declared MULTISITE_AMOEBA_PREPARATION_PLAN
and MULTISITE_NATIVE_PARAMETERS_PLAN. ALL FIVE real frameworks nowpass:
PARV4CPVCD/EF(1611atoms,onebackgroundCa,onewater,frameworkQ-5) and AEQ1SL8
EF1/EF3/EF4(2866atoms,twobackgroundCa,threewaters,Q-6). Source ACE0-ALA1bond,
allheavy/Hcoordinates/protonation/ions/waters retained. Noenergy/DFT/MACE/response.
Artifacts multisite_amoeba_capability_v3/result.json, config_v2/config.json.
v1failedbeforetyping becauseoldselectedmetalalias; v2matchedtemplatesbutOpenMM
GKdefaultBondiradiustablehasnoCa. Bothfailurespreserved. v3usesOpenMMpermanent/
polarizationparameterswithGKexplicitlyunavailablethere, thenactualnativeTinker
SOLUTEpreflight. Samequalifiedbackend/AMOEBA2018, noinstalledenvchange/radiusguess.
NativeCa(type358)Q2,alpha.55A3,radius1.82485A; maxpermanentparamdiff2.22e-16.
Scopedoptionalflagsinmace_amoeba_capability.py; nativevalidatorrecognizesCa.
Newscriptmace_multisite_amoeba.py. Legacydefaultsremainold. New4tests pass3.216s; oldAMOEBA/nativeframework11regressions pass11.094s. Nativeprep83.401715
wall/79.253806parentCPU +3.239334nativechildCPU; noGPU/cluster. Report
MULTISITE_NATIVE_PARAMETERS_REPORT.md. No fullhybrid scores on these sites yet.

NEXT: scopedcommitownfiles/SESSIONS. Implement
source-backed quantum-core bridge for multisite cases without droppingbackground
Ca/ACE/waters; selectedmetalrealIDs currentlydifferfromsinglemetalalias, account
explicitly. Existingrepairedgeneric sitecores undernonpqq_direct_site_benchmark;
use same declared peptide policy, supportedwaterHmapping. Needexplicitstate
accounting beforedeclaring quantum/scoringmanifest. Preserveparvalbumincrossstudy
and aequorin orderedvector(noinventedlabels). Think about a specific environmental
accuracy test fromthe2FW0 failure; don't justretuneorrepeatallfailedmodels.
Goalactive, baseline/defaultunchanged. Noexpansionjoblive. Checkactualqueue;
PLM1200796,H2001200809,otheragentMopB1201160untouched. No projectcompute/timebudget.

## Latest: responsive trial-density candidate passes development gates

Jobs 1201135 and 1201137 completed. All four alpha-minus-GGR contrasts pass:
10.127062/10.668296 kcal/mol for 1F6S, 9.252283/9.793516 for 6IP9
(extended/connected GGR). GGR connected-minus-extended = -0.541233142,
passing unchanged absolute <=2 criterion. All 32 numerical/identity/rigid/
radius checks pass. These are two consumed biological groups, not blind or
broad validation. No absolute reference/bands or combined gradient exists.
Baseline remains unchanged; active goal is not complete.

Protocol: saved_responsive_trial_density_AMOEBA2018_GK_proxy_POLAR_short_hybrid_v1.
Same full physical geometry/cavity/permanent environment/MACE as the frozen
parent. All density-dependent terms replaced together using eight preassigned
saved responsive wavefunctions; their old generating-field energy is subtracted
once. Trial density is not self-consistent with the current AMOEBA/GK model.
No labels, waters, thresholds, geometry, radius or per-case density selection
changed to obtain this result. All 16 boundary initializations/12 checks pass.

Actual outputs: trial_density_gk_fields_v1/report_job_1201135/result.json,
trial_density_gk_hybrid_v1/collection_job_1201137.json, and
trial_density_gk_comparison_v1.json under workspaces/mace_omol_20260917.
Frozen implementation/config: trial_density_gk_pipeline_v1. Reports and exact
commands: diagnostics/mace_omol_20260917/TRIAL_DENSITY_GK_{REPORT,RESULT,COMMANDS}.
Twelve distinct new tests plus 13 old regressions pass across focused runs;
no final scientific integration test remains skipped. No new DFT/MACE/fit.
Eight actual density queries; native 34 energies/32 field queries/40 responses,
37 verified environment-only reuses. No failed scientific attempt/retry.
Total new jobs: 42,256 allocated core-s, 4,846.341587 actual CPU-s, zero GPU.
Fresh matched production-pair cost still unmeasured; archived costs additional.

Matched component audit: small differential intrinsic response energies mask
larger reaction-field and induction changes. MACE short term is unchanged;
removing it algebraically leaves GGR partition -5.656808 versus full -0.541233.
This is a development component audit, not independent validation or proof
that MACE alone explains the ordering improvement.

NEXT: freeze this passing candidate and prepare wider real-case evaluation.
Read existing evidence curation and compatible physical preparations first.
Preserve biological groups and evidence strata; no fitted test-set threshold.
No new trial-density job remains live. Check live queue before using resources;
PLM1200796, H2001200809 and other agents' MopB jobs remain untouched.

## Completed: full density/GK/MACE candidate; accuracy gate fails

Job1201074 completed all51static energies,48field queries and60response solves,
no retries. All32numerical/identity/rigid/radius checks pass. Ordering2/4fails
all-case gate; GGRconnected-minus-extended=-3.783707453869283kcal failsabs<=2.
Alpha1F6S minusGGRextended/connected=-0.4377995353/+3.3459079186;
alpha6IP9=-3.5462673929/+0.2374400610. Do not select favorable core/structure.
These remain two consumed biological groups, not four independent observations.
Maximum convergence correction2.264684e-7kcal; rigid2.842171e-12;
+/-5%commonmetal-radius relative effect<=0.000197334; vacuumidentity0.
Primarymaxatomicmu0.208933--0.336394eA, no large-response flag.

Native solver cost570wall x64CPUs=36480allocatedcore-s,3659.722079actualCPU-s,
zeroGPU; nativekernels124.133763s, children210.915941wall/3316.199736CPU.
Preparation41.400926wall/38.763210CPU; build1.775239wall/1.713074CPU.
All6distinct hybrid tests nowpass (4before,5aftercalls,remaining1aftercollection).
Fullsource/assumptions/results/operations in DENSITY_GK_HYBRID_REPORT.md,
DENSITY_GK_HYBRID_RESULT.json,DENSITY_GK_COMMANDS.md. Actual products under
workspaces/mace_omol_20260917/density_gk_hybrid_v1; frozen software_v1.
No newDFT/MACE in this pilot; reused upstream costs additional. Numerical
implementation credible, partition/predictive gatefail, nativecomponentcheap,
complete matchedproduction costunmeasured. Baseline/default unchanged.

Prerequisites complete: density_multipole_coupling_v1 job1201063 ran8saved-density
orca_vpot utilities, allnumericchecks pass.37points/site; maxQrefinement
0.000206909396kcal; no nuclearHessian/newDFT.3native moment exports; actual
roundedQtensortrace retained, not zeroed. Five tests pass10.237s. Cost846wall
x8=6768core-s,2951.569CPU-s,zeroGPU. DENSITY_MULTIPOLE_COUPLING_REPORT/RESULT.
density_gk_boundary_v1 ran16native initialization-only preparations,12paired
reference/cavity checks pass. Realmetal retained in nativeindex; fullphysical
atoms4698/1932/1898, correct AMOEBAresiduecharge ledger/source mask, common
Ca2018GKmetalradius1.82485A (not LaFF). Three tests pass39.789s; nativewall
16.248823s/CPU16.064919s. DENSITY_GK_BOUNDARY_REPORT/RESULT. Allinputsimmutable.

Historical next step, now completed: TRIAL_DENSITY_GK_PLAN.md was declared before
new responsive-density query or
fullmodelscore. Reuse eight1200983protein-field responsive quantum densities in
exactly the same current functional, removing their old generating-field energy
once and recomputing every density-dependent term. No newDFT/MACE/fit needed.
Explicit trialdensity, not currentAMOEBA/GKselfconsistency; no endpointdensity
selection, changedthreshold or oldtotalGBreuse. Declared8utilities49points/site,
16nativeinitializations, same51static/48field/60response logicalinventory;
compatible unchanged environment outputs may be reused with exactreceipts.
Newprotocol saved_responsive_trial_density_AMOEBA2018_GK_proxy_POLAR_short_hybrid_v1.

Source audit implemented scripts/mace_trial_density.py audit and completed:
trial_density_gk_source_audit_v1.json. All8actualwavefunctions available, exact
geometry/microstate/source mapping; old-field subtraction reproduces intrinsic
core energies within1.165e-10kcal. Actual intrinsic polarizationcost37.355376--
41.323737kcal, allnonnegative. Audit15.358869wall/14.253791CPU, no newsciencecalls.
Initial audit wrongly treated validator status dict as manifest; fixedAPIread,
no scientific rerun. The49point preparation, responsiveboundary and dependentexecution are now
implemented as detailed in the live section above. Do not restartvacuumjob.

Goal remainsactive; preserve all unrelated dirtyfiles. Reports/current/vault/
SESSIONS recorded in scopedcommits. PLM1200796andH2001200809untouched;
rechecklivejobs. No per-analysis approval gate, projectbudget or promotion.

## Previous: native energy identity and supplied-field response both pass

Goal remains active. Two completed component experiments; no new metal score.
Job1201055:12native direct-field queries+6no-response energies, combined with
12old mutual controls. All15checks pass. Maximum native energy-contraction
residual3.183231456205249e−11kcal, rigid field error3.986394547794703e−15e/Å².
Exact formula I_GK=−.5*(electric/dielec)*sum(mu_solv,d*F_solv,p) reproduces
native ep+es−es_static. Keep d/p fields distinct; source-frozen arrays arezero.
Six rotated outputs were initially rejected only by regenerated-coordinate
validation. Logged parameters/coordinates exactly match frozen parent and XYZ;
reparse recovered all withzero native reruns. Local regeneration agrees; no
numerical size/hardware cause of the compute-node replay difference is claimed.
Use collection_recovery_pinned_v1.json and its frozen recovery implementation,
not the initial incomplete collection. Five actual tests pass19.316s.
Cost38wall×64=2432allocatedcore-s,106.172actualjobCPU-s,0GPU; nativekernels6.304s.
Manifestc86f03720e7f9aa6eca9f775fc96253b8184cb6a7f56ace011281d069cc586e9.

Job1201061: six supplied-field native mutual/GK replays, all9checks pass. Actual
native fields are inputs, not new quantum fields. All16iterations at1e−7Debye;
max dipole replay1.8214596497756474e−16eÅ, energy3.092281986027956e−11kcal.
Rigid energymax6.99565134709701e−8kcal, vector4.742733983320591e−15eÅ.
Only nativeinduce0c name/input declaration/dfieldassignment changed; exact inverse
source test proves all solver/kernel/convergence lines unchanged. Pinned library
and original source untouched. Four tests pass7.062s. Cost25wall×64=1600core-s,
199.544actualjobCPU-s,0GPU; kernels4.325634061s. BuildV1failedonTinker'smaxval
symbol shadowing intrinsic; V2uses equivalentanycomparison, samecriteria. Both
builds retained; failedprepare rejected unavailable frontend before science.
Manifest60e6da85420d4942093f4f66fc1888877bb53bf2459aaeede33cdefb4cd32682.

Productsunderworkspaces/mace_omol_20260917/native_field_{accounting,input}_v1;
accounting software_v1, input software_v2. Reports, exactcosts, runbooks and
plans NATIVE_FIELD_ACCOUNTING_* / NATIVE_FIELD_INPUT_* in matchingdiagnostics.
No newDFT,MACE,geometry,training or biological classification in eitherexperiment.
Old failures stay; baseline/default unchanged. No nativecontrol jobstilllive.
PLM1200796running and H2001200809pending are unrelated/untouched; recheckqueue.

NEXT DENSITY_MULTIPOLE_COUPLING_PLAN.md declared before outputs:3native
initialization/rotation-only moment exports +8saved-density orca_vpot queries,
37observation points per real environment atom to check spatial potential
Hessians (NOT nuclear DFT Hessians). Exactsource,steps,tolerances,resourcesinplan.
No code/preparation/submission yet for this next experiment. Quadrupolecoupling
must use actual internal native Qtensor (kmpole already dividesby3), fulltensor
contraction, not guessed XML conventions. Use unchanged exterior moments only
ascomponentdiagnostic until charge-closed source/boundary preparation exists.

Plan also states a concrete candidate expression:vacuumDFT +exactdensitydirect
coupling +GKpermanentproxytransfer +environmentinductiondifference +MACEshort
full−core. Exactdensity drivesinduction; projectedCHELPG isonlyaGKproxy, not a
rescued fieldfit. Fullphysicalmetal/cavity and AMOEBAresiduechargeledger still
need implementation/explicitparameters before any fullmetalmanifest. NoCPCM,
missingcomponentzero, oldbands or unique non-electrostaticMACEclaim. Continue
this concrete route toward a real accuracy test, not another literature survey.

## Latest: distributed source trials closed; direct-field accounting next

Both declared fits completed on the eight real normalized vacuum endpoints.
Job1201022 (V1) failed6/8endpoint field screens; job1201026 (two-center V2)
failed4/8, all alpha, and the alpha6IP9 paired bare-U0 flag (−1.231046kcal).
All potential, numerical and differential-vector screens pass. These are
spatial density-representation tests, not biological predictions; U0 is not
an environmental correction. V1/V2 validation locations differ. Neither
source is qualified; stop this bounded fit trial, no parameter/order sweep.
See DISTRIBUTED_SOURCE_{REPORT,RESULT,COMMANDS} and
DISTRIBUTED_SOURCE_SAMPLING_{PLAN,REPORT,RESULT} in diagnostics/mace_omol_20260917.

Actual16fits+16native saved-density potential utilities, zero newSCF/MACE/FF
energies. Combined3904allocatedcore-s,1623.268actualjobCPU-s,zeroGPU. Eleven
real parser/algebra/integration tests pass11.623s. InitialV2preparation failed
before scientific calls due to shadowed function; retained failure and recovery
with identical scientific settings. Its timing is unavailable, not zero.
Workspaces distributed_source_v1 and distributed_source_v2_recovery_v1 under
workspaces/mace_omol_20260917 retain frozen code, native outputs and receipts.
Manifest hashes98820483e24705f632a993400f01aaabbd55ffe6551005c8a7db998040e67282
andf07f44bd2c41d3887a8ebeb6894a1eed9d438dc216b64cf73e55a35552a1108a.

NEXT NATIVE_FIELD_ACCOUNTING_PLAN.md:12native direct-field queries and6actual
no-response energies on the already consumed ion-excluded framework controls.
Reuse12old mutual energies/dipoles. Test vacuum/GK contractions against actual
native component differences at frozen1e−7kcal tolerance before per-site quantum
field injection. No new fitted source, DFT, MACE, gradients or hybrid score.
This is an energy-accounting prerequisite, not full model validation. Plan
exists; implementation/preparation/execution not yet begun at this checkpoint.
Goal active; baseline/default unchanged. Live1200796PLM running and1200809H200
pending untouched. Per-analysis permission rule already durably removed.

## Latest: fast native polarization; saved-density fields expose a source weakness

Goal active. Commit72481f1 implements/tests the native solver and field runner.
Tinker1201015 completed12framework energies: all30checks pass, frozen induced
variables exactlyzero, standard11/tight16iterations. Refinement max.0013258445
kcal, rigid max6.99547e-8kcal. Kernel17.122408376wall-s. Allocation37s×64=
2368core-s, actual467.119CPU-s, zeroGPU. Five tests pass8.005s. This is an
ion-excluded framework control, NOT a La/Ca score or full-system validation.
See TINKER_FRAMEWORK_SOLVER_REPORT/RESULT/COMMANDS and workspace
workspaces/mace_omol_20260917/tinker_framework_solver_v1 (manifestSHA
 ee7cadb5d9be64a70b0faf920f4f48213d4a9830982431e686afc3962347e2ed).
Native GK bulk78.3hardcoded, internal1; esincludesnonpolar, notpureRF.

Electric-field job1201017 COMPLETE:8orca_vpot calls on saved normalized VACUUM
r2SCAN3c densities, zero newSCF/fit/MACE/FF energy. Alloutside-QM-support atoms
probed (GGRext4640/conn4587,alpha1880/1843),12offsetpoints/site for two central
spacings. These are observation-point derivatives, NOT nuclear gradients.
Numerical pass: max1.45403e-8au vector change, pairedU0 change1.13235e-6kcal.
Four actual tests pass5.788s (previously explicitintegration skip).
All8endpoint field screens FAIL: fittedweightedrelative errors21.4–33.5%,
projected27.9–40.7%. All4Ca-La vector screens pass3.28–5.83%, but quadratic
response also uses common endpoint fields. Paired projected−exact bare U0
errors: GGRext−1.792883744,conn−.852195257,alpha1F6S−3.141361990,
alpha6IP9−2.728425113kcal. 3/4fail1kcal flag. >=3A stratumfields2–8%, U0errors
<.088kcal; this is predeclared diagnostic, NOT rescue/omissionofnear sites.
U0=-.5sum(alpha E²) is diagonal undamped diagnostic, NOT environmental energy.
It omits native covalent scaling/damping/mutual/GK/protein permanent fields.
No unique claim that these errors cause prior scorer failures.

Job203s×8=1624core-s, actual711.484CPU-s,0GPU; utility summedwall702.264582902s,
699.04CPU-s,maxindividualRSS318076KiB. Preservedreport qm_electric_field_v1/
report_job_1201017/result.json, cost.json. ManifestSHA
0b639b8502fdade3b8dd7727563908cbdf1d2237ae4008b8cfe329ce872d0673.
Compact QM_ELECTRIC_FIELD_REPORT/RESULT/COMMANDS underdiagnostics/mace_omol_20260917.

NEXT DISTRIBUTED_SOURCE_PLAN.md declared BEFOREfits/validation: ONE uniform
charge+dipole density representation on existing physical support;8fits using
savedactual centerpotential+fields, chargeconstrained, fixedSVDregularization.
8newnativevpot calls at separate fixedtranslated spatial probes for validation;
no newSCF/MACE/FF energies. Exactfitformula,steps,gates,resourcesinplan. No code,
preparation, fits, validationoutputs or submissions for this next step yet.
Reuse/pin original explicit_field_short_v1 potentials only after verifying
wavefunctions/physicalIDs/order with fieldsource; both are normalizedVACUUM.
Do not refit/tune on heldout probe outputs or biological labels. Historical
MBIS monopoles already tried; native extraction575–748s/endpoint and some
worsepotential/coupling errors. Multipole-capabilitydocs are not efficacyproof.
Full source damping, covalent boundary charges, subtraction and predictive
partition checks remain unresolved. Do not launch a complete hybrid scoreyet.
Baseline/concurrentPLM1200794/95/96 and pendingH2001200809 untouched; recheckjobs.
No approvalneeded within goal; no production promotion/push. Continueautonomously.

## Latest: native Tinker backend built; three parameter round-trips complete

Current goal turn made progress: pinned CPU/GPU source, retrieved both primary
La supplements, built CPU library/parameter adapter and executed3real native
framework initializations. NO scientific energy/force/DFT/MACE/GPU calls yet.
Commit for preceding AMOEBA milestone is d76b93e; this section supersedes it.

CPU Tinker26.2 commit87050685eff8840d312e2a332cc82c33f63c7c3d; source pristine.
Build1201010 failed because upstream CMake omitted existing uatom.f; one
deferred target_sources include repairs build integration, not scientific
source. Recovery1201011 completes, same64CPU node. Total2560allocated core-s,
285.312actualCPU-s,zeroGPU. Isolated workspace tinker_software_v1 has manifest_v2,
driver, compat include, failed/successful receipts and executable SHA
aee994b1b15e5bf5d093d80c3106f6b767d72d68f46115caafdf46700356c738.

New scripts mace_tinker_capability.py, mace_tinker_probe.f90 and4passing tests
(2.501s). Exactly3framework exports/native parses: GGR4697/alpha1931/1897 atoms;
coords,bonds,charges,axes exact; local multipole difference max1.071e-10 MD units,
polarizability/damping roundoff2.22e-16. Explicit metal remains unparameterized.
Source masks57/51/54 from existing projection support. Native POLARIZABLE is
actually overwritten later in kpolar: all57/51/54 became allowed again. Adapter
restores douind after mechanic and preserves damping/polarity. Subsequent SCF
mask behavior is NOT YET tested. No claim of a working complete hybrid.

Workspace tinker_capability_v1 manifest SHA
f992467d4cb9ab69b24112df268fe93b38a82123a80238f5b722305c7e4aedaf.
Native logs/results and mappings all saved. Export4.314806wall/4.063931CPU-s;
native parses3.220401wall/2.049645CPU-s. Reports TINKER_CAPABILITY_{REPORT,RESULT,
COMMANDS}; source/parameter acquisition pinned under tinker_sources_v1.
Supplement s001 TXT SHAeaa7da82c749c482d2c0b5110aa3cf1f1098e59600b016ed8e155471a32e06cd;
s002 PDF SHA5f112c7aa4976e1c365d77711678b762ab672c368e2cf6e6b0b6ae8e118f2100.
Actual La09 SI amide parameters differ from paperTable1; keep distinct. SI
TableS3 reports~49kcal monodentate-acetate error, so no universal accuracyclaim.

GPU now tinker-gpu commit44bcd7c898e68f6827c112df79b8361e443730d2, pins separate
Fortran33dcaf084a3e85730559d4e3a3188d9dd87a0f10. Inspected dispatch has no GK
calculation; SOLV enum/comments are not support. GPU copies polarity/pdamp but
not douind. No GPU build. CPU coupled GK is the next affordable candidate.

NEXT TINKER_FRAMEWORK_SOLVER_PLAN.md declared BEFORE energies:12framework-only
CPU energy evaluations (3structures×all-polarizable/frozen/frozen-tight/rigid
frozen-tight), native electrostatics+GK, exact same atoms/charges. No La/Ca
score from ion-excluded framework controls. Verify frozen induced dipoles,
convergence,0.01kcal numerical/rigid gates and cost before full-source hybrid.
Need implement native energy frontend/runner using existing pinned library,
record all actual radii/flags; source current energy/prmkey/induce examined.
induce0c already zeroes all vacuum/solvent dipoles at each solve, and uses
douind throughout. Pin solver settings before output. No solver tasks prepared
or submitted yet. Active goal continues; baseline and all unrelated jobs intact.

## Latest: AMOEBA framework preparations pass; source-damping backend gap isolated

Completed `AMOEBA_CAPABILITY_PLAN.md` with exactly three real preparations:
GGR1GLG 4697 framework atoms, alpha1F6S 1931, alpha6IP9 1897. All protein and
retained-water coordinates, identities, hydrogens and disulfides unchanged;
each metal remains explicitly unparameterized in the full physical ledger.
Framework charges −6/−7/−7; waters 0/2/3. No energy, force, DFT, ML, optimization
or GPU calls. Installed AMOEBA2018/GK, OpenMM8.5.1. Six tests pass (3.659s),
manifest pins/dry-run pass. Preparations 69.9076944924891 wall-s,
69.65033699600001 process CPU-s, peak RSS480652KiB, no Slurm allocation.

New opt-in tool `scripts/mace_amoeba_capability.py`, tests of same name;
workspace `workspaces/mace_omol_20260917/amoeba_capability_v1/` with serialized
frameworks, all physical mappings/multipoles/covalent maps and receipts.
Manifest SHA4bf2e37dad96c9eba1fb9be45bffae2f45aad318a0e8efa8678d5c6e838c9339.
Compact `AMOEBA_CAPABILITY_REPORT.md`, `_RESULT.json`, `_COMMANDS.md` and
`AMOEBA_ACCOUNTING.md` under diagnostics/mace_omol_20260917.

Installed API/header plus matching upstream8.5.1 source verified: GK induction
is coupled and evaluated inside multipole kernel; separate GK force-group
energy returns zero in Reference. Do not build an additive scalar from it.
OpenMM exposes per-atom Thole minimum, not published La POLPAIR. Paper Table1
amide pair .349/.299 exceeds La default .250; can't implement with unchanged
default/min rule. La09/La22 are not AMOEBA2018. Frozen distributed QM source
also is not the paper's formal +3 polarizable ion. No complete hybrid qualified.
Source pins/verification `amoeba_sources_v1/`; installed/upstream headers byte
identical. Downloaded PMC file is CAPTCHA; paper was read through web tool.

NEXT `TINKER_CAPABILITY_PLAN.md` declared: inspect/pin maintained Tinker and
Tinker9 support for POLPAIR, frozen induced variables and GK; retrieve exact
primary parameters; one isolated CPU/reference build if needed; at most3real
topology exports, NO scientific energies/forces/DFT/MACE/optimization. Do not
run an energy analyzer as a format check. Resolve boundary/source damping and
accounting before declaring an energy pilot. PATH/selected bin directories
have no Tinker executable (only tinker_to_amber converter); broader installation
inventory not done. No installation or Tinker run has started. Goal active.
Latest live jobs still unrelated PLM1200794 running, 1200795/96 dependent,
H2001200809 pending; recheck. Baseline/default unchanged; no new approvals.

The sections below are historical checkpoints; their NEXT entries are superseded
by this top entry.

## Latest: responsive density/solvent candidate complete;3/4 ordering, boundary fails

NativeDFT1200983 (8),utilities1200986 (8CHELPG+8vpot),GB1200989 (44) complete.
All8variational,5near-boundary representation and50solver/componentchecks pass.
Orderingimprovesfromfrozen1/4to3/4: alpha1F6SminusGGRext/conn +2.795344719/
+7.233492066; alpha6IP9 -2.107829947/+2.330317400. GGRpartition-4.438147347
failsunchanged2. No favorablecore/structureselection orproductionpromotion.
Responsepaired<.18kcal;updatedGB shiftsalphaR+2.15/+2.44 andGGRR-3.18/-3.27.
Thus solventupdateddensitymatters despitealmostcancelledelectronicresponse.
Notbroadvalidation:2consumedbiologicalgroups,qualifiedcross-studyalpha evidence.

Artifactsunderworkspaces/mace_omol_20260917:
responsive_quantum_v2/manifest.json SHA
c50a571b0945e89f239194d1d8311c607a39a76dc16b550683e0c406fc4e2c78;
responsive_quantum_result_v2.json andpreservedcollectorresponsive_quantum_collection_v2;
responsive_charges_v3/manifest.json SHA
b20e5c45abfbe895da1d5776ed79c6db7ea24ab696b6ad0ecbdf9f6704e4e46b;
responsive_charge_report_v2/result.json;
responsive_GB_v1/manifest.json SHA
4c89ed5c9af1b64c80eaace9d5cfff466c3144686b9bfb7bf69d574c16caae7c;
responsive_GB_result_v1.json, responsive_cost_v1.json. Compactreport/results/
commandsRESPONSIVE_FIELD_*. Nativeoutputs/wavefunctions/failedattemptsretained.

Recovery:initialcollectorwronglymatchedconditionalnumerical-gradientwarning;
preciserepairtestedagainstactualanalyticoutputs+corruptedrealheader,noDFTrerun.
Firstchargepreplocalvariable-shadowingfailurebeforemanifest;initialutility
preflight1200984failedbeforecalculationsdueexactcaplambdaequality. Same-node
read-only1200985proved1.11e-16floatingdifference;V3permits1e-12replay withrecorded
weightsunchanged. Scientificgatesunchanged. All21distincttests pass. Costs45529
allocatedcore-s,37362.387actualCPU-s,66GPU-s inclfailedpreflight/replay. Cached
sourcecostsadditional;fullproductionpaircostunmeasured. NoMACE/trainingcalls.

NEXT AMOEBA_CAPABILITY_PLAN.md declared: readmaintainedsolver/source+primary
parameters, atmost3realproteinparameterizationpreps (samephysicalfullGGR/alpha),
NOenergies/forces/DFT/ML/fit/optimization. Determine coherentpermanent/induced/GK
andcore-subtraction/boundaryaccounting beforedeclaring any energypilot. Notnew
method:JacobrequestedAMOEBAandoldMACEplanlinksLa-parameterpaper. Installed
read-onlyinventoryamoeba_installed_inventory_v1.json:OpenMM8.5.1 AMOEBA2018/GK
files/APIs available;La/modelcoverageUNVERIFIED. No score/backendpromotion.
Researchgoalactive;thiscandidatefailedrobustnessbutproducedusefulinformation.
Preservebaseline andconcurrentPLM1200794/95/96,pendingH2001200809;inspectlivejobs.

## Latest completed: explicit electrostatics + MACE short component

Jobs1200980/1200981 complete:8vpot+8shortcore calls,0newDFT/GB/chargefit/
wholeMACE. All4algebra and5near-boundary charge-quality checks pass. Projected
Ca/La coupling errors[-.099381,-.541941,-.455469,-.686763]kcal; GGRpartition
error-.442560. NearestQM/cap probe1.167–1.202A. Charge-fit refinement alone
is not supported as the principal fix. CandidateFAILS: partition-4.139186004
vs2; alpha-minus-GGR[-2.523559,+1.615627,-7.739381,-3.600195],1/4pass.

Sources manifestexplicit_field_short_v1 SHA
90d553d2480a1efb3911fdc9169774610c9b124d40db7e17e6b4ad45af32d439;
result explicit_field_short_result_v1.json, cost explicit_field_short_cost_v1.json.
All underworkspaces/mace_omol_20260917. CompactEXPLICIT_FIELD_SHORT_REPORT.md,
RESULT.json,COMMANDS.md. Eightdistincttestspass (4prep10.525s,3legacy1.138s,
1actual.031s afterinitialskip). Frozenenvdrypass. Costs63GPU-s,1240core-s,
177.707actualCPU-s;7.561124short-model-s;peakGPU480389120bytes. Utility88.160987
summedwall-s. AllcachedDFT/whole/fit/GB costs separate. No productionpromotion.

NEXT RESPONSIVE_FIELD_PLAN.md declared BEFORE newembeddedoutputs. Use same
permanentcharges/corestates and nativeDFT method, letcoreelectrons respond;
recomputeCHELPG/GB fromnewdensity, preserve shortfull/core reuses. This is an
explicit approximate one-way fixed-field+reaction model, notselfconsistentGB.
PriorPQQpermanentresponse helpedonly~.974partition, notnewmethod/expectedwin.
8newnativeDFT,8CHELPG,8vpot(combinedexterior+actualenvironmentprobes),44GB;
4environment-onlyGB reuses,0MACE. Oldgatesunchanged +variationalresponsecheck.
NOTIMPLEMENTED/PREPARED/SUBMITTEDYET. Verifyinstallednativeenergyaccounting,
reuseexistingrunner/helpers, no newworkflow oradaptiverefitting. Fullplanhas
exactcostevidence/inventory. Goalactive, baselineprotected, continuedauthorized.

## Latest: full-boundary GB and direct-coupling audit completed

Job1200975 COMPLETED48solvercalls,46/46 numericalchecks pass. CandidateFAILS:
GGRpartition+4.482882260 vs2; allfouralpha-minus-GGR differences -120to-126kcal.
No parameter/label/geometry rescue. Report FULL_BOUNDARY_GB_REPORT.md/RESULT.
Fullreportfull_boundary_GB_report_v1/result.json; manifestfull_boundary_GB_v1 SHA
69f88af3655ea6d9784d682e528f7cf3a780410b93dda999671f7446e425e446.
Cost57GPU-s/912core-s/117.953CPU-s;4.379188012solver-s,167908KiBpeakRSS.
GPUmemoryunmeasured. Four newtests+sixnative+twolegacyGB pass;actualintegration
0.004s afterinitialskip. No newDFT/MACE/chargefit/training. Oldbaselineunchanged.

Separate saved-state auditQMFF_coupling_audit_v2 passes20checks; no newscientific
calls. Alpha directR+132.187/+123.936 nearlycancelsGBcross -135.462/-128.192;
GGRextended direct-0.233/cross-10.058,connected+9.656/-15.357. OMOLcontext is
entangled and hasnoverifiedmatchingdirectterm. Do notjustaddC toOMOL:possible
doublecount. Beyond36Aalpha noatoms,GGRtail~-13.2 cannotrepair~120kcalfailure.
AuditV1preflight comparedlive/frozenimplementationpaths despiteidenticalhash;
V2verifieshashesandallscientificfields, no physics/gate change. Resources retained.

NEXT EXPLICIT_FIELD_SHORT_PLAN.md declared BEFOREnewpotential/shortoutputs.
Reuse earlier qualified MEDIUM MACE-POLAR shortcomponent (notnew); its old
PQQfixedfield andalpha-only failuresremain. Candidate expression:
DFTvaccore + EXACTsaved-density directcoupling + unchangedfullGB + short(full-core).
No OMOLtotal, MACEelectrostatic/electronterm orcoreCPCM. Eightnewcore shortcalls
+eightnativevpotcalls atactualnonzeroFFsites;0newDFT/GB/wholeMACE/chargefit.
Existingcouplingpointapproximation isdiagnostic; uniformexactdirect chosenbefore
outputs. Fixedpairedpoint-vs-exacterror<=1kcalpercase/partition, same2partition
andfour>0.02orderingcriteria. Forcespartialonly; no relaxation/entropy/threshold.

field_short_reuse_audit_v1.json verifies6exactwhole mediumshortreuses from
mace_global_benchmark_20260916/mace_v1/medium/collection_job_1200701.json;
536taskgeometrycomparisonsin3namedarchives found0exactnormalizedcorematches.
Qualifiedshortadapter polar_scale_shift_exact_readout_v1;
workspaces/mace_short_engine_20260916/pilot_v1/collection_job_1200736.json.
New16lowlevelcalls NOTimplemented/prepared/submittedyet. Implementminimalexisting
runnerdispatch/sourceadapter, finitepreps/tests/dryrun, then launchdocumented
GPU/CPUallocations autonomously. Allournewjobs terminal; preserveconcurrentjobs.
No productionpromotion,push,permissiongate orprojectcompute/timebudget. Goalactive.

## Latest: normalized charges complete; next full-boundary solvent test declared

Job1200970 COMPLETED eight CHELPG+eight vpot utilities, zero failures. All
endpoint and paired fit/projection quality gates pass. Differential field
relativeRMS0.3493–0.5532%; after geometric cap projection0.3715–0.5852%.
Absolute endpoint potential errors are larger; alpha1F6SLa projected10.0724%
passes the predeclared absolute branch. This is NOT a solvent-energy/accuracy
pass. No full environment charge model or GB scalar exists yet.

Report ../mace_omol_20260917/NORMALIZED_CHARGE_REPORT.md and compactRESULT;
fullnormalized_charge_report_v1/result.json SHA
2689bcdfa9b3f2d0358c5e21fa8b5b38533b36502c7e5c3f9826ee7487041d10.
Manifestnormalized_charge_v1 SHA
31caf15249cbfab1ff381f1185c7243e77c69a003c735f52ca80e5ad63bc3fc6.
Four distinct real-fixture tests pass;3preparation5.094s,actual0.027s after
initialskip. Cost242wall-s/8CPU,1936allocatedcore-s,897.082CPU-s,1523816KiB
peakRSS. Utilities891.622784summedwall-s. No newDFT/MACE/GB/GPU/training.
Cost normalized_charge_cost_v1; utility cost additional to prior subtotals.
All source copies/receipts/matrices/potentials preserved, sourceGBWunchanged.

NEXT: FULL_BOUNDARY_GB_PLAN.md declared AFTER charge results and BEFORE full
charge assembly/solver energies. New protocol
normalized_QM_projected_ff19SB_full_OBC2_vacuum_hybrid_v1. Assemble Pq_QM plus
identical ff19SB exterior, no FFcharge on projection support, exact residue
formal-charge ledger and local bond-neighbor redistribution. Preserve full
physical boundary, all waters, assembly, H coordinates. Add ONLY reaction-field
G to vacuum hybrid; no bare Coulomb or coreCPCM. Fixed existingOBC2parameters.
48specified solver calls, no newDFT/MACE; see fullplan for identity,rigid,
component,Reference-CUDA and same partition/ordering gates. Not implemented,
prepared or submitted yet. Do not execute unmanifested solver work. Next actions:
implement minimal source/boundary adapter and existingrunner dispatch; prepare,
realfixturetests/dryrun; submit documented A5000 allocation; collect/compare.
No new permission required. No project CPU/time budget or defaultpromotion.

Matched-H phase c63055a. Charge phase commit recover fromgitlog. Allournewjobs
terminal. Inspect concurrentPLM/H200 state before acting; do not cancel/duplicate.

## Latest result: matched H complete; consistency improves, accuracy fails

Jobs1200950/1200951 COMPLETED eightDFT/eightMACE, zero failures; six whole reuses.
13/13 numerical checks pass. GGR partition+0.118496483 (old-H+1.829806242)
passes2. Alpha1F6S minus GGRextended/connected +4.329631810/+4.211135326 passes;
alpha6IP9 -0.805098984/-0.923595467 fails. Fixedall-fourgateFAIL2/4. No aqueous
score, inherited threshold, relaxation or broad accuracy claim. Baseline unchanged.
[Report](../mace_omol_20260917/MATCHED_H_REPORT.md), MATCHED_H_RESULT.json.

Source-H normalization is the pre-existing ff19SB/water rule; all heavy atoms,
caps, protonation and inventories unchanged. Exact prepared/quantum/MACE pins
remain below. Full result matched_H_report_v1/result.json; cost matched_H_cost_v1.
Five distinct new tests and six native regressions pass; actual integration3.973s
following initial explicit skip. No new whole inference/solver/training.

DFT751wall-s/64CPU,48064allocatedcore-s,43768reportedCPU-s. MACE67GPU-s,
1072core-s,83.861CPU-s;model3.856576160s,802472448GPUbytespeak. Model subtotalV16
451success/4failedcalls,12589GPU-s,264656core-s,30942.989CPU-s. Both recent DFT
phases are additional/separate. Do not repeat completed jobs/tests. Goalactive.
Next scientific step not yet declared at this checkpoint; investigate a coherent
solvent extension rather than claiming this vacuum candidate is validated.

Earlier submitted/pending notes below are historical and superseded here.

## Latest execution checkpoint — intact gradients qualified

This supersedes the earlier active-gradient notes. Full job1200886 COMPLETED;
all six calls and11/11 numerical checks pass. WholeGGR has4698 mapped atoms;
Ca/La gradient evaluations take40.699184/40.789196seconds and12,011,144,704bytes
peak GPU (11.19GiB). Both center energies equal the archived scalar. Along the
fixed metal-to-GLN140/O direction, grad(R)=+8.836946204modelkcal/A; actual odd
0.01A change0.088339618 versus predicted0.088369462. No physical-force,
relaxation, entropy or biological-improvement claim. Production unchanged.

Final adapter `omol_exact_checkpointed_edge_product_autograd_v3` requalified
on the exact73-atom1H4I core (job1200885,23/23checks). The actual fullV2 CUDA
OOM was resolved by replacing index_add with equivalent scatter_add: saved
edge-message tensors no longer consume GPU memory outside checkpoints.
Core/native max gradient error9.02e-14eV/A. Earlier coordinate/report roundoff,
TorchScript early-stop, collector-serialization and OOM attempts remain visible.
No coordinates, scientific tolerances or criteria changed during recovery.

Authoritative artifacts under workspaces/mace_omol_20260917:
- core: masked_gradient_core_v4/manifest.json, core_report_v5 is named in full:
  `masked_gradient_core_report_v5/result.json` (all checks pass).
- whole: `masked_gradient_full_v3/manifest.json`,
  `masked_gradient_full_report_v1/result.json`, plus all-atom TSV/paired array.
- costs: `masked_gradient_cost_v1.json`, seven terminal jobs1200863/64/65,
  1200878/84/85/86.27successful calls,2failed model calls,410GPUallocation-s,
  6560allocatedcore-s,470.555reportedCPU-s. Two preflight failures made0calls.
- cumulative: `intact_engineering_status_v12.json`:389successes,4failedcalls,
  12022GPU-s,255584allocatedcore-s,30231.464reportedCPU-s.
- seven final real-fixture tests pass:6in12.883s plus mapped-export test1.531s;
  none skipped. Earlier native regressions and all resource receipts retained.

[Gradient report](../mace_omol_20260917/MASKED_GRADIENT_REPORT.md) includes replay
commands and immutable attempt details. Source files: mace_omol_gradients.py,
mace_omol_gradient_worker.py, mace_omol_gradient_run.py; minimal dispatch in
mace_omol.py. No production score/energy-only adapter changes.

## Latest scientific results after the gradient milestone

Response job1200888 COMPLETED40/40calls, zero failures. Exact archived DFT inputs
and physical Jacobians were used. Own analytic derivatives pass24/24; direct
DFT deformation-energy comparisons pass6/24; curvature passes19/24; DFT-gradient
anchored predictions pass24/24 under their separate0.02 absolute floor. The
anchored pass does not override five curvature failures. Maximum direct error
0.320106; anchored0.019487kcal-scale. No relaxed score is supported. Report:
`masked_response_report_v1/result.json`; compact MASKED_RESPONSE_REPORT.md.

Cost320GPUallocation-s,5120allocatedcore-s,416.572reportedCPU-s;20.148532summed
model-s,1,031,137,792bytespeakGPU. Three distinct real-fixture tests pass; actual
integration was explicitly skipped before outputs then passed after completion.
CumulativeengineeringV13:429successful calls,4failedcalls,12342GPU-s,
260704allocatedcore-s,30648.036reportedCPU-s. No newDFT/solver/training.

A separate static candidate was then declared: DFT/CPCMcore + masked(full-core).
Its cached partition prerequisite FAILS: connected-minus-extended GGR shifts
DFT−7.343500873, maskedcore−27.503809106, hybrid+20.160308232 versus2kcal gate.
All original-H shared source/cap mappings pass (max4.99e-11A). Zero new inference;
DO NOT run its conditional six whole-chain calls. Full result:
`masked_context_partition_v1/result.json`; compact
MASKED_SUBTRACTIVE_CONTEXT_REPORT.md. Two actual-fixture tests pass4.384s.
No adaptive change, learned multiplier or omitted failure. Baseline unchanged.

## Shared learned-neutral feature complete: partition failure

Eight real center gradient calls completed1200901;16/16 numerical checks pass.
All four charge0La native energies agree exactly;maxgradienterror9.77e-15eV/A.
Connected-minus-extended GGR:DFT/CPCM−7.343500873,learned−18.255784798,
hybrid+10.912283925kcal-scale;declared2gateFAIL. DO NOT run conditional six
whole calls or sweep charge categories. Physical charge/spin/coordinates unchanged.
Report ../mace_omol_20260917/SHARED_NEUTRAL_FEATURE_REPORT.md.

Manifest shared_neutral_core_v2 SHA4b3420ef72fc742b25cc407b1992073e3a55da2d28f25f15c9055781bb7db814;
result shared_neutral_core_report_v1/result.json. Eight scientific tasks unchanged
fromV1, which failed before inference at auxiliary native report host roundoff.
V2 retains exact raw receipts/arrays and all original thresholds/decisions.
Recovery shared_neutral_reference_recovery_v2.json. Three new tests and six
native regressions pass;actual integration executed after initial explicit skip.

Both jobs1200893/1200901:82GPUallocation-s,1312allocatedcore-s,93.244CPU-s;
modelsum6.595379632s,peakGPU1,031,137,792bytes. NoDFT/solver/training/optimization.
Costs shared_neutral_cost_v1 and cumulativeintact_engineering_status_v14:
437successfulcalls,4failedmodelcalls,12424GPU-s,262016core-s,30741.280CPU-s.
Local resources and preflight failures retained. Baseline unchanged;goal active.

## Matched vacuum complete: solvent explains most of the discrepancy

Eight endpoints completed1200905, no failures. Connected-minus-extended GGR
R_vacuum=-25.674002864 versus R_CPCM=-7.343500873. Their solvent difference
is+18.330501991kcal/mol. Raw-zero learned core shift=-27.503809106, leaving
vacuum hybrid+1.829806242, PASS at fixed2 tolerance. Shared-neutral residual
is-7.418218065, FAIL. Component identities close; no old candidate is relabeled.
Report ../mace_omol_20260917/MATCHED_VACUUM_REPORT.md.

All eight actual gradients pass native-state/ECP/component and coordinate checks.
Four distinct real-fixture tests pass; integration initially skipped then passed
1.812s. No new MACE/solver/training, optimization or numerical DFT gradients.
Job cost618s wall/64CPU,39552allocatedcore-s,35369reportedCPU-s. Separate cost
record preserves MACE forward counts rather than counting DFT as model calls.

Manifest matched_vacuum_v1 SHA
f43b66d7dd708082e643ef42f9ed939424787afabdcbd85f319cf7b237a98960;
matched_vacuum_report_v1/result.json,matched_vacuum_cost_v1.json.
New script mace_omol_vacuum.py, frozen source and existing ORCA runner/locks.
Bounded archiveaudit:165inputs/25geometrymatches/zero vacuum counterparts.
Baseline/default unchanged. No aqueous score or force/curvature accuracy claim.

## Whole vacuum hybrid complete: consistency passes, ordering fails

Job1200924 completed6/6 original-H whole calls.11/11numerical checks pass;
partition+1.829806242 passes2. Ordering only1/4passes:alpha1F6S versus
GGR[extended,connected]=[+0.558714448,-1.271091793];alpha6IP9=
[-0.265848708,-2.095654950]. Fixedall-fourgateFAIL. No threshold/label/geometry
change or aqueous score. Report ../mace_omol_20260917/VACUUM_HYBRID_REPORT.md.
Manifest vacuum_hybrid_v1 SHA
8d054667b104384a49f3cbaa7b97ddd712b8c620c4660f2907f79afa4bc72d7e;
fullvacuum_hybrid_report_v1/result.json. Three new and six native tests pass;native regression12.283s. Actual
integration3.895s afterinitialskip. CostV2 adds final regression resource pins. Initialdry-run rejected execute-only CLI
argument beforecompute; corrected frozen dryrunpassedwithout sourcechanges.

Cost98GPUallocation-s,1568allocatedcore-s,117.848CPU-s;47.057306305model-s,
6186040832bytespeakGPU. No newDFT/solver/training/gradient/optimization.
Model engineering subtotalV15:443success,4failedmodelcalls,12522GPU-s,
263584allocatedcore-s,30859.128CPU-s. MatchedvacuumDFT39552core-s/35369CPU-s
is ADDITIONAL and separately recorded. No double count or zero substitution.
Newmace_omol_vacuum_hybrid.py, smallmace_omol dispatch; existingworkersunchanged.

## Matched H normalization: implemented, jobs submitted

The declared MATCHED_H_NORMALIZATION_PLAN.md is now implemented in
scripts/mace_omol_matched_h.py, with small existing MACE dispatch. Four real
preparation tests pass5.157s; actual integration explicitly pending. Frozen
quantum and MACE dry-runs pass. Both reuse existing executors, locks and workers.

Prepared: matched_H_prepared_v1/preparation.json SHA
72952416d3394b2a23a2663a29e78ff244d17f55a505b6e861577c46a4f5cd3c.
Quantum: matched_H_quantum_v1/manifest.json SHA
d63801217ed7b6092cf67f4224b34b28574199aefc1ca291abde34c07fee6275;
job1200950,64CPU/four16-rank tasks, largest cores first for utilization.
MACE: matched_H_mace_v1/manifest.json SHA
5db20a83f79a332e63c3e399bc8c33e3ffc8de823cb4f945c6b5c7cef8fc96d2;
job1200951,A5000/16CPU/64474MiB,eight energy-only core evaluations.
Inspect live state; do not infer termination from this checkpoint.

Every source H comes from the previously frozen ff19SB/water normalization.
Heavy coordinates and sigma caps remain exactly unchanged; charge, donor,
protonation, water and element inventories match. Source-H changes21/47/16/18.
Full normalization arithmetic replays within1e-12A; derived max-error summaries
may vary within that same tolerance between hosts, not change scientific inputs.
Six actual normalized whole energies reused. Alpha duplicates agree within0.01;
lexically first task selected, all alternatives retained. No new whole inference.

After both jobs finish, use frozen matched_H_quantum_v1/implementation/
mace_omol_matched_h.py report --quantum[quantum manifest] --mace[MACE manifest]
--output workspaces/mace_omol_20260917/matched_H_report_v1. Then actual integration
and native regression, costs/report/vault and scoped commit. Old original-H
failure remains immutable. Same fixed2kcal partition and four>0.02 ordering
criteria; no aqueous/broad accuracy/response claim. Goal active.

Vacuum DFTphasecommitted28a86e0. Initialfindingemail18:24UTC acceptedlocalrelay;
wholehybridfollow-up accepted localrelay;failedaccuracyresult reported. Allnewjobsabove
terminal. PendingnativeH2001200809 and otheragents'PLM1200794–96 untouched.

Shared-neutral commit296f809;response/context215c794;gradientc22a881. Do not repeat passed
engineering gates. Other agents'PLM1200794–1200796 and pendingnativeH2001200809
preserved;recheck live state. No push/default promotion.

## Actual result: masked MACE candidate works on the supported PQQ panel

Protocol `mace_omol_intact_charge_feature_ablation_descriptor_v1`.
The same pinned OMOL checkpoint and whole-chain physical preparation are used,
with only the raw total-charge embedding zeroed before the joint projection.
Actual physical charges, spin and weights remain recorded and unchanged.
This is a learned descriptor, not a quantum electronic-state energy or binding
free energy. It lacks responsive long-range electrostatics; whole proteins are
outside published model training sizes.

- Development1200828:42/42newforwards, all numerical/three-order gates pass.
  XoxF−MxaF44.329179modelkcal; alpha1F6S−GGR12.741688;
  alpha6IP9−GGR4.112398. The two alpha structures are one qualified affinity
  comparison. Sodium score change−1.07e-8 is architectural consistency.
- Canonical1200830: **COMPLETE100/100newforwards**,8exact crystal reuses.
  All25calibration proteins separate, gap9.108591modelkcal. MxaF1H4I17.579956
  is Ca; XoxF4MAE61.909135 is La. **2/3transfers valid and correct;1KB0
  unsupported. The declared all-three transfer gate remains false.** All cases
  consumed; composition already separates the canonical panel. No broad
  affinity accuracy or incremental value beyond composition established.
  [Canonical report](../mace_omol_20260917/CHARGE_ABLATION_CANONICAL_REPORT.md).

Actual canonical cost3546GPU/allocation-s(59m06s)oneA5000,56736allocatedcore-s,
3925reportedCPU-s(seconds precision),2657.578770summed model-evaluation-s,
13751386112bytes peakGPU,1769204KiBpeakRSS. Local work separately measured.
Existing bound pairs measured49.85685–60.772997forward-s(median53.068669),
53.744777–65.495802worker-wall-s(median57.267661), excluding prep/controller/report.

## Two-call reference and runnable research classifier

The exact descriptor factorization is
R_mask=[T_bound,Ca−T_bound,La−(C_Ca−C_La)]*23.06054783061903.
Fixed1H4I native detached-node-plus-embedding terms are Ca−18430.794927644074
and La−850.2720512362149 model eV. These are not aquo/quantum-ion energies.
Factorization development60checks and full-panel189checks pass; full-panel
maxerror2.357415596e-8 vs declared0.01tolerance. Zero additional inference.

Use `scripts/mace_omol_prepared.py`: audit,prepare(optional--factorization),
existing mace_hybrid dry-run/execute/collect,report(optional--calibration).
Source-backed whole-chain preparation required; arbitraryXYZ is unsupported.
Strict peptide connectivity, atom/source replay, charge/spin/paired coordinates,
cofactor/water state and actual endpoint receipts enforced.

IMPORTANT driver for prepared-input operations:
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
(OpenMM). GPU workers use the separately pinned existing MACE venv. Never
install into either to hide a mismatch. MACEvenv works for pure panel reports.

Matching two-call PQQ bands:Ca<=42.346663801320496;La>=51.455254788976355;
otherwiseinconclusive. PQQ functional association only; no non-PQQ affinity band.
Original four-call bands stay immutable and numerically separately typed.
Reference: workspaces/mace_omol_20260917/masked_calibration_v1/reference.json
SHA217a127f7ba8f6cccb59d1c63b2d389e9022449a5feb09f58eba2e3b04cba54d.
Frozen current reporter: masked_calibration_source_v1/implementation/.
[Commands](../mace_omol_20260917/PREPARED_INPUT_COMMANDS.md),
[reference report](../mace_omol_20260917/MASKED_CALIBRATION_REPORT.md).

Executed no-new-inference interface reports:
prepared_interface_pqq_calibrated_report_v1:17.57995558550866,Ca;
prepared_interface_ggr_calibrated_report_v1:23.98742145552895,classificationnull,
outside_canonical_PQQ_calibration_scope. Reference and direct algebra replay.
All five actual-artifact calibration guard tests PASSED in675.322seconds,
none skipped; masked_calibration_validation_v1.json pins tests/log/resources.
Standalone inspected PDF/SVG/PNG figure inmasked_calibration_figure_v1. Earlier4prepared,3factorization,
8ablation/interface and6native regressions passed. No synthetic science fixtures.

## Artifacts and live work

All candidate products under workspaces/mace_omol_20260917/ unless specified.
- charge_ablation_canonical_v1/manifest.json SHA
  2415b4dacd755530dec888017d76afe8895708213a36ef09497c7264e341476f.
- collection_job_1200830.json SHA
  144d7e312716672345358b319823da8cd4b61a02782a849ea1b013c921b660de.
- charge_ablation_canonical_report_v1/result.json SHA
  9c12a55ac09aab4c68f2399a1ff784ebd078d6f6a8ea020b85fccb577d017f25.
- factorization_report_v2/result.json SHA
  35ba8bb53b209f650b8392ff03f0c426d0f95455b39b18c082e25aabb5434fbe.
  V1auxiliary NumPy reduction differences preserved; v2uses math.fsum.
- prepared_interface_pqq_two_call_v2/manifest.json has0newtasks/2realreuses.
- prepared_interface_ggr_two_call_v1/manifest.json has0newtasks/2realreuses.
- prepared_interface_ggr_fresh_v1 contains4UNEXECUTEDtasks. Do not duplicateGGR.
- Engineering statusV8:348successfulforwards,2OOMs,11282GPUallocation-s,
  243744allocatedcore-s,29420.954reportedCPU-s. Includes4invalid native1KB0
  diagnostic calls and prior failures; local work separate.
- Native H200job1200809 remains pending unchanged14-taskmanifest. Inspect live
  state; do not cancel/duplicate/alter its manifest or any executor lock.

Earlier own commits2fcf2f3,287d3d0,30cb9b4,410390a. This checkpoint accompanies
the scoped canonical-calibration implementation/results commit; recover its
hash fromgitlog. Many unrelated dirty files are concurrent/historical.
Never blanket-stage/reset/stash. SESSIONS own append
must be staged separately from another agent's working additions.
Email accepted by local relay20260917T145316Z_masked_canonical_pass. No delivery
confirmation claimed. Vault note now includes the canonical milestone, costs, limitations and source finding.

## Important source finding: 1KB0 is not a simple terminal-cap repair

[Raw-source diagnosis](../mace_omol_20260917/ONE_KB0_RAW_CHEMISTRY.md): deposited
TRO512 is present, but fixed-core-oriented normalization removed it and
whole-chain parsing connectedLEU511directlytoGLU513(4.750089A). A separate
missing574–578loop produces573to579(18.969263A). Raw HEC802 is covalently joined
toCYS604/CYS607 and coordinatesHIS608/MET647. Heme was explicitly excluded from
the old normalized source. Do not fabricate a whole-chain pass by converting
TROtoTRP, omittingheme, guessingironstates, or bridgingthegap. New full-system
support needs actual modified-residue/heme/loop chemistry under a new policy.
Old fixed-core1KB0 remains valid under its own recorded scope and unchanged.

## Earlier failures remain results

Native whole-chain OMOL canonical gap−306.006296; no bands. DisconnectedNa
shifted four scores4–15kcal/mol and reversedXoxF/MxaF. Checkpoint charge features
alias−100..−6and+11..+100; verified behavior, no inferredtraininghistory.
Masked representation is a separate development response, not erased native
history. Read INTACT_CANONICAL/SPECTATOR/LOCALITY reports, chargeembeddingaudit.
Prior POLAR+GB reversed allthree directions; local/hybrid partition and paired
mechanics failed. Relaxation/entropy remain unavailable. Do not rediscover them.

## Next useful work

Calibration tests, scoped implementation/results and vault update are complete. The operational
candidate milestone merits reporting, but goal remainsactive: all-three transfer
support and broad/composition-challenging predictive usefulness remain open.
Existing CHALLENGE_PANEL_CURATION.md records eligible versus unresolved affinity
controls. AQUALYSIN_EVIDENCE_FOLLOWUP.md records read-only follow-up: adding/-char/en
recovers the primary2002PDF. The2019Ca-1assignment remains an inference from
structural analogy, not direct assay-site mapping. No new label/preparation/score.
Do not spend effort merely recounting provenance: prioritize chemically sound
new inputs and tests that can reveal predictive weaknesses. No new user approval
is needed within this goal. No automatic production/default promotion.


## Expanded GGR robustness: completed failure, 2026-09-17

Jobs1200845/1200846 COMPLETE four forwards, all numerical checks pass.
GGR ordered [1GLG,2FW0,2FVY] = [23.987421458,44.884127357,45.975500410]
model kcal. Both new structures reverse both alpha-minus-GGR comparisons:
only2/6margins pass; all-case robustness FALSE. No label/threshold/input rescue.
The separate PQQ result stands; broad affinity improvement is unestablished.
See ../mace_omol_20260917/GGR_STRUCTURE_ROBUSTNESS_REPORT.md.

Actual cost114GPUallocation-s,1824allocatedcore-s,118.469reportedCPU-s,
52.124891summed model-s,6108277248bytespeakGPU. Cumulative engineeringV9:
352successfulforwards,2OOMs,11396GPU-s,245568allocatedcore-s,29539.423CPU-s.
All costs include prior failures; local work separate. Four source/actual-result
regressions pass6.713s,none skipped. No newDFT/solver/training/forces.

Source bridge: scripts/mace_omol_source_prepare.py with explicit raw inventory,
exclusions, source metal/paired-coordinate and peptide checks; original1GLG
physical coordinates exactly reproduced, raw1KB0TRO512 explicitly rejected.
Policy omol_source_backed_chain_A_ff19sb_H_explicit_exclusions_v1.
New actualpreps ggr_source_bridge_{2fw0,2fvy}_v1/source_preparation.json;
source-row/evidence/exclusions JSONs in ggr_structure_sources_v1/.
Source tests and initial failed comment-byte assertion preserved; no scientific
coordinate discrepancy. Full result ggr_structure_report_v1/result.json.
Saved readout diagnosis places changes across metal and several donor residues,
not a unique causal term. Same donor identities; several distance shifts0.1–0.2A.

## Multisite family coverage: complete, supporting gate fails

PlanMULTISITE_PANEL_PLAN.md; resultMULTISITE_REPORT.md/MULTISITE_RESULT.json.
Jobs1200851–1200855 COMPLETE tenforwards, allnumericchecks pass. Identical
all-Caenergies acrossselectedsiteatomorderings (reportederror0).
Orderedparvalbumin[CD,EF]=[27.431183345093064,64.34764681174364]modelkcal.
CDexceedsonly1GLG;EFexceedsallthreeGGR. Supportingcross-studycontrasts4/6pass,
all-casegateFALSE. No sub-kcalaffinitysiteorderclaim. Aequorin[EF1,EF3,EF4]=
[13.04767291865186,33.17334887989371,51.5641437611785]. Siteunresolvedassay,
no bestsiteselection,physicalzero,PQQband orbinaryassayreproductionclaim.

New scripts/mace_omol_multisite.py, multisite_report.py. Policy
omol_intact_multisite_fixed_background_Ca_source_acetyl_ff19sb_v1.
Rawsourceheavyatoms/covalentlinksreplay;4CPVactualACE0cap,H,bondretained.
4CPV1611atoms,1backgroundCa,1water,QCa−3/QLa−2.1SL82866atoms,
2backgroundCa,3waters,QCa−4/QLa−3;its10missingN-terminalresiduesexplicit.
Onlyselectedmetalchanges.Background-Caindicesmustbeexplicit;default
single-metalcheckandallpreviouspreparationpathsremainunchanged.
V1prepfailedbeforeinferenceatoldsingle-metalguard;V2fixedthespecifictechnical
limitation.Noalteredscientificselection.Unknownrawchemistry/gapsstillfail.

Artifactsunderworkspaces/mace_omol_20260917/:
- multisite_recipes_v1/{PARV_4CPV,AEQ_1SL8}.json
- multisite_prepared_v2/CASE/preparation.json (all5prepared)
- multisite_scoring_v1/CASE/manifest.json and actualcollections/receipts
- multisite_reports_v1/CASE/result.json; multisite_comparison_v1/result.json
- multisite_cost_v2.json, multisite_sacct_v1.tsv. CostV1hadoneauxiliarytest
  resourcehashcapturedinflight;V2correctsit,scores/allocationcostunchanged.
- Sevennew/source/resulttestspass14.973s,fourlegacyPQQ/GGRtestspass28.313s,
  onenewactualmultisitereporttestpasses9.923s;12distincttests,noneskipped.

Actualmultisitecost216GPUallocation-s,3456allocatedcore-s,221.486reported
actualCPU-s;63.758972model-s,80.695251workerwall-s;3936891392bytespeakGPU.
Localpreparation/dryrun/report/testreceiptsseparate. CumulativeV11through1200855:
362successfulforwards,2OOMs,11612GPU-s,249024allocatedcore-s,29760.909CPU-s.
NoDFT/solver/training/forces/relaxationcalls. Goalactive;productionunchanged.

## Next declared work: analytic descriptor gradients

MASKED_GRADIENT_PLAN.md declaredbeforeimplementation/inference. Add a separate
exactbackward-capableedge/productadapterwithcheckpointing;currentqualified
energy-onlyadaptersmuststayunchanged. Sameweights,mask,spin,float64andalgebra.
StageA10forwards onreal73-atom1H4Icores from ablationdevelopmentV2:2nativegrad,
2batchedgrad,2rotatedgrad,4signedmetaldisplacements. StageBconditional6forwards
onreal4698-atomGGR1GLG:2gradcenters,4signedmetaldisplacements. Frozen numerical
criteria and exactsources inplan. No gradientimplementation/callsdone yet.
Outputsarederivativesofmaskeddescriptor,notvalidatedphysicalforces. Do not
addrelaxation/entropy merelybecauseagradientorcurvaturematrixcanbecomputed.

GGRcommit1038b55;nextownscopedcommitcoversmultisitecode/result,nextplan.
Vaultandagentguideupdated. LatestemailacceptedrelaynotifiedGGRfailure;do not
spamrepeats. NativeH2001200809stillpendingunchanged;inspectlivejobs. Allnew
A5000jobsabovearecomplete.Noautomaticpromotion/push/per-analysispermissiongate.
Preserveotheragents'PLMjobs,watchers,dirtyfilesandlocks.


## Gradient implementation history

The original in-flight notes are superseded by the latest checkpoint above.
All recovery versions and failure receipts remain under workspaces; see the
completed gradient report for the full sequence. The next scientific task is
the declared response screen, not another rerun of passed engineering checks.
