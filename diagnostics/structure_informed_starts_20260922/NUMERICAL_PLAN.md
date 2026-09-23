# Q88JH5 near-identical geometry numerical check

Approved by root under Jacob's contained parallel-pilot authorization: exactly
8 newLa nativeGFN2 endpoints, noMACE/DFT, geometry optimization or broaderrescore.
The primary alternative-start experiment is already complete and stays immutable.

Question: why does itsLa vacuumGFN2 energy change by0.007657590594Ha between
near-identical accepted geometries, while nativeMACE andALPB barely change?
The old andnew runtimeinput hashes agree; both containNoAutostart,MaxIter500,
SmearTemp300,UseXTBMixertrue,charge−2/singlet/438electrons and8MPIranks.
Hosts differ. Both reportSCFconverged, but oldprintedMAX/RMS densitychanges exceed
their displayed tolerances. This is an unresolved numerical/electronic-state
sensitivity, not evidence that the new geometry improves discrimination.

Freeze two actualXYZfiles: originaladaptiveLa and templateLa-searchLa forQ88JH5.
Copy bytes without coordinate rewriting. For each usevacuum andALPBwater, under
(1) unchanged primary input and(2) same input plus`Convergence Tight` inside%scf.
Keep explicit`UseXTBMixer true`; this is the previously demonstrated native-tight
recipe family, not ordinarySCF, MORead or a restart. NoexistingGBW/xtbw is copied.
Onefresh directory/start per task, all8 on one64CPU/128GiB CPUallocation with8×8
concurrency using theexistingrunner. MaxIter remains500, includingtighttasks.

Record all8 raw energies, actual mixer, electron/charge/spin, temperatures,
convergencecycles, printedenergy/MAXdensity/RMSdensityresiduals andtheir printed
tolerances, fractionaloccupations, charges, runtimeinputs/parameterexports and
receipts. A normaltermination doesnot erase densitydiagnostic failures. The
requirednative mixer andTight energytolerance≤1e−8 must appear in actual outputs;
otherwise markthat numerical variant unsupported, with no silent replacement.

Compare primary exact-XYZ repeat to eachoriginalreceipt; compare primary/tight
and old/new geometry within each numericalpolicy for bothendpoint energies and
matchedALPB-minus-vacuum transfers. Retain the prior qualification scales
0.1kcal/mol perendpointtransfer and0.2kcal/mol forcontrast diagnostics. These
are numerical acceptance scales, not bandpadding or a thresholdrefit. This La-only
check cannot supply a newly qualifiedfullCa/Laprotocol or a new absolute reference.
No favorable-branch selection, origin-historyoverwrite or automaticpromotion.
