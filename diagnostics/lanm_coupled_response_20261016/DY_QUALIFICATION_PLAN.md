# La/Dy native-GFN2 import qualification on one real peptide core

Declared16October2026 before new backend calls. Root owns execution. This is
a small prerequisite for a coherent cheap solvent calculation in the monomer
response experiment, not a new affinity model or a rerun of failed whole-LanM
vacuum subtraction. No production configuration changes.

## Fixed question and inputs

Does a full native parameter profile with only the already-corrected La arrays
and the three explicitly corrected Dy arrays reproduce author-parameter
standalone GFN2 energies and populations on the same real48-atom LBT3 core?

Use preparation12469's exact7CCO matched_core_La and matched_core_Dy XYZs,
charge−1. All coordinates, atom membership and solvent conventions are fixed.
GFN2 is the original spin-independent f-in-core Hamiltonian, effective
multiplicity1/uhf0; retain physical Dy(III) sextet as separate metadata. This
test cannot qualify open-shell response, DFT or MACE by association.

The eight endpoints are La/Dy × vacuum/ALPB(Water) × nativeORCA/standalonexTB.
Use a new task-local full JSON derived from corrected_la.xtb.json, changing only
Dy kcn/shpoly/lgam to the values in DY_PARAMETER_READINESS.md. Declare La and
Dy explicitly in the existing parameter guard policy. Preserve every other
field and declared shell order exactly. Import via the established diagnostic
guard, preflight and postflight, with a task-local full file. No native-default
fallback, global reordering or silent adoption by the La-only runtime registry.

## Settings and acceptance, fixed before execution

NativeORCA6.1.1 uses Native-GFN2-xTB, NoAutostart, UseXTBMixer true, SmearTemp300K,
MaxIter500, TolE1e−12, and its unchanged printed density criteria. One native
process per endpoint; derive MaxCore from actual allocation with25%headroom.
Use the existing pinned native runtime and serial invocation. StandalonexTB6.7.1
uses its original installed author parameters, --gfn2, --uhf0, --etemp300,
--acc0.0001, --iterations500, --norestart, JSON/full populations, and ALPB water
gsolv where specified. The exact command writes these as separate proper CLI
arguments. Solvent temperature298.15K and electronic smearing300K are distinct.

Each endpoint must have successful process termination, actual parameter/state
evidence, finite total energy and full atomic populations, matching atom/electron
counts and charge closure≤1e−5e. Native admission also requires all original
printed SCF residuals within their printed criteria and property/stdout energy
agreement≤1e−9Eh. No successful-SCF claim from normal termination alone.

For each metal/medium require absolute native–standalone energy discrepancy
≤1e−6Eh and maximum all-atom population discrepancy≤1e−5e. These are the existing
runtime-control tolerances, not fitted to this molecule. Direct corrected-native
versus standalone222-atom La controls previously passed them; the separate
standalone-versus-tblite three-body-dispersion failure is not an excuse to assume
this pair differs by the same term. Also report each metal's ALPB-minus-vacuum
component discrepancy at the same1e−6Eh resolution, without replacing a failed
absolute gate. No offsets, tolerance changes, hidden retries or automatic wider
qualification. Record incomplete/failed cells in the full eight-cell denominator.

The full profile, policy, inputs, source preparation, executable/helper/runtime
pins and script are bound before calls. This small real-core test qualifies only
its declared implementation/state, not transferable lanthanide affinity. Any
later registry explicitly limits elements, state and media and cites real results.

## Execution and follow-up

Use shared standard Slurm, at most four independent one-thread endpoint workers;
no memory/time override and no computation on login. Count all eight attempted
calls, startup failures, wall time, allocatedCPU-s and observed memory. Completion
wakes are mandatory. Reuse successful cells if only reporting needs repair.

Passing permits a separately declared whole-peptide solvent/force scout. A failed
gate requires a concrete numerical/representation diagnosis before new chemistry;
the independent nativeMACE representation result can still be interpreted. Do
not resume the old14-ion LBT3 campaign, change PQQ, or label this an affinity test.
