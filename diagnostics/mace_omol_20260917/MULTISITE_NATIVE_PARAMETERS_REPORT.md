# Multisite preparation: actual ACE, waters and background Ca supported

**All five native framework parameter checks pass.** This is preparation
support, with no new metal score, quantum endpoint, MACE inference, energy,
force or induction solve. The failed GGR transfer remains failed.

| Prepared source/site | Physical atoms | Background Ca | Retained waters | Framework charge excluding selected metal |
|---|---:|---:|---:|---:|
|4CPV CD|1611|1|1|−5|
|4CPV EF|1611|1|1|−5|
|1SL8 EF1|2866|2|3|−6|
|1SL8 EF3|2866|2|3|−6|
|1SL8 EF4|2866|2|3|−6|

Actual source coordinates/protonation, all source atoms, the ACE0-C–ALA1-N
bond and the declared site-water union are retained. The selected QM metal
remains in the physical ledger; only its force-field parameterization is
pending. Native framework/background atoms retain real identities, permanent
multipoles and polarizability. Background Ca uses the installed type358:
charge+2, polarizability0.55A³, GK radius1.82485A. No La parameter was guessed.
Maximum native-versus-OpenMM charge/polarizability discrepancy2.22045e−16.

## Compatibility fixes and preserved failures

v1 rejected the original selected-metal source identity because the old
single-metal interface expected the alias `metal`; multisite preparations
retain actual source IDs. The opt-in now accepts that actual identity.
All five v2 preparations matched templates but failed in OpenMM's unused GK
builder: its default Bondi radius table lacksCa. No radii were invented and
no installed package was changed. v3 explicitly leaves OpenMM GK unavailable,
uses its qualified AMOEBA2018 permanent/polarization parameters, and obtains
actual GK parameters from the existing native Tinker SOLUTE backend. These
are the parameters already used by the hybrid's scientific energy solver.

The legacy single-metal helpers remain the default. New explicit flags admit
background calcium and native GK preparation; real non-Ca background atoms
still fail. Native parameter validation now recognizes declaredCa atomic
number20. Source support/core mappings, quantum densities and full hybrid
scores on these five sites remain unavailable; preparation is not accuracy.

Protocol: `source_multisite_background_Ca_AMOEBA2018_native_GK_capability_v2`.
Plan: [native recovery](MULTISITE_NATIVE_PARAMETERS_PLAN.md), following the
[initial capability scope](MULTISITE_AMOEBA_PREPARATION_PLAN.md).
Products: `workspaces/mace_omol_20260917/multisite_amoeba_capability_v3`;
explicit frozen configuration: `multisite_amoeba_config_v2/config.json`.
Both failed versions and their per-case reasons/receipts remain available.

v3 preparation83.401715wall/79.253806parentCPU seconds; five native children
3.283880wall/3.239334CPU (childwall already included in preparationwall).
v1 preflight0.947469wall/0.830171CPU; v2 failures8.529910wall/8.375095CPU.
No cluster/GPU allocation. Four real/corrupted-fixture tests pass3.216s;11legacy AMOEBA/native framework
regressions pass11.094s. No energy integration test is claimed.

Next: map the real site cores to these unchanged whole physical systems,
including selected-site versus exterior waters and backgroundCa. Declare
full score runs only after paired state and core/source accounting work.
Parvalbumin remains supporting cross-study evidence; aequorin stays the
orderedEF1/EF3/EF4 vector with no invented site-resolved assay labels.
