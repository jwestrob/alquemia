# Exact two-forward factorization of the charge-masked descriptor

Declared2026-09-17 before this audit. Engineering/algebra under the active goal,
using already computed development outputs; no new molecular inference. The
running100-call canonical experiment, its criteria and its inputs are untouched.

For the masked short-range model, disconnected components receive the same
fixed spin feature and zero charge feature. Installed ScaleShiftMACE sums node
readouts plus the explicitly captured embedding readout. There is no additional
long-range charge term in this model. Thus, if the separated metal has no edges:
T_detached,M = T_protein + C_M,
where C_M is the detached metal node plus its embedding readout. Paired protein
coordinates and species are identical, and the global charge feature is masked.
Consequently the existing four-forward descriptor should equal
R_mask = [T_bound,Ca - T_bound,La - (C_Ca - C_La)] * 23.06054783061903.

C_M is a model readout term, not a measured/computed quantum isolated-ion energy,
aquo correction, fitted intercept or biological calibration. No physical affinity
zero follows. This would be an exact computational factorization of the existing
modified descriptor, not a new scientific model or repaired physical Hamiltonian.

Use the detached primary1H4I Ca/La metal readouts from completed job1200828 as
the fixed reference; that anchor is selected before inspecting this audit.
Verify actual receipt, readout closure, zero metal edges, identical nonmetal
coordinates/state, unchanged mask/weights and charge/spin records. Test every
completed detached variant in that42-task run, including alpha repeats/rotation/
farther displacement and GGR sodium states. Compare each detached metal term
with the fixed anchor, and compare paired nonmetal readouts. Test all five
primary scores and the sodium score against direct four-forward algebra.

Freeze acceptance at the existing0.01 model-kcal energy/contrast tolerance.
Retain unrounded component terms and errors. No favorable reference selection,
averaging, fitted offset or threshold change. Failure means no two-call path.
If it passes, implement an explicit opt-in factorization setting in the prepared
input interface. Validate by reusing actual bound endpoints without any new
forward, retain the four-call path, and type the reference/cache identity.
After canonical job1200830 completes, verify the same algebra across its25
calibration cases and two supported transfers before claiming panel-wide
qualification. This is numerical equivalence, not another accuracy observation.

Outputs: source-pinned factorization report/reference, precise algebra and
component accounting, interface tests and measured local cost. New model/DFT/
solver/gradient/training calls: zero. Do not cancel or shorten existing jobs.
