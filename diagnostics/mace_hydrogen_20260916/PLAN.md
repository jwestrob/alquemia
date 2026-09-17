# Correct stretched protein hydrogen bonds, preserving chemical states

Approved under the active MACE discriminator goal. The read-only bond audit of
the existing 1H4I physical state found mean deviations from the already used
Amber ff19SB equilibrium lengths of +0.1033 A (3459 C-H), +0.1767 A (916 N-H)
and +0.2247 A (92 O-H). Protein heavy-bond RMS deviations are much smaller.
Both models' traced charge-normalization denominators are well conditioned on
the observed La endpoints; no near-singular denominator has been identified.

Question: does an independently defined, uniform protein-hydrogen preparation
correction reduce the full-system model disagreement and excessive force
response? This is a consumed development test, not a new blind validation.

For every protein H, use actual bonded connectivity from the existing ff19SB
topology, require exactly one bonded heavy atom, and project its bond to the
force-field equilibrium length along the existing heavy-to-H direction. Apply
the same rule everywhere. Preserve all heavy atoms, PQQ atoms, metal position,
protonation/charge, assembly, water inventory and ordering. Do not minimize,
change donor membership or repair heavy-atom deviations in this experiment.
This independent standard-protein bond rule is not fitted to a score or label.

Prepare one common physical geometry and paired La/Ca XYZs. Evaluate the same
analytic medium and large checkpoints, each on the full La/Ca pair: **four MACE
calls, zero new DFT**. Retain read-only charge tracing. Compare raw contrasts,
charges, forces and update stages to the archived uncorrected geometry. There
are no matching DFT endpoints for the changed core hydrogens: hybrid scores,
aquo-normalized S and classifications remain unavailable.

Check paired coordinates/charge parity, exact preservation of heavy/cofactor
atoms, single attachment for each H, target bond lengths, orientation, source
mapping and rigid-transform covariance of the coordinate operation. Inference
uses the already verified analytic kernel without altered weights, cutoff,
precision or electronic model. Charge closure still requires 1e-5 e. Report all
effects rather than changing criteria after seeing them; no claimed success
from moving a score toward a desired label. Geometry changes define new IDs.

Existing A5000 allocation policy: one GPU/16 CPUs/64474 MiB host RAM. Matched
expected model time is approximately 2+4 minutes, plus preparation/startup/trace
overhead. No artificial stopping budget. Preserve the original baseline and
all source and experimental records. This is an input-preparation experiment,
not a validated relaxation or thermodynamic correction.
