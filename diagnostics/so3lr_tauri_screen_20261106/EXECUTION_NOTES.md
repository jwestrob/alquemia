# Execution notes

Before any model prediction, GPU scout49306 remained PENDING with
QOSMaxGRESPerUser; other-session dnv1_fold45665_1 remains running. Preserved the
scheduler evidence and cancelled only our unstarted49306. Static CPU evaluation
on standard_expand uses identical models, precision and cutoffs; it cannot supply
a measured H200/bulk-solvent throughput. Separate CPU output directories prevent
mixed-backend artifacts. No resource limit or other user's allocation was changed.
The native no-cutoff API requires dispersion damping=None when lr_cutoff=None;
source-level validation corrected that argument before any Tauri evaluation.
Six reference49309–49314 inputs/settings are frozen in REFERENCE_MANIFEST.json;
no duplicates/retries/root restart. OLD_REFERENCE_AUDIT.json reproduces+18.51456
and−9.03750 from native gas/sc_gas outputs; those methods/geometries/gradients do
not provide the physical-spin comparator at the current24 fixtures.

CPU scout49325 evaluated its first real frame, then failed on a Python Path/string
concatenation in the NPZ writer, before saving any predicted energy/force. The
failed allocation/logs remain; those values are unavailable, not inferred. The
writer now builds the path in one expression and snapshots executed source/panel/
environment at startup. Recovery uses a fresh v2 directory, identical scientific
settings. This is an output failure, not a potential failure.

Saved reference check49393 parsed the first completed LaW8 energy/gradient,
physical counts, occupations and positive stability eigenvalues. Its generic
1e-8-A coordinate comparison was finer than ORCA's seven-decimal Bohr coordinate
export. The collector now derives half a printed coordinate unit plus2e-9 A for
conversion convention; exact source/input coordinates are still asserted before
this export check. This fixes a file-precision check, not a scientific gate.
The snapshot remains preserved. LaW8 still fails three printed Tight residuals.

Full CPU panel49387 completed:86 saved rows, including the two native scout rows
reused by hash; all48 primary L/M predictions and all six implementation checks
pass. Analysis49392 uses saved arrays only. H200 throughput remains unmeasured.

Reference export check49416 revealed a second unit detail: ORCA uses the legacy
Bohr-to-Angstrom constant0.5291772083, documented in its6.0 manual and confirmed
by6.1.1's twelve-decimal property coordinates. The collector now verifies those
native property coordinates against the archived XYZ to1e-11 A, then verifies
the lower-precision engrad tail against the property file at its printed rounding.
It uses the same native constant for force conversion. This resolves both
format/conversion artifacts without changing geometries or scientific criteria.
Source: https://www.faccts.de/docs/orca/6.0/manual/contents/general.html#units-and-conversion-factors


Second writer failure49347 saved both native W8 predictions and calculated the
ASE comparison, then failed on the equivalent Path/string expression in its JSON
filename. That expression was fixed. Recovery49384 reused the two native NPZs by
hash, recomputed the ASE check and completed in a fresh v3 directory. No potential
settings changed. The subsequent full panel reused those two native predictions.

The original water REFERENCES/FORCES/CYCLES manifests and the saved donor
materializer/original radial-rotation generator were reviewed. Historical targets
remain the original uncorrected ECP-PBE0-D4 forces/energies, with native states
explicit; no CP/correlated fragment labels or force-field operators substitute.
AcW7 readouts apply the same original acetate axes and rigid OHH shell directions
to its unchanged groups; this is saved-force analysis, not another calculation.

Root declared SCIENCE_STOP before cancelling only49311/49312/49314. All three
completed La references had already failed the mandatory frozen Tight residuals;
both balanced contrasts require LaW8, so unfinished Dy could not repair screen
qualification. All six were invoked, three completed, three stopped after2686s
each. Their intermediate SCF energies are not final energies or fitted targets.
The1-kcal contrast budget was not applied without qualified references. No retry,
solver/method/threshold change, replacement endpoint, MD or additional chemistry.

Final reference audit49450 and comparison49464 are saved-only. Reproduction
49465/49466/49467 preserves originals and reproduces all five numerical tables,
reference rows and decisions. All23 jobs, including failures and cancellations,
are in COMPUTE_COST.json. The two completion-watch receipts have been consumed;
delayed queue messages are acknowledgements only, with no job restart.

Final review normalized committed TSV line endings to LF and set the writer
accordingly. The verified original/replay CRLF exports are preserved and their
canonical LF contents match exactly. No numeric values or raw arrays changed.
