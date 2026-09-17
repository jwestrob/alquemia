# One bounded source-sampling follow-up

Declared after DISTRIBUTED_SOURCE_REPORT.md and before any new fit or spatial
validation output. Active-goal authorization applies. V1 failed6/8endpoint
field screens; its coefficients, output and tests remain immutable.

Question: does sampling the actual density at two nearby sets of positions
reduce the gap between atom-center fitting and nearby-field accuracy, without
adding multipole orders or tuning the fit? The observed training error1.5–2.4%
versus separate-probe8–14% motivates a sampling test; no biological labels enter.

Reuse all8normalized vacuum endpoints and their original physical source sites,
formal charges and CHELPG priors. Keep the exact charge/dipole basis, ell=1Å,
unweighted potential/field rows, charge constraint and
lambda=1e−4*sigma_max(A_nullspace). No quadrupoles, regularization sweep,
per-case model selection, source geometry, protonation or water change.

Training now concatenates the actual atom-center potential/field observations
and the native observations at (+0.17,−0.11,+0.13)Å from1201022. Use the native
values, never the previous fitted predictions. Those former validation points
are explicitly consumed development data in V2. Duplicate each atom's alpha
only for reporting the two training observations; fitting remains unweighted.

For new separate spatial validation, use the fixed opposite displacement
(−0.17,+0.11,−0.13)Å from each original environmental atom. This samples the
other side of the atom-center neighborhood; it has not been evaluated before.
Do not change it based on outputs. Pack its center plus the same12derivative
offsets (0.001/0.0005bohr) into8newnative orca_vpot calls, one per endpoint.
No refitting to these new values. All nuclear/source coordinates remain fixed.

All V1acceptance criteria remain: per-endpoint and Ca−La field screen weighted
RMS≤1e−4au OR relative≤10%, potential RMS≤0.005au OR relative≤10%, numerical
field≤1e−6au, paired U0 numerical difference≤0.01kcal, paired U0error≤1kcal,
charge closure≤1e−9e and rigid replay≤1e−8au. Report training versus new spatial
validation, coefficient magnitudes and both old monopole comparisons on the
same new probes. U0 is still not an environmental energy or score.

Exactly8newfits and8native potential utilities; zeroSCF/DFT, MACE, force-field
energy, new population analysis or optimization. Same8CPU/16GB shared allocation,
one thread per worker. V1cost244job seconds/819.184CPU seconds; training matrices
double in row count but retain their columns. Count all costs/failures. No
project time/CPU limit. Preserve pinned implementation, explicit source receipts,
task locks and restart history under a new workspace/protocol ID.

This is the one denser-sampling follow-up, not an open-ended parameter search.
A fail remains a fail and does not trigger an automatic regularization sweep,
quadrupole expansion or geometry change. A pass still qualifies only the
source representation on these spatial probes. Native source damping/covalent
scaling, boundary charge closure, full-cavity energy subtraction, partition
behavior and biological utility still need coherent independent tests.
