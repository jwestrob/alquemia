# One-cell canonical CC DIIS recovery — source only

Question: does a larger DIIS history converge the original failed canonical
Dy +0.005 Å-equivalent rotation cell within its unchanged 100-iteration budget?
Only `cc__Dy__rotation__plus_005__complex` (failed job18866) is eligible.
Its failed raw record remains immutable and unavailable as a CCSD(T) energy.

The sole intended native algorithm change is adding `MaxDIIS 12` to `%mdci`;
the original native output printed7. Keep MaxIter100, STol1e−8, the actual
LShift0.2, RHF/CCSD(T)/NoFrozenCore, VeryTightSCF/PModel/MaxIter300, the exact
eight-center geometry, charge+2/effective singlet, Dy ECP55/40 explicit and
correlated electrons,189 orbital functions, original basis and executable.
No old orbitals/amplitudes are copied; retain NoAutostart. MaxCore alone derives
from the actual default allocation: floor(0.75×memory_MiB/4). Four single-thread
MPI ranks on one node, no memory/time override or extra solver fallback.

One new native invocation at one fixed new output path; no automatic retry,
other-cell release, collection replacement, fitting or fresh validation. A root
release records exact source/manifest/preflight/input/output scope; this is
machine-readable provenance, not another user permission request. No release
is supplied by preparation. No watchers are launched.

The current v2 canonical deck/parser and its JSON-canonical raw comparison are
reused. All original native convergence, state, basis, coordinate, final-energy
and triples-closure checks remain. Also require exact native printed MaxDIIS12,
MaxIter100, STol1e−8 and LShift0.2. Failed native output retains available raw
SCF/CC trace and solver settings, with no substituted unconverged energy.
Original SCF flags, T1/single-reference/frozen-f and finite-basis limitations
remain; even successful one-cell recovery cannot qualify the full eight-cell
rotation panel or physical method accuracy.

Zero-call preparation verifies original failure/raw pins, source/panel/cell joins,
exact one-line deck delta, protected output path and one-call cap. Focused
negative checks cover unauthorized cell/output and changed/absent readback;
no broad prior suite or molecular call is repeated. Report actual preparation
cost; eventual root review must account native allocation and failure costs.
