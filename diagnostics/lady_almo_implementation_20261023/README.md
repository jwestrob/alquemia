# Molecular ALMO constraint

This implements the real, closed-shell AO-block ALMO variational constraint on
PySCF's ordinary full-system RHF/RKS energy. It allows polarization inside each
fragment while preventing occupied orbitals from expanding into other fragments'
atomic basis functions. It does not replace the Hamiltonian or fit a correction.

**Core numerical tests pass; cross-code and La/Dy validation are incomplete.**
Jobs 22332/22333 verified orbital derivatives, density/support invariants,
constrained stationarity, one-fragment recovery and release to ordinary SCF.
Read [CURRENT.md](CURRENT.md) for actual jobs and unresolved comparisons.
Both bounded native CP2K attempts stopped in isolated-water SCF before ALMO;
the strict comparison rejected them. This validation round is terminal.
No La/Dy charge-transfer conclusion is admitted yet.

## Files

- [almo.py](almo.py): density, analytic orbital derivative, constrained optimizer,
  numerical invariants, and ordinary-SCF constraint release.
- [THEORY.md](THEORY.md): equations, parameterization, convergence and limitations.
- [archive.py](archive.py): exact archived La/Dy inputs and basis/ECP conversion.
- [run_system.py](run_system.py): baseline, isolated-fragment and ALMO execution.
- [tests/run_tests.py](tests/run_tests.py): independent numerical correctness tests.
- [tests/run_cp2k.py](tests/run_cp2k.py) and
  [tests/compare_cp2k.py](tests/compare_cp2k.py): native implementation comparison.
- [submit.py](submit.py): immutable source snapshots, Slurm and completion receipts.

## Library interface

```python
from almo import solve_almo

# mf is a configured real, closed-shell PySCF RKS or RHF backend.
# fragment_aos partitions every AO index exactly once.
# fragment_nocc contains the number of DOUBLY occupied orbitals per fragment.
# C0 has shape (number_of_AOs, sum(fragment_nocc)); forbidden entries are zero.
result = solve_almo(
    mf, fragment_aos, fragment_nocc,
    initial_coefficients=C0,
    gradient_tolerance=1e-7,
    energy_tolerance=1e-10,
)
if not result.converged:
    raise RuntimeError(result.status)
```

Use converged isolated-fragment orbitals for the frozen state in an energy
decomposition. A core-Hamiltonian initial guess is available for solver use but
does not define that same frozen reference. Fragment occupied overlaps are
retained; the resulting density generally has off-diagonal fragment blocks.

## Reproducible project execution

The isolated environment is
`workspaces/lady_almo_implementation_20261023/venv/bin/python`.
Its `requirements.lock`, `PIP_REPORT.json` and `ENVIRONMENT.json` retain versions
and downloaded wheel provenance. Every task has its own `CONFIG.json`, source
snapshot, Slurm output, `SUBMISSION.json` and completion receipt.

All numerical execution in this project uses `standard_expand`. For example,
from the repository root, after checking CURRENT for existing work:

```bash
python diagnostics/lady_almo_implementation_20261023/submit.py \
  --name UNIQUE_RUN_NAME --cpus 8 --memory-gb 32 --hours 2 -- \
  @PYTHON@ @SOURCE@/run_system.py --repo @ROOT@ \
  --system aquo8_central_La --output @TASK@/result \
  --threads 8 --max-memory-mb 28000 --mode baseline
```

Use `--mode fragment --fragment 0` for the metal or `--fragment 1` for the whole
shell. The ALMO mode requires `--baseline-dir` and two `--fragment-dirs`, in
metal/shell order, pointing to admitted saved results. It verifies their inputs,
runtime and numerical settings before reuse. It does not repeat their SCFs.
Existing completed work must be reused instead of rerunning this example.

## Limits

This is a first-generation finite-AO-space constraint. A successful optimization
establishes constrained stationarity, not a global minimum or a physical Dy
electronic-state assignment. The quartet retains the archived effective singlets,
including Dy ECP55; it is not an all-electron Dy sextet calculation. CT values
depend on this partition and finite basis. Neither those values nor the fixed
cluster substitution contrast are binding free energies.

The implementation currently provides combined frozen, polarization and CT
terms. It does not invent separate Pauli/electrostatic components. Archived D4
is added identically to the full and constrained state at each fixed geometry,
once, and is reported separately.
