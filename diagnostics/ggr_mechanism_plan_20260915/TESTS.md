# Executed verification

These checks use pinned real structures, archived ORCA outputs and, when
available, the actual new gradient files. Deliberately malformed inputs are
identified corrupted copies of those artifacts. No scientific output is mocked.

Python: `/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python`.

| Command (from repository root) | Result | Scope |
|---|---|---|
| `python -m unittest discover -s tests -p 'test_ggr_*.py' -v` | 31 passed, 19.497 s | Preparation, exact coordinates/charge/waters, source mapping, unsupported chemistry, PQQ pins, physical Jacobians, task packaging; four actual-gradient integration checks |
| `PYTHONPATH=scripts python -m unittest discover -s tests -p 'test_affordable_benchmark.py' -v` | 4 passed, 0.494 s | Archived reference algebra, no threshold transfer, pinned preparation, corrupt/partial receipt rejection |
| `PYTHONPATH=scripts:tests python -m unittest test_affordable_development.ArchivedRegression -v` | 5 passed, 0.716 s | Released baseline energies/bands, endpoint correction signs, corrupt-output rejection, full configuration cache identity |
| `python -m unittest discover -s tests -p 'test_ggr_benchmark_export.py' -v` | 5 passed, 6.958 s | All 68 inherited rows/fields/groups/scores intact; nine appended development pairs; malformed collection and missing-endpoint handling; actual C records retained separately with null corrections |
| `PYTHONPATH=scripts:tests python -m unittest test_ggr_sensitivity.RealGradientIntegration -v` | 4 passed, 8.822 s | Repeated on all four completed real gradient centers, including the 111-atom connected model; actual atomic numbers, ECP, coordinates, components and paired mappings |

The independent preparation audit also checked all five Stage A coordinate
and hash mappings, added donor distances and nonbonded cap contacts. All 18
Stage A/B endpoint outputs passed normal termination, SCF, provenance and score
collection. Printed SCF+D4+gCP closes to final energy within the archived parser's
2×10⁻⁹ Hartree printing tolerance; no component was inserted as zero.

The initial real analytic-gradient checks are documented in
[GRADIENT_INTEGRATION.md](GRADIENT_INTEGRATION.md). Directional energy consistency
is a separate scientific integration result, reported in REPORT.md when complete;
software-test success does not substitute for it. No dense or numerical DFT
Hessian, numerical DFT gradient, relaxation, entropy or environmental calculation
is part of these tests.

The complete `rebuild_report.sh` operation also ran successfully against all
38 finished endpoints into `workspaces/ggr_mechanism_20260915/report_v1`. Its
derived ledger preserves all 68 original records, adds nine development pairs,
retains 20 separate sensitivity endpoints and links the exact 194,048 allocated
core-seconds. Four final figures were visually inspected; all 12 derivative
checks have correct gradient/amplitude units and no clipped labels. The
conditional half-step scientific branch was not triggered or executed.
