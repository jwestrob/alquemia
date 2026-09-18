# GMRES solves the equations, but does not earn a runtime expansion

Both real coarse2FW0 endpoints completed with the unchanged native PCM
operators. The independently reevaluated equation residuals pass. Eighteen
of21frozen checks pass; the complete qualification is **FAIL**. No tolerance
was changed, no wider GMRES run is launched, and no biological score follows.

| Quantity | Ca | La |
|---|---:|---:|
| Energy, kcal/mol | -83.93538012794836 | -40.850965390325555 |
| Difference from native result, kcal/mol | -1.30901978e-6 | +1.23540573e-8 |
| Dielectric true relative residual | 3.77545952e-11 | 4.48092325e-11 |
| Single-layer true relative residual | 7.08353058e-11 | 8.18481513e-11 |
| Dielectric iterations | 604 | 506 |
| Single-layer iterations | 78 | 76 |
| Dielectric wall seconds | 490.249037 | 411.107083 |
| Single-layer wall seconds | 4.350987 | 4.048632 |

The three failures are the saved native Ca composed-equation residual
(2.59814206e-8versus1e-9), Ca energy agreement(1.30902e-6versus1e-6kcal), and
paired contrast agreement(1.32137e-6versus1e-6kcal). These tiny energy differences
are far below the physical numerical-convergence scales, but the declared
engineering gate remains failed. The earlier coarse reciprocity error0.122802kcal
also remains; changing the linear solver does not repair it.

The useful finding is the cost split:901.356120seconds in the finite-dielectric
system versus8.399619seconds in the single-layer system. More than99%of the
two-system solve time is spent in the former. This motivates the separately
declared conductor-like approximation in DDX_CPCM_PLAN.md, chosen without any
new affinity classifications. It is a different physical approximation, not
a dropped term reported as exact PCM.

## Executions and limitations

Kernel build1201295 completed in89allocated wall seconds:5696allocated core-s,
76.976reported CPU-s,82012KiB sampledRSS. All19native numericalFortranmodules
remain byte-identical; a thin access adapter was added in an isolated build.

Qualification1201296 completed2forward endpoints/4linear solves. Slurm945wall
seconds,64CPUs,60480allocated core-s,58042reported CPU-s(whole-second precision),
478588KiB sampledRSS. Driver921.696707wall/58020.4547CPU seconds and455476KiB
peakRSS. No GPU, DFT, MACE, forces or native iterative rerun. Local preparation,
imports and tests are separate;3actual-source/preparation tests pass29.197s.
This is not a matched end-to-end production benchmark.

Manifest: `workspaces/mace_omol_20260917/ddx_krylov_v1/manifest.json`, SHA256
`69a7e063bd5abe511ac7d266fdaae0972c7a4faaf6a4a029642eb06d67d740f1`.
Actual collection: same directory,`collection_job_1201296.json`. Native vectors,
solutions, true residuals and callback histories remain under`attempt_0001`.
Read-only summary and accounting: `ddx_krylov_report_v1/result.json`.
Original native recovery jobs remain undisturbed. Baseline/default unchanged.
