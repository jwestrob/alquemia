# Uniform native CHELPG check: promising, not promoted

All four standalone ORCA6.1.1 utilities completed at the [predeclared native
defaults](CHARGE_FIT_PLAN.md), with verified charge sums and atom order. No new
SCF calculation or biological comparison was performed. Against independent
quantum potentials, CHELPG improves exterior RMS and direct coupling agreement
on these four consumed development states. No scheme was selected per case.

| State | CHELPG wall s, one CPU | MBIS exterior relative RMS | CHELPG exterior relative RMS | MBIS − exact coupling, kcal/mol | CHELPG − exact coupling, kcal/mol |
|---|---:|---:|---:|---:|---:|
| 1h4i_qm33_La | 58.889 | 0.053225 | 0.025069 | +2.241769 | +1.069303 |
| 1h4i_qm33_Ca | 58.593 | 0.028505 | 0.014112 | -6.616972 | +0.949101 |
| 1h4i_qm36_La | 73.329 | 0.036323 | 0.013510 | -5.295331 | +0.991340 |
| 1h4i_qm36_Ca | 73.130 | 0.018564 | 0.009234 | -6.093414 | +0.463666 |

All four CHELPG fits pass the historical exterior-potential rule. Their paired
Ca-minus-La coupling errors are −0.120202428 and −0.527673470 kcal/mol for
qm33/qm36. The corresponding MBIS errors are −8.858740778 and −0.798083137.
This independently checked interaction error matters more than a fitting score
alone. Residual error and orientation sensitivity still need attention before
using any fitted-charge model in a global score. None has been promoted.

Prior native MBIS population sections took 574.688–747.686 seconds per vacuum
endpoint, versus standalone CHELPG58.593–73.329 seconds here. Different nodes,
MPI allocation and timer categories prohibit a controlled speedup ratio. Actual
Slurm allocation and per-utility CPU/RSS records are preserved; this is a useful
engineering candidate, not proof of production cost or improved affinity.

Artifacts: `workspaces/density_embedding_20260916/chelpg_v1/` contains the
manifest, copied wavefunctions, native utility logs, execution receipts, and
`comparison_v1.json`. `utility_accounting_v1.json` in the parent workspace also
includes the original failed density utility calls and their same-input retry.
Six real-artifact/parser checks pass. Baseline/default and old MBIS results stay
unchanged. The separate four-field-endpoint response test is still running.
