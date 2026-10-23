# Saved six-vacuum CHOLMOD comparison

**All six sufficient factorization tests and all six inherited dense criteria pass.** The declared execution completed once with six factors, six dense eigensystems, six solves, and no missing rows. No source or threshold changed after release.

| Saved system | Dimension | Dense min / max | Sufficient U | Shifted dense margin |
|---|---:|---:|---:|---:|
| aquo9_central_La__dimer_01 | 12 | 0.619830643 / 1.51295724 | 1.83091255 | 0.619828812 |
| aquo9_central_Dy__dimer_01 | 12 | 0.61984503 / 1.50612937 | 1.81342809 | 0.619843216 |
| aquo9_central_La__near__metal_two_water | 21 | 0.617141581 / 1.61079402 | 2.09454874 | 0.617139486 |
| aquo9_central_Dy__near__metal_two_water | 21 | 0.61723319 / 1.60136011 | 2.06331967 | 0.617231126 |
| aquo9_central_La | 84 | 0.611100629 / 1.81559327 | 3.08699421 | 0.611097542 |
| aquo9_central_Dy | 84 | 0.611302696 / 1.81219188 | 2.69623154 | 0.6113 |

Maximum factor/solve backward errors: 4.00970621e-16 / 2.09322858e-16, below every unchanged100*n*machine-epsilon limit. Maximum saved-spectrum replay difference: 3.33066907e-15, below the inherited arithmetic bounds. Raw normalized antisymmetry is5.551115123125783e-17 for each case; reconstruction is checked against the original full matrix, without averaging or diagonal/eigenvalue repair.

The factorization retains the full couplings that the earlier block-norm compression discarded. On these admitted vacuum operators, even the conservative infinity-norm upper bound leaves the shifted matrix comfortably positive. This supports the numerical backend path on these six fixed matrices through84 coordinates. It does not establish sparse fill, cost, definiteness of the isolated periodic real-space part, full-PME parity, or stability of either approximately35,874-coordinate box.

No insufficient case occurred here. In general, failure to factor B−1e−6 U I would remain a failure of that sufficient test, not proof that the full relative-spectral criterion fails. The separate manufactured indefinite-downdate behavior remains relevant: return success alone cannot admit positive definiteness.

Cost: successful computation body0.811192s wall/0.713304s CPU; full child process1.551031s outer wall/1.237267s child CPU, maxRSS89,120KiB. One numerical thread, zero molecular/native/PME calls, no saved-box calculations. All stdout/stderr and the process receipt are retained.

## Evidence

- Result: `5e92e23a92208a6e33ae4b2abb5730df3ed108cc81048584439e6625a037d431`
- Process receipt: `a5f6f4bf2c6bc8593620abb1a97ab70d5a2146953be0d4f426532a3cfce4a44e`
- Frozen configuration: `38e18214d8fe385816379cb0dbe32bf14694a7956104cdc25420ff8a03e747f7`
- Root release: `47fab4b25d5ca7b6d22ce2d225c8420c5ccd03ba698c415cf81b2fc790438be5`
