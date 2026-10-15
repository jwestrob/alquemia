# Canonical rebased geometry package

Metadata freeze is complete. Python AST and shell syntax checks passed; Heisenberg independently cleared the frozen source with no material blocker. Numerical geometry, Jacobian checks and the new exact A0AC rebase regression are **UNRUN**. No job was submitted by this agent.

Root's one-CPU geometry-only command, from the repository root:

```bash
sbatch workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_relaxation_v1/implementation/RUN.sbatch \
  "$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_relaxation_v1/manifest.json" \
  4d7a83d44c0cf3b7631adfa334d273bc0db17389062b5908936e3a12b829415c
```

The wrapper requests one CPU on standard, constrains numerical libraries to one thread, and has no host, memory or walltime override. Root owns submission and completion wake. This command cannot perform an energy evaluation or search.

| Artifact | SHA256 |
|---|---|
| Frozen manifest | `4d7a83d44c0cf3b7631adfa334d273bc0db17389062b5908936e3a12b829415c` |
| PREPARE.py | `be3356652accc7ce0d413180cd05eb465e5a5b716106b1adcaac107728ade0ec` |
| RUN.sbatch | `156e5772ba0b8673612ddff031e9c0f15c44d46715b57b8e3266985b44b904ae` |
| PLAN.md | `1359be5d59ac288fd2b3cdb0fe002814a88097df7c40ff77ba211fdd5d638861` |
| Actual 4936 fixed-flip preparation | `5f69e8c7d57a71d6601d7ced65e8dcc419c8ad2ac154f755141b6faaf395c968` |
| Actual prior A0AC rebased preparation | `6405de39b9d3ffa88aea7f99f895a685ddcd8ec5b1b37529766b72c365d9d57d` |
| Completed A0AC search collection | `2ed05f5a241f4f1eccbb2a680428a4b032a51449df56ad323e138b83cf32fd82` |

## Scope and reuse

All 25 canonical sources remain in their original order. Context counts, charges, cap inventories and complete source-union sizes are read from the exact inputs. Original four-mode membership/bounds remain unchanged, including sources whose terminal χ2 is not selected. No optimizer, q0 energy or force gradient is fabricated.

The metadata join confirmed the exact A0AC La-conditioned sample 3 source: original tasks/maps/state, six baseline MACE receipts/XYZ pins, and both flipped XYZ hashes match the executed branch. The job must then pass a numerical geometry-only regression against its existing rebased context/checks, cap residuals and full-source union. It exports the actual prepared mapping/contact/task/XYZ and completed search-result aliases; repeated A0AC energy/search work is disallowed. The other 24 sources remain explicitly unqualified for energy/force q0.

Outputs are `execution/PREPARATION.json` and `EXECUTION.json` beneath the frozen workspace, with all 25 rows and explicit failures. Each admitted row records q0 geometry, Jacobian/bond/cap/source admission and original mode membership. These checks do not establish a clash-free or equilibrium conformer, validate the full search domain, or produce scores. Subsequent energy work remains a separate decision based on the ongoing solvent result.
