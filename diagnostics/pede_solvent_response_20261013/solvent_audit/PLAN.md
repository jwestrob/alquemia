# Exact saved solvent-response decomposition

Read all nine already-analyzed PedE sources (five two-Ca, four compatible La),
both metals, vacuum and ALPB. No new molecular calls, charge analysis, optimizer,
scorer or calibration. Root submits this finite one-CPU Slurm analysis.

Question: does the sample3→sample1 solvent double difference arise from the
printed electrostatic, surface, hydrogen-bond, reference term, or the remainder
outside Gsolv? The same-extrema comparison is already specified by the completed
energy analysis; it does not filter the full population.

Use actual xTB summaries and the identity, for each metal/source:

`U = E_ALPB − E_vac = Gelec + Gsasa + Ghb + Gshift + I`

`I = E_ALPB − Gsolv − E_vac`.

Gelec and Ghb use the converged ALPB charge distribution. I is the complete
remaining intrinsic-energy change, including electronic-temperature contributions
and density-dependent terms. Do not call it pure polarization relaxation or
discard it to construct a new score. The summary components already belong to
the total; no second Gsolv is added. Keep printed SCC/repulsion/dispersion terms
for accounting, not a second additive decomposition of nested quantities.

Before execution, closure tolerance is 1e-9 Eh (printed entries have twelve
decimal places). Verify exact saved identities/pins, matching medium geometries,
charges −3/−2, 588 explicit valence electrons, and the unchanged complete source
population. Missing required components fail visibly; absent separately printed
entropy is unavailable, never zero. Retain actual temperature/reference/warning
lines and the original state/cap/charge records. No admission rule is relaxed.

## Sources checked

- [Official ALPB documentation](https://xtb-docs.readthedocs.io/en/latest/gbsa.html):
  Gsolv is already included in total energies; the printout splits electrostatic,
  surface, hydrogen-bond and reference terms.
- [xTB6.7.1 Born implementation](https://raw.githubusercontent.com/grimme-lab/xtb/v6.7.1/src/solv/gbsa.f90):
  `getEnergyParts` uses the current atomic charges; Ghb is charge-dependent.
- [xTB6.7.1 SCC implementation](https://raw.githubusercontent.com/grimme-lab/xtb/v6.7.1/src/scc_core.f90):
  self-consistent electronic energy includes solvation and Fermi-temperature terms.

Exact downloaded source hashes are under the owned workspace `upstream/` and
pinned by the execution configuration. Source inspection does not establish that
a model is accurate for these charged metal/cofactor fragments. The independent
DFT reference branch will qualify the response before choosing a correction.

## Invocation and output

`DECOMPOSE.py --config <frozen CONFIG.json>` requires real Slurm before analysis.
Root uses one standard CPU; no memory, walltime, account or QOS overrides.
Outputs are `RESULT.json`, `COMPONENTS.tsv`, measured receipt. New molecular
evaluations and population recalculations are both zero. Do not rerun a completed
output directory. A separately named failure/recovery would be required for an
actual parser defect; no scientific settings change follows automatically.
