# Failed all-electron Dy–water: bounded recovery assessment

**17260 failed electronic convergence.** The missing `endpoint.engrad` is a consequence of the native abort, not a collection-only defect. One solver-only recovery is defensible; this assessment executes none.

## Actual evidence

The original four-atom, charge +3, multiplicity-6 calculation used unrestricted PBE0-D4, scalar ZORA(MP), SARC-ZORA-TZVP on Dy, ZORA-def2-TZVP on water, SARC/J, RIJCOSX, DefGrid3 and VeryTightSCF. Native output confirms 73 electrons, 164 basis functions and the automatic one-center ZORA treatment required by this gradient calculation. The RHF-to-UHF warning is the expected open-shell switch, not evidence that a restricted wavefunction was calculated.

DIIS/SOSCF approached a small residual, then AutoTRAH activated after iteration 53. The saved output contains 300 macro rows (0–299), including 109 Newton–Raphson (NR) rows and 16 explicitly rejected TRAH rows. Its smallest macro orbital-gradient norm was **9.681655e−6**, already at the first TRAH row, versus the required **2e−6**. The final norm was **9.875006e−5**; the last 20 norms span 6.788516e−5–1.982739e−4. Thus additional iterations did not show sustained convergence. The final NR micro-solver residual of 1.650990e−6 is not the macro orbital residual and does not establish SCF convergence.

The final NR update predicts a positive energy change of about 1.875025e−6 Eh; late macro rows repeatedly switch between NR below the configured **1e−4** entry threshold and TRAH above it. This motivates testing that switching decision. It does not prove it is the only cause: even the first TRAH step was rejected after predicted/actual energy changes disagreed.

Native output explicitly aborts after 5291 cycles, with no final energy, orbital populations, spin expectation or gradient. The output exactly matches the health-review snapshot. Overlap eigenvalues show no discarded linear dependency (minimum 6.303e−4, cutoff 1e−7); no explicit OOM is reported. Scheduler cost was **10,041 s × 4 CPUs = 40,164 allocated CPU-s**, measured TotalCPU 38,512 s; batch MaxRSS was 360,896 KiB (a scheduler field, not a summed-rank measurement). The adapter's returncode 0 and `failed_collection` do not override the native failure.

## What this says about the electronic state

The intended scalar-relativistic sextet sector was attempted, but no stationary electronic solution was obtained. Multiplicity 6 alone does not establish the Dy-localized 4f occupation, low-energy determinant, spin purity, or correct ground state. No converged spin/orbital analysis exists. The failure neither invalidates this spin sector nor validates the frozen-f reference. Spin–orbit levels, thermal populations and basis/core-treatment uncertainty remain unresolved. Even a successful AE/ECP comparison also changes basis and relativistic treatment, so it cannot isolate an f-electron contribution.

## One proposed recovery, for a separate root release

Run **one fresh versioned Dy–water EnGrad calculation** from the identical PModel guess and geometry. Retain every method, charge, multiplicity, basis, grid, tolerance, resource and iteration setting; add only this nested setting inside the existing `%scf` block:

```text
trah
  QuadRegionStart 1e-6
end
```

This delays the NR switch until below the unchanged TolG=2e−6 convergence target, keeping the trust-region solver active in the observed oscillatory region. ORCA documents this parameter as the NR-entry gradient threshold. Its value changes the optimization algorithm, not the Hamiltonian, occupations or convergence tolerance. This is a targeted hypothesis with no guarantee of convergence. [Official ORCA 6.1 SCF documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html#input-parameters).

Retain MaxIter300, four MPI ranks, one thread/rank, default scheduler memory/time and the existing actual-memory MaxCore rule (the original grant yielded 3000 MiB/rank). Preserve the original failed directory and do not seed from its unconverged GBW. One endpoint only; no automatic extension or second fallback. Input admission must prove the sole scientific-deck delta above and native startup must confirm it.

Admission still requires native SCF convergence with its active orbital-gradient criterion, normal termination, complete finite four-atom gradient, original energy/gradient closure and state/basis checks. Retain the actual spin expectation, occupations and spin populations if printed; missing analysis cannot become a physical-state pass. A different converged determinant remains possible even with unchanged settings. Only after the other matched endpoints are admitted can the balanced water-to-donor comparison be assembled. This endpoint alone cannot qualify forces or physical selectivity.

Changing multiplicity, core treatment, basis, functional, spin–orbit Hamiltonian, or using nonzero smearing would be a different state/model test. Increasing the iteration cap alone does not address the observed cycling. None is included here.

## Provenance and scope

`ASSESSMENT.json` records full source pins, native line locations, the macro summary and saved-only cost. The active plan's Stage 2 requires SCF and state qualification separately. No native input changed; zero molecular calls, submissions, restarts or new reserved-validation outcome reads.
