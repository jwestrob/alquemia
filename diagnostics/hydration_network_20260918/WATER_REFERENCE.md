# Bulk-water reference: implemented and computed

This is an independently specified reference, not an occupancy model or a new
production aquo reference. The site calculations still need bound-state free
energy terms. Native r2SCAN-3c gas-water energy is tied to a converged gas-water
geometry; a separate fixed-geometry CPCM water calculation is retained only as
a comparator. Both use TightSCF/DefGrid3 and singlet neutral H2O.

## Energy convention

At 298.15 K and pure liquid water at 1 bar, approximate:

`mu_liquid = E_gas,min + ZPE_gas + [Hgas(T)−Hgas(0)] − T*Sgas(T,1bar)`
`            + RT*ln(psat/1bar) + V_liquid*(1bar−psat)`.

This uses ideal vapor and an incompressible liquid pressure correction. Tabulated
harmonic frequencies supply an explicit harmonic ZPE approximation, not an exact
anharmonic ZPE or a new quantum Hessian. Source values and units are pinned in
[WATER_THERMOCHEMISTRY.json](WATER_THERMOCHEMISTRY.json). No terms were fitted to
protein classifications. No extra 55.34 M correction is appended: the liquid
chemical potential is already specified through phase equilibrium.

| Contribution | kcal/mol |
|---|---:|
| Harmonic gas ZPE | +13.472284 |
| Hgas(298)−Hgas(0) | +2.367113 |
| −T Sgas at 1 bar | -13.456228 |
| RT ln(psat / 1 bar) | -2.045533 |
| Liquid pressure correction | +0.000418 |

Their sum is **+0.338054 kcal/mol** relative
to the gas electronic minimum. Saturation pressure is
0.031667487 bar. Gas E=-76.418935329701 Eh;
mu_liquid=-76.418396606095 Eh within the stated approximations.
All unrounded terms/receipts: [result index](REFERENCE_RESULT.json).

## How to combine it with site states

For adding a water, `DeltaG_add = G_site(n+1)−G_site(n)−mu_liquid`.
The G terms must include the bound-water vibrational and configurational
contributions. Subtracting mu_liquid directly from bare electronic site energies
would omit these contributions. In particular the reference includes roughly
13.47 kcal/mol gas ZPE per water: its corresponding bound term cannot be silently
omitted. If a diagnostic assumes that gas-like internal ZPE carries into the
bound water, record that assumption explicitly; the corresponding electronic
reference offset is -13.134230
kcal/mol, with bound vibrational shifts and basin terms still unknown.

Missing bound-state entropy and non-electrostatic solvent terms remain null.
They are not zeros, fitted penalties, or occupancy probabilities. Site-state
energy rankings and the correction required to reverse a ranking can be reported
without claiming a complete free energy. The old aquo gauge/bands are unchanged.

## Actual execution and recovery

Job1201830 requested one CPU but Slurm's exclusive node allocation exposed64CPUs.
The existing runner inferred16MPI ranks from that allocation; MPI rejected the
launch because the request supplied one slot. Both attempts stopped during
startup, before a scientific SCF result. Original outputs/receipts are retained.
Added an optional explicit MPI-rank/concurrency setting to the existing runner;
default behavior remains unchanged. Identical scientific inputs were prepared
in a fresh directory and successfully rerun with one rank as job1201831.

Actual scheduler cost including failure:25s summed allocation wall,1,600
allocated core-seconds (64 allocated CPUs on both jobs), zeroGPU. Successful
reference work used5 gas optimization
energy evaluations and1
CPCM energy evaluation. No quantum Hessians or numerical gradients. This is a
shared reference cost, not a per-protein calculation.

Sources: [JANAF gas thermochemistry](https://janaf.nist.gov/tables/H-064.html),
[CCCBDB harmonic frequencies](https://cccbdb.nist.gov/exp2x.asp?casno=7732185&charge=0),
[NIST vapor pressure](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4),
[published water concentration](https://pmc.ncbi.nlm.nih.gov/articles/PMC2700946/).
Plain CPCM lacks non-electrostatic solvation terms; see the
[ORCA manual](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html).
