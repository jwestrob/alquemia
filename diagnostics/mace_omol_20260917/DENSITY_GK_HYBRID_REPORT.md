# Frozen density / AMOEBA-GK / MACE: stable numerics, accuracy gate still fails

All111native tasks completed first try as job1201074. The candidate passes all
32numerical, identity and prescribed cavity-sensitivity checks. It passes only
**2/4development ordering comparisons** and fails the unchanged GGR partition
gate: **-3.783707454kcal/mol** connected-minus-extended versus abs<=2.
Keep the production baseline. This candidate does not earn calibration or
benchmark expansion; the broader MACE research goal remains active.

## Actual predictive and partition results

Larger `R=E_Ca-E_La` is more La-like. Compare relative site differences; these
raw contrasts have no aquo reference, inherited zero or calibrated bands.

| Alpha-minus-GGR comparison | kcal/mol | Frozen direction gate |
|---|---:|---|
| 1F6S minus extended | -0.437799535 | fail |
| 1F6S minus connected | +3.345907919 | pass |
| 6IP9 minus extended | -3.546267393 | fail |
| 6IP9 minus connected | +0.237440061 | pass |

Both alpha structures must exceed both GGR representations by0.02kcal/mol.
Do not select the favorable core or structure. These are two already consumed
biological groups: GGR direct same-assay direction and qualified cross-study
alpha evidence. Alpha's differing explicit water inventories remain fixed;
these are not four independent affinity observations or prospective validation.

GGR connected-minus-extended component audit:

| Contribution | kcal/mol |
|---|---:|
| Vacuum native DFT | -23.951208420 |
| Exact density-to-permanent-multipole coupling | +12.034985832 |
| GK permanent reaction-field transfer | -1.043006589 |
| Environmental induction transfer | +4.059946620 |
| MACE short full-minus-core | +5.115575103 |
| **Total** | **-3.783707454** |

Induction and learned context make material contributions, but their sum does
not close the representation error. This decomposition does not uniquely
identify a physical cause or establish the experimental interpretation. Compared
with the older frozen ff19SB/OBC2 candidate, the new model moves1/4ordering to2/4
and partition-4.139186to-3.783707; the older responsive ff19SB/OBC2 candidate had
3/4ordering and-4.438147partition. Multiple model components differ between these
protocols. Neither the new calculations nor the previous ones support claiming
robust accuracy improvement or superiority to the production baseline.

## What the checks establish

All51static native energies,48direct-field queries and60coupled response solves
succeeded. Maximum primary endpoint correction change on1e-7to1e-9Debye solver
refinement is2.264684e-7kcal/mol. Maximum full rigid-replay component error is
2.842171e-12kcal/mol; both vacuum identity corrections are exactly0 in the stored
outputs. The largest prescribed +/-5%common-metal-radius effect on a relative
comparison is0.000197334kcal/mol. These effects are far smaller than the remaining
ordering/partition errors; changing these numerical settings is not a supported
remedy. This radius check covers only the prescribed metal radius, not every
possible cavity or dielectric approximation.

Primary maximum atomic induced dipoles range0.208933to0.336394eA, below the
declared1eA large-response flag. All source induced components stay exactlyzero.
Correct source charges, unchanged exterior moments, paired physical boundaries,
and unrounded residuals are verified in real native outputs. No solver divergence
or enormous dipole explains this candidate's failure.

## Energy expression and limits

Protocol `vacuum_density_direct_AMOEBA2018_GK_proxy_POLAR_short_hybrid_v1`:

```
A_M = E_DFT,vac(core,M) + V_density(Q,E)
      + [es_static(Qproxy,E) - es_static(0,E)]
      + [I(density+proxy_RF,E) - I(E)]
      + T_short(full,M) - T_short(core,M)
I = -0.5*(electric/dielec)*sum(mu_solv,d * F_solv,p)
```

The exact quantum density supplies permanent direct coupling and induction
driving fields; projected endpoint CHELPG charges supply only the GK proxy.
Native environment d/p channels and mutual/GK response are retained. Native
permanent vacuum energy is audited but excluded to avoid adding Q-Q and Q-E
terms already represented by DFT and density coupling. Same physical cavity
cancels nonpolar/environment-only terms. No CPCM or second full solvation is
added. Source density remains frozen in vacuum, the GK density representation
is approximate, and MACE's jointly trained short readout is not uniquely
nonelectrostatic. Nuclear geometry/protonation/assembly/waters stay fixed.
Combined gradients, mechanical correction and absolute calibrated scores remain
unavailable. Full expression and pre-output rules: [plan](DENSITY_GK_HYBRID_PLAN.md).

## Execution, tests and cost

Job1201074: **570wall seconds on64CPUs**, **36480allocated core-seconds**,
**3659.722079actual CPU-seconds**,2094872KiB sampled batch peak RSS,zeroGPU.
Native kernels sum124.133763045s wall; complete native child processes sum
210.915940717s wall/3316.199736CPU. Build1.775239s wall/1.713074CPU;
preparation41.400926s wall/38.763210CPU. No failed build, preparation, scientific
attempt or retry in this pilot. Exact Slurm JSON and execution receipts retained.

This is incremental development cost with repeated numerical controls.
Reused DFT, charge fitting, quantum field queries and MACE costs are additional;
a matched end-to-end production-pair cost is not yet measured. The new native
solver component is practical at this size; a claim of complete production
affordability would require that full accounting. The original vacuum8endpoint
DFT job alone cost48064allocatedcore-s; the new density-multipole job added6768.

Six distinct source/receipt/algebra/cache/error/native integration tests pass.
Before execution:4pass/2explicit skips in5.139s. After native calls:5pass/1skip
in22.325s. The remaining collection/native parser test passes16.558s once the
real report exists. No final integration test remains skipped or replaced by
fabricated scientific success. Eight additional density/boundary tests pass in
their separate reported runs. Commands: [operations](DENSITY_GK_COMMANDS.md).

Unrounded components, all variant contrasts, source hashes and costs:
[DENSITY_GK_HYBRID_RESULT.json](DENSITY_GK_HYBRID_RESULT.json).
Full immutable products: `workspaces/mace_omol_20260917/density_gk_hybrid_v1`.

**Assessment:** numerical implementation credible on these controls; physical
partition robustness and predictive gate fail; native component affordable,
complete production cost unmeasured. Recommendation: retain baseline and pursue
the specific missing quantum-density response question, using saved responsive
wavefunctions before spending more DFT. Do not tune this candidate's radius,
solver tolerance, threshold, labels or core selection to force a pass.
