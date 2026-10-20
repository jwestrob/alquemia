# Energy lowering survives; reference forces remain substantial

20 October2026. Native15545 and saved collector15838 complete the exact
step92-to725 whole-peptide La comparison. PBE0 plus the same fixed-charge ALPB
supports the energy lowering, but the MACE-converged endpoint is not stationary
under the admitted PBE0-plus-ALPB forces. Neither result establishes that PBE0
is the physically correct force field.

## Same deformation, different electronic models

All energies below are **step725 minus step92**, in kcal/mol. These are native
endpoint differences; the separate optimizer's172-to725 work is not substituted.

| Component | MACE model | PBE0 reference |
| --- | ---: | ---: |
| Whole-peptide vacuum energy | +16.362172 | +15.133560 |
| Same fixed-FF-charge ALPB | −45.442067 | −45.442067 |
| **Composite work** | **−29.079895** | **−30.308506** |

The model-minus-reference discrepancy is **+1.228612kcal/mol**, outside the
unchanged1kcal reporting resolution. Both models nevertheless support lowering
for this same deformation. That resolution is not a native admission gate,
confidence interval or physical uncertainty estimate.

The lowering is driven by the shared solvent term; both vacuum terms increase.
This comparison therefore checks the electronic cost of the motion but does not
independently validate the solvent benefit or solvent polarization.

## Model convergence does not transfer to reference stationarity

Force values are kcal/mol/Å. RMS means atom-vector RMS over all256 atoms.

| Quantity | Step92 | Model-converged step725 |
| --- | ---: | ---: |
| MACE+ALPB force RMS | 1.086248 | 0.029611 |
| MACE+ALPB maximum atom force | 5.674648 | 0.075055 |
| PBE0+ALPB force RMS | 3.639505 | 3.587757 |
| PBE0+ALPB maximum atom force | 7.904442 | 7.868978 |
| MACE−PBE0 force-error RMS | 3.579797 | 3.589135 |

The original model stopping limits were maximum0.10 and RMS0.03. Step725
satisfies them for MACE+ALPB, as previously reported. The reference residuals
remain substantial despite the successful model minimization. Continuing the
same model optimization is not a remedy for this cross-model discrepancy.

Along the frozen normalized step92-to725 chord, final composite force projection
is−0.030674 for MACE and−1.754370 for PBE0. This is a descriptive projection of
the saved forces, not an independently qualified derivative or search direction.
Full force vectors and both original endpoint reports are preserved unchanged.

## Native admission and remaining qualifications

The new endpoint completed one ORCA6.1.1 vacuum PBE0-D4/def2-TZVP/RIJCOSX
analytic-gradient call with the archived La basis/ECP46 assets, charge−2,
physical/effective singlet, DefGrid3, VeryTightSCF, NoAutostart and PModel guess.
All256 atoms are physical; no external point charges or synthetic caps occur.
The actual12MPI allocation supplied48000MiB, giving MaxCore3100MB/rank.

Normal termination, native SCF convergence and the full finite gradient pass
the original admission. ORCA reports convergence by its gradient check after18
cycles. Orbital gradient2.1020e−7 and orbital rotation5.8796e−7 pass their printed
2e−6 thresholds. Energy change7.1741e−9, maximum/RMS density changes and DIIS
error do not pass all their printed thresholds. These observations are retained;
the admission rule was not changed after seeing this outcome. An early small-gap
warning remains evidence, not proof of the final ground state.

Whole-system PBE0 finite-difference consistency has not yet been demonstrated.
The native result and collector both retain physical-state and force-consistency
qualification as false. Conditional terminal chemistry, the vacuum dianion,
finite basis, ECP/no-SOC and post-added solvent assumptions remain.

PBE0-D4 is an independent comparator, not MACE's training Hamiltonian. The
pinned MACE-OMOL target is ωB97M-V/def2-TZVPD; the dataset also uses different
lanthanide spin treatment. Exact metal basis/ECP equivalence is not established.
Thus the force discrepancy cannot be uniquely assigned to learned-model error.
See the [OMol25 calculation details](https://arxiv.org/html/2505.08762v1#S2.SS6)
and the [local method evidence](../lbt_7ccn_reference_direction_v1/HAMILTONIAN_CONTEXT.md).

## Next bounded physics check

Root declares [one directional reference consistency test](../lbt_7ccn_reference_direction_v1/SPEC.md):
two fixed singlepoint energies around exact step725, along the normalized full
MACE−PBE0 vacuum force-difference vector. The existing central gradient is reused.
Step size0.01Å in the full Cartesian norm and the numerical comparison rule are
frozen before either displaced SCF. Settings and physical state stay fixed.

This tests whether PBE0's own energy changes support its reported derivative in
the direction relevant to the discrepancy. It does not assume that PBE0 must be
right, change the geometry by optimization, reopen the closed Mex branch, or
authorize a functional/grid sweep. One sampled derivative does not qualify the
whole force field. No metal-preference or natural-protein transfer claim follows.

## Execution, review and preservation

Native15545 used25,985scheduler seconds on12CPUs: **311,820 allocated CPU-s**,
zero GPUs. Collector15838 used7seconds on1CPU: **7 allocated CPU-s**. The earlier
preflight15543 used5CPU-s. This declared phase totals311,832allocated CPU-s;
the reused15321 reference and previous optimizer costs are recorded separately.
Only one new DFT call occurred; the collector made zero new molecular calls,
optimizer steps or geometry evaluations. Allocation is not measured utilization.

- Native RESULT: `a53125b98e597ea0473e06d58316ac2b32c914e9647da595f163cf0aa76a8c19`.
- Collection RESULT: `407bdf88e42c43cb2ae48ed1493cbbc008f5c0eb957628787a6aa1c1b795bb0d`.
- New force report: `e6c0d039`; old force report: `4f96db5d`.
- Source: `f251ca1e`; run config: `e995bda3`; collection config: `cd2189ab`.
- [Native archive](native_result_v1/ARCHIVE.json): `5a9fffa1`.
- [Collection archive](collection_result_v1/ARCHIVE.json): `62909dbb`.
- [Independent native review](ACTUAL_NATIVE_REVIEW.md) and
  [collection review](ACTUAL_COLLECTION_REVIEW.md).

Archives preserve native text/gradient/input/runtime/receipts, full force tables,
source bindings and actual accounting. The larger density and restart files
remain hash-pinned in their original workspace; exclusions are explicit in the
native archive. Original results, prior failures and production are unchanged.
