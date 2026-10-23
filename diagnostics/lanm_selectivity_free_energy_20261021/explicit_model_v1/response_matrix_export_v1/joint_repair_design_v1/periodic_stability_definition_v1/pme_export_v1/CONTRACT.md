# Initialized native PME snapshot: source-only exporter

## Exact hook and scope

The source authority is the frozen `periodic_native_prototype_v1/source_v1` manifest c4feea15. The isolated patch adds `call pmenormexport` **inside udirect1, after qfac(0) and its nonperiodic branch, immediately before fftback**, then appends the helper to induce.f. No existing source tree is edited. At this point the current grid, moduli, spline coefficients, igrid and chunk table were initialized by getchunk/moduli/bspline_fill/table_fill, and the entire qfac has been stored. The original qgrid, induced dipoles and all solver/native arrays remain read-only to the helper.

Opt in with exactly `PME-NORM-EXPORT 1`: export the first initialized udirect1 state in that process only, to `pme_norm_export_001.dat`. Absence returns before allocation/output or accessing native PME arrays. Further initialized states are not dumped and inherit **no admission or cached bound**. This is a static diagnostic interface, not an MD logging or stability implementation. A separate caller must bind source, executable, complete input/key/geometry and output hash; a future run at different geometry requires a new snapshot and new admission.

The helper accepts orthorhombic periodic TINFOIL only. It exports active `douind` sites and rejects nominal alpha<1e−8 on those sites, making the native PCG floor irrelevant in this initial contract. Inactive sites contribute no response columns and retain their native permanent-field role. Export capacity is at most1,048,576 grid points,20,000 active sites and spline order2..8; unsupported dimensions/aliased stencils are rejected. These are bounded reporting capacities, not physical thresholds. It checks the actual pmetable owner for every stencil point, compatible chunk partitions and the native igrid+2 stencil offset. No periodic execution guard is removed: the inherited backend/stability hooks remain fatal.

## Records and replay

The file contains version/header, initialization sequence, grid dimensions/order and counts, FFT backend, chunk layout, aewald/volume, the Cartesian-to-fractional matrix A[b,a]=nfft_b*recip(a,b), qfac min/max/maxabs/zero/evenness statistics, and scatter norm statistics. For every active site it records native atom index/type, igrid, alpha and Cartesian coordinates plus the first two stored theta entries for each axis/stencil slot. It then writes **every actual qfac value** in native Fortran grid order, including masked and zero entries. This is O(M+Np) records, not an approximation to native multipliers or a dense response matrix. END carries active/theta/qfac/stencil counts. New-file open refuses collisions; failed writes leave incomplete evidence, and the reader rejects missing, duplicate, out-of-order or trailing records.

`SNAPSHOT.py` replays the norms from those stored quantities. For each active site and wrapped stencil point,

```
d=(theta1' theta2 theta3, theta1 theta2' theta3,
   theta1 theta2 theta3')
B[g,ia]=sqrt(alpha_i) sum_b d_b A[b,a].
c1=max_column sum_grid |B|; cinf=max_grid sum_column |B|
cF=sum_grid,column B²
U_R,nominal=M*max|qfac|*min(c1*cinf,cF).
```

No grid aliases occur with order≤every dimension, so cF is formed from actual entries. Different sites are different columns. Computing these quantities needs O(Np³+M) arithmetic and one O(M) row-sum array; raw replay also stores qfac and O(Np) splines. The M factor is required by the native unnormalized forward/inverse FFT pair. For nonnegative, conjugate-even stored weights the intended discrete reciprocal stiffness is PSD and U_R,nominal bounds its norm in exact arithmetic. The report records negative counts and max conjugate difference rather than modifying qfac or claiming evenness from a formula. Half the conjugate difference times the same multiplier-free bound gives a nominal two-channel asymmetry perturbation estimate.

Native/helper statistics and independent saved replay are ordinary-double arithmetic. They are not outward certificates, and this task establishes no FFT/scatter/gather application error allowance. Subsequent numerical qualification must resolve those effects using the inherited conventions; it need not introduce a new formal interval-factorization prerequisite. U_R alone neither supplies a positive real-space lower bound nor admits the full relative spectral criterion. No solver convergence, static snapshot, or valid file makes that inference.

## Prospective finite controls

One thread, artificial arrays only: dense scatter versus streamed norms and Fourier operator; correct alpha/Cartesian normalization; conjugate-pair handling and odd perturbation; wrapping and no-stencil-alias rejection; parser complete/duplicate/truncated/trailing-record checks; source hook placement and absent-path early return. Fixed ordinary-double comparison allowance1e−12. No native/compiler/physical box/candidate evaluation is performed. These controls qualify the source/replay design only; compilation, actual initialized-state parity and resource measurement remain unperformed.
