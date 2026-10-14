# PedE: independent DFT retains the response direction, with a smaller magnitude

The eight-endpoint reference is complete. Native r2SCAN-3c/CPCM gives a source-dependent Ca−La solvent contrast of **+10.0342163383 kcal/mol**, versus **+21.2566283014** from corrected GFN2/ALPB. Their difference is **−11.2224119631 kcal/mol**. This exceeds the predeclared2 kcal/mol mechanistic triage scale. The existing solvent response deserves improvement, but setting it to zero is not supported: the independent model retains a substantial response in the same direction.

This is a fixed-geometry diagnostic, not demonstrated classifier improvement or proof that either result is the ground-state response. No score, calibration or production setting changed.

## Exact comparison

Both sources retain the same192 atoms, fourteen caps, dry inventory and PQQ-state hypothesis. Ca−3 and La−2 are physical singlets; paired coordinates and maps match. Each method uses its own complete endpoint totals:

`U_M(s) = E_water,M(s) − E_vacuum,M(s)`

`D = [U_Ca(1) − U_La(1)] − [U_Ca(3) − U_La(3)]`.

All displayed values below are kcal/mol. Unrounded Hartree endpoints, coefficients and source pins remain in the result; conversion used the unchanged627.509474 convention once. The printed CPCM dielectric contribution was not substituted for U.

| Quantity | DFT/CPCM | GFN2/ALPB |
|---|---:|---:|
| U_Ca(sample1) | −358.7366442021 | −424.3822138993 |
| U_La(sample1) | −255.4089128006 | −326.7216079226 |
| U_Ca(sample3) | −353.5027489481 | −367.6529149582 |
| U_La(sample3) | −240.1408012083 | −248.7356806800 |
| U_Ca−U_La, sample1 | −103.3277314014 | −97.6606059767 |
| U_Ca−U_La, sample3 | −113.3619477398 | −118.9172342782 |
| **D, sample1−sample3** | **+10.0342163383** | **+21.2566283014** |

The individual metal responses reveal more than the partly cancelling contrast:

| Sample1−sample3 response | DFT Ca | DFT La | GFN2 Ca | GFN2 La |
|---|---:|---:|---:|---:|
| Water endpoint | +298.2581258589 | +291.9326759148 | +237.0238688962 | +216.1370387419 |
| Vacuum endpoint | +303.4920211128 | +307.2007875071 | +293.7531678373 | +294.1229659844 |
| **Water−vacuum transfer** | **−5.2338952540** | **−15.2681115923** | **−56.7292989411** | **−77.9859272425** |

Both methods stabilize sample1 more strongly through solvent, especially for La, but GFN2's individual source sensitivities are much larger. The large common source-energy changes concern unrelaxed saved structures; they are not protein affinities or equilibrium source populations.

## Population response is distributed beyond the metal

Both independently reported DFT population schemes show stronger solvent-induced redistribution in sample1. Below is water−vacuum charge change in elementary-charge units; negative means increased assigned electron population. PQQ includes all27 atoms. Residue entries contain only the atoms retained in the local context, with synthetic caps separate.

| Scheme / source / metal | PQQ | Asp220 | Pro309 | Arg386 |
|---|---:|---:|---:|---:|
| Mulliken /1 /Ca | −0.602917 | −0.419517 | +0.483997 | +0.373421 |
| Mulliken /1 /La | −0.728765 | −0.431965 | +0.548104 | +0.393811 |
| Mulliken /3 /Ca | −0.185112 | −0.220944 | +0.081605 | +0.008871 |
| Mulliken /3 /La | −0.321198 | −0.235066 | +0.063430 | +0.011504 |
| Loewdin /1 /Ca | −0.587243 | −0.411357 | +0.452160 | +0.376463 |
| Loewdin /1 /La | −0.699186 | −0.423860 | +0.517389 | +0.395307 |
| Loewdin /3 /Ca | −0.161319 | −0.222353 | +0.064360 | +0.018373 |
| Loewdin /3 /La | −0.285615 | −0.226830 | +0.048738 | +0.018917 |

Metal charge changes are smaller: Mulliken Ca−0.017563/+0.001537 and La+0.075765/+0.121352 for samples1/3; Loewdin Ca−0.012045/−0.005792 and La+0.035668/+0.050937. Other residues and caps also respond. The shared qualitative pattern makes a metal-only explanation insufficient, but it neither establishes a particular electron-transfer pathway nor assigns integer oxidation states. The schemes are not interchangeable measurements, and agreement between them does not validate the underlying density.

## What this narrows down

The [saved GFN2 decomposition](../solvent_audit/REPORT.md) attributes its D to electrostatics+38.355221, intrinsic remainder−18.604272, hydrogen-bond+1.702550, surface−0.196870 and reference0 kcal/mol. Those are compensating parts of the self-consistent response. Dropping the intrinsic remainder or using printed solvation energy alone would change the defined model.

The DFT result independently reproduces a substantial source-sensitive response and accompanying cofactor/context redistribution, while disagreeing on its magnitude. It strengthens the case for examining coupled local electronic/solvent response. It does not isolate the error to ALPB, GFN2, the cavity, a particular fragment or an SCF solution: DFT/CPCM changes the electronic and solvent treatment together. The [cofactor audit](../cofactor_review/REPORT.md) found no demonstrated atom/H/cap preparation defect, but did not establish an electronic minimum or redox population.

**Next decision:** qualify the electronic-response discrepancy before choosing an affordable replacement correction. Retain both source geometries and both metals; include PQQ and the responding context groups in that diagnosis. No fitted scaling, source removal, redox change or new class boundary follows from these two inspected development structures. The eight endpoints do not calibrate a replacement scorer, establish biological transfer or show that accommodated PedE failures are repaired.

## Admission, artifacts and cost

All8/8 cells passed the existing diagnostic collection rules, raw final-SP/population replay, state/coordinate and native-receipt/hash checks. Five used TRAH and passed the actually printed1e−5 orbital-gradient criterion; the other three completed their native SCF route. Normal termination and those criteria do not establish a unique lowest electronic solution, stricter density convergence, solvent gradients or a binding free energy. The original failed collection and successful collection-only recovery remain preserved.

[EVIDENCE.json](EVIDENCE.json) pins the full result, config, receipts and reused component report. Authoritative result: `workspaces/pede_solvent_response_20261013/comparison_v1/execution/RESULT.json`, SHA256 `0df3271eb35b617bac1c79caf66be24cf5fe52f3da9c631f21dabd7bdf01bba0`.

Collector3377 completed0:0: 0.357787s analysis wall,0.284590s process CPU,25236KiB peak RSS; scheduler1s on one CPU. It reused eight molecular endpoints and performed zero new molecular calls. Molecular execution costs are retained individually in the native receipts; this small collection cost is not a new-site scoring runtime. Root owns their campaign accounting and next experiment. Production and the primary failed classifications remain unchanged.
