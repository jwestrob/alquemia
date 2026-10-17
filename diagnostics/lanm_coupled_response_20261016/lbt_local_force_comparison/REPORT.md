# Local MACE–PBE0 force comparison

Saved comparison **12938 completed all four declared 48-atom cases**. The
[generated tables](GENERATED_REPORT.md) and [exact result](RESULT.json) retain
every force vector, native gradient, physical/cap map and SCF diagnostic.

## Meaning

MACE and PBE0 give the same Cartesian component signs for the target-metal force
in each of the four fixed structures. Their target-vector differences are
0.130118, 0.193107, 0.099720 and 0.185623 eV/Å for 7CCO/La, 7CCN/La,
7CCO/Dy and 7CCN/Dy respectively. These are reported vector differences, with
no force-accuracy threshold or qualification claim.

The differences also occur on represented real atoms: their RMS vector
differences are 0.164303, 0.194809, 0.171673 and 0.173151 eV/Å in the same
order. Synthetic-cap RMS differences are smaller in every case
(0.073341, 0.075533, 0.072577 and 0.073339 eV/Å). The comparison therefore
cannot attribute the method difference solely to the seven artificial caps.
Each core contains one target, 40 other real atoms and seven caps.

Together with the separately admitted same-metal energy comparison, these
results give no basis for replacing the pronounced local source penalty with
the flatter GFN response merely to reduce source dependence. This comparison
does not identify which aspect of whole-protein mechanics limits discrimination.
Jacob subsequently stopped the proposed clamped-core structural diagnosis.
The next planned work is the [coupled mechanics comparison](../../lanm_coupled_mechanics_20261017/PLAN.md),
separating restricted motion from inaccurate energetic costs.

## Limits

All exact local GFN receipts are energy-only, so GFN force comparisons remain
unavailable. PBE0 analytical derivatives have not passed a native finite
derivative check, and its retained SCF residual/state warnings remain relevant.
The frozen-f effective singlet and MACE physical Dy sextet metadata are different
electronic treatments; spin–orbit accuracy is unqualified. These are raw cluster
forces, not whole-peptide forces assembled through a cap chain rule. No affinity,
source-robust protein classification or default-method change follows.

## Cost and provenance

This saved-only job used one CPU for **10 scheduler seconds**, with recorded
TotalCPU **8.288 seconds**, application wall time **9.201953770010732 seconds**
and application CPU time **8.1563738 seconds**. It made **zero molecular calls**
and created no geometries. Parent MACE12937 supplied four new forces with exact
zero energy replay deltas; parent molecular costs remain separate.

[ACCOUNTING.json](ACCOUNTING.json) pins the actual scheduler receipt, generated
report, copied result, implementation, configuration and all parent results.
The immutable execution result SHA-256 is
`947433d3266386acc2a7305af3a4e3a70c1258cd808bde776ac238afc5a51da7`.
