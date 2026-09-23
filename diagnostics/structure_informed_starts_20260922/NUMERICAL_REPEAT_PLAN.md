# Q88 exact-primary repeats: executed scope replaces unsubmitted Tight proposal

Before numerical submission, root verified the
[ORCA6.1 native-xTB SCF documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/semiempirical.html#scf-with-native-xtb):
native special SCF settings override other ORCA SCF settings. A printed generic
Tight tolerance is therefore insufficient to establish actual tighter convergence.
The earlier native-tight experiments and their outputs remain archived, but their
label is not evidence that all native-mixer criteria were tightened.

Use root's explicitly authorized unchanged-input repeat alternative. The prepared
`numerical_v1/manifest.json` (two policies includingConvergenceTight) is preserved
**unsubmitted**. No numerical output was inspected before this change.

Executed request: **oldadaptiveLa andnewtemplateLa exactQ88XYZ × vacuum/ALPBwater
× two separate fresh primary-input repeats =8LaGFN2 calls**, zeroMACE/DFT. Same
nativeHamiltonian,NoAutostart,300K,MaxIter500,UseXTBMixertrue,charge−2,singlet,
8MPIranks, fresh workingdirectories withoutGBW/xtbw. Bothrepeats run on one64CPU,
128GiB allocation at8×8concurrency. Do not select a favorable repeat or average
this diagnostic into historical scores. All original outputs stay immutable.

Retain all diagnostics specified in NUMERICAL_PLAN.md, including every printed
energy/MAXdensity/RMSdensity residual and tolerance, whether each passes,
SCFcyclehistory, occupations/charges, source/input/runtime hashes and receipts.
Compute each exact-XYZ repeat difference, within-host repeat spread, old/new
coordinate difference and matchedLaALPB-minus-vacuum contrast. Keep0.1kcal/mol
endpointtransfer /0.2kcal/mol contrast diagnostic scales unchanged. No bandpadding,
new threshold, ordinarySCF retry or numerical score correction. This tests
reproducibility, **not successful native tightening or a complete newCa/Lascore**.
