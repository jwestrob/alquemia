# Mex electronic response: differential signs agree, magnitudes differ

**All six predeclared Dy-minus-La electronic gradient projections have the same
sign in the approximation and archived PBE0 reference.** The largest differential
discrepancy is 2.934454 kcal/mol/Å along target z. Weak own-metal and common-load
sign disagreements remain. This supports a contained response experiment, not an
accuracy pass or an affinity prediction.

The comparison uses identical repaired Mex8FNS coordinates, caps and external
field: 209 local atoms, 1840 field sites and 2049 physical atoms, with local
charge 0 and three fixed Nd(III) spectators. The approximation is
`MACE_vacuum + GFN2_embedded − GFN2_vacuum`; the reference is embedded PBE0-D4.
Both cap-anchor derivatives and the field gradients are mapped to physical atoms.

## Saved result

These are signed **energy-gradient** projections; force loads have opposite
signs. Each fixed direction is normalized over all physical Cartesian coordinates,
including hydrogens. The donor directions rotate the complete sidechain about
CA→CB and preserve the carboxylate group. Units are kcal/mol/Å.

| Direction | PBE0 Dy−La | Approximation Dy−La | Approximation−PBE0 |
|---|---:|---:|---:|
| Target x | −2.144963 | −1.790690 | +0.354272 |
| Target y | −4.908377 | −6.821325 | −1.912948 |
| Target z | +5.434087 | +2.499633 | −2.934454 |
| Asp84 CA–CB | −9.698416 | −12.055650 | −2.357235 |
| Asp86 CA–CB | +0.466401 | +0.658475 | +0.192074 |
| Glu95 CA–CB | −1.342031 | −1.526020 | −0.183988 |

Asp86 has weak electronic loads with opposite signs for La
(PBE0 +0.428455, approximation −0.454356) and the common mean
(+0.661656 versus −0.125119). Adding identical MM and fixed-charge ALPB terms
also exposes a total La Asp84 sign difference: PBE0 plus shared terms gives
−0.583277, versus +1.264644 for the approximation. Thus matching differential
signs does not guarantee the same own-metal descent direction. The full saved
La, Dy, common and differential projections, including separate MM/ALPB loads,
remain in [RESULT.json](RESULT.json); no locally normalized loads are summed into
an energy or label.

## Scientific decision and limits

Proceed with the already declared [two-direction physical derivative check](../MEX_FORCE_CHECK_PLAN.md).
It tests whether the combined force is the derivative of the implemented energy,
including external-site and cap-anchor motion. If admitted, proceed to the bounded
own-total-energy protein response with both metals cross-scoring the same candidate
pool. Preserve the discrepancies above; do not tune directions, scales or thresholds
to make the methods agree. This comparison does not replace force consistency or
justify a larger optimization/DFT campaign.

The reference is ORCA 6.1.1 PBE0-D4, light-atom def2-TZVP and explicit lcecp1-TZ
metal bases/AuxJ, RIJCOSX/DefGrid3/VeryTightSCF/DoEQ false. La ECP46 and Dy ECP55
both use 832 explicit electrons and effective restricted-singlet valence. MACE
uses physical multiplicities 1/6; GFN2 has its own f-in-core effective singlet.
These are different electronic representations. Archived physical-state and
force-consistency qualification remain false; small-gap warnings and absent
spin–orbit/stable-state evidence remain. There is no matching isolated PBE0 pair
to separate MACE error from the field correction. The shared FF-charge solvent
is not independently validated by electronic agreement.

The result provides local response evidence at one fixed occupied source. It
establishes neither finite relaxation energies, source robustness, full-protein
binding thermodynamics nor a natural-protein versus peptide classifier.

## Cost and provenance

Job **12686 COMPLETED, exit 0:0**, in **4 scheduler seconds on 1 CPU**; no GPU.
Application wall time was 3.770516 seconds and scheduler measured CPU time was
3.366 seconds. **Zero molecular calls, displaced geometries or optimizer steps.**
This report copies saved projections and reads terminal accounting; it performs
no new projection calculations.

- [Actual Slurm result](../../../workspaces/lanm_coupled_response_20261016/mex_reference_response_v1/execution/RESULT.json), SHA256 `1f3174a2849787c2878512114db578de398cd545f76b41db3f04118370f53091`.
- [Compact result](RESULT.json), SHA256 `287222ad9fb7b747e76302d1d3f4b8f1cd22d07995da624104b4adc3d7a53f9a`.
- [Terminal accounting](ACCOUNTING.txt) and [frozen comparison plan](../MEX_REFERENCE_RESPONSE_PLAN.md).
- The actual result pins admitted static 12685 (`ff7b27d4`), both archived PBE0
  rows, raw gradients, source identities, six Jacobians and all replay residuals.

No production/default, shared status or experimental label changed.
