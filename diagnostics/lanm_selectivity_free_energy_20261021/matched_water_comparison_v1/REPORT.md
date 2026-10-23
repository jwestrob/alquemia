# Matched water-removal comparison

All **18 new endpoints** are admitted: six ECP-PBE0, six MACE, and six transferred-AMOEBA results. Surviving metal/ligand coordinates, ordering and charge are exactly preserved; AMOEBA types, remaining bonds, keys and parameter files match their hydrated parents. Zero molecular calls were made by the collector or review.

## Main result

This is the metal-dependent energy of adding the water shell at fixed surviving geometry: `(Dy wet − Dy dry) − (La wet − La dry)`. Negative values favor water addition more for Dy. All values are kcal/mol; these are fixed-configuration electronic/model energies, not free energies.

| Donor | ECP-PBE0 | MACE | Transferred AMOEBA |
|---|---:|---:|---:|
| Monodentate acetate | −38.768 | +43.977 | −39.255 |
| Bidentate acetate | −39.094 | −27.892 | −32.659 |
| Acetamide | −45.364 | −26.785 | −46.547 |

MACE minus reference is **+82.745, +11.202, +18.579** kcal/mol. The large first discrepancy survives independent raw-energy arithmetic and actual MACE charge/spin/graph checks; it is not a cycle-sign or missing-reference artifact. This identifies a strong model dependence of environmental response. It does not by itself identify a specific learned term, establish the DFT reference as physical truth, or establish which endpoint causes the discrepancy.

Transferred AMOEBA minus reference is **−0.488, +6.435, −1.183** kcal/mol. Its good energy-increment agreement in two fixtures does not validate its donor interactions. Its corresponding hydrated exchange contrasts still differ by −6.811, −7.838, +7.388 kcal/mol, and the water-induced paired ligand-radial derivatives are −10.115, −10.963, +11.705 versus reference +10.749, +9.800, +4.907 kcal/mol/Å. Thus scalar energy cancellation coexists with force-response differences.

## Interpretation and separate distance control

The wet and matched-dry calculations share the same coordinates for every retained atom. This removes the old donor-distance confound from their difference. The result includes the total response to the added water shell: metal–water, ligand–water, polarization/electronic redistribution and other interactions. It is not a unique decomposition into any one pair or a protein-affinity prediction.

The collector separately recomputes original-dry and matched-dry exchange contrasts using the **same aquo8 reference**. For bidentate acetate, the placement-only change is +11.795 kcal/mol in ECP-PBE0, −0.758 in MACE and +13.023 in transferred AMOEBA. The snapshot has no original-dry monodentate/amide DFT endpoints; those two placement comparisons remain explicitly unavailable. They are not zero and do not block the complete new hydration comparison.

## Verification and retained limits

The collector replays raw ORCA energy/gradient units and coordinates, raw Tinker energies/gradients, exact geometries, chemical states, native admission and 242 artifact pins. A separate Decimal calculation and sum-before-projection force audit passed **443 checks**: maximum energy/contrast arithmetic residual 1.17×10⁻¹⁰ kcal/mol and radial residual 2.85×10⁻¹⁴ kcal/mol/Å. Eight saved-only controls test sign, cancellation, missing values, charge/atom imbalance, force sign/dimensions, nonfinite values and changed artifacts.

Physical La/Dy multiplicities remain 1/6; the frozen-f ECP reference uses its declared effective singlet. Native printed SCF flags and unqualified physical-state/method-accuracy flags remain in the machine-readable result. No model is promoted, no parameter is fit, and no experimental selectivity is claimed.

## Artifacts

- Collector/config: `COLLECT.py` and matching workspace `CONFIG_v1.json`.
- Immutable result: matching workspace `collection_v1/RESULT.json`, SHA256 `211911fea05f84992d921ccbb05cd71f9caa07313aa7435dffb34a47c499370c`.
- Independent audit: matching workspace `AUDIT_v1.json`; test receipt `TESTS_v1.json`.
- Source declaration: `SPEC.md`; `AUDIT.py` and `TEST_COLLECT.py` are saved-only.

The output contains all three reference contexts, endpoint energies/forces, exact result pins, explicit missing statuses and full paired radial responses. Collection took 1.13 seconds on one local thread.
