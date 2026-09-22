# Separate adaptive shared-pool comparison

The existing `inspect` and `calibrate` commands now recognize the separately
named adaptive protocol from its actual collection and pinned manifest. Base
shared-pool behavior and the frozen `f96b2425…` reference remain unchanged.
`compare` still accepts only the existing base `primary225` population.

## Scientific boundaries

- Adaptive protocol: `nikasha_adaptive_angular_common_geometry_native_OMOL_GFN2_ALPB_v1`.
- Reference: `Nikasha_adaptive_angular_common_geometry_canonical25_reference_v1`,
  with distinct mathematical and operational suffixes.
- Require exactly `adaptive4` plus `adaptive26` and the complete original 30-source
  union. The same 25 canonical identities, labels, class-extrema rule and strict
  >0.02 model kcal/mol minimum gap apply. Crystals, PLM sources and folds never
  calibrate either variant. An unavailable canonical member leaves bands null.
- Base/adaptive collections or references cannot mix. Model, software, executable
  and protocol-specific selection settings must match.
- Native GFN2 iteration ceilings 125 and 500 represent unchanged convergence
  tolerances and parameters only under the actual pinned MaxIter qualification.
  A 500-ceiling manifest or recovered pilot requires that successful diagnostic:
  both controls pass, the formerly failed cell completes, all three rows retain
  parameter identity and charge checks, and the ORCA pin matches. These numerical
  differences remain separate from the Hamiltonian compatibility signature.
- Inspection/reference outputs preserve the primary failure, explicit recovery
  source, qualification pins, ceiling/policy and additional attempts. The recovered
  pilot must retain the original per-case primary pools. No missing cell is zeroed
  or replaced by a baseline score.

## Actual available operation

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/nikasha_pool_compare.py inspect \
  --collections workspaces/adaptive_accommodation_20260922/common_pool_v1/recovered_collection_v1.json \
  --output workspaces/adaptive_accommodation_20260922/common_pool_v1/comparison_extension_v2.json
```

The executed version of this command wrote `comparison_extension_v1.json` and its
Markdown companion: four available cases, two unknown PLM labels, one retained
primary failed case and three explicitly recorded diagnostic attempts. No adaptive
reference has been fitted by this extension. Calibration must wait for the actual
remaining 26 result collection; its command interface is unchanged.

## Checks

Twelve tests pass with zero skips. They use the actual recovered pilot and archived
base results, including a full replay of the frozen base canonical reference.
They reject incomplete adaptive calibration, base/adaptive mixing, use of the base
reference for adaptive scores and an explicitly corrupted real recovery with its
qualification removed. No adaptive energy or successful scientific fixture was
fabricated. Full adaptive calibration remains unrun until its real inputs exist.
