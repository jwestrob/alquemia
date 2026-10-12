# Ca-conditioned coverage: 107 new prepared pools would still be needed

**Inventory only, 12 October 2026.** No coordinates, proposals or molecular
energies were generated. Any extension waits for interpretation of the current
corrected transfer result and its measured throughput.

The known benchmark has **125 Ca-conditioned folds from 25 proteins**. Original
source preparation admitted 110 (60 Ca-class, 50 La-class), from 22 proteins.
The other 15 remain unavailable: every Ca sample 0–4 for C5AXV8, Q88JH5 and
Q9Z4J7 failed with `exact selected PQQ required`. This inventory neither repairs
nor reinterprets those exclusions. All 125 original source-coordinate pins were
verified; no new folds are needed.

## Exact reuse versus missing work

| Category | Sources | Status |
|---|---:|---|
| Existing expanded own-three pools | 3 | A0A3F2YLY8 Ca samples 1/3; A0ACD6B9F2 Ca sample 4 |
| Supported source, expanded preparation still needed | 107 | New graph/state preparation must pass before scoring |
| Original PQQ-preparation exclusion | 15 | Retained in the 125-source denominator |

The three existing pools have verified exact coordinates, state/maps and all
18 MACE cells. Their 36 corrected GFN2 inputs per backend belong to root's
current transfer run; reuse requires actual successful component receipts.
No additional compatible expanded Ca preparations or proposals were found.
The older ten-fold archive contains 110 prepared Ca sources, but **none has the
same fragment list as its protein's frozen canonical 4.3 Å envelope**. Its
energies/proposals cannot simply stand in for the expanded model.

## Same-policy extension, if justified

Apply each protein's existing canonical **La-three-source 4.3 Å fragment union**
to its Ca-conditioned sources, exactly as the three declared Ca probes did.
Keep the archived core/proton/water state, physical cap mappings, bounded
optimizer, common origin/Ca/La pool, corrected solvent qualification and actual
canonical-only bands. Do not fit bands to these Ca folds.

This is a broader test of the existing fixed-membership probe rule. Discovering
a new Ca-only or mixed-metal three-source union changes membership selection
and is a separately versioned preparation experiment. The current inventory
creates no new triple aggregation rule and does not replace the existing
100 La-source triples or their exclusions.

For the 107 missing pools, conditional on successful preparation:

- 107 expanded source/context preparations and 214 bounded MACE searches.
- 214 MACE origin evaluations plus 214 cross-proposal evaluations; the 214
  final proposal evaluations are obtained within the searches. The resulting
  matrix has at most 642 cells before exact-coordinate deduplication.
- At most 1,284 new GFN2 inputs per backend: 2,568 primary-plus-independent
  calls. Optimizer trial count/cost cannot be inferred from endpoint count.

These are prospective upper counts for the declared finite source set, not
executed calculations or an automatic launch instruction. New preparation can
add explicit unavailable cases; it must not silently substitute an older carve.

## Reproducible inventory

`INVENTORY.json` contains every source, conditioning/class stratum, original
status, target union pin and exact reuse status. The stdlib-only audit is:

```bash
python diagnostics/orca_parameter_audit_20261012/ca_conditioned_inventory/AUDIT.py
```

A separate read-only archive scan confirmed the three Ca probes in the search,
pool and source-selection manifests and zero Ca sources in the 104-pair primary
transfer. No reserved/unknown PLM labels were opened. Production is unchanged.

**Recommendation:** finish and interpret the current corrected transfer first.
If its benefit and throughput justify extension, prepare the 107 missing
source/context requests under the existing probe rule and retain all 125 in
the evaluation denominator.
