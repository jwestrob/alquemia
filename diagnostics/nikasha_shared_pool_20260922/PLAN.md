# Nikasha: common geometry competition

## Authorization and purpose

Jacob's 22 September handoff explicitly requests execution: recover completed
work, implement a shared candidate pool, then pursue force-selected local
accommodation and supported PLM delivery. Existing discretionary contained
MACE/GFN2 authorization remains. No old project budget or per-analysis approval
loop is reinstated. The separate unqualified native-candidate DFT campaign stays
dry-run only. No new DFT, optimization, protonation, water occupancy or folding is
part of this first experiment. Root owns shared scoring and integration; Khoury
recovers old outcomes, second_shell handles current naming/PLM input inventory,
water_basins examines archived physical forces and eligible modes.

The existing proposal scorer gives each metal only its own proposal and origin.
Both metals should compete over the same finite geometries before interpreting
their energy difference. No existing scorer, protocol ID or result is overwritten.

## Energy and fixed candidate pool

For each single source context and fixed chemical state, collect its actual q0,
Ca proposal and La proposal. With the established native OMOL100M and primary
native ORCA GFN2 recipe, evaluate

E_M(q) = E_OMOL,vac,M(q) + E_GFN2,ALPBwater,M(q) - E_GFN2,vac,M(q).

Use existing component units and conversion constants exactly once. The new
mathematical descriptor is R_pool = min_q E_Ca(q) - min_q E_La(q). Retain each
individual work relative to its own origin, every matrix cell and selected
geometry. Row minima can only lower each metal's energy; no sign is prescribed
for the change in R.

Also retain the existing operational origin tolerance: use the row-minimum
candidate only if it lowers energy from q0 by more than 0.10 kcal/mol; otherwise
select q0. Report both mathematical and operational results, plus the old own-
proposal/origin score. Ties use deterministic candidate order (origin, Ca proposal,
La proposal); report actual energy gaps rather than candidate populations.

Require common atom order/membership, source preparation, physical cap mapping,
PQQ/proton/water state and exact paired origin coordinates, apart from metal and
its corresponding charge. This first input set is dry PQQ with common origins.
Different wet-origin coordinates are explicitly unsupported here; later wet pools
must retain both origins and account for their additional evaluations, never
average them. No minimization across folds or changed context compositions.

Deduplicate numerical coordinate copies at the established 1e-12 Angstrom maximum
Cartesian difference, with identical atom identities after canonicalizing only
the substituted metal. Retain every source coordinate pin and the representative
chosen in deterministic order; no coordinate averaging. Reuse an existing cell
only with matching geometry/state/method provenance. A failed required old or new
cell makes the pooled result unavailable. Keep the separately named valid old
score visible; do not retry old failures or substitute it as pooled success.

## Finite execution sequence

1. Real integration pilot: 1H4I, 4MAE and the two original sample0 PLM compression
   contexts from the completed 30-source proposal archive. At most eight new MACE
   cells and sixteen GFN2 singlepoints, all optimization and four own cells per
   context reused. No new reference label; both PLM labels remain unknown.
2. After technical end-to-end success, complete the remaining 26 original
   canonical/crystal contexts, reusing the pilot. Total across these 30 sources:
   at most60 new MACE cells and120 new GFN2 endpoints, before deduplication.
3. Apply the same implementation to all225 declared primary noncanonical source
   folds. Preserve17 preparation failures, the existing missing q0 and any actual
   failed proposal cells. Do not spend new cells on a source already known to have
   an unavailable required cell. Upper bound before those exclusions:416 MACE and
   832 GFN2 calls. No new optimization, selected failure subset or score-driven
   retries. All cases remain in the report denominator.

Use one warm H200/32CPU/200000MiB worker per batch and existing native ORCA task
runners, eight concurrent eight-rank tasks per64CPU/128GiB shard. At most four
disjoint solvent shards for the full fold panel. Actual finite manifests and
dry-runs precede submission; count preparation, failures and allocated resources.
No project-time/compute cap. Keep input/output paths explicit and immutable.

## Calibration, evaluation and next decision

Preserve old-band transfer. A separate pool reference may use exactly the original
25 canonical training members and the same extrema rule/minimum-gap condition.
Freeze it before inspecting pooled fold outcomes; report mathematical and
operational pool variants separately under their declared references. Crystals,
PLM predictions and noncanonical folds never determine bands. No fitted mixture,
threshold search or post-result threshold movement.

Report complete known-class denominators, common-coverage decisions/transitions,
raw paired changes, structural spread, strict La4/Ca5 and all100 previously declared
La triples. Compare to the strongest supported existing scorer after recovery.
These consumed structures are development evidence and replicas are not independent
biological observations. A near-zero effect closes this finite-pool question while
retaining the scoring primitive. Improvements on unknown PLM scores alone are not
accuracy evidence. Subsequent force-selected searches reuse this common-pool rule.
