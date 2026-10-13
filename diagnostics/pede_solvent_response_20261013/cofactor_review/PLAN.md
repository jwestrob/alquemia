# PedE PQQ chemical representation: saved-source audit

2026-10-13. This bounded audit asks whether the exact PedE preparations contain a demonstrated cofactor or cap defect, and whether their inherited PQQ geometry has independent experimental support. It runs only under a real Slurm allocation. No energies, forces, geometry generation, protonation, optimization or folds.

## Fixed scope

- Reuse the completed geometry3340/3344 result: five two-Ca PedE raw sources, five La-conditioned raw sources, and nineteen saved local geometries. Unsupported La sample1 remains raw-only. Do not rerun their geometric feature extraction.
- Audit source-to-local PQQ membership/coordinates for the nine admitted origins, using the exact context preparation, core preparation, maps and Ca/La endpoints. Require all24 PQQ heavy atoms and exactly the three declared invariant H; record source chemical-state metadata. Check all14 cap records for PQQ involvement, and report the retained/omitted identities rather than infer missing bonds from distance alone.
- Compare the existing bond metrics with every complete PQQ instance in three already-consumed archived crystal sources:1H4I,4MAE and1FLG. Preserve every instance and its source identity; instances are not independent proteins or redox observations. These reference geometries are descriptive, not target minima or new biological labels.
- Pin the current official PQQ CCD as a **chemical graph/bond-order reference**. Its coordinate model comes from1CRU and its ideal coordinates are explicitly unavailable; do not manufacture ideal bond lengths from it. Match generic names to CCD names with the fixed full graph correspondence, verified against both existing naming-schema bond lists and CCD heavy connectivity. Carboxyl oxygen pairs are equivalent under the scored deprotonated state; report their unordered distances as well as exact source labels. The two quinone-region C–O roles have unique correspondence.
- Compare the exact old Ca, new two-Ca and canonical La AF3 input PQQ SMILES and source CIF descriptors. They declare neutral carboxylic acids at fold time; scoring uses the fixed oxidized PQQ3−, C14H3N2O8, singlet state. Record this existing hypothesis change explicitly. A graph/charge declaration does not constrain the self-consistent electronic density to a unique local oxidation state.

## Measurement and admission policy

Use the same0.90–1.85Å broad connectivity interval already used by preparation, solely as a replay check. No new fine bond-length, planarity or class threshold. Exact source/local coordinate tolerance1e−6Å, consistent with pinned prepared-coordinate serialization. Paired endpoint coordinates must agree within1e−12Å. PQQ hydrogen lengths are checked against the already-declared1.01Å(N–H) and1.08Å(C–H) construction with2e−6Å serialization tolerance. These are implementation checks, not force-field or redox validation.

Compute each experimental PQQ's26 heavy graph-bond lengths and ring-plane deviation with the same ring atom membership as the saved analysis. Do not align proteins, copy crystal coordinates into folds, or select the reference closest to a desired label. Output all entries, including missing/ambiguous references. Capture nearby explicitly represented source heterogens descriptively only if present in the archived metadata; no new interaction-energy interpretation.

Interpret source/geometry/state separately. In particular, 1H4I's primary report interprets its PQQ as semiquinone; 4MAE contains a metal-bound PEG species omitted from the historic dry scoring core. Neither is automatically an oxidized dry-pocket minimum. The source CCD's chemical component name alone is not experimental oxidation-state evidence.

## Delivery

One one-CPU Slurm analysis, pin-checked fresh output. `ANALYZE.py freeze --workspace ...` performs metadata preparation only; `execute --manifest ...` guards Slurm before importing NumPy/Gemmi. Root owns submission/wake. Report geometry comparisons, PQQ/H/cap checks, explicit limitations and a remedy only if supported. Do not alter production or histories, introduce new proton/redox states, shorten a bond by hand, or launch a new quantum campaign in this branch.
