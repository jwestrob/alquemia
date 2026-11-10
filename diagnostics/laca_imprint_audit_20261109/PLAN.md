# La/Ca imprint audit — frozen plan

Frozen 2026-11-09, before any structure file was downloaded and before any outcome
(coordination number, chelate count, distance or rule call) was computed for a paired site.
Execution: biotite login node, single-threaded HTTP + light parsing; no Slurm, no QM/ML.

## Question

Imprint hypothesis: CN, chelating-carboxylate count and mean M–O record which ion is bound,
so the rule `CN ≥ 9 or ≥ 2 chelating carboxylates → La` answers "which ion is in this
structure". Preference hypothesis: the features reflect the site and persist across ions.
Test: the same native-Ca protein site crystallised once with Ca and once with a lanthanide.

## Frozen definitions (reused unchanged)

Geometry: `cod_reference_geometry_2026-10-13/scripts/protein_shells.py::process` (md5
ca1e68fd0cc5bc6c6d0e0b357f068318) imported unchanged, with gemmi from the
`conda_envs/lanm_qmmm` environment it was written for. Two in-memory adaptations, neither a cutoff:

1. `TARGET` extended from (La, Ca) to Ca + {La Ce Pr Nd Sm Eu Gd Tb Dy Ho Er Tm Yb Lu}.
   `metal_o_cut` already treats all of these as large cations (3.1 Å).
2. Shape and arrangement calls (`cshm.all_shapes`, `go.arrangement`) are stubbed to `{}` and
   `BV_R0` returns the La value for any lanthanide. These outputs do not enter CN, chelate
   count, mean M–O or the rule.

Features: `features(s, protein=True)` copied by AST extraction from `la_ca_diagnostics.py`
(md5 9204a0d9613df6541cda945276caad3b) with its WATER/NEUTRAL sets; PRIMARY = 3.0 Å.
Rule call: La iff `cn >= 9 or n_chel >= 2`.

## Candidate discovery (RCSB search + data APIs)

1. Ln entries: X-ray, resolution ≤ 2.5 Å, containing a non-polymer component in
   {LA CE PR ND SM EU GD TB DY HO ER TM YB LU}.
2. Their protein polymer entities' RCSB 95 % sequence-identity cluster IDs.
3. Ca partners: X-ray, ≤ 2.5 Å, containing component CA, with a protein entity in one of
   those clusters. Pair identity is re-checked by global alignment of the two chains
   (identity over the shorter chain ≥ 0.95).
4. Download cap ≈ 1,000 mmCIF files. If clusters need more, clusters are shuffled with
   `random.Random(42)` and taken in order; per cluster, at most 8 Ln and 8 Ca entries
   (lowest resolution value first, ties by ID). The cap and any sampling are reported.

## Site matching

For each Ln ion (first model, altloc '' or 'A'): coordinating residues = protein residues
of the asymmetric unit with an O/N donor ≤ 3.0 Å (symmetry-image donors are ignored for
superposition). Superposition atoms = backbone N, CA, C, O of each coordinating residue and
its sequence neighbours ±2, mapped to the Ca-entry chain(s) of the same cluster through the
chain alignment; ≥ 8 matched atoms required. Kabsch superposition of the Ca entry onto the
Ln entry; a Ca ion of the partner entry within 1.5 Å of the Ln ion is the same site. All
chain assignments of the partner are tried; the closest Ca ion is kept per partner entry.
Superposition RMSD is recorded, unfiltered. A Ln site with no coordinating protein residue
is skipped.

## Native-Ca curation (named claim)

A pair is NATIVE-CA when ≥ 2 of the coordinating residues of the Ca-deposit site map
(RCSB UniProt alignment) to UniProt positions carrying a "Binding site" feature whose ligand
is Ca(2+), where the UniProt entry is Swiss-Prot reviewed or the feature's evidence is
non-automatic (anything other than ECO:0000256/0000259/0007829/0000313/0008006). Primary-citation reading is out of scope for
this pass. Pairs without that statement form the UNVERIFIED stratum.
Excluded from NATIVE-CA (reported as a count): UniProt protein name, entity description or
cofactor comment naming a lanthanide or lanmodulin, "lanthanide-binding tag"/"LBT", XoxF/
ExaF/PedH-type Ln-dependent dehydrogenases, de novo / designed / synthetic-construct entities.
Reverse pairs: UniProt Binding-site ligand is a lanthanide ion at the Ln site's residues.

## Statistics

Unit of call: (entry, site) — chain copies of one site in one entry are averaged.
Protein = UniProt accession (fallback: 95 % cluster). Within a protein all its paired Ln
deposit-sites and all matched Ca deposit-sites are averaged.
- p = protein-weighted mean fraction of Ln deposit-sites called La (native-Ca, primary
  stratum = Ca-deposit site has ≥ 4 protein donors, matching the benchmark's Ca inclusion).
- q = same for the matched Ca deposit-sites.
- Bootstrap: 2,000 resamples of proteins with replacement, `numpy.random.default_rng(42)`;
  95 % percentile interval on p − q.
- Paired changes (Ln − Ca), per Ln deposit-site minus mean of its matched Ca deposit-sites,
  averaged within protein, bootstrapped as above: ΔCN, Δmean M–O, Δchelates.
Frozen readings (primary stratum, ≥ 10 proteins required for either verdict):
- IMPRINT SUPPORTED: p − q ≥ 0.40 and interval excludes 0.
- IMPRINT REFUTED: p − q ≤ 0.10.
- otherwise INCONCLUSIVE (report n).
Secondary, descriptive only: all native-Ca regardless of donor count; UNVERIFIED stratum;
light (La–Eu) vs heavy (Gd–Lu) lanthanides; deposit-pooled (unweighted) p, q.

## Part 3 (power)

Smallest accuracy difference between two discriminators detectable at 80 % power,
α = 0.05 two-sided, La side with family (n = 3) as unit; Ca side with site (n = 125) and
deposit (n = 42) as units; paired (same sites) and unpaired variants; exact tests.
