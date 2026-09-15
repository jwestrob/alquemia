# La/Ca discriminator external benchmark expansion

**Frozen:** 2026-09-15, before any energy in this expansion was inspected
**Objective:** determine whether the energetic method carries La/Ca information beyond the canonical PQQ-MDH D+2-Asp/core-charge split.
**Parent PQQ protocol:** `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`

## Question and evidence boundary

The existing 25-protein PQQ-MDH calibration and its 1H4I/4MAE crystal transfer
are complete. They are not rerun here. Their perfect separation is confounded:
the D+2 Asp motif, fixed-core atom count, and fixed-core charge all track the
biological label exactly. This expansion asks whether frozen calculations recover
independent experimental La/Ca directionality and transfer to new proteins and
architectures.

Thermodynamic affinity, functional cofactor use, native occupancy, and regulatory
response remain separate labels. Unresolved cases are predictions, never silently
converted into truth labels. The PQQ-MDH threshold of 19.1586 kcal/mol applies only
to exact `fixed_core_v3` canonical eight-bladed PQQ-MDH inputs. Non-PQQ and
six-bladed scores are reported on their own frozen protocols and are never pooled
with that threshold.

## Lanes

### 1. Canonical eight-bladed PQQ-MDH transfer

Use exact `fixed_core_v3`: complete PQQ(3-), anchor Glu/Asn, catalytic Asp and its
H-bonded Arg/Lys, plus D+2 only when it is an acidic direct donor; dry vertical
La/Ca swap at fixed coordinates, r2SCAN-3c/CPCM(water), native ORCA basis/ECP.

- Immediate new crystal control: Q46444/QH-ADH, PDB 1KB0, Ca biological class.
- New structural replicates: 5GB1C MxaF; LW13 XoxF/MxaF; five Huang et al. strict
  Ln-positive enzymes.
- Boundary predictions, excluded from binary accuracy: CBE67239.1,
  WP_012592127.1, and WP_245258612.1.

Every modeled target uses the same preregistered seed/sample policy. Structural
acceptance and core-role selection are completed without energetic output. All
valid replicates are scored; no best-energy structure selection is allowed.

### 2. Direct non-PQQ site-direction benchmark

Use a separately frozen explicit-shell vertical-swap protocol. A site contains
the complete, preregistered directly coordinating residue fragments and only the
preregistered crystallographic inner-shell waters. La and Ca arms have identical
nonmetal coordinates. No PQQ threshold is used. The centered exchange score is
reported against the frozen symmetric CN8 aquo gauge; positive means La-favoring
relative to aquo exchange and negative means Ca-favoring.

Primary falsifiers:

- GGR/MglB P0AEE5, PDB 1GLG, single site: direct same-assay result is about
  29-fold Ca over La; expected score direction is negative.
- Calbindin D9k, two sites in the same protein: canonical C-terminal site is
  La-favoring and pseudo-EF N-terminal site is Ca-favoring. The required common
  holo scaffold/site mapping is frozen before execution.

Locked secondary/site-ensemble challenges:

- Aequorin 1SL8: all three structural sites are scored and reported as a fixed
  ensemble because the weak protein-level Ca preference cannot be assigned to
  one site. No favorable-site postselection.
- Carp parvalbumin 4CPV: both sites are scored separately as supporting
  La-positive evidence; the Ca comparison is cross-study.
- Aqualysin I 4DZT: run only after the experimentally assayed low-affinity Ca
  site is mapped unambiguously.

### 3. Architecture and coupled-cofactor challenges

These are included in the campaign but cannot be forced through either lane above.
Each receives a frozen architecture-specific state before execution:

- PqqT 9B1U: PQQ-bound, no unique crystallographic metal state; site placement
  must be preregistered.
- Six-bladed PQQ proteins: B2, C5, Ca-sGDH P13650, and CcPDH A8P0V4.
- HRP-C: retain both Ca sites and the heme/axial environment.
- ConA: retain the obligatory adjacent S1 Mn while testing S2.
- Mex-LanM: report its three high-affinity sites as an endpoint range only;
  do not claim a local frozen carve reproduces folding/cooperativity.
- Alpha-lactalbumin: two native metal-conditioned geometries crossed with both
  ions; conflicting affinity literature means no truth label.

Lanpepsy is recorded but not treated as an independent-site benchmark because its
multimetal channel and occupancy coupling require an incremental multimetal
protocol. Pike beta-parvalbumin is reserved for within-lanthanide work.

## Frozen decisions and failure logic

1. Existing completed PQQ targets are referenced, not recomputed.
2. No label-dependent structure choice, failed-run deletion, site choice, water
   addition, threshold change, or fragment-radius adjustment is permitted after
   energies are visible.
3. A numerical failure is rerun only unchanged. A chemistry failure creates a new
   protocol ID and repeats the entire affected lane.
4. A direct-shell method that cannot recover GGR's Ca direction and the opposite
   directions of the two calbindin sites has not escaped the anionic-pocket/core-
   charge confound, regardless of PQQ-MDH AUROC.
5. Near-zero controls are judged for uncertainty/calibration, not forced into a
   crisp classification.
6. Results are released by evidence stratum with every invalid or unresolved case
   retained in the ledger.

## Execution policy

Prepared immutable inputs are consumed by one isolated, finite, restartable
SLURM queue on `standard,memory`, excluding the canonical bad nodes. It uses the
full allocated node with allocation-aware MPI, checkpoints each completed leg,
and exits when its durable queue is empty. The new allocation is placed above the
owner's other pending jobs without cancelling running work. A cluster-side watcher
records terminal state first, then wakes Codex thread
`019f6eec-2399-70a3-bf9f-76e17e864b34` with deduplicated retryable events.

The master ledger in `benchmark_manifest.tsv` is authoritative for target status;
lane-specific immutable manifests remain authoritative for atomic coordinates,
charges, selectors, hashes, execution receipts, and results.
