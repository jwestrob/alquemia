# Intact monomeric peptide response reduces—but does not resolve—source sensitivity

The complete four-search/shared-pool pilot succeeds numerically and preserves
covalent geometry. Source-local shared-pool selection, using **Dy minus La**,
reduces the absolute difference between the two LBT3 source contrasts from
**72.966739 to57.895436kcal/mol**. This is an improvement in structural sensitivity
within one model, not robust lanthanide classification or experimental transfer.
The global shared menu continues to select a7CCO-derived geometry for La and a
7CCN-derived geometry for Dy. A single global menu removes source choice by
construction; the source-local comparison is the actual robustness diagnostic.

## Executed comparison

All256 atoms remain represented in native MACE plus author-GFN2(ALPBwater−vacuum).
Both original sources and all four admitted terminal geometries form six distinct
candidates. The complete12-cell paired matrix reuses eight cells and adds only
four opposite-metal composite evaluations. Exact state/model/geometry pins and
all component values are in POOL_RESULT.json; no alignment or approximate
coordinate deduplication was used. Experimental terminal chemistry and absent
explicit NMR waters remain conditional assumptions.

| Source | La own accommodation | Dy own accommodation |
|---|---:|---:|
| La-conditioned7CCO | −2.880329 | −5.182784 |
| Lu-conditioned7CCN | −20.626547 | −7.857699 |

Works are kcal/mol on each metal's own complete energy. The expanded global pool
changes its Dy-minus-La contrast by−4.977370kcal/mol relative to the original
shared two-source menu. This is a conditional energy change; no hydrated-ion
reference, affinity, population or new experimental preference is inferred.

## What limits the response

All four searches stop because the next allowed line-search steps exceed the
same0.20Å displacement envelope. The limiting atoms are Asp7OD1/OD2, rather than
newly damaged bonds or contacts. No search reached a stationary minimum.
The Lu-conditioned La case retains large active forces and omitted donor/backbone
loads, including Asp3CB–CG and the Trp9/Ile10 backbone junction. Thus the current
small motion set and envelope leave genuine model strain unresolved.

Full residual forces also include internal carboxylate responses. For example,
7CCO Asp7 has large opposing bond-parallel force descriptors even though its
C–O lengths remain near the force-field reference. These descriptors do not prove
bad deposited hydrogens, incorrect experimental structures or faulty MACE.
They argue against explaining all residual force by a frozen protein backbone,
and against allowing unrestricted bond deformation without local qualification.
Exact limiting atoms, donor distances, all64 mode loads and bonded-force tables
are preserved by the saved-only diagnostic12542, referenced in DIAGNOSTIC_SOURCE.json.

## Costs and next decision

Search12533:107s; remaining-three12538:225s; cross-pool12541:53s. Each used one
H200 and12CPUs. Saved diagnostic12542 used1s/1CPU. Together this accommodation
round used32 new MACE evaluations,64 standaloneGFN2 evaluations,385 allocated
GPU-seconds and4621 allocatedCPU-seconds, including the diagnostic. All new
molecular calls succeeded. Whole monomer-branch cumulative allocation through
this point is5942CPU-s/433GPU-s; prior failures/preparation are retained. Allocation
is distinct from observed CPU work or backend latency. Exact scheduler records
and per-backend costs are in the adjoining accounting/results files.

The next experiment should distinguish insufficient allowed motion from missing
mechanical coordinates, retaining both sources/metals and these completed results.
Do not adjust a classifier threshold to this two-source peptide or repeat the
completed searches. No PQQ change, DFT campaign, biological classifier or natural
LanM validation is implied. The active protein-level preference goal remains unmet.
