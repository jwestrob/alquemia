# Consumed experimental target for LanM La/Dy work

16 October 2026. Read-only evidence recovery plus this note; no numerical
calculation, molecular call, job, or reserved-outcome access.

## Target and smallest comparison set

The defensible primary ordering is **stronger apparent La-over-Dy selectivity
for wild-type Hans-LanM than for wild-type Mex-LanM at pH 5**. Define the
experimental selectivity descriptor symbolically as
`S(P) = Kd,app(P,Dy) / Kd,app(P,La)`; the qualitative target is
`S(Hans) > S(Mex)`. No ratio or free energy is calculated here.

This is one consumed protein-level comparison. It does not imply that Hans has
tighter absolute La binding than Mex, assign separate labels to EF1/2/3, or
validate a single-EF3 score. The measured response includes coupled occupation,
folding and, for Hans, assembly. The completed fixed-source electronic comparison
remains a conditional mechanism diagnostic.

| Construct | Role in a small comparison | Available structural basis |
|---|---|---|
| Wild-type Hans | Primary selectivity target | Retain both consumed 8DQ2 and 8FNR sources; these are correlated states of one protein |
| Wild-type Mex | Primary reference | Consumed 8FNS source; declare its crystallographic construct and Nd conditioning |
| Hans R100K | Optional, more controlled mechanistic comparison with wild-type Hans | Same-assay experimental reference below; no prepared mutant source was found in the scoped inventory, so a separate versioned preparation is needed |

R100K tests whether a model can retain the near-unchanged Dy response while
weakening the La response. It is an auxiliary mechanistic check, not another
independent wild-type protein. No larger library is needed to define these
questions, and no experiment is launched by this note.

## Traceable primary measurements

| pH-5 CD apparent affinity | La, pM | Dy, pM |
|---|---:|---:|
| Hans wild type, Table S2 | 68 ± 7 | 2,600 ± 700 |
| Hans R100K, Table S7 | 120 ± 10 | 2,700 ± 400 |

Uncertainties are SEM: three independent titrations for each WT metal and
R100K La; two for R100K Dy. Conditions: 15 µM monomer,
25 °C, 20 mM acetate, 100 mM KCl; chelator-buffered CD at 222 nm, Hill fits,
EDTA for La and EGTA for Dy. Dy's main response is incomplete, with an additional
weak response beyond the accessible range. These are apparent conformational
affinities. [Primary Supplementary Methods/Tables S2 and S7](https://media.springernature.com/original/springer-static/esm/art:10.1038%2Fs41586-023-05945-5/MediaObjects/41586_2023_5945_MOESM1_ESM.pdf)

The [2023 primary article, Fig. 1e](https://www.nature.com/articles/s41586-023-05945-5)
compares earlier Mex pH-5 measurements from
[Deblonde et al., 2020](https://doi.org/10.1021/acs.inorgchem.0c01303).
This is not one contemporaneous assay batch. Exact Mex La/Dy values,
uncertainties and full matching conditions remain unverified here: the official
[Fig. 1 source workbook](https://media.springernature.com/original/springer-static/esm/art:10.1038%2Fs41586-023-05945-5/MediaObjects/41586_2023_5945_MOESM11_ESM.xlsx)
contains Hans/R100K results only. Do not substitute pH-7.2 Mex values or fit a
quantitative Hans-minus-Mex target from a plotted estimate.
The unavailable matched Mex numbers do not block the qualitative ordering above.

## Keep observables and physical states separate

- **Assembly:** Hans La dimers are much tighter than Dy dimers; the reported
  bounds/values are `<0.4 µM` and approximately `60 µM`. The ITC experiment used
  pH 7, 30 °C, MOPS/KCl and two metal equivalents. Mex is monomeric in the
  reported comparisons. These are dimer dissociation quantities, not metal Kd.
- **Occupancy:** solution loading, the two-equivalent CD transition and crystal
  occupancy are different observations. EF4 is weak; four deposited sites do
  not establish a universal four-ion assay state.
- **Separation/capture:** the binary-column result concerns Nd/Dy, immobilized
  constructs and competing solution. Hans loading changes upon immobilization.
  Neither separation factor nor captured equivalents is a La/Dy solution Kd.
- **Fluorescence:** intrinsic-fluorescence displacement and metal-sensitized
  luminescence are assay-dependent signals. A signal change alone is not an
  independent metal-affinity measurement. [Primary Methods and Extended Data](https://www.nature.com/articles/s41586-023-05945-5)

The current model represents conditional monomers. Its actual source states are
Hans8DQ2: QM/total waters 2/74, spectators La/La/Na; Hans8FNR: 3/179,
Dy/Dy/Dy; Mex8FNS: 7/167, Nd/Nd/Nd. Electronic spectators are formal charges;
the matched classical ledger also contains ion-specific pure 12–6 terms.
Common protein identities between the Hans sources do not make their hydration
or occupancy states interchangeable. No source population follows from their
raw total energies.

## Assembly evidence and the monomer model's limit

The experimental mechanism couples coordination, dimerization and the
heterogeneous Dy conformational response. It does **not** identify dimerization
as the sole or quantitatively dominant source of the full La/Dy affinity gap:
R100K strongly increases monomerization while changing the La response much less
than that full gap and leaving Dy nearly unchanged. Partner Arg100 contacts
EF3 Asp85/Glu91/Asp93; that interaction network is absent from the current
isolated-chain model. [Primary structural and mutagenesis evidence](https://www.nature.com/articles/s41586-023-05945-5)

| Deposited structure | Biological assemblies in the actual source metadata | Implication |
|---|---|---|
| 8DQ2, La, pH 7 | Assembly 1: A+C; assembly 2: B+D; both dimers | Chain A alone is an extracted subunit |
| 8FNR, Dy, pH 7 | Assembly 1: A+B; assembly 2: C+D; both dimers | No deposited matched dissociated-monomer alternative |
| 8FNS, Nd | Assembly 1: monomer | Reference construct is structurally distinct |

These assignments are present in local `8DQ2.cif` assembly categories and
`8FNR.pdb` REMARK 350, and agree with
[RCSB 8DQ2](https://www.rcsb.org/structure/8DQ2) and
[RCSB 8FNR](https://www.rcsb.org/structure/8FNR).
In 8FNR, EF1 occupancy/order differs across chains; alternative dimers therefore
cannot be treated as matched-state replicates. The article's deposited structure
set is 8DQ2/8FNR/8FNS, all wild type. No R100K structure was found in that set,
the scoped local inventory, or a bounded RCSB search; this is not an exhaustive
claim about every repository.

**Bounded implication:** a matched monomer/dimer comparison could isolate the
partner's conditional contribution, but extracting a chain is not an observed
equilibrated monomer. The current monomer EF3 model can test local accommodation;
it cannot determine dimer populations or the complete multi-site folding/occupancy
response. Adding a partner is a specific missing-physics hypothesis, not proof
that it will fix the ordering.

## Local evidence and boundaries

Repository-relative paths:

- Archived original article: `workspaces/benchmark_set_20260915/evidence/lanthanide/hans_primary.{xml,txt}`;
  capture index: `workspaces/benchmark_set_20260915/evidence/lanthanide/source_captures.json`.
- Prior consumed La/Dy interpretation: `diagnostics/lanm_series_followup_20260923/REPORT.md`
  and `DY_TRANSFER_REPORT.md`; completed six-origin state/pins:
  `diagnostics/claude_lanm_20261012/03_comparison_root_recovery/{FINAL_REPORT.md,FINAL_RESULT.json,FINAL_PROVENANCE.json}`.
- Exact prepared contracts:
  `workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/{Hans8DQ2,Hans8FNR,Mex8FNS}/INPUTS.json`.
- Actual crystallographic source files:
  `workspaces/benchmark_set_20260915/sources/pdb/8DQ2.cif`,
  `hans_lanm_protenix_benchmark/references/8FNR.pdb`,
  `workspaces/benchmark_set_20260915/evidence/lanthanide/8FNS.cif`.
- Construct and assay inventory: `diagnostics/spicy_lams_inventory_20260923/{REPORT.md,EVIDENCE.md}`;
  metadata-only structure index:
  `workspaces/spicy_lams_inventory_20260923/inventory_v2/STRUCTURES.json`.

The SpyCI-LAMBS inventory describes a distinct pH-3, immobilized competitive
capture/ICP-MS assay. Its numerical ion profiles and PLM holdouts were not
opened. Their labels cannot be imported into this pH-5 solution-affinity target.
The inventory's previously documented limited qualitative exposure remains
recorded; this note makes no new blanket-blindness claim.

**Interpretation rule:** a coupled-response pilot can test a specified physical
contribution to the consumed ordering. A protein-level prediction requires a
declared treatment of sites, occupancies, hydration, protonation, assembly and
their uncertainties; neither a favorable EF3 displacement nor summing isolated
site energies supplies that treatment automatically.
