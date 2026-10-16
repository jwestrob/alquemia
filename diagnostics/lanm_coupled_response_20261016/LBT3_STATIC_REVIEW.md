# LBT3 finite scout: independent static review

16 October 2026. Read-only source/evidence review plus this note. No scientific
execution, imports of the molecular stack, numerical calculations, or submissions.

## Current disposition

**No blocking issue found in the corrected preparation, executed preparation
receipt, or bound eight-cell MACE runner/configuration.** The original
preparation/plan hydrogen mismatch is resolved. This review does not interpret
new MACE energies or replace their scientific collector.

Reviewed runner `LBT3_MACE.py` SHA256:
`a71e19843396ac1c3d487d347c58fedc78bbac1bd8a7b2a931be6dd772df4002`.
The final runner adds explicit runtime-version and Python identity checks to
the previously reviewed `8e196a04` source.
Reviewed wrapper `LBT3_MACE.sbatch` SHA256:
`63d84ead84c0827a9f0f37a17aac85ad5b856d37aa02ac20a4d40596d576cd43`.

Reviewed corrected preparation `LBT3_PREPARE.py` SHA256:
`f0dbe4436523985cecd5607983ca634710d3df76a4b1c0751feeb2ac77cf3cc6`.
Reviewed corrected plan `LBT3_PLAN.md` SHA256:
`d2263867bcfe6e2282f69a20e37751b68586301fd9213952708e55a84364dd0d`.

Actual preparation job 12469:
`workspaces/lanm_coupled_response_20261016/lbt3_preparation_v1/execution/RESULT.json`,
SHA256 `649a12d4c78a72058c7446184a67cf252d3b37536c76beb700685713d3e881ec`.
Both sources pass; eight tasks are prepared. Read-only verification of 37 file
pins found no mismatch. Saved reports preserve all inherited H/heavy atoms,
16 peptide bonds, and exact shared real-atom maps. Two native missing-H calls
completed; internal iteration counts remain unknown, with zero MACE/QM calls.

Bound MACE configuration:
`workspaces/lanm_coupled_response_20261016/lbt3_mace_v1/CONFIG.json`, SHA256
`3e47b08f6de893783e72c24e90809b5ccdef653c27b111c892dc12e7d28bd17e`.
All 314 referenced file pins match, including the checkpoint. Its model object
equals the native precedent. All eight tasks equal actual prepared tasks apart
from added cache-key metadata; the saved Mex receipt is pinned at
`a5c6946a33887456fb3303396ee8c5388dab30caa219cde6fd416db0d6c10dd3`.

## Findings and resolution

1. Initial `LBT3_PREPARE.py` SHA256
   `b6109f5dd5c19bfeddd2cf9c2756354f0b91263e8138c87d676828c40151d254`
   called `mace_global_prepare.protein()`, which invokes radial hydrogen
   normalization. The plan instead required every inherited hydrogen coordinate
   to remain unchanged. Its final equality gate covered only heavy atoms.
   Corrected source removes `protein()` and exports the native topology directly.
   It checks all 252 inherited H/heavy atoms after missing-H completion and again
   after serialization. Only OXT/H2/H3 may be added; hydrogen bond lengths are
   diagnostics. The exact archived `build_cluster()` uses these shared physical
   coordinates, with separate synthetic-cap mapping and four tasks per source.
   This resolves the static mismatch; actual preparation 12469 passes its gates.
2. The primary literature does not establish the exact terminal chemistry.
   Main MD methods use holo charge −2 and negatively charged Asp/Glu, whereas
   SI Fig. S5 states peptide charge −4 at pH 6. Neither the experimental
   Materials/ITC methods nor SI sequence figure documents amidation or another
   terminal modification. Thus free NH3+/COO− termini with peptide −5/holo −2
   are a declared preparation hypothesis, matching the main MD total, with the
   conflicting SI statement retained. Missing OXT does not establish amidation.
   The corrected plan explicitly retains this discrepancy and preparation
   assumption. This limits interpretation; it does not block a conditional vacuum diagnostic.
   [Primary article](https://doi.org/10.1038/s41598-020-76527-y),
   [Supplement, Figs. S1/S5](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41598-020-76527-y/MediaObjects/41598_2020_76527_MOESM1_ESM.pdf).

## Native model and state compatibility

- Exact precedent:
  `workspaces/lanm_series_followup_20260923/prepared_v2/manifest.json`, SHA256
  `6eab0d244f99b07073e4d2c00e75793f2dc280bdb9ca0985a57bd359a081eeb4`.
- Current `scripts/lanm_series_followup.py` SHA256
  `433fdf85448182841a967a6855d2941b55976d4f429beacf77d0da802bb299c9`
  and `scripts/mace_omol.py` SHA256
  `1194c14bb0ba1130212f8a365228808475694e126dbaebefc370a24dc6d981ef`
  match that manifest's implementation pins.
- Checkpoint SHA256
  `9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`,
  `ScaleShiftMACE`, native `omol` head, float64, nonperiodic vacuum,
  no additional dispersion. Its charge/spin categorical embeddings admit the
  declared states; `input_batch()` checks actual batch values.
- Existing convention is **spin input = multiplicity**: La 1, Dy 6.
  The new runner preserves it and checks all-electron parity; this does not
  validate Dy's spin-orbit or f-electron physics.
- New runner checks the exact two-source/two-representation/two-metal denominator,
  paired atom/coordinate equality, permitted elements, charges −2/−1, actual
  loaded model/head/units, native state batch, finite energy and absence of forces.
  It uses one warm model and stops subsequent calls after a failure.
- The wrapper uses the original MACE Python environment, one GPU and 12 CPUs,
  without memory/time overrides. Runtime versions are recorded. Before reuse,
  compare them with existing Mex receipts: torch 2.8.0, mace-torch 0.3.16,
  ASE 3.29.0; do not silently reuse if the method differs.

The actual reusable Mex vacuum receipt is
`workspaces/lanm_series_followup_20260923/prepared_v2/mace_remaining/result.json`:
all six EF1–3/La–Dy rows are complete, 44 atoms each, charge −1, with the native
1/6 multiplicity checks. These are the small carved MACE endpoints, distinct
from the later hydrated/embedded EF3 QM preparation with charge 0. Reuse their
saved vacuum energies only; no native GFN2 solvent component is qualified for Dy
by the corrected-La registry.

## Scientific interpretation

The consumed LBT3 ITC table directly measures La and Dy and supports Dy-over-La
binding at pH 6 and 10 °C. Deposits 7CCO/7CCN are La/Lu-conditioned NMR ensembles,
with coordination restraints; model 1 is an explicit fixed-source choice.
Absence of deposited waters does not establish a dry experimental site.

The whole-versus-fragment contrast of Dy-minus-La energies is a balanced
double difference: the common elemental offsets cancel despite the different
representation charges. It measures the complete representation change,
including omitted D13 charge and changed capping, not an isolated mechanical
scaffold effect. Likewise the specified Mex-minus-LBT3 exchange sign is correct:
positive indicates greater conditional relative light-ion preference in Mex.
Neither construction is an experimental binding free energy.

Retaining individual-source results and labeling any source-pool minimum as a
mathematical summary avoids assigning populations to two NMR conformers.
The full-protein Mex observable has no independent EF1/2/3 labels; the proposed
vector must remain unsummed conditional site diagnostics. Different assay
conditions, hydration and folding/occupancy remain explicit. These limitations
are already appropriate to a finite MACE-vacuum representation scout.
