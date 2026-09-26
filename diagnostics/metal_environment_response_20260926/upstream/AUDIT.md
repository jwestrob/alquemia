# MACEPOL-EF / MLMM4AMBER capability audit

2026-09-26. Code inspection only; no model inference, installation, molecular
calculation, cluster submission, or author contact occurred.

## Decision

**Continue the embedded-reference scout; the paper-matched ML comparison is
blocked by an unidentified/unavailable trained checkpoint and a missing public
paper-compatible engine interface.** Do not substitute upstream MACE-POLAR.
The public model code does contain the important atomwise potential coupling;
the README's mean-field description is outdated. This is a more promising
starting point than the README alone suggested, but not an executable validated
metal model.

## Pins and access

| Component | Actual inspected pin |
|---|---|
| ClickFF/MACEPOL-EF, main | `e2b0aeed27c2822790a6994c9a64327e5f131158` |
| ClickFF/MLMM4AMBER, feat/boundary | `4af44eb9a9ce91428dad142cbaa8bdd598e7a09f` |
| MLMM master inspected interface | `140fcc6f7470169ac070441ae95d35b752694650` |
| MLMM sander23 inspected interface | `5d83e48096c957fb597f445f40434691c5c1fb60` |
| Fine-tuned checkpoint SHA256 | **Unavailable: no checkpoint obtained** |

Shell GitHub DNS failed; GitHub connector fetches supplied API trees and pinned
raw source. Local copies live in
`workspaces/metal_environment_response_20260926/upstream/`.
`pins.json` records commit metadata; `source_files.json` records downloaded-file
hashes. API tree responses were not truncated. Both repositories have no release
assets. The model main/evolved trees contain legacy MACE foundation binaries,
not an identified MACEPOL-EF checkpoint. Other named model branches remain
unqualified, not assumed to be the paper release.

## Paper, checked against Methods 5.1–5.3

The author-hosted full-text mirror describes closed-shell training on H, B, C,
N, O, F, Si, P, S, Cl, As, Se, Br and I. Reference calculations used ORCA 6.1.0,
ωB97M-V/def2-TZVPD, def2/J, RIJCOSX, DEFGRID3 and self-consistent VV10 for
gradients. The architecture adds atomwise charge–potential coupling and a
field-responsive head to MACE-POLAR. Its electrostatic ML energy already
includes interaction with MM charges; another classical Coulomb sum would
double-count it. MM-internal and cross Lennard–Jones terms remain classical.
Forces require derivatives through both scalar potential and vector field,
plus intrinsic coordinate derivatives. The paper describes Gaussian-smeared MM
charges, bonded exclusions/1–4 scaling, and link sites in the electronic region. These
claims do not qualify metal response or the publicly accessible engine.

Source: [paper full-text mirror](https://www.researchgate.net/publication/414451620_High-throughput_physics-based_enzyme_engineering),
DOI [10.64898/2026.09.15.751901](https://doi.org/10.64898/2026.09.15.751901).
Direct bioRxiv access returned HTTP403; no inaccessible figures were interpreted.

## Actual model interfaces: consequential details

Code links below use the model pin above. Local line numbers match the downloaded
source after transport-newline normalization checked against Git blob hashes.

- `mace/modules/extensions.py:680`: `PolarMACE.forward` accepts a data dictionary,
  including `external_field`; line754 reads `external_potential`. Both graph and
  per-atom potential shapes are accepted, but atomwise qV is used in the
  **per-atom field branch** (lines1085–1129). Supply `[N,3]` fields and `[N,1]`
  potentials explicitly. Merely giving graph-shaped fields is not equivalent.
- Line1100 computes `q_atom * external_potential_per_atom_scalar`. Missing
  potential at this low level permits a graph-mean-field fallback. Do not use
  that fallback for a finite charged protein environment.
- `mace/data/atomic_data.py:389–404` negates physical fields at input. Public
  units are V/Å for field, V for potential, eV for energy, Å for positions and
  eV/Å for forces. Potential is not negated. The environment variables
  `FORCE_NO_NEGATE` and `FORCE_DISABLE_QV` alter the physics and must be recorded
  and forbidden in the candidate runner unless explicitly part of its protocol.
- **Spin trap:** `extensions.py:900–901` uses `Q ± (total_spin - 1)` to constrain
  the two spin channels. This implementation expects multiplicity, not 2S.
  A closed-shell input is1. The README and exporter smoke fixture's0 examples
  cannot establish a correct physical spin. Actual checkpoint training metadata
  remains necessary, especially for Dy.
- **Potential trap:** ASE calculator default `info_keys` (line139) maps charge,
  spin and field, but omits potential. AtomicData passes extra properties through
  (line450), then **fills absent external_potential with zeros** (line468).
  Consequently the default ASE path can silently use qV with zero potential.
  A future adapter must demand explicit phi, supply a custom property mapping or
  a direct batch, and test the actual received arrays.
- **Cache trap:** `MACECalculator.check_state` (lines340–365) ignores ndarray
  values in `atoms.info`; field/potential-only changes may return cached results.
  Clear results on every embedding evaluation or use a separately keyed adapter.
  Scalar charge/spin changes are inspected, but no numerical cache test was run.
- `extensions.py:1157` differentiates energy against local positions. With
  separately supplied phi/field tensors this supplies intrinsic local forces,
  not the complete moving-environment derivative. The adapter must additionally
  differentiate with respect to phi and field and propagate both to local and
  MM coordinates. No MM atoms enter the model's default coordinate list.
- `scripts/convert_polar_to_pt.py:224` defines `PolarMACETraceWrapperEF`, with
  explicit potential/field inputs and differentiable returned energy. Use
  `--polar-schema ef`; the standard schema is field-only. Export defaults to
  float32, so precision must be chosen and qualified. Traced energy/density
  outputs require caller-side autograd; returned force arrays are not supplied.
- The architecture's `atomic_numbers` and learned embeddings determine actual
  executable elemental support. Without the fine-tuned weights, Ca/La/Dy
  inclusion cannot be verified. Organic training restrictions do not prove
  exclusion from the executable, nor does upstream elemental coverage prove
  accurate field response.

Pinned source: [extensions.py](https://github.com/ClickFF/MACEPOL-EF/blob/e2b0aeed27c2822790a6994c9a64327e5f131158/mace/modules/extensions.py),
[AtomicData](https://github.com/ClickFF/MACEPOL-EF/blob/e2b0aeed27c2822790a6994c9a64327e5f131158/mace/data/atomic_data.py),
[calculator](https://github.com/ClickFF/MACEPOL-EF/blob/e2b0aeed27c2822790a6994c9a64327e5f131158/mace/calculators/mace.py),
[exporter](https://github.com/ClickFF/MACEPOL-EF/blob/e2b0aeed27c2822790a6994c9a64327e5f131158/scripts/convert_polar_to_pt.py).

## Checkpoint identity and engine blocker

`scripts/eval_macepol_0428.py:25–44` lists private `/data2/xuw74/...` checkpoint
paths. Its production tag `macepol_ef` identifies an April26 best_ep1 checkpoint
with **no field negation**, while later qV candidate checkpoints use negation.
This contradicts assuming every mentioned model follows the README's present
convention. No public URL/hash identifies the paper's selected checkpoint.

The inspected engine's `sander/mlp_models.hpp` registers ANI, MACE-OFF,
MACE-OMOL, Egret, AIMNet, SpookyNet and EANN, but no MACEPOL-EF/POLAR backend.
`mlp_potential_forward` in `sander/mlp_potential.cpp:96` takes only species,
coordinates and atom count. Its separate charge/multiplicity globals are for
AIMNet. It has no phi/field arguments or their response outputs. Master and
sander23 interfaces likewise lack the paper coupling. Existing `mlp_ewald.F90`
and boundary utilities are useful code, not evidence of the new model's complete
derivative implementation.

Pinned engine: [model registry](https://github.com/ClickFF/MLMM4AMBER/blob/4af44eb9a9ce91428dad142cbaa8bdd598e7a09f/sander/mlp_models.hpp),
[forward API](https://github.com/ClickFF/MLMM4AMBER/blob/4af44eb9a9ce91428dad142cbaa8bdd598e7a09f/sander/mlp_potential.cpp).

## Licensing and tests

Model source LICENSE.md is MIT. Engine README claims CC BY4.0, but the actual
`sander/LICENSE` is BSD3-Clause: retain both records and clarify the scope before
redistributing a bundled engine. Base AMBER distribution rights are separate.
No checkpoint license was supplied; source licensing cannot be assumed to cover
unreleased weights. Paper mirror states CC BY-NC4.0.

Downloaded Python sources passed AST parsing. File hashes and pinned Git blob
checks are recorded. These are source integrity checks, not molecular validation.
Model loading, elemental support, response/gauge/rigid-transform/force checks
and performance are all **unrun/unavailable** pending the exact weights.

## Unsent author question

Could you share the checkpoint and SHA256 used for the MACEPOL-EF results in
10.64898/2026.09.15.751901, its redistribution/use license, and the matching
MLMM4AMBER commit or branch implementing potential/field sensitivities? We find
qV support in MACEPOL-EF e2b0aeed, but public engine4af44eb9 has no POLAR backend
or phi/field API. Which field sign and spin convention applies to the released
checkpoint, and does it retain Ca, La and Dy embeddings despite the nonmetal,
closed-shell field-response training set? A minimal finite charged-system input
and expected energy/force output would help qualify the coupling independently.

No message was sent. Next scientific action remains the independent embedded
Ca/La reference scout; no large ML/MM integration is justified yet.
