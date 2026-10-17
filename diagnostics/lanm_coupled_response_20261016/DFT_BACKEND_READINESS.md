# DFT replacement: existing whole-scaffold readiness

16 October 2026. Bounded source/artifact inspection only; no scientific computation,
new molecular input, submission or production change. Paths below are repository
relative. `W` means `workspaces/lanm_coupled_response_20261016`; `E` means
`workspaces/metal_environment_response_20260926`.

## Decision

**Mex can reuse its two exact embedded PBE0 origins in the current whole-scaffold
ledger. The LBT48 carve cannot yet be substituted into an equivalent prepared
QM/MM partition.** The ongoing isolated 48-atom reference tests local electronic source
response; it does not supply embedded peptide energies or establish a whole-protein
DFT preference model.

## Mex: two origins already available

The compatible replacement is

`E = PBE0_embedded(Q209; shifted1840PC) + existing_MM_ledger + fixedFF_ALPB(R2049)`.

Replace **both** `MACE_vacuum` and `GFN2_embedded − GFN2_isolated`. The embedded
DFT term already contains interaction with the point charges; add no separate
QM–MM Coulomb, GFN field correction or full-protein FF vacuum energy.

| Exact evidence | Path and SHA256 prefix |
|---|---|
| Completed La_A / Dy_A PBE0 cells | `E/lady_repaired_origins_Mex8FNS_v1/FINAL_COLLECTION.json`, `b35f2d3647c4769d` |
| Their native manifest | same directory, `manifest.json`, `1da171c5fe98b041` |
| Repaired physical/core/field source | `E/lanm_ef3_CboundH_repaired_v2/Mex8FNS/INPUTS.json`, `4b3b98983af0cb5f` |
| Existing successful DFT/current-model join, job 12686 | `W/mex_reference_response_v1/execution/RESULT.json`, `1f3174a2849787c2` |
| Current MM/ALPB assembly, job 12685 | `W/mex_static_v2/execution/RESULT.json`, `ff7b27d490de058d` |

The native recipe is PBE0-D4/def2-TZVP/def2-J, RIJCOSX, DefGrid3,
VeryTightSCF, EnGrad, with the archived embedded `DoEQ false` setting. Both core
states have charge 0, effective multiplicity 1 and 832 explicit electrons; physical
La 1/Dy 6 and ECP 46/55 remain explicit frozen-f distinctions. The exact La/Dy core
XYZ hashes start `7d379007` / `c7e2823a`; both 1,840-point fields are `ff8bc9bf`.
These match current 12569 mappings. Preserve the three fixed Nd spectators, all
waters and 2,049 physical atoms; this is an EF3 substitution, not protein affinity.

The reusable classical rows are
`E/scaffold_origin_transfer_v1/Mex8FNS/{La,Dy}_A.json`
(`1c894f1c` / `d8cdcb0c`). They contain retained bonded terms, MM–MM Coulomb/LJ
and real-QM–MM LJ, explicitly excluding QM–MM Coulomb and C4. Reuse the exact
saved origin ALPB rows `W/mex_static_v2/execution/{La,Dy}/solvent/RESULT.json`
(`a72b423d` / `26ab5c68`). Their common charge vector, `12b03bb0`, is the unshifted
2045 native FF charges plus four formal +3 ions, total +2. This global solvent
vector and the shifted embedding field are distinct existing conventions.

`MEX_REFERENCE_RESPONSE.py::gradients` already replays both origins and assembles
`PBE0_plus_shared_MM_ALPB`. Its pinned `metal_environment_force_assembly`
`prepare_mapping/scatter/assemble` maps real-QM, both cap-anchor Jacobians and
all 1,840 point-charge gradients onto 2,049 physical IDs. Native QM/PC derivatives
are Hartree/bohr; conversion is 627.509474/0.529177210903, with force equal to
negative gradient. Use the DFT `metal_environment_reference.read_pcgrad` parser,
not the xTB-specific external-gradient layout.

**Missing work:** a saved-only origin energy/gradient assembly needs zero new
molecular calls. A changed target, donor, cap or nonzero point charge needs a new
embedded DFT evaluation for each metal. Even the O97 MM-only fixture changes the
DFT field. Existing `MEX_FORCE_CHECK.py` / `MEX_RESPONSE.py` can update physical
coordinates, caps, fields and maps, but no displaced DFT energies on those current
response geometries were established here. The 12686 projections are saved slopes,
not endpoint energies. DFT rows retain `force_consistency_qualified=false` and
`physical_state_qualified=false`; the approximate-model 12709 finite-difference
pass does not qualify this replacement total surface.

## LBT: native assets exist; the 48-atom boundary is not ready

For each source 7CCO/7CCN, `W/lbt3_preparation_v1/execution/<source>/` contains
`PHYSICAL.json`, `PEPTIDE_ATOMS.json`, `peptide_system.xml`, `BONDS.json` and
`MATCHED_CORE_MAP.json`. Preparation result SHA256 is
`649a12d4c78a72058c7446184a67cf252d3b37536c76beb700685713d3e881ec`.
The native ff19SB system covers 255 peptide atoms; the physical system adds target
index 255. Both XMLs hash `cddaa7004e795399`. Global solvent charges are already
defined as peptide −5 plus target +3, holo −2, with no cap-population projection.
The four origin FF-ALPB cells are in job 12715,
`W/lbt_matched_model_v2/execution/RESULT.json` (`2644f6d8bf3dff9b`).

The legacy48 map contains 41 real atoms and seven synthetic H caps, formal
holo charge −1. Its exact source coordinates are linked to full 256, but it is **not**
an existing215-external-atom embedding/ledger:

- Five CB–CA caps and the W9 C–CA cap have ordinary two-anchor definitions.
  Their metadata still needs conversion to the existing physical scatter contract.
- W9 cap 47 (`H_link_C_N`) is placed along the negative normalized sum of C→O
  and C→CA unit vectors. Its `removed_neighbor` is an idealized missing amide-N
  direction, not a physical N ID. Its force depends on C, O and CA. The Mex mapper
  requires retained/omitted real IDs and two-anchor Jacobians; silently treating
  this as an ordinary C–N cap would give an incorrect whole-scaffold derivative.
- No prepared shifted external charge field, boundary charge-closure policy,
  selected bonded terms/exceptions or real-QM–MM LJ ledger was found for this
  carve. The native XML supplies peptide parameters, not a complete metal/MM
  partition. Metal LJ and cross exclusions require an explicit existing-parameter
  binding; adding the full peptide FF energy would double-count internal QM terms.

Thus an embedded 48-atom route needs a declared boundary/charge/interaction construction
and its derivatives before new electronic calls. The isolated 48-atom reference cannot
replace those steps. Keeping Q=whole 256-atom avoids these boundaries and makes the
external/MM terms zero, as in the existing matched model, but whole 256-atom DFT would
be new work; no compatible completed whole 256-atom DFT cells were identified.

## Smallest conditional next step

After the local-reference evidence is available, reuse the Mex origin pair for a
saved-only DFT+MM+FF-ALPB assembly. If pursuing motion, declare a finite derivative
check of that exact replacement before using its forces. Do not rerun solved
origins or import the compact-vacuum SCF recipe into the archived embedded path.
For a balanced Mex/LBT comparison, choose and declare the LBT electronic region
first; mixing embedded-Q209 Mex with isolated 48-atom LBT plus an undefined outer energy
is not a complete comparison. No optimization, new default or larger campaign is
implied by this readiness note.
