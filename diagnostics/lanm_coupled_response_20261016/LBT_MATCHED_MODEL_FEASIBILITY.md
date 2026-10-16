# LBT3 comparison under the Mex energy rule: feasibility

16 October 2026. Existing-file/source inspection only. No numerical analysis,
molecular call, submission, parameter change or source-artifact mutation was
performed. This note prepares a possible test; it does not launch a campaign or
replace the separately pending Mex force qualification.

## Finding

**The four original whole-LBT3 states support a small comparison using the same
partitioned energy rule and the same physical FF-charge solvent convention as
Mex.** The simplest explicit partition puts the entire 256-atom peptide–metal
system in the electronic region Q. Its external region is empty:

`E = MACE_vac(Q) + [GFN2_emb(Q;empty) − GFN2_vac(Q)] + MM(empty) + G_ALPB(qFF,R)`

Thus the field correction and separate MM ledger are algebraically zero, and
`E_LBT = MACE_vac(whole) + G_ALPB(qFF,whole)`.

These zeros must be recorded as consequences of the declared partition, not as
completed numerical evaluations. No GFN2 call is needed to subtract a method
from itself. Do not add the full peptide FF vacuum energy: its internal mechanics
already belong to whole-peptide MACE. FF parameters supply solvent charges only.
This uses the same energy construction as Mex, but different region boundaries;
it does not establish partition independence or equivalence to whole-protein MACE.

## Existing evidence and exact sources

Repository-relative paths below are immutable existing artifacts.

| Evidence | Existing path / identity |
|---|---|
| Whole preparation, job12469 | `workspaces/lanm_coupled_response_20261016/lbt3_preparation_v1/execution/RESULT.json`, SHA256 `649a12d4c78a72058c7446184a67cf252d3b37536c76beb700685713d3e881ec` |
| Four whole MACE energies and forces, job12500 | `workspaces/lanm_coupled_response_20261016/lbt3_forces_v1/execution/RESULT.json`, SHA256 `4ea82417d3a621b5b267e024c2c6229da5d83d1c9d969ffaf12fb816bdb6024d` |
| Fixed-charge native solvent interface qualification, job12566 | `workspaces/lanm_coupled_response_20261016/frozen_solvent_force_v1/execution/RESULT.json`, SHA256 `1301d3b555098574f96c75ee6beb162c5932adee480662f55ac3ae5213d12c0c` |
| Completed Mex static assembly, job12685 | `workspaces/lanm_coupled_response_20261016/mex_static_v2/execution/RESULT.json`, SHA256 `ff7b27d490de058d6fbc494ff524a144f533e63e754f6dd3282d58c8fc904d0e` |

For each source `s = 7CCO, 7CCN`, the preparation directory
`workspaces/lanm_coupled_response_20261016/lbt3_preparation_v1/execution/s/`
already contains `PREPARATION.json`, `PHYSICAL.json`, `PEPTIDE_ATOMS.json`,
`peptide_system.xml`, `BONDS.json`, `HYDROGEN_AUDIT.json`, `completed_peptide.pdb`
and exact `whole_La.xyz` / `whole_Dy.xyz`.

- Both are model1 of the consumed 17-residue sequence `FIDTNNDGWIEGDELLA`.
  7CCO is La-conditioned; 7CCN is Lu-conditioned, not an observed Dy structure.
  They have the same atom identities/topology and modeled protonation, but are
  different conformations. Retain their separate coordinate frames; the scalar
  exchange comparison needs no alignment or coordinate interpolation.
- Each physical map has 255 peptide atoms followed by target index255, ID
  `metal`, placeholder element `M`. Replace only that placeholder with the actual
  La or Dy species. No caps, external atoms, waters or spectators occur.
- `PEPTIDE_ATOMS.json` already stores native `charge_e`, ID, element and the exact
  prepared coordinates for all255 peptide atoms. The preparation source obtains
  these from the unconstrained ff19SB `NonbondedForce` and checks peptide charge−5.
  Both `peptide_system.xml` files have SHA256
  `cddaa7004e7953998c739ebcdf2fd01071c92a82f908bc353530ed6a416da8e7`;
  the ff19SB XML pin is `086bfdb1c05e7d5cb1d330e6c39fab5056a95e89492cdbd3044c746db13e9d09`.
- Use those unmodified255 FF charges plus formal target charge+3, giving holo−2.
  Join by physical ID and assert the existing ordering; never use archived GFN2
  populations, infer charges from `PHYSICAL.json`, or renormalize the vector.
  The La/Dy charge vectors must be identical, while their cavity element differs.
- Terminal completion added only H2/H3 at Phe1 and OXT at Ala17. The declared
  free NH3+/COO− termini and deprotonated D3/D7/D13/E11/E14 remain assumptions.
  `LBT3_PLAN.md` retains the primary SI peptide−4 discrepancy; this note does not
  resolve experimental termini or authorize a new protonation assignment.
- Job12500 contains four **flat rows**, with `task_id`, `xyz`, `physical_source`,
  `energy_eV`, `forces_eV_A` and `input_state_check` directly on each row, not
  under `row.task`. All are complete, whole256/q−2/target255, physical La1/Dy6.
  Its native OMOL checkpoint/head/float64 lineage is the one reused by Mex12685.

## Minimal missing preparation and possible later cost

1. Freeze a four-cell manifest: 7CCO/7CCN × La/Dy, original256 coordinates only.
   Verify source/FF/model/runtime pins, all physical IDs, exact paired XYZs and
   the existing MACE state receipts. Preserve all four cells and both sources.
2. Export the256-entry `qFF` vector with17-digit round-trip precision, plus its
   per-ID provenance. Check closure−2, native FF/physical ordering, identical
   La/Dy vectors and consistent peptide charges across the two source frames.
   Record whole-system Q, empty external/MM sets and zero cross-boundary terms.
   No new topology, hydrogen preparation, caps, ion LJ parameters or force field
   are required for this partition.
3. Prepare a narrow fixed-charge native ALPB adapter using the already reviewed
   256-atom wrapper `FROZEN_SOLVENT_FORCE.f90` and12566's pinned library/build.
   Use the same water/P16/gsolv/230-point/zero-salt settings as Mex; FF charges
   are exact declared model values, so do not label them as rounded GFN charges.
   The current generic wrapper's archived-charge error-bound bookkeeping needs
   that explicit adjustment, as in `MEX_STATIC_SOLVENT.f90`; its energy functional
   and gradient routine do not change.
4. Re-admit/reuse all four12500 MACE results. The missing fixed-origin panel would
   require **four new ALPB evaluations, zero new MACE/GFN2/MM/DFT calls** if exact
   reuse passes. This is a prospective ceiling, not an execution receipt.

Job12566 evaluated saved GFN charge vectors, and12510 evaluated self-consistent
GFN2 water-minus-vacuum corrections. Neither supplies the required FF-charge
LBT solvent energies. Old accommodated or one-water composite scores also cannot
be substituted. No modeled water or relaxed endpoint is added by this proposal.
The legacy48-atom donor carve has charge−1 and synthetic caps; using it instead
would require a new explicit partition/field/MM ledger and is not the minimal route.

## Balanced exchange and interpretation

For each original LBT source `s`, define

`R_s = [E_Mex,Dy − E_Mex,La] − [E_LBT,s,Dy − E_LBT,s,La]`.

This is the energy of
`Mex–La + LBT–Dy → Mex–Dy + LBT–La`.
Both sides retain one La, one Dy, all scaffold/cap atoms, the same Mex waters and
three Nd spectators, and total charge(+2)+(−2)=0. Positive R favors the left
side within this conditional model: relatively greater La preference in Mex.
Take balanced differences separately for MACE, field correction, MM and ALPB
before conversion/summation; the common elemental offsets cancel. Retain R_7CCO
and R_7CCN separately, including disagreement or unavailable cells, without
selecting the favorable source. A future common-source menu needs its own rule.

This is a fixed-origin, partition-dependent exchange diagnostic. Mex is one EF3
state at fixed three-Nd occupancy; LBT is a dry represented peptide with uncertain
experimental termini. Mex boundary mechanics, missing solvent configurational
sampling, unequal hydration inventories, MACE physical-spin versus GFN2 frozen-f
representations and unqualified total derivatives remain limits. No ion hydration
reference, apo folding free energy, water chemical potential, site summation or
experimental binding affinity is obtained. Old whole-LBT self-consistent-solvent
offsets and calibration thresholds do not transfer.

Both source families and their qualitative preference evidence are already
consumed development material. LBT is an engineered peptide, not a second natural
LanM; the two assays do not provide a matched quantitative exchange target. No
reserved numerical outcomes, new labels or reserved structures were inspected
for this feasibility note. It adds no blind-validation claim.
