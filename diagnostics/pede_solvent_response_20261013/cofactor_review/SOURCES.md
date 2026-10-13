# Chemical-state evidence and reference limits

Read 2026-10-13. These sources identify chemical conventions and reference limitations; they do not validate a PedE score or assign its in-protein redox state.

## Actual folding and scoring records

The original PedE Ca-fold input and the new two-Ca input use

`C1=C(C2=C(C(=O)C(=O)C3=C2NC(=C3)C(=O)O)N=C1C(=O)O)C(=O)O`.

This source string explicitly contains two quinone-region carbonyls and three neutral carboxylic-acid groups. The raw AF3 CIF also retains that descriptor in `_chem_comp.pdbx_smiles`. The canonical La input is included in the frozen comparison. The scoring preparation instead declares `pqq_ox_3minus_v1`, C14H3N2O8, formal−3, singlet; it supplies invariant pyrrole/ring hydrogens only and preserves source heavy coordinates. This is a known fold-time/score-time protonation difference, not evidence of a reduced-PQQ substitution. No protonation or redox population is inferred here.

The exact inputs, source preparations, schemas and geometry archives are pinned by `workspaces/pede_solvent_response_20261013/cofactor_review_v1/manifest.json`. The name correspondence is explicit: generic C5–O1 corresponds to CCD C5–O5; generic C6–O2 corresponds to CCD C4–O4. The complete correspondence is checked as a heavy-atom graph on Slurm. Numeric suffixes alone are not atom equivalence.

## Primary chemical and structural sources

- [RCSB PQQ chemical component](https://www.rcsb.org/ligand/PQQ) describes neutral PQQ as C14H6N2O8 with two keto groups. The downloaded [CCD record](https://files.rcsb.org/ligands/download/PQQ.cif), SHA256 `fe241f87f7f3d3dec74ad761f1d2cfeb75cc5305353b33d075c8a61bd863c4d2`, assigns CCD C4–O4 and C5–O5 double bonds. Its model coordinates are attributed to1CRU, and ideal coordinates are explicitly missing. The CCD supplies chemical connectivity; its stored coordinates are not idealized oxidized-PQQ minima. Carboxyl proton removal changes the formal charge without requiring a change in the quinone reduction state.
- [1H4I and its primary crystallographic paper](https://www.rcsb.org/structure/1H4I), Ghosh et al., *Structure*1995, [DOI10.1016/S0969-2126(01)00148-4](https://doi.org/10.1016/S0969-2126(01)00148-4), explicitly interprets PQQ's nonplanar carbonyl geometry as supporting semiquinone. This is a source-author interpretation, not a new redox assignment by this audit. Consequently1H4I is an experimental structural comparator with a state caveat, not a certified oxidized target geometry.
- [4MAE structural deposition](https://www.rcsb.org/structure/4MAE) identifies PQQ, Ce(III), and PEG. The historical preparation independently records exclusion of the nearby PEG atom15P603/OXT. A dry, metal-substituted local calculation is therefore conditional on a state different from the intact crystal environment. This audit does not claim that the deposition's PQQ component name proves a particular electronic redox state.
- The locally archived1FLG source comes from the [ExaA deposition](https://www.rcsb.org/structure/1FLG) already consumed by the structural-Ca audit. It is an additional saved structural comparator, not PedE itself or a new experimentally assigned La/Ca/redox label.

## Existing research records

`diagnostics/pqq_state_hypothesis_20260922/PROPOSED.md` already distinguishes oxidized PQQ3− from a proposed fixed-proton one-electron-reduced PQQ4−/doublet test. That test remained proposed. `diagnostics/chemical_states_20260919/` studied alpha-lactalbumin internal proton transfer; it did not qualify PQQ bond geometry. The separate electronic-accuracy report left larger PQQ checks unsubmitted.

A total charge and multiplicity input does not impose integer local fragment charges on a self-consistent electronic calculation. The formal preparation state and actual computed charge/spin distribution must therefore remain distinct. The parallel charge/solvent audit addresses density behavior; this branch supplies structural/chemical consistency evidence only.
