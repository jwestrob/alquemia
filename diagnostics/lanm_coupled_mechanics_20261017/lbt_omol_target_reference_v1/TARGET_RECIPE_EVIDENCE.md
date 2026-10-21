# OMol25 target recipe: whole 256-atom La peptide

20 October 2026. Source, documentation, JSON and hash inspection only. No scientific imports, model loading, input generation or molecular calculations were performed.

## Decision supported by this evidence

The published method and its external orbital-basis coefficients are available. A separately declared **ORCA 6.1.1 reproduction of the published OMol25 recipe** is feasible and would narrow the Hamiltonian mismatch behind the MACE–PBE0 discrepancy. It must not be described as an exact reproduction of the original ORCA 6.0.0 backend or of a particular training label. The installed executable, its internal auxiliary/ECP libraries, defaults and converged electronic solution have not been shown equivalent to the historical backend.

The original PBE0 directional result remains failed under its frozen gate. This document neither changes that result nor authorizes a retry, optimizer, grid sweep or molecular run.

## Primary authority and reproducible settings

The [OMol25 paper, §2.6 and Appendix A](https://arxiv.org/html/2505.08762v1#S2.SS6) specifies ORCA 6.0.0, ωB97M-V/def2-TZVPD, RI-J/COSX, tight SCF and DEFGRID3. Singlet lanthanide systems use UKS with a 20° beta HOMO–LUMO rotation. Its Appendix A describes dataset filtering of force outliers, spin contamination, density-integral electron counts, COSX errors and negative gaps. Those filters are not proof of a stable electronic state or a numerical derivative for this peptide.

An [OMol25 author explicitly directs reproduction to FAIRchem's calculator source](https://github.com/facebookresearch/fairchem/issues/1934#issuecomment-4155539511). The following are literal settings in [FAIRchem `calc.py` at commit `d4e10d0f4e87593922fda1d1ffb7dd28e504cadd`](https://github.com/facebookresearch/fairchem/blob/d4e10d0f4e87593922fda1d1ffb7dd28e504cadd/src/fairchem/data/omol/orca/calc.py#L123), not reconstructed from a method nickname:

| Item | Source declaration |
| --- | --- |
| Functional and orbital basis | `wB97M-V`, `def2-TZVPD` |
| Requested property | `EnGrad` |
| Coulomb/exchange approximation | `RIJCOSX`, `def2/J` |
| Symmetry and solver | `NoUseSym`, `DIIS`, `NOSOSCF`, `NormalConv` |
| SCF thresholds/limit | `%scf Convergence Tight maxiter 300 end`; `%scf THRESH 1e-12 TCUT 1e-13 end` |
| Integration | `DEFGRID3`; no explicit grid replacement in the generator |
| External orbital basis | `%basis GTOName "def2-tzvpd.bas" end` |
| Spin for this singlet metal case | Explicit `UKS`, plus `get_symm_break_block` |
| Other output | `ALLPOP`; dipole/quadrupole; reduced Lowdin/Mulliken orbital populations and bond orders, Fock matrix and orbital energies |

The [single-point wrapper](https://github.com/facebookresearch/fairchem/blob/d4e10d0f4e87593922fda1d1ffb7dd28e504cadd/src/fairchem/data/omol/orca/recipes.py#L26) defaults to 12 processes and `nbo=False` (`NONBO NONPA`). The standalone input writer instead defaults to NBO enabled. This is a real interface difference; energy/gradient reproduction should declare which wrapper convention is followed. Neither wrapper requires NBO to define the electronic method. The caller can override settings, so these public defaults alone do not identify every historical dataset job.

### La basis, ECP and state

The [bundled basis file](https://github.com/facebookresearch/fairchem/blob/d4e10d0f4e87593922fda1d1ffb7dd28e504cadd/src/fairchem/data/omol/orca/basis/def2-tzvpd.bas) contains the actual H/C/N/O/La orbital coefficients. Its La section is at lines2467–2503, with the additional diffuse P exponent `0.02057322150`; the calculator's basis inventory specifies 43 La functions. Its La ECP section at lines4732–4747 explicitly reads `LA-ECP GEN 46 3`. The calculator independently assigns La an ECP size of 46. These are concrete assets; the old comparator's custom La orbital basis must not be substituted.

The historical [`om-data` basis addition](https://github.com/Open-Catalyst-Project/om-data/commit/c5bfc3a1a7e560983ff7ba67923bec6ff323d9f5) predates the paper and is specifically labeled as an updated lanthanide basis. The current FAIRchem and historical `om-data` basis files have different byte hashes but **identical whitespace-separated tokens**, checked as text only. Preserve one complete file by hash rather than retyping coefficients. The generator does not explicitly set `ECPName`, `NewECP` or an external auxiliary file: its global def2 basis choice supplies the library ECP, while `GTOName` supplies the orbital basis and `def2/J` selects the built-in auxiliary. Actual printed ECP coefficients/core count and orbital/auxiliary basis assignments therefore still need native verification. [ORCA basis/file/ECP documentation](https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/basisset.html#reading-orbital-and-auxiliary-basis-sets-from-a-file).

For this unchanged fixture, the admitted metadata already records charge−2, multiplicity1, 1030 explicit electrons and 515 electrons in each spin channel. Applying the published indexing rule to those saved counts implies **beta rotation `{514, 515, 20, 1, 1}`**; this is a source-level inference, not a generated or executed deck. ORCA counts orbitals from zero and operator1 is beta. UKS permits spin symmetry to break; the rotation does not require a broken-symmetry final solution, identify spin localization or prove a pure singlet. [Rotation semantics](https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/initguess.html#changing-the-order-of-initial-guess-mos-and-breaking-the-initial-guess-symmetry).

The writer does **not** infer metal treatment from atomic numbers: the caller must select `Vertical.MetalOrganics` (or the input-writer's `Oss`) for a singlet. Default vertical would omit explicit UKS/rotation. The generator leaves the underlying guess unspecified. ORCA6.0 documents `PAtom` as default and automatic GBW restart if a matching file exists. Thus inheriting the PBE0 `Guess PModel`, or copying its GBW, would change this recipe. [Initial-guess documentation](https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/initguess.html#huckel-and-patom-guesses).

### VV10 and the two integration systems

`wB97M-V` automatically invokes **post-SCF, non-self-consistent VV10** in ORCA6.0. The generator contains neither `SCNL` nor D3/D4. Adding D4 or switching to self-consistent VV10 changes the target. [ORCA6.0 B97-V implementation](https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/model.html#the-b97-v-family). ORCA6.1 documents the same default, with `b=6.0`, `C=0.01`, unit NL scaling and analytic gradients available. [ORCA6.1 VV10 documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/dispersioncorrections.html#basic-usage).

The documented ORCA6.0 DEFGRID3 SCF grids are XC angular6 and COSX angular2,3,4. The corresponding published radial IntAcc table gives XC4.959 and COSX4.020,4.338,4.871. These rounded table values should not be substituted for internal defaults. The separate VV10 defaults are `VdwAngularGrid 2`, `VdWIntAcc 5.0`, adaptive pruning and `VdwDistTCut 10` Å. The generator overrides none of them. Preserve DEFGRID3 and capture actual XC/COSX/NL headers; equal keywords across releases are not proof of byte-identical grids. [ORCA6.0 integration documentation, Tables7.2–7.3 and nonlocal grids](https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/numint.html).

## Remaining gaps and bounded practical next step

Only ORCA6.1.1 was present in the inspected `/groups/banfield/users/jwestrob/bin/ORCA/` installation. There is no pinned local ORCA6.0.0 executable or archived OMol25 La job in this evidence. The [official raw-output release](https://github.com/facebookresearch/fairchem/blob/main/docs/molecules/datasets/omol25_elec.md) offers original input/output/gradient files through Argonne, initially for the4M subset; access needs group approval. None were downloaded here. Current generator source is a maintainer-endorsed reproduction reference, not a per-training-example build or wrapper receipt.

Recommendation: separately declare **one fixed central energy/gradient scout on6.1.1 using the published recipe**, keeping exact step725 physical IDs/coordinates, charge and multiplicity. Preserve the recipe's external basis, UKS rotation, Tight/300/DIIS/NOSOSCF, DEFGRID3 and post-SCF VV10. Admit actual method/basis/ECP/UKS/counts, convergence and finite gradient; retain printed residuals, density integrals, gap and S². Do not reuse the old RHF/PBE0/D4-specific parser as its scientific admission. The generic pinned runner/resource renderer can be reused with a separately reviewed method parser and the existing12MPI allocation policy. Native runtime and cost remain unmeasured for the diffuse, range-separated UKS calculation.

After actual central admission and cost review, the same two already-frozen displaced singlepoints could test this recipe's derivative, with a predeclared unchanged rule and no automatic retry. Until then its forces are descriptive. A smaller MACE discrepancy would support Hamiltonian mismatch as part of the old difference; a persistent discrepancy would concern transfer to this **published-method port**, with backend-version and state uncertainty retained. Neither result uniquely isolates functional versus basis versus spin treatment. It does not validate the physical ground state, shared ALPB, an equilibrium structure or metal preference. No numerical result from this new method exists yet.

## Exact evidence pins

SHA256 values below were read or computed from source/metadata bytes only. Remote files were read in memory; this document does not stage an executable input or new basis asset.

| Authority | SHA256 |
| --- | --- |
| FAIRchem commit `d4e10d0f4e87593922fda1d1ffb7dd28e504cadd`, `src/fairchem/data/omol/orca/calc.py` | `6cbf04e71579e9021ba31628f04463a59b28bc063834dc395581f579133c9843` |
| Same commit, `orca/recipes.py` | `15d2e4d2c3e204998f29a07433aa56ad942d273d88e95998b4e1258d18a9d842` |
| Same commit, `orca/basis/def2-tzvpd.bas` | `7389c56f0cff9e76098b0757a34d1abdd79c48e9249597cea074a0106a3462a6` |
| `om-data` commit `32f10cd9dfe6740ef9d4ff1ed9d2956e45081b07`, `omdata/orca/calc.py` | `0210b962a66bc3acb5cfe242ddd32ab3db8cdd58b7e4cff323d4d7253b82067b` |
| Same historical repository commit, `omdata/orca/basis/def2-tzvpd.bas` | `cb5719b24b233553456390c97c443ada5d37f2cb0e48154fffc50688542844ca` |
| `workspaces/mace_omol_20260917/software_v1/release.json` | `1855c8105133f34f98cb387c980c42abfacebb07a59505d8d09b0edaaf7d6fa5` |
| `workspaces/mace_omol_20260917/software_v1/software_manifest.json` | `7d2448b279b5191a529ef329a60442635601576794f9737efc64343c9a6e0168` |
| `workspaces/mace_hybrid_20260916/software_v1/venv/lib/python3.11/site-packages/mace_torch-0.3.16.dist-info/METADATA` | `ce3a831df793896eb140235a93dc0281e161db56bb89a0790b9d1cae3154935d` |
| `workspaces/lanm_coupled_mechanics_20261017/lbt_7ccn_converged_reference_v1/run/execution/RESULT.json` | `a53125b98e597ea0473e06d58316ac2b32c914e9647da595f163cf0aa76a8c19` |
| Same execution, `new_La/endpoint.inp` (old PBE0 comparator) | `b94c22fc46ff47bd8c678f2e49634948fa99195638f6f69823d0b99f075b102d` |
| `/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca` | `38b5f057452fef275c0a1b98d270ad03d0411dab70453820d1c678bd732f6c83` |

The saved release/manifest bind the official100M MACE-OMOL checkpoint `MACE-omol-0-extra-large-1024.model`, SHA256 `9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`; it was not loaded or rehashed here. Installed metadata identifies its DFT target as `wB97M-VV10`. That shorthand is resolved by the paper and generator above, not interpreted as an additional dispersion correction.
