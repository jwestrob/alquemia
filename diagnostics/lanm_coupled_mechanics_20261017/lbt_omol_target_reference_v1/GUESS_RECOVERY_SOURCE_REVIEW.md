# Initial-guess recovery source review

21 October 2026 UTC. Source, native text and metadata inspection only; no scientific execution.

## Finding and recommendation

Recommend one separately declared **PAtom → PModel** initial-guess recovery. Preserve the failed job and every molecular setting, including the external def2-TZVPD basis/ECP46, wB97M-V with actual SCNL, UKS charge−2/singlet, beta rotation `{514,515,20,1,1}`, grids, DIIS/NOSOSCF, thresholds and 300-iteration limit. This is a derived-recipe variant with a changed initialization, not an exact historical-default reproduction. Convergence and the selected SCF solution remain uncertain.

Job 16602 terminated during guess construction, before molecular SCF: `endpoint.out:2154` reports `SCF orbitals not available for element 57(La)` from `orca_guess/guess_eht.cpp:712`; line2225 identifies polarized atomic densities, and line2645 identifies GUESS error termination. This is evidence of unavailable PAtom atomic-orbital data on this path, not failed SCF convergence or exhausted memory.

## Primary documentation and known La coverage

The [ORCA 6.1 initial-guess manual](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html#pmodel-guess) describes PModel as diagonalization of a Kohn–Sham matrix built from superposed neutral-atom densities, available for HF and DFT and particularly suitable for heavy elements. [HCore](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html#one-electron-matrix-guess) diagonalizes the one-electron matrix but commonly gives excessively compact orbitals. [Hueckel and PAtom](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html#huckel-and-patom-guesses) use minimal atomic bases; changing to Hueckel alone has no demonstrated La coverage here. These facts favor PModel over an untested HCore/Hueckel fallback. The [rotation interface](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html#changing-the-order-of-initial-guess-mos-and-breaking-the-initial-guess-symmetry) operates on initial MOs; retaining the same indices and angle does not make the new orbitals identical to the old guess.

There is direct local La coverage evidence. Completed job15545, on the same 256-atom endpoint and ORCA6.1.1 executable, prints at `new_La/endpoint.out:2307–2333`:

```text
INITIAL GUESS: MODEL POTENTIAL
Loading Hartree-Fock densities                     ... done
Atomic calculation to generate a density for atom 255 (La with ECP)
  Multiplicity assumed   :    2
...
Returned from the atomic calculation with res=0
...
                      INITIAL GUESS DONE ( 157.3 sec)
```

That run used the prior PBE0 orbital basis and converged normally. It proves this build can construct a PModel La/ECP46 guess, **not** that the diffuse UKS/SCNL recovery will converge. No equivalent success is claimed for HCore or Hueckel.

**Scope/cost detail:** PModel can launch an internal atomic calculation to generate the La density. The old output explicitly invokes `endpoint.runtime_atom57.inp`. A recovery would request one whole-molecule endpoint through one runner invocation, while allowing this native initialization subcalculation; it must not be reported as one total internal solver invocation. Preserve the generated atomic input/output and include its time in total allocated cost. The atomic density's printed multiplicity is an initialization detail, not a change to the molecular singlet declaration.

## Published-recipe distinction and bounded recovery

Pinned [FAIRchem calc.py](https://github.com/facebookresearch/fairchem/blob/d4e10d0f4e87593922fda1d1ffb7dd28e504cadd/src/fairchem/data/omol/orca/calc.py#L123) and recipes.py do not select PAtom, PModel, HCore or Hueckel. Our failed adapter explicitly added PAtom based on the manual's default description. Keep that historical reasoning and failure intact. A PModel-only replacement changes initialization, not the declared final energy functional, but can converge to a different electronic stationary point; no stable-state or ground-state claim follows.

Use fresh native files, NoAutostart and no prior GBW. Admit the actual model-potential guess, unchanged beta rotation, method/state/basis and final native evidence. Retain pending method review, including actual SCNL/gradient treatment. One declared attempt only: no automatic alternate guess, SCF extension, grid change, SP pair or restart. Root subsequently authorized preparation of a separate variant; no calculation was performed for this review.

## Exact local pins

Paths start under `workspaces/lanm_coupled_mechanics_20261017/` unless marked diagnostic.

| Evidence | SHA256 |
| --- | --- |
| `lbt_omol_target_reference_v1/run/execution/RESULT.json` | `881852372e543e2f2c721e1b162997cb348c2bc3ddbbbb7d6e910560f7b5c55f` |
| Same execution, `central_La/endpoint.out` | `770979dbf59e185ea50ead6687278985f3bc0c4333e0c81604a0ed4756f990c0` |
| Same execution, `central_La/endpoint.inp` | `d482826aa2384c0e00f8c8fb0cfa5dc965e22a3c2a99f17a4e608602a45200c6` |
| `lbt_7ccn_converged_reference_v1/run/execution/RESULT.json` | `a53125b98e597ea0473e06d58316ac2b32c914e9647da595f163cf0aa76a8c19` |
| Same execution, `new_La/endpoint.out` | `6cddbc036f4a10ff99ff27e42ef3c289f4ba2f252ce42c36102d6f55495fbf29` |
| Diagnostic `lbt_omol_target_reference_v1/assets/calc.py` | `6cbf04e71579e9021ba31628f04463a59b28bc063834dc395581f579133c9843` |
| Diagnostic `lbt_omol_target_reference_v1/assets/recipes.py` | `15d2e4d2c3e204998f29a07433aa56ad942d273d88e95998b4e1258d18a9d842` |

No existing file, input, job or result was modified for this review.
