# Native VV10 method discovery: job 16602

21 October 2026 UTC. Read-only native-header, source and primary-documentation inspection; no scientific imports, calculations, job changes or edits to frozen inputs.

## Correction and decision

**ORCA 6.1.1 has automatically selected self-consistent VV10 (SCNL) for the published `wB97M-V EnGrad` input.** The earlier claim that this EnGrad calculation would retain post-SCF VV10 was incorrect. This note supersedes that interpretation in `TARGET_RECIPE_EVIDENCE.md`, the declaration's anticipated method description and the adapter's anticipated post-SCF evidence label; their original bytes remain preserved. The actual method must govern subsequent interpretation.

Continuing the unchanged single central calculation is scientifically useful as a **published-input ORCA 6.1.1 port with native SCNL treatment**. This is an observed implementation behavior, not an operator retuning of the input. It is not yet an admitted final result, an exact ORCA 6.0.0 reproduction or proof of the original training-label treatment. Root owns the continuation decision. The existing pending-method-review status remains necessary.

## Actual native evidence

Source: `workspaces/lanm_coupled_mechanics_20261017/lbt_omol_target_reference_v1/run/execution/central_La/endpoint.out`, lines 218–227, observed at `2026-10-21T03:16:33.670980+00:00`:

```text
WARNING: B97M-V or wB97M-V functional requested.
  ===> : libXC variant has now been activated.
         Any manually given DFX/HFX/RangeSep scaling parameters are IGNORED!


WARNING: Analytic gradients for DFT-NL dispersion can only be done using the self-consistent method,
  ===>   Switching to SCNL

WARNING: Self-consistent DFT-NL dispersion not compatible with AutoTRAH.
  ===>   Switching AutoTRAH off
```

The input echo retains `wB97M-V def2-TZVPD EnGrad`, RIJCOSX, DEFGRID3, UKS and the declared beta rotation. No explicit SCNL keyword was added. The native header therefore establishes three actual behaviors: libXC activation, self-consistent NL selection and AutoTRAH deactivation. It does not establish completed convergence, electronic stability, successful final symmetry breaking or gradient accuracy.

The output was still changing. The pins below identify the exact observed prefix and excerpt, **not a terminal hash of the whole output**. Hashes include original newline bytes.

| Evidence | SHA256 |
| --- | --- |
| `endpoint.out`, first 227 lines, 14,749 bytes | `4df2b57a925796d9430aa212b54838ad269938343d971bda53023e5679f6fd3f` |
| `endpoint.out`, lines 218–227 only | `69fc16f94b62820c8c30f4d0127e52aeabaf53f82684e8f393f78ddca8326622` |
| `central_La/endpoint.inp` | `d482826aa2384c0e00f8c8fb0cfa5dc965e22a3c2a99f17a4e608602a45200c6` |
| `central_La/endpoint.runtime.inp` | `d028084243c29e768fa5ab609614de7eb4ad3504ba0552800f51e55a93becd23` |
| `central_La/endpoint.runtime.json` | `f7daac838334126e20573a6df5df16bd680785ae6ad30a40cf424122fc14c134` |
| `execution/MANIFEST.json` | `1b543ed7c470cd9a8df7b0a7cfa3329cf6c8f30be862f896ec15049c9a44d43e` |
| `run/ROOT_NATIVE_METHOD_OBSERVATION.json` | `16e1c029e1ec17c6d28a761da162a22e5f1634dc7a5c59ffb4760f1aae2ca686` |

Paths in the table are relative to the same run/execution or run directory as indicated. Runtime evidence records job 16602, 24 MPI ranks and the unchanged allocation-derived memory policy.

## What the official documentation establishes

The [ORCA 6.0 manual, §7.4.2.9.1](https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/model.html#the-b97-v-family), describes post-SCF NL as the B97-V default. However, its [§7.4.2.9.3, self-consistent DFT-NL](https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/model.html#self-consistent-computations-with-the-dft-nl-dispersion-correction), specifically places analytic-gradient availability within the self-consistent subsection. That supports interpreting the published EnGrad request as potentially invoking SCNL; it does **not** explicitly document the automatic switch. No original ORCA 6.0.0 output or executable was inspected, so **automatic SCNL for the original 6.0.0 EnGrad workflow remains unverified**.

The [ORCA 6.1 manual, §3.4.2.1](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/dispersioncorrections.html#basic-usage), likewise distinguishes default non-self-consistent NL from explicit SCNL and lists analytic-gradient availability. Its general default description alone does not determine the effective gradient calculation. For this run, the native warning supplies the decisive evidence. SCNL optimizes orbitals and density with the full DFT plus VV10 potential, whereas post-SCF NL evaluates the correction on the preceding DFT density.

The [maintainer-endorsed FAIRchem calculator](https://github.com/facebookresearch/fairchem/blob/d4e10d0f4e87593922fda1d1ffb7dd28e504cadd/src/fairchem/data/omol/orca/calc.py#L123) requests `EnGrad` and contains no explicit SCNL selection. Thus matching its input keywords does not prove a post-SCF target. The [OMol25 paper](https://arxiv.org/html/2505.08762v1#S2.SS6) identifies ORCA 6.0.0; an original raw EnGrad output or explicit author clarification of the effective NL treatment would resolve the remaining historical ambiguity. No author message or additional calculation was sent or run.

## Consequences for interpretation and any later derivative test

The final receipt must be reviewed as SCNL/libXC with AutoTRAH disabled, retaining actual NL energy/gradient evidence and all native state/convergence limitations. Absence of a post-SCF correction heading is not evidence that VV10 was omitted. No parser-driven rerun is warranted.

`TARGET.py` uses Cartesian-component RMS, `sqrt(sum(F_ij**2)/(3*N))`, for both force comparisons. Some earlier summaries use RMS atom-vector norm instead. Compare the same definitions or the full saved vectors; a change of RMS convention is not a reduction of force error.

Simply removing `EnGrad` from a future single-point deck risks reverting to the documented post-SCF default and comparing different energy definitions. **No single-point calculation is automatically released.** Any later derivative test needs a separate declaration and native verification of matching self-consistent NL treatment at both displaced points, alongside the already frozen direction and numerical gate. This note authorizes no such execution.

## Frozen provenance retained

| Existing authority | SHA256 |
| --- | --- |
| `TARGET_RECIPE_EVIDENCE.md` (superseded interpretation preserved) | `cb20795b46d2d413b1d494cd6462d68ca0754638942a8e4993267ef39d8b30bb` |
| `TARGET.py` | `aa62db4f37de3c50765ad81115dd75cb4aef2a99d92b4465d9bbe59df9a4a74f` |
| `SPEC.md` | `0f11bb6a187391a5284fd8962fd9738672401b1c7cc5c370c1418effb1bc0ae4` |
| `CONFIG_RUN.json` | `1c9f22df1c9ac096568dccac7689a5e6c9f8af4d8f8ae2bd0a3d589aa9ca0c83` |
| Public `assets/calc.py` | `6cbf04e71579e9021ba31628f04463a59b28bc063834dc395581f579133c9843` |
| Public `assets/recipes.py` | `15d2e4d2c3e204998f29a07433aa56ad942d273d88e95998b4e1258d18a9d842` |
| Public `assets/def2-tzvpd.bas` | `7389c56f0cff9e76098b0757a34d1abdd79c48e9249597cea074a0106a3462a6` |

Only this new discovery note was written. The old failed PBE0 directional result, current native job, source, settings and inputs remain unchanged.
