# Same-deformation whole-peptide model comparison

17 October 2026 UTC. Complete: jobs 13140, 13141 and 13142. Root reviewed actual
native results and saved replay; independent integration review agrees.

## Result and decision

**The exact same La deformation lowers the hybrid energy by 28.244859 kcal/mol
but raises whole-peptide MACE plus fixed-charge ALPB by 11.290076 kcal/mol.**
This is a disagreement in the cost of motion under two Hamiltonians. It does
not establish which model is accurate, nor isolate the discrepancy uniquely to
MM stiffness: the electronic region and its treatment also differ.

Both metals use the same before coordinates and the same after coordinates.
The latter are the actual 80-step, nonconverged La endpoint from job 13017.
Work means `E(after) - E(before)`; positive values mean an energy cost.

| Whole-peptide model | MACE vacuum | Fixed-FF ALPB | Total work |
| --- | ---: | ---: | ---: |
| La | +5.908035 | +5.382042 | +11.290076 |
| Dy | +5.134286 | +5.948114 | +11.082401 |
| Dy minus La work | -0.773749 | +0.566073 | -0.207676 |

All table values are kcal/mol. The Dy-minus-La work difference is below the
predeclared 1 kcal/mol reporting resolution. That resolution is neither an
uncertainty estimate nor proof that the difference is zero. These deformation
works are not metal-binding energies or a preference prediction.

The previously collected hybrid La work on this exact displacement is:

| Hybrid component | Work, kcal/mol |
| --- | ---: |
| Local MACE vacuum | -10.371217 |
| GFN2 embedding difference | +0.025876 |
| MM interaction ledger | -23.281560 |
| Same fixed-FF ALPB | +5.382042 |
| Total | -28.244859 |

The solvent contribution is shared and therefore does not explain this
between-model La disagreement. No absolute energies are subtracted between
Hamiltonians. Hybrid after-Dy remains unavailable; no hybrid Dy work is inferred.

The next accuracy evidence comes from the already-running independent DFT
scouts. Whole-peptide job 13051 evaluates **before-La only**, so its immediate
test is the starting force and computational cost. It cannot by itself validate
the finite before-to-after work. Further reference cells remain subject to the
existing admission and cost review. No parameter adjustment, renewed failed
hybrid optimization, new geometry search or production change follows here.

## Forces and coverage

All four whole-model energies are available. Three whole-MACE force sets are
available; the archived before-Dy request was energy-only and remains so.
The full available force arrays and component projections are retained in
[the collected result](result/COLLECTION_RESULT.json).

| Cell | Total force along frozen deformation chord | Total force along initial hybrid Phe1 CA direction |
| --- | ---: | ---: |
| Before La | +2.466204 | +0.798547 |
| After La | -1.076003 | -44.464198 |
| Before Dy | unavailable | unavailable |
| After Dy | -1.120943 | -45.113217 |

These are force projections in kcal/mol/Å on the two directions frozen by
preflight 13049. Negative means force opposing the named direction. They are
endpoint observations, not integrated work or evidence of a stationary state.

## Execution and qualification

The [specification](SPEC.md) was frozen before execution. The finite ceiling was
two whole-MACE requests and one fixed-FF ALPB request. Actual calls match it.
There were zero DFT, GFN, MM, optimization or new-geometry requests in this round.

- Preflight 13140 passed source/hash joins, coordinates and atom order, charge
  and spin inputs, old native energies, available forces, three old ALPB raw
  outputs and hybrid component closure. It made zero molecular calls.
- Run 13141 completed the two after-state MACE energies/full forces and the
  missing after-Dy ALPB energy/gradient. All native admission checks passed.
- Collector 13142 independently replayed the saved native outputs, physical
  inputs, requests, force arrays and comparison with zero molecular calls.
- No executed attempt failed in this finite adjunct. Earlier failures elsewhere
  remain preserved; the closed hybrid SCC recovery was not reopened.

The 256 atoms, target index 255, charge -2, free NH3+/COO- termini and five
deprotonated acidic groups are unchanged. MACE uses physical La/Dy
multiplicities 1/6. The effective singlet convention of the independent frozen-f
DFT reference is distinct. Terminal chemistry remains conditional; the endpoint
is nonstationary and the electronic terms are vacuum models with an added
fixed-charge continuum contribution. This is not an equilibrium or affinity test.

MACE uses the original OMOL checkpoint `9b64b4fd…`, `omol` head, float64 and
unchanged native helper. Versions are torch 2.8.0, mace-torch 0.3.16 and ASE
3.29.0. The native fixed-FF ALPB binary and charge vector are pinned in the
actual configurations and collected result; no native GFN default is used.

| Job | Work | Allocation | Elapsed | Allocated CPU-seconds | Allocated GPU-seconds |
| --- | --- | --- | ---: | ---: | ---: |
| 13140 | Saved preflight | 1 CPU | 5 s | 5 | 0 |
| 13141 | Two MACE and one ALPB calls | 12 CPUs, 1 H200 | 11 s | 132 | 11 |
| 13142 | Saved collector | 1 CPU | 3 s | 3 | 0 |
| Total | | | | **140** | **11** |

All numerical work ran through Slurm. No memory, time or priority overrides
were requested. Scheduler accounting is retained in [SACCT.txt](result/SACCT.txt);
application times are recorded separately in the native results.

## Provenance and restart

| Artifact | SHA-256 |
| --- | --- |
| SPEC.md | `e46f387c4646e4406b273a5752faadaa37279f9587bf0ef83348dbc45f2f40ee` |
| MODEL_COMPLETION.py | `f35e2225473b1df10288c532a930bdd49186eaaf2e3bab40c65026c44024d8bc` |
| Actual preflight result | `66bd7fa0cf4aa8789c68c1512d8fb8195e60d9dca594479a9f9af22f503e162c` |
| Actual run result | `12c85eb939cbde298c896a80701c91b3eab8b7fb0b1beeb35d1acacdfce5088f` |
| Actual collection result | `64f27e852d2403600ef7a73616c7e52cf50cfbace75689a5b932fde926a06ef0` |

Execution workspace:
`workspaces/lanm_coupled_mechanics_20261017/lbt_whole_model_completion_v1/`.
The exact executed commands, Python environment, configurations and watcher
records are preserved in the three submission JSON files under `result/`.
The adapter exposes `--stage preflight|run|collect --config FILE --output FRESH`;
use [PREFLIGHT.sbatch](PREFLIGHT.sbatch) or [RUN.sbatch](RUN.sbatch) with those
recorded arguments on Slurm. Completed chemistry must not be resubmitted from
a delayed notification. Root's [review](result/ROOT_REVIEW.json) records proactive
inspection; an actual delivered callback is acknowledged separately.
