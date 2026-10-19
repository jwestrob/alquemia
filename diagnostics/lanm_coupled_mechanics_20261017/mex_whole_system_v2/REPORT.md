# Mex whole-system mechanics: fixed-motion result

19 October 2026 UTC. The complete eight-cell test shows that **changing the
internal-energy model changes which metal benefits more from the same motion**.
It establishes feasibility and model sensitivity; independent Mex accuracy
remains untested. The candidate is not ready for a unique-state optimization or
production use.

## What was held fixed

Both metals use the exact start and terminal coordinates of the previous La
search, job 12970: all 2,049 physical atoms, including 167 waters and three Nd
spectators. The terminal was admitted but did not converge. The new energy is
whole-system native MACE-OMOL vacuum plus the same fixed-charge ALPB. Four saved
ALPB evaluations were reused. No MM, GFN field, extra dispersion, caps, new
geometry or optimization was added.

The declared denominator is two geometries × two metals × two aggregate-spin
inputs. All eight cells completed, with one scout forward and seven further
forwards. No successful cell was repeated. The manifest, source, raw energies,
full force arrays and exact component assembly passed collection and independent
[review](COLLECTION_REVIEW.md).

## Deformation work

After minus before, kcal/mol. Negative means the model lowers the energy along
this particular motion.

| Model / aggregate-spin input | La | Dy | Dy minus La |
| --- | ---: | ---: | ---: |
| Whole MACE + ALPB, La 10 / Dy 15 | −1937.990241 | −1942.176256 | −4.186016 |
| Whole MACE + ALPB, La 4 / Dy 9 | −1944.455151 | −1948.078345 | −3.623194 |
| Previous hybrid comparator | −2945.113073 | −2934.389151 | — |

The old hybrid favors La more along this motion; both whole-system inputs favor
Dy more. The whole-system work remains strongly favorable, with ALPB opposing
the vacuum lowering by +552.925885 for La and +551.835625 for Dy. The previous
hybrid is a comparator, not an independent reference. These are within-model
work differences, not comparisons of absolute energies from different models.

Changing the aggregate-spin input changes La work by −6.464910 and Dy work by
−5.902089 kcal/mol. Both exceed the frozen 1 kcal/mol reporting resolution.
The change in their difference is only +0.562821, but this cancellation does not
establish a physical spin state. The aggregate categories cannot identify local
Nd moments, the reversed ion or exchange coupling. Neither input is selected
after inspecting its energy.

## Interpretation and next decision

This replacement changes scaffold, partition-boundary, long-range and spectator
treatment together. It cannot identify one defective term. The saved force
changes describe after minus before, not errors against DFT or convergence
norms. No affinity, binding free energy or robust La/Dy preference follows.

The finite Mex panel is complete. Its material spin sensitivity blocks treating
one selected conditioning input as a unique prediction. No new spin sweep,
optimizer or reopening of the failed embedded-DFT derivative branch follows.

The independent 7CCN whole-peptide DFT challenge is now running its first
before-La scout, job 14795, on 12 MPI ranks. That test asks whether whole-MACE
deformation work transfers to another fixed peptide pair. Its eventual result
will not establish the mixed-ion Mex electronic state. Root will review native
admission and actual cost before the remaining three declared reference cells.

## Evidence and costs

Final collector 14794 made zero molecular calls. Result SHA-256:
`4c0c8b7379a2210f64267f0df94b1df72b54902d7233ec4b9614d4e941349d27`.
Its original path is
`workspaces/lanm_coupled_mechanics_20261017/mex_whole_collection_v1/execution/RESULT.json`.
The executed source SHA is `dd611ca0…`; manifest `acdaa462…`;
collection config `1536e026…`. Full hashes are in the independent review and
archived execution evidence.

The first H200 forward took 1.076935 s after 6.613859 s model loading, with
81,684,463,616 bytes peak allocated GPU memory on a 150,109,880,320-byte device.
Both frozen continuation limits passed. These are measured native/application
costs; allocation costs are recorded separately from actual scheduler receipts.
The five Mex jobs used **345 allocated CPU-seconds and 28 allocated GPU-seconds**,
including the failed preflight. The independent LBT saved preflight adds
3 CPU-seconds; the running DFT scout is excluded from these completed totals.
Exact original-to-copy hashes, full arrays, configs and raw accounting are in
[the execution archive](execution_evidence_v1/ARCHIVE.json), with
[scheduler costs](execution_evidence_v1/ACCOUNTING.json).
The original zero-forward preflight 14788 failed on a source-name assertion;
v2 corrected only that literal. The failed attempt remains preserved.

The broader mechanics goal remains active and unmet. Production is unchanged.
