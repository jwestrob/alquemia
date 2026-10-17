# Whole peptide reference challenges the hybrid deformation benefit

17 October 2026 UTC. All four fixed before/after La/Dy reference cells are
complete and admitted under the original native collection rules. **Whole-peptide
PBE0 and whole-peptide MACE assign an energy cost to the tested motion, while the
hybrid model assigns a large La benefit.** This challenges the hybrid mechanics for
this specific deformation. The small La/Dy difference does not establish a metal
preference.

## Same coordinates and solvent term

The comparison uses the exact 256-atom 7CCO starting geometry and the saved
La-proposed endpoint from the baseline hybrid response. Both metals are evaluated
on both geometries. The endpoint is the original nonconverged 80-step result;
neither geometry is a DFT minimum. No refolding or geometry repair enters this
test.

Whole vacuum PBE0 replaces all internal electronic and classical terms. Each
model then receives the same saved fixed-force-field-charge ALPB solvent term.
Define deformation work as `W_M = E_M(after) - E_M(before)`: positive work means
that the endpoint costs energy relative to the start. Define the metal contrast
as `W_Dy - W_La`. These quantities compare one deformation, not binding affinities.

## Deformation work

All entries are kcal/mol. Values come from the unchanged full four-cell collector
and the pinned earlier whole-model completion.

| Model or contribution | La | Dy |
| --- | ---: | ---: |
| Whole PBE0 vacuum | +6.897552 | +5.774542 |
| Shared fixed-charge ALPB | +5.382042 | +5.948114 |
| **Whole PBE0 plus ALPB** | **+12.279593** | **+11.722656** |
| Whole MACE vacuum | +5.908035 | +5.134286 |
| **Whole MACE plus ALPB** | **+11.290076** | **+11.082401** |
| **Hybrid plus ALPB** | **−28.244859** | unavailable |

The reference metal contrast is −0.556937 kcal/mol, compared with −0.207676 for
whole MACE. Both are below the original 1 kcal/mol reporting resolution. That
resolution is a reporting convention, not an uncertainty estimate. The hybrid
Dy work remains unavailable from the closed prior response branch.

For La, the hybrid work comprises −10.371217 from its MACE region, +0.025876 from
the GFN2 field difference, −23.281560 from its classical interaction ledger, and
+5.382042 from the shared solvent term. This decomposition identifies large
terms in the hybrid benefit; it does not uniquely identify a faulty parameter
or an incorrect interaction count.

## Forces support the finite deformation discrepancy

The projections below use the direction frozen before the DFT calculation: the
initial hybrid force on terminal Phe1 CA. Units are kcal/mol/Å.

| Cell | Whole PBE0 plus ALPB | Whole MACE plus ALPB | Hybrid |
| --- | ---: | ---: | ---: |
| Before La | +3.341027 | +0.798547 | +51.063006 |
| After La | −40.851355 | −44.464198 | −1.895426 |
| Before Dy | +3.345638 | unavailable | +51.064870 |
| After Dy | −40.850094 | −45.113217 | unavailable |

The hybrid initially drives a much larger terminal force and leaves a much
smaller opposing force at its endpoint. Both whole-peptide methods strongly
oppose further motion along that terminal direction at the endpoint. Initial
reference projections are positive, so these results do not imply that every
infinitesimal displacement along the proposed motion is unfavorable.

For La, the global force difference RMS relative to PBE0 plus ALPB is 10.141916
before and 8.850041 after for the hybrid, versus 3.511987 and 3.516438 for whole
MACE. MACE therefore agrees more closely on this fixture but retains substantial
force differences. This is not a general validation of MACE. The full 256-atom
forces, normalized deformation-chord projections and component records are in
the [archived result](full_result_v1/RESULT.json) and its four force reports.

## Admission and interpretation limits

All four native calculations terminated normally and supplied full gradients,
with the original native SCF-convergence and gradient-admission rules preserved. Printed SCF residuals
remain nonuniform: before-La passes orbital-gradient and rotation criteria;
the other three pass only the orbital-gradient criterion. Every cell retains
`all_printed_residuals_pass=false`. Early SCF small/negative-gap warnings are
retained; final orbital tables have positive gaps. Native
admission establishes neither electronic-state stability nor numerical accuracy
of the analytic gradient; both qualification flags remain false.

The calculation is conditional on free NH3+/COO− termini and five deprotonated
acidic side chains. Experimental terminal chemistry is unresolved. The vacuum
dianion has charge −2; adding fixed-charge ALPB afterward does not self-consistently
polarize its density. La uses ECP46 and Dy uses ECP55 with frozen 4f electrons,
effective singlets and no spin-orbit coupling. Physical multiplicities remain
La 1 and Dy 6. Both effective electron counts are 1030. These assumptions and the
nonstationary endpoints limit the scope of the reference.

## Decision and next check

Do not promote or continue this challenged LBT hybrid response. The present
evidence points to an inadequate force balance and deformation cost in the
hybrid on this fixture. It does not distinguish a displaced preferred local
structure from stiffness or partition effects.

The next finite check compares already-saved baseline and expanded electronic
regions at the exact starting coordinates, for both metals, against these
whole-peptide reference forces. It preserves the existing projections and
component definitions. No new molecular calls are needed. The earlier classical
term audit already closed all 52 Phe1 CA terms without an inventory defect;
repeating that count would not resolve the discrepancy. Expanded-region
after-state evidence remains unavailable.

Mex central consistent-grid calculation 13281 remains separate and active at
this checkpoint. Its actual native admission and cost still gate at most two
displaced calls. The broader mechanics goal remains unmet, including independent
physics qualification and gated natural-LanM transfer. Production is unchanged.

## Execution and preservation

| Job | Work | Elapsed seconds | Allocated CPUs | Allocated CPU-seconds |
| --- | --- | ---: | ---: | ---: |
| 13051 | Before La | 23741 | 12 | 284892 |
| 13290 | After La | 24919 | 12 | 299028 |
| 13291 | Before Dy | 23058 | 12 | 276696 |
| 13292 | After Dy | 24449 | 12 | 293388 |
| 13424 | Full saved collection | 3 | 1 | 3 |

The four molecular allocations total 1,154,004 CPU-seconds. Including original
preflight 13049, scout collection 13283, continuation preflight 13288, controls
13289 and full collection 13424 gives **1,154,051 allocated CPU-seconds**. This
excludes the earlier whole-model completion and other mechanics campaigns.
Allocation accounting is not measured utilization. Exact scheduler records and
scope are in [COSTS.json](full_result_v1/COSTS.json).

The full collector ran with zero new DFT, MACE, GFN2, MM or ALPB calls. Its
`RESULT.json` SHA256 is
`9229f7bebb533dd44a024bb80eca10eded2b3a64cc3972b245d65d82b5ecda97`.
The archive includes exact config, declaration, submission, native result pins,
four force reports, scheduler costs, root review and independent review.
The collector and after-La admission passed independent review. Original
failures and missing comparators are retained. Proactive review and delivered
notification acknowledgements remain separate records.
