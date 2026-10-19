# Independent scout result review

19 October 2026. **Before-La scout 14795 satisfies the unchanged native
admission. No blocker to the declared remaining-three continuation was found,
subject to root's exact scout and measured-cost decision bindings.** Review was
limited to source, saved fields, raw output, identities and hashes; no numerical
replay, scientific imports or new molecular calls were performed.

## Admission and execution

Actual RESULT is `complete`; only `before_La` is `available`. It records one
DFT attempt and one native runner invocation, with zero new MACE, ALPB, GFN,
MM, geometry or optimization calls. The other three cells remain
`not_requested`; `full_four_available` is false and the comparison unavailable.
The exact scientific binding matches admitted zero-call preflight 14793.
Task, singleton manifest and execution receipt identities agree. Thirteen
direct source/result/native artifact pins were independently hash-verified.

ORCA 6.1.1 terminates normally with return code zero and a complete engrad
receipt. The saved gradient has all 256 three-component rows, admitted by the
unchanged finite-gradient parser. Raw output shows the specified PBE0-D4,
def2-TZVP/def2-J RIJCOSX, DefGrid3, VeryTightSCF EnGrad recipe, charge −2,
effective singlet, 1030 explicit electrons, La ECP46 at target 255 and MaxIter
125. ECP, COSX, XC and dispersion gradient stages finish. The atomic-guess
multiplicity 2 is distinct from the printed molecular multiplicity 1.

Twelve MPI processes and one thread per rank are bound in the receipt.
Allocation-derived memory is 48,000 MiB with 25% reserve and MaxCore 3100 MB
per rank. The terminal scheduler record is COMPLETED, exit 0:0, 23,440 seconds
on 12 CPUs. That duration lies within the SPEC's recorded 7CCO endpoint cost
range; it supports an unchanged finite continuation after the required factual
cost receipt and decision are pinned. Application wall time is not utilization.

## Convergence and scientific limits

The native route is explicitly **“Energy Check signals convergence” after
15 cycles**, with `ConvCheckMode = Total+1el-Energy` and `ConvForced = 0`.
Only the printed energy change meets its tolerance. MAX density, RMS density,
DIIS error, orbital gradient and orbital rotation flags remain false:
`all_printed_residuals_pass = false`. A small HOMO/LUMO-gap warning occurs during
SCF. These observations are preserved under the original admission, without
promoting native convergence to all-criteria convergence or physical stability.

Final gridx3 exchange integration and nuclear gradients finish. Its printed
exchange-energy correction is not a force-error estimate; no independent
whole-system finite-difference qualification follows. Both
`force_consistency_qualified` and `physical_state_qualified` remain false.
The vacuum dianion, additive fixed-charge solvent, conditional protonation and
physical-versus-effective spin assumptions remain. This is one endpoint of a
fixed 7CCN deformation of the same peptide, not a completed work comparison or
independent-protein accuracy validation. Both Dy MACE force comparators remain
unavailable by declaration. No settings change or scout repetition is indicated.

## Exact pins

| Artifact | SHA256 |
| --- | --- |
| Scout RESULT | `5fe74c9bcd8ca90a616117387f8960a2fdcacfd62cee92f4cfe30ddd424b9f39` |
| REFERENCE.py | `c0cca6dd4f44619fc2e3201f024b7913f6b7e9d8e8c4fca694fe17bd942d0ecf` |
| CONFIG | `8f0dffb2c05425cf1ff4ea7149cebf08ee42b452f452c3cd4bd096118a59fb2c` |
| SPEC | `a8cb688c8d5f8f1f109a6c1436da24e68410cb92528f503d123cc83650307447` |
| Preflight RESULT | `64db262eadb942246cff51a8bccca5b5b52520bf802cf93ff6f241610691be39` |
| Singleton MANIFEST | `ba4e8059b5907de9dfeab2b1e3109dc9eaaa805615eabe1cc0f81bc7d4295f81` |
| Native execution receipt | `01571ff9eff767a89523c2c3608696eaa2c16a61a28a58300b7182e0e1695b0f` |
| endpoint.out | `457c1c6d15caae482ee83b63c2f2500bacc388508d1a358e5b98050591899a55` |
| endpoint.engrad | `685cc58cdfb98aa2e02048dbb00bf86c160ff5eb95cf90da9b7c9b0a4e747c04` |
| Runtime sidecar | `e9785f03a2da7d53b69b9634f937325273686a26ddf12b0c0436be489039767f` |

Actual artifacts reside under
`workspaces/lanm_coupled_mechanics_20261017/lbt_whole_7ccn_before_La_v1/execution/`.
The continuation must retain exactly `after_La`, `before_Dy`, `after_Dy`, with
one invocation each and the original four-cell comparison denominator.
