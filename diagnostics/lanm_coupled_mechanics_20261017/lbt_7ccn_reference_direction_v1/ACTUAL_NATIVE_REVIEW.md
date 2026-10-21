# Independent completed native review

**Minus 15847 and plus 15848 both meet the unchanged singlepoint admission rules.** No native blocker was found to the declared saved collector. Native convergence does not mean all printed residuals pass or that the reference is physically qualified.

## Exact evidence

| Artifact | Minus SHA256 | Plus SHA256 |
| --- | --- | --- |
| RESULT.json | `8e73f55a229e89bd2cec3f7b9e1ec4553ff82d46417b9936e18cc2454294c5db` | `9cedffa6dcd0ce1a53dd8cea5b29a78b11477f8c4b8a3053e33d5dcd5f37ae32` |
| Native output | `fffa1676b73aa11a2d215cab422ecc895f214a63aa29a8ba5b4280ecbe1aebdf` | `2130829d5881a130d1ad6a2c911aa96d735824b7f4cae2806c4a9c259a735444` |
| Native receipt | `28f3cac75c8a4b490a3dcd83688b0717364360595d192d08c169fef5fd9e19ba` | `995826b7468720f95c46d08b4b1fccb2a190a70565dffa4ca496271fd8151ae5` |
| MANIFEST.json | `a10fe4f3b7d29443407fb46f6a9567c2a2ba6f63c466c5759a87e9c6c3c868ce` | `0dcd90f6bacb0689d67bb4eef3759f52a985aa53dc22215559e28ee5704c7029` |

Both results bind source `cd4396c02cf50905307acbf90772d062ae4a6a12aaf9914e01562167dca1ab32`, CONFIG_RUN `dd4875fabd6d17a12d8c73ad7254068ba83a7070b33d079c62b3b7c2dcd74ff7`, release `62cd2ead`, preparation `779412e9` and the identical frozen DIRECTION `f20be914`. Twenty-seven selected artifact pins were independently verified, including inputs, physical/position records, runtime artifacts, executable and native runner/renderer. Both task signatures exactly match their prepared singleton.

Each result is complete with its requested row available, one runner invocation and one new DFT call. All central, other-metal, MACE, ALPB, GFN, MM and optimizer counts are zero. The unused sibling row in each singleton result remains not attempted; the complete pair consists of the two distinct selected rows.

## Native settings and stopping route

Both raw outputs and receipts show ORCA 6.1.1 normal termination, return code 0, no launch error, and 12 MPI ranks with one thread/rank. Actual molecular headers confirm RHF, charge −2, multiplicity 1, 1030 explicit electrons and La ECP46. All 256 atoms and their prepared physical-ID order are retained. PModel executes explicitly; the input retains NoAutostart and the original PBE0-D4/def2-TZVP/def2-J/RIJCOSX/DefGrid3/VeryTightSCF recipe, without a field or implicit DFT solvent.

These are singlepoints: EnGrad is absent, `engrad_path=null`, and neither receipt contains an engrad artifact. The runner's `expected_engrad_complete=true` means no requested gradient is missing; it does not report a newly calculated gradient. Runtime memory remains 48000 MiB allocation, 25% reserve and MaxCore 3100 MB/rank. MaxIter 125, ConvCheckMode Total+1el-Energy and ConvForced 0 are unchanged.

Both outputs explicitly stop by **Gradient check signals convergence**, after 18 cycles.

| Printed residual | Minus | Plus | Threshold | Both pass? |
| --- | ---: | ---: | ---: | --- |
| Energy change | 6.5920e-9 | 9.3787e-9 | 1e-9 | No |
| Maximum density change | 4.1204e-5 | 5.7703e-5 | 1e-8 | No |
| RMS density change | 2.7305e-7 | 3.5022e-7 | 1e-9 | No |
| DIIS error | 1.2496e-3 | 1.2446e-3 | 1e-8 | No |
| Orbital gradient | 2.0366e-7 | 2.1443e-7 | 2e-6 | Yes |
| Orbital rotation | 6.5258e-7 | 5.5211e-7 | 2e-6 | Yes |

`all_printed_residuals_pass=false` remains explicit. Both preserve the early small-gap warning −0.061, without treating it as a final gap or stability result. Final gridx3 integration completed, with printed exchange corrections −0.000002460 Eh (minus) and −0.000002232 Eh (plus). These scalar corrections are not force-error estimates.

## Energies, costs and scope

The unique final native energy strings match the saved rows exactly:

- Minus: **−6765.683736715161 Eh**.
- Plus: **−6765.681619165231 Eh**.

Saved scheduler events report COMPLETED/exit 0:0, with 21623 seconds on 12 CPUs for minus and 21616 seconds on 12 CPUs for plus. Native receipt elapsed times are 21615.944008 and 21608.010095 seconds. Allocated CPU capacity is not measured utilization; root owns final accounting.

The two energies are ready for the frozen decimal-string collector. No derivative was recomputed during this review. One-direction consistency, full-gradient accuracy, electronic-state stability, solvent, a physical minimum and preference are not established by native admission. The original no-retry rule and all state/SCF limitations remain.

Review used saved JSON, literal metadata equality, native text and hashes only. No scientific imports, molecular calls, source changes or submissions were made.
