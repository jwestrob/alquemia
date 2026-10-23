# Electrical error in the scale of the consumed references

This supplement preserves the five-step preliminary report. It summarizes the frozen consumed targets with the **same weights as the original objective**, without fitting or a new physical threshold. Its calculated candidate ratios use only the preserved step5 snapshot. The root-reported step6 values below are progress context, not terminal admission.

| Observable | Reference RMS with original weights | Origin error RMS | Step5 error RMS | Root-reported step6 error RMS |
|---|---:|---:|---:|---:|
| Pair Cartesian component, D | 2.271140 | 0.614949 | 0.634584 | 0.63842 |
| Full-shell C Cartesian component, D | 0.895190 | 0.072309 | 0.083650 | 0.086877 |
| Four conditional susceptibility states, a₀³ | 17.975770 | 0.034161 | 0.448326 | 0.536302 |

Step5 error/reference-RMS ratios are 27.94%, 9.34% and 2.49%; the root-reported step6 values correspond approximately to 28.1%, 9.70% and 2.98%. These are descriptive aggregate scales, not tolerances or per-observable relative errors. In particular, the large susceptibility increase relative to the tiny origin residual does not mean susceptibility itself changed by that factor.

## Individual reference scales and cancellation

The pair objective gives each of four parent contexts equal weight and averages its three Cartesian components over every water. All 34 literal pair targets are retained. Individual pair-vector magnitudes are 3.8936–3.9844 D; their component RMS by parent is 2.2480 (W8 La), 2.2833 (W8 Dy), 2.2575 (W9 La), and 2.2954 D (W9 Dy). Thus the roughly 0.64 D Cartesian residual is substantial on the pair scale even though its increase from the already imperfect origin is only about 0.023 D. Component-relative percentages are unsuitable for near-zero transverse coordinates.

The full-shell C targets contain strong geometric cancellation:

| Parent | Reference C z component, D | Largest absolute reference x/y component, D |
|---|---:|---:|
| W8 La | 0.000965548 | 8.55×10⁻⁷ |
| W8 Dy | 0.000832664 | 2.14×10⁻⁷ |
| W9 La | 2.244444077 | 5.88×10⁻⁷ |
| W9 Dy | 2.139826124 | 2.43×10⁻⁷ |

A single C percentage mixes nearly cancelling W8 moments with order-2 D W9 moments. The terminal synthesis should retain all four signed three-component errors; aggregate RMS alone cannot show which cancellation is preserved. A percentage of a near-zero W8 or transverse component is not a useful accuracy measure.

The four fixed susceptibility targets are energy-curvature estimates at h=0.001 atomic units:

| State | Reference, a₀³ |
|---|---:|
| La dimer | 15.472821701 |
| La trimer | 22.016524156 |
| Dy dimer | 13.158106697 |
| Dy trimer | 19.880690161 |

Their RMS scale is about 18 a₀³. The terminal report should list each model value and signed error; an aggregate 0.54 a₀³ error does not imply the same percentage for every state. Derived trimer-minus-dimer reference increments are 6.543702455 (La) and 6.722583464 (Dy), with Dy−La only 0.178881010 a₀³. A small error relative to the individual states can matter much more for this cancelling difference. These derived differences were not four additional fitted observations, and their errors cannot be inferred from the aggregate state RMS alone.

The response references retain their original finite-field estimator, frozen electronic-state/ECP/PBE0 and numerical limitations. The additional field rows remain descriptive; no new estimator is chosen. Neither these scale ratios nor surviving energy/force gates establish physical electrical accuracy, transfer, or a global limitation of the Hamiltonian.

Evidence: `ELECTRICAL_CONTEXT_VALUES.json` pins TARGETS `3f4f90c8`, OBJECTIVE `af7d26f6`, all34 pair vectors, four C vectors, four susceptibility references, original weighting and frozen step5 VALUES. Reference-scale extraction used one thread, 0.0058 s wall/0.0043 CPU-s, zero molecular/optimizer calls and no live fit read in that extraction. The accepted step6 artifact had been read once during the preceding context inspection; it is not re-extracted or numerically re-admitted here. No further growing-fit reads are needed before the root supplies terminal pins.
