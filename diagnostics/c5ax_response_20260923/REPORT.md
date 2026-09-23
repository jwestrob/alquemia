# C5AXV8: one persistent replica ambiguity, not a solver failure

Completed 2026-09-23 using archived results only. No molecular evaluations,
new coordinates, thresholds, optimization or production changes.

## Main finding

La-conditioned sample 3 is **correct under released static scoring and DFT**.
Its first new inconclusive call occurs when the union-context method applies
its own reference bands to the **same numerical source score**. Subsequent
accommodation changes this score little. The strict native solver leaves it
inconclusive: R = −405457.848334666 kcal/mol in model energy units, 1.468395
below the frozen strict La boundary. This is not evidence that a numerical
failure still needs rescue, or that this particular structure must be repaired.

| Sample 3 method | R, model kcal/mol | Own-reference call |
|---|---:|---|
| DFT | −405390.913164 | La |
| Released static | −405457.712287 | La |
| Union static | −405457.712287 | Inconclusive |
| Local adaptive | −405457.877547 | Inconclusive |
| Union adaptive | −405457.877547 | Inconclusive |
| Union precision | −405457.848121 | Inconclusive |
| Union precision, strict native scalar | −405457.848335 | Inconclusive |

Released La_min is −405459.111382; union-static La_min is −405456.684469;
strict-adaptive La_min is −405456.379940. These are different calibrated
protocols. Changing the band to recover this case would be outcome-based tuning.

## All sources retained

The canonical source is **La sample 1**, not sample 0. All five La sources
have the same 172-atom union chemistry, Ca/La charges −3/−2, singlet endpoints
and ten caps. Their state-signature SHA256 is
`7df1549c6afe8a77b6635b80bce2fb141ae415bcb88d7507b913f35b88c8df19`.
All five Ca-conditioned sources retain their existing preparation exclusions:
the selected ion is 25.02–25.21 Å from the intact PQQ. No ion was moved or
source silently substituted. That is a placement problem, not absent PQQ atoms.

The following origin/pool values all use the same **strict scalar recipe**;
the native and solvent columns sum to the origin contrast.

| La sample | Native origin R | Solvent contribution | Composite origin R | Pooled R | Accommodation ΔR |
|---|---:|---:|---:|---:|---:|
| 0 | −405335.685819 | −117.426047 | −405453.111866 | −405452.236296 | +0.875569 |
| 1, canonical | −405341.000313 | −114.778053 | −405455.778366 | −405456.379940 | −0.601574 |
| 2 | −405337.932800 | −108.994221 | −405446.927021 | −405445.107240 | +1.819781 |
| 3 | −405337.820629 | −119.898715 | −405457.719344 | −405457.848335 | −0.128991 |
| 4 | −405340.887333 | −92.692985 | −405433.580317 | −405435.358279 | −1.777961 |

The source spread is much larger than the accommodation effect. The solvent
contribution spans 27.205730 kcal/mol versus 5.314494 for native MACE.
Composition is matched across these union sources, so this is an actual
geometry-dependent model response, not different atom counts. Repeatable
electronic energies do not establish physical accuracy of that response.
The strict noncanonical La4 median is −405448.671768, correctly La-supported;
Ca5 and the balanced group remain unavailable. These replicas are one biological
protein, not five independent labels.

## What the motions actually do

Source roles are Glu198, Asn275, catalytic Asp317 and extra-acid Asp319.
Sample 3 starts with the shortest Glu198 OE1–metal distance (2.455 Å) and
longest Asn275 OD1 distance (2.619 Å) among its five La sources. Its paired-origin
force selector spends three of four coordinates on Glu198 (χ1/χ2/χ3), plus
Asp319 χ1. Asn275 is fixed. Other replicas select different common subspaces.

The Ca proposal reaches the unchanged 0.8 Å displacement boundary. Glu OE1
moves to 2.315 Å and extra-Asp OD2 to 2.394 Å. The La proposal is interior
(0.669 Å extent), with Glu OE1 at 2.563 Å and extra-Asp OD2 at 2.553 Å.
Each metal chooses its own proposal in the common pool. At strict scalar
settings, Ca work is −4.753122 kcal/mol and La work is −4.624130: substantial
individual relaxation nearly cancels in discrimination (ΔR = −0.128991).
Native differential work is −1.214491, offset by +1.085500 solvent work.

At the La proposal, selected normalized vacuum-MACE gradients are at most
0.000103 kcal/mol/Å, while omitted Asn275 χ1 remains +5.896, Asp319 χ2 −4.849
and Asp317 χ1 −3.897. These are derivatives in physically normalized torsion
directions, not composite solvent-aware forces or estimated energy gains.
Other replicas also retain large omitted loads, so this is a potential
**generic subspace limitation**, not proof of a C5AX-specific repair.
Metal translation also remains loaded, but is outside the proposed angular test.

## Numerical and physical interpretation

Changing SLSQP stopping precision moved the sample-3 Ca proposal by at most
4.76×10⁻⁵ Å, changing its native energy by 6.11×10⁻⁸ kcal/mol. The pooled
score shifted +0.029426. Tightening native scalar stopping then shifts it only
−0.000214, with the same choices. Neither explains the remaining 1.468395 gap.

The supported next question is whether a small, uniform extension of the
paired-force-selected angular subspace accesses useful **differential**
accommodation. Adding Asn by name would be a case-specific choice. A six-mode
extension should instead keep the first four selected modes and use the next
two physically independent directions under the existing alternating
differential/individual-load rule. Its feasibility is documented separately;
no such optimization has run here. Improvement is not guaranteed: lowering
both endpoint energies can leave R unchanged or move it in either direction.

## Artifacts and verification

`workspaces/c5ax_response_20260923/AUDIT_v2.json` contains all ten sources,
exact original/strict pool pins, force arrays, mode projections, coordinates,
source roles and state signatures. SHA256:
`7bee1f61c3ba4ae06f295b15268176a5eb0224654fcc53b7f82db593dc82804c`.
The earlier v1 snapshot remains intact. The script replays actual coordinates
and analytic Cartesian-force projections, verifies paired maps and saved active
gradients, and requires all ten distinct source identities. Strict225 comparison
SHA256: `6ed4c1f4acf4359010922ff3f088c5b36105caf52dc92c734e27dd7e35d0df46`.
Only local file reading and algebra ran; molecular calls and cluster allocations
were zero. Baseline, labels and calibration artifacts remain unchanged.
