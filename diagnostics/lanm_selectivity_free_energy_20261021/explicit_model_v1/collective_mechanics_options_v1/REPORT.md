# Next mechanics direction: separate pair driving fields from collective response

**The next repair should address electrical response, not add an energy offset to one coordination contrast.** The admitted W9 dipoles show two opposing errors: too little response in isolated metal–water pairs and too little suppression when those pairs form a shell. A small common-basis metal–two-water comparison can establish whether the collective error is already present at three-body order. Its exact finite design is in `TRIMER_DECLARATION.md`; preparation performs no molecular calculations.

## What the completed evidence establishes

| W9 z dipole, Debye | La quantum / model | Dy quantum / model |
| --- | ---: | ---: |
| Full interaction C | 2.244444 / 2.107200 | 2.139823 / 2.017552 |
| Sum of pair interactions P | 3.929194 / 2.916725 | 3.959697 / 2.729250 |
| Collective remainder B=C−P | −1.684749 / −0.809525 | −1.819874 / −0.711698 |

The approximately 6% full-interaction errors conceal much larger pair and collective errors. In B, Dy−La is −0.135124893 Debye for the reference and +0.097827 for the model. Independent saved arithmetic reconstructed all vectors and printing bounds from the admitted rows: 80 checks pass, zero new calls. W9 differences exceed the combined printing bounds (at most 0.000132 Debye per component); those bounds do **not** measure reference accuracy. W8 cancellation by symmetry is not a physical success criterion.

Common-parent-basis electronic decomposition plus separately retained original D4 gives the paired energy B of +9.036179 kcal/mol versus model +2.797568. Full C remains wrong by −4.663007 kcal/mol under that alternative reference. Whole-shell CP and equalizing fragment bases did not remove the discrepancy. Both fixed positive EXCHPOL amplitudes reduced model B further and worsened the contrast. Metal-alpha probes also supplied no admissible improvement. The failed consumed donor development is separate evidence against simply extending the same fit; it does not establish the same microscopic cause.

Saved native component sums confirm that **the model's paired B** is polarization: +2.79756793 kcal/mol. The vdW contribution cancels to about 1.2e−14; permanent multipoles and reported bonded terms cancel to zero. Full-energy/component closure differs only at printed precision. This identity does not make the quantum remainder a uniquely defined polarization energy. Quantum charge populations also disagree on the small differential charge-transfer sign; functional, frozen-f and retained SCF stopping limitations remain.

## Concrete repair architecture to assess

The smallest useful response extension is **metal-pair-scoped separation of the permanent field that drives induced dipoles from the induced–induced interaction tensor**. Existing radius/epsilon terms cannot do this; the present original-Thole pair coefficient is used in both direct and mutual screening. A stronger isolated-pair response and stronger collective suppression must be constrained separately, using P and the three-body energy/dipole, then checked against full C/B. No sign or coefficient is selected here.

This must come from one stationary energy functional, schematically

```
U_response(mu; R) = 1/2 mu^T A(R) mu − mu^T E_direct(R)
A(R) = alpha_inverse − T_mutual(R)
```

with the exact native polarization-group exclusions and dual-field conventions retained or consistently rederived. This schematic is an implementation requirement, not a claim that replacing a single native array realizes it. A symmetric positive-definite response Hessian is required throughout the finite panel; forces must include both field and tensor coordinate derivatives. A parameterized field edit without its matching energy/gradient is unacceptable.

Only interactions involving the declared metal may change. Every water–water, water–protein and metal-absent term must retain the original route, including inside a metal-containing system. Smooth radial kernels, intact local frames, and permutation/rotation covariance preserve differentiability. No discrete coordination count is needed. The eventual implementation must reproduce the baseline when disabled, preserve pure-water energies/forces, and pass energy/force derivatives and induction convergence before any fitting. GPU and periodic support would require separate work.

### What the installed engine actually supports

| Source capability | Consequence for this task |
| --- | --- |
| `damptholed` / `dampthole`, `induce`, `epolar1` | Direct and mutual paths and derivatives already exist. A **pair-scoped** extension needs source changes across field, energy/analysis and derivative routes; it is not currently a safe keyword-only feature. |
| Alternate direct Thole (`use_tholed`) | Global flag is activated by atomic direct damping. Zero pair-direct coefficients then return undamped factors, rather than falling back to original Thole. A metal-only trigger can therefore alter water–water direct fields in its cluster even if isolated pure water is unchanged. This rules out a simple global switch as the proposed repair. |
| `CHGPEN` | Class parameters exist, but any penetration sites switch the global polarization route away from original Thole. Permanent penetration is pair additive; only its self-consistent field consequences can change B. A metal-only key does not provide a complete, water-preserving density model. |
| `CHGTRN` / `echgtrn*` | Installed CT is an explicit exponential **pair energy**, not variational transferred charges. Its direct contribution cancels in matched B/E3 and cannot repair these dipole differences by itself. |
| Bond/angle charge flux | Depends on covalent geometry. It supplies no new metal-coordination response at these fixed intact-water geometries; adding metal bonds would change the model topology. |
| `extra`, `extra1`, `extra3` | User hooks are empty. A smooth metal–water–water energy could be implemented with analytic derivatives, but it cannot repair failed zero-field dipoles unless coupled to an electrical response model. |

These distinctions follow the pinned local implementation. AMOEBA+ literature provides precedent for distinct short-range electrostatic/response terms, not validated La/Dy coefficients or permission to replace the current water model ([Liu, Piquemal and Ren, 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6615954/)).

## One finite discriminator before implementation choices

Measure electronic E3 and common-origin dipole mu3 for the nearest and farthest water pair in each consumed W8/W9 parent, both metals. Use the same parent basis for all seven inclusion–exclusion terms. Only the trimer and water-pair terms are new: **16 reference single points and at most 16 model analyze-EM calls**, no gradients, fitting, new geometry, or amplitude sweep. Both selected pair results remain in the denominator regardless of whether they help interpretation.

This probes angular/overlap dependence while every M–O distance remains 2.50 Å. A far pair may have a symmetry-suppressed dipole; it remains in the panel and is never replaced using an energy or dipole result. If the pair deficit and wrong collective dipole response already appear at three-body order, the direct/mutual response separation has a concrete finite target beyond one full-shell contrast. If E3 fails but mu3 agrees, an additional short-range energy term becomes plausible locally, while the established full-shell dipole failure still needs a response repair. If both three-body observables agree, a three-body-only patch is poorly motivated: higher-order response or reference uncertainty remains. None of these cases uniquely identifies CT, penetration, Pauli repulsion or polarization.

The eight sampled E3 values cannot reconstruct full B. An exhaustive pair inventory would need 256 new reference fragments at this two-metal/two-parent scope, and even that omits metal-containing orders above three. It is not proposed. The pending correlated-reference uncertainty work remains complementary; no result here changes original development targets, passes the failed candidate, releases donors, or qualifies bulk sampling.

## Status

Saved/source evidence only, one local thread. `EVIDENCE.json` records exact source/result pins and component arithmetic. `PREPARE.py` renders the finite literal masks/templates and model inputs from admitted parents; no execution method is exposed. Root owns review and any later native release.
