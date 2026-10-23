# From intrinsic exchange to the LBT3 assay

Saved-evidence derivation only. No new literature, validation outcomes, constants, molecular calls or fitting. This supplements the existing [observable contract](../OBSERVABLE_CONTRACT.md) and preserves its unresolved metadata.

## 1. The exact one-site bookkeeping

Use the reaction **PLa + Dy(aq) → PDy + La(aq)**: negative exchange free energy favors Dy. Choose a chemically specified apo reference state P₀, a free hydrated M³⁺ reference species M₀, and a bound reference state H₀,M. La and Dy reference complexes must have the same peptide proton inventory. Equilibrate conformations and hydration within every ensemble claimed by the calculation; a fixed geometry is not a binding free energy.

First use a **concentration convention in one fixed medium and temperature**. Define

\[
K_{d,0,M}=[P_0][M_0]/[H_{0,M}],\qquad
Q_P=[P]_{\rm unbound,total}/[P_0],
\]
\[
\alpha_M=[M]_{\rm unbound,soluble,total}/[M_0],\qquad
Q_{H,M}=[PM]_{\rm bound,total}/[H_{0,M}].
\]

The numerator of α excludes peptide-bound metal. QH includes all states counted as bound by the assay, including relevant protonation, hydration, conformation and peptide–metal–buffer states. These factors are equilibrium ratios at the specified reservoir conditions, not fitted offsets.

Then

\[
K_{d,\mathrm{cond},M}=K_{d,0,M}\,Q_P\alpha_M/Q_{H,M},
\]
\[
\boxed{\Delta G_{\mathrm{ex,cond}}=\Delta G_{\mathrm{ex},0}
+RT\ln(\alpha_{\rm Dy}/\alpha_{\rm La})
-RT\ln(Q_{H,\rm Dy}/Q_{H,\rm La})},
\qquad
\Delta G_{\mathrm{ex},0}=RT\ln(K_{d,0,\rm Dy}/K_{d,0,\rm La}).
\]

The **same apo polynomial QP cancels exactly**. Holo populations generally do not. Stronger unbound Dy sequestration adds a positive correction; greater Dy-bound statistical weight adds a negative correction. This formula assumes a reversible, soluble, one-site equilibrium with fixed free reservoir composition. Precipitation, kinetic capture, oligomerization, ligand depletion that changes speciation, or a signal selective for only some bound states requires the corresponding larger equilibrium/measurement model.

## 2. Protonation and side reactions

A state polynomial can be written schematically as

\[
Q_X=\sum_s\exp[-\beta\{G^\circ_{Xs}-G^\circ_{X0}
-\Delta n_{H,s}\mu_H-\sum_j\nu_{j,s}\mu_j\}],
\]

using one consistent standard-state convention. The additional species and proton reservoirs must be balanced. Protonated peptide states belong in QP or QH; hydrolyzed/anion-complexed unbound metal belongs in α. A bound hydroxo or ternary MES complex belongs in QH, so correcting α alone does not account for it. A proton-release MES constant cannot be inserted as an ordinary deprotonated-ligand formation constant without its proton factor.

At otherwise fixed medium, the proton-linkage sign is

\[
\frac{\partial\Delta G_{\rm ex,cond}}{\partial\mathrm{pH}}
=RT\ln10\,[\langle\nu_H\rangle_{\mathrm{binding,Dy}}
-\langle\nu_H\rangle_{\mathrm{binding,La}}].
\]

Positive νH denotes protons taken up from solution during binding, including both bound and unbound speciation. Equal nominal peptide charges do not establish equal mean proton uptake. Buffer blank-heat subtraction in ITC does not remove these equilibrium weights; proton-buffer contributions to measured enthalpy are not a standalone correction to binding free energy.

## 3. What cancels, and what requires evidence

| Contribution | Paired La/Dy exchange |
|---|---|
| Same unbound apo protonation/conformation ensemble | QP cancels algebraically under the same construct, medium and equilibrium model |
| Common concentration standard c° | Cancels for matched 1:1 association conventions |
| Free-metal MES, hydrolysis or anion sequestration | Only the common part cancels; αDy/αLa remains unknown where not measured |
| Bound protonation, folding, hydration, electronic or ternary states | QH,Dy/QH,La remains; include them in sampling or an explicit state sum |
| Leading charge-only activity terms | Equal charges support cancellation in that approximation; specific ion/complex interactions need not cancel |
| Terminal chemistry or a different sequence | Does not become irrelevant by taking a metal ratio; it changes the two binding ensembles differently |
| Shared assay calibration/active fraction | Not a general cancellation theorem for a nonlinear fitting model |

For thermodynamic activities, define dimensionless Kd° and activity coefficients γ. The microscopic concentration constant obeys

\[
K_{d,0,M}^{c}=c^\circ K_{d,0,M}^{\circ}\Gamma_M,
\qquad \Gamma_M=\gamma_{H_{0,M}}/(\gamma_{P_0}\gamma_{M_0}).
\]

Thus conversion from a thermodynamic exchange to the above concentration convention adds **RT ln(ΓDy/ΓLa)**. This term is not added if the computed/empirical quantity already uses that same fixed-medium convention. Equal +3 ions and equal complex charges alone do not establish ΓDy=ΓLa.

## 4. Match what the simulation actually includes

An aqueous leg that samples chloride association already includes that contribution to its unbound ensemble. A bound leg that samples hydration/conformations already includes those bound weights. **Do not apply α or QH a second time.** For a shared reference convention, the omitted-state correction from a simulation ensemble to an assay ensemble is

\[
RT\ln\!\frac{(\alpha_{\rm Dy}/\alpha_{\rm La})_{\rm assay}}
{(\alpha_{\rm Dy}/\alpha_{\rm La})_{\rm sim}}
-RT\ln\!\frac{(Q_{H,\rm Dy}/Q_{H,\rm La})_{\rm assay}}
{(Q_{H,\rm Dy}/Q_{H,\rm La})_{\rm sim}}.
\]

If temperature, solvent composition or reference chemical potentials also change, retain their differential transfer/standard-state contribution separately. The ratios alone do not convert a NaCl box into a MES-buffered assay. Equivalently, transfer the actual bound and aqueous ensembles: δGex = δμH,Dy−δμH,La−δμM,Dy+δμM,La. Finite box, counterion and restraint conventions must be reconciled too.

## 5. Keep the two assays separate

**2020 ITC:** compare with its retained conditional ΔG difference at **283.15 K, pH 6, 50 mM MES and 100 mM NaCl**, after the above mapping. The assay's terminal chemistry is unresolved. An amidated calculation is therefore a construct-conditional comparison; no statistical error bar repairs that ambiguity.

**2025 competition:** first predict the competition signal from conditional constants in the same medium. Under the existing one-site/linear-signal assumptions, with p the total unbound peptide ensemble concentration,

\[
PT=\frac{pT_{\rm tot}}{K_T+p},\quad
PI=\frac{pI_{\rm tot}}{K_I+p},\quad
P_{\rm tot}=p+PT+PI.
\]

At half of the zero-competitor signal let y=PT₀/2, p=yKT/(Ttot−y), and A=Ptot−p−y. Then **IC50=A(1+KI/p)**. A shared Cheng–Prusoff denominator cancels from the *reported* Ki ratio, but the exact depletion offset A does not:

\[
K_{I,\rm Dy}/K_{I,\rm La}
=(IC50_{\rm Dy}-A)/(IC50_{\rm La}-A).
\]

The already consumed, same-buffer conditional KT must not receive a second MES correction. Existing saved arithmetic found a small depletion sensitivity for the reported concentrations; this was not a raw-curve refit. The affinity-assay temperature, numerical uncertainty, I10/Y10 text discrepancy and actual fitting implementation remain unresolved. Without temperature, a conditional **ln(Kd,Dy/Kd,La)** can be reported; a kcal/mol value requires an explicitly labeled temperature scenario.

## 6. Honest conditional reporting

Report: (1) the computed exchange with construct, protonation/state ensemble, temperature, medium and standard-state conventions; (2) sampling uncertainty; (3) the explicit omitted-state/activity correction formula; (4) separately labeled sensitivity scenarios and unresolved metadata. If the correction is not bounded by primary evidence, report **no numerical matched-assay prediction or empirical mapping confidence interval**. Do not set the missing correction to zero or infer it by forcing agreement with consumed LBT3 data.

The existing evidence admits neither a matched La/Dy MES correction nor a finite empirical bound. The Eu illustration and simple hydrolysis screen are not substitutes. Mechanistic qualification and intrinsic exchange work can proceed, but an experimental-agreement claim remains conditional on this mapping. No assay term explains a fixed-geometry force or energy failure merely by being unknown.

Provenance: [PLAN](../../PLAN.md), [target report](../REPORT.md), [observable contract](../OBSERVABLE_CONTRACT.md), [saved mapping results](../ACTUAL_MAPPING_AND_COLLECTOR_REVIEW.md), and [bounded MES follow-up](../mes_primary_followup_v1/REPORT.md). Exact read-file pins are in `PROVENANCE.json`; no held-out outcomes were opened.
