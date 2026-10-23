# Saved baseline oxygen-stiffness leverage

Question: can the proposed oxygen-only stiffness directions materially change
the four observed lab-z curvatures within the declared small parameter domain,
or are they locally too weak? This calculation uses only the admitted original
model19 field outputs and complete model/reference comparison. No native code,
new molecular call, fitted coefficient, optimization or extra field is used.
Root executes this saved calculation once after source review.

## Fixed derivative and units

The frozen design is `../INTERPRETATION.md`: a drives the permanent field;
`sigma_i=f(r_i)[b_M+c_M eta_i]` changes oxygen radial inverse susceptibility,
with `eta_i=sum(j!=i) f(r_j)`. f is the original S2U radial expression with native
exponents 4 Å⁻¹, normalized at 2.50 Å. At zero parameters,

```
dH_i/db_M = f(r_i) n_i n_iᵀ / alpha_O
dH_i/dc_M = eta_i dH_i/db_M;       d chi/d a_M = 0
```

Use actual key-over-parameter POLARIZE resolution (water-O type36, 0.837 Å³)
and all ordered induced sites. Keep the complete native site/name/type/parent
mapping and exact real-water graphs. The complete OHH polarization groups have
zero intragroup d/p factors and identical intergroup factors; Ed=Ep here and
the same external field is added to both. This equality justifies the single
response-column formula for these four fixtures, not for arbitrary proteins.

Recover every site's three-component response separately at h=.001 and .0005:

```
x_i = [u_i(+h) − u_i(−h)] / [2h * DEBYE * BOHR]
g = K / (HARTREE * BOHR)
d chi/db_M = −g sum_i [f(r_i)/alpha_O,au] (n_i·x_i)^2
alpha_O,au = alpha_O,A³ / BOHR³
```

Equivalently the prefactor with alpha in Å³ is `K*BOHR²/HARTREE`; retain the
actual native constants rather than setting g to one. x is in a0³ and each
Jacobian entry is a0³ per dimensionless parameter. Permanent moments have zero
field derivative at these fixed coordinates. The site-z sum must reproduce
the already admitted moment-slope estimator. Also report metal/O/H site-z sums
as a descriptive partition of this model, without quantum site attribution.

Native saved site components have six decimals in Debye. Propagate
`eps_x=0.5e−6/(h*DEBYE*BOHR)`, `eps_s=||n||1 eps_x`, and
`eps_(s²)=2|s|eps_s+eps_s²` through both derivatives. Derived increment/paired
bounds sum absolute component bounds. Compare both step estimates and retain
any discrepancy exceeding this printing allowance; this is not a solver-error
bound. The original finite induction residual cannot bound solution error
without conditioning evidence.

## Outputs and interpretation

Retain four state rows and six columns ordered a_La,b_La,c_La,a_Dy,b_Dy,c_Dy.
Derive La/Dy trimer−dimer and Dy−La directions without new observations. Report
singular values, stiffness-column conditioning and rank at relative 1e−6 for
both six-amplitude and c=0 nested four-amplitude families. The two drive columns
are exactly zero; maximum curvature-only rank is four. Repeating steps or
estimators cannot create new physical rank. Dipole/energy/force constraints
remain necessary to identify the drive amplitudes.

For each existing comparison row report a conservative **linearized** absolute
correction envelope: the smaller of the coefficient-box bound
`.2 sum |J_b,c|` and the bound implied by every `|sigma_i|<=.2`.
Apply absolute coefficient sums to derived contrasts. Retain the nested-four
envelope, site-print allowance and target printing bound. These bounds do not
solve for coefficients or assert simultaneous reachability of the targets.

One saved column x=A⁻¹Bz gives the exact baseline derivative under the stated
linear/equal-dual assumptions. It does not supply A⁻¹ on all oxygen-radial
directions, the full spectrum, or a finite-parameter stability/feasibility
certificate. An observed target shift outside the envelope exceeds local
first-order leverage only. A shift inside the envelope is not a successful
model or proof of attainability. No physical tolerance is introduced.

`CALCULATE.py --config CONFIG.json --out FRESH_DIRECTORY` writes immutable
RESULT.json and REPORT.md with one CPU/thread. Source preparation performs only
syntax validation; no repeated native-admission suite or synthetic framework.
