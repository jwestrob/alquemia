# Possible compact solvation transfer — inspected, not executed

This is a candidate after the pending native DFT results, not a second submitted
experiment or a default change. Parent asked to wait for those results before
choosing a new physics experiment.

## Existing observations

For the unchanged 1H4I/4MAE neutral expansions, native DFT context−core ΔR is
−2.0327768101303505 / −4.425466673350552 kcal/mol. Its included CPCM dielectric
term changes contribute diagnostic Δcontrasts +32.47643850359059 /
+30.834958996581587 kcal/mol. Native OMOL vacuum changes are −22.85019933769945 /
−35.247174426855054 model-kcal. These terms belong to different methods/densities;
subtracting the CPCM diagnostic does not produce an actual vacuum DFT endpoint.
The comparison only motivates separating solvent response from charge-feature
sensitivity. Charged contexts also add real Asp groups about4.5–6.8Å from the
metal. Their field can physically favor La, even without an implementation error.

## Coherent approximate energy expression

At one fixed geometry/state, a possible new score would use:

`E_composite,M = E_OMOL,vac,M + E_GFN2,ALPB,M − E_GFN2,vac,M`.

Then `R_composite = R_OMOL,vac + Delta_solv,Ca − Delta_solv,La`.
The low-level difference includes the solvent-induced GFN2 electronic response
and native ALPB terms at identical nuclei. It replaces no CPCM term because the
chosen high-level learned endpoint is vacuum. It is an approximation to the
solvation of an OMOL-represented state, not a self-consistent hybrid density or
protein dielectric boundary. The entire compact cluster sees the solvent
boundary. No term is to be added to native CPCM DFT. No thermodynamic water
occupancy, nuclear entropy or binding free-energy claim follows.

## Software evidence and limitations

The [ORCA6.1 native xTB manual](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/semiempirical.html#native-gfn-xtb-and-gfn2-xtb)
confirms native GFN2, ALPB support without an external executable, the default
xTB mixer, parameter export and300K electronic smearing. Its restart warning
means both states would start independently in empty directories, using
NoAutostart and no old xtbw. The same page recommends native ALPB over more
expensive continuum routes for xTB. This documents implementation availability,
not La/Ca accuracy or a validated hybrid.

The already archived whole-protein run exported native Ca/La parameters but
produced zero converged endpoint energies. Its9141-atom startup/memory failure
must not be described as a failed chemical test of compact GFN2. Existing actual
exports include Ca4s/4p/3d and La5d/6s/6p shells, reference valences2 and3. They
contain no per-solvent parameter export. Any new compact pilot must therefore
verify actual ALPB settings/output and converged states, not infer a successful
solvation calculation from this old export. No parameter, radius or charge is
to be fitted to classification labels. Native defaults would be frozen globally.

Potential numerical checks would compare actual vacuum/solvent convergence,
state identity, sane charges, retained parameters and a repeated tighter SCF
pair before interpreting small signals. Native GFN2 alone is a byproduct cheap
comparator, separate from whether its solvent transfer helps OMOL. A fixed panel
must include PQQ and the complete alpha/GGR structural matrix; no per-case backend
choice. Exact finite task count/resources await the current DFT result and a
separate declaration. No compact GFN2 calculation has been launched here.

## Equal total charge is informative, with limits

The archived original core charge also perfectly separates PQQ families:
Ca-endpoint Q=−2 for14Ca-family references and Q=−3 for11La-family references.
The added Asp makes the three Ca-family outliers Q=−3, identical to ten neutral-
expanded La-family Ca endpoints. All13 expanded cases still separate at the
2.315459model-kcal gap. This is a harder equal-total-charge structural comparison,
so narrowing alone does not establish lost predictive accuracy. The smaller
margin limits robustness claims and motivates verification. All are still
consumed calibration structures; donor composition, geometry and grouping remain
confounders. This subset is a descriptive diagnostic, never a new favorable
calibration or conditional routing rule.
