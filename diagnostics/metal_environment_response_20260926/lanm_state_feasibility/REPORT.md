# Native r2SCAN-3c can represent Dy explicitly; the research adapter needs extension

**Recommendation:** retain native ORCA6.1.1 r2SCAN-3c for the first local La/Dy reference. Use physical La(III) singlet and Dy(III) sextet in the one-metal electronic region. No basis substitution or effective-singlet workaround is required by the installed library. Actual embedded Dy gradient execution/convergence is still untested: this read-only audit establishes capability and state consistency, not a successful molecular endpoint.

## Installed evidence

The installed `orca_exportbasis` utility exported the default `def2-mTZVPP` La/Dy orbital bases and ECPs, saved in `native_mTZVPP_La_Dy.bas`; binary/export hashes are in PINS.json. This reads the installed library and performs no SCF or molecular calculation.

| Element | Native ECP | Replaced electrons | Explicit neutral-atom electrons | f-electron treatment |
|---|---|---:|---:|---|
| La | def2-ECP | 46 | 11 | La(III) is the closed-shell, empty-4f hypothesis |
| Dy | def2-ECP | 28 | 38 | 4f remains explicit; Dy(III) has the 4f9, five-unpaired-electron hypothesis |

The export prints Dy f basis functions and the 28-core ECP, contradicting any assumption that this DFT recipe shares native GFN2's f-in-core treatment. Dy's ECP has lmax5; La's lmax3. The exported citations identify Dy's Ce–Yb28-core ECP as Dolg/Stoll/Preuss,JCP90,1730–1734(1989), and La46-core as Dolg/Stoll/Savin/Preuss,TCA75,173–194(1989).

The [official basis table](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basisset.html) assigns def2-mTZVPP to r2SCAN-3c with default def2-ECP through Rn. The [composite-method documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/3cmethods.html) supports elements through Lr and describes the native r2SCAN/D4/gCP combination. Neither establishes lanthanide-response accuracy; it establishes that no elemental substitution is needed.

## Exact archived EF3 state counts

Recomputed directly from all six pinned XYZ files and charge−1, not copied from prose. STATE_AUDIT.json retains their paths/hashes.

| Source | Atoms | La all/explicit electrons | La alpha/beta | Dy all/explicit electrons | Dy alpha/beta |
|---|---:|---:|---:|---:|---:|
| Hans8DQ2 | 50 | 254/208 | 104/104 | 263/235 | 120/115 |
| Hans8FNR | 50 | 254/208 | 104/104 | 263/235 | 120/115 |
| Mex8FNS | 44 | 234/188 | 94/94 | 243/215 | 110/105 |

Every multiplicity/electron parity passes. ORCA expects multiplicity2S+1:1 for La,6 for Dy. The Dy sextet means alpha−beta5, not spin input6 or an effective singlet. UKS follows automatically for multiplicity>1; make the research state explicit and parse executed UKS. [ORCA wavefunction conventions](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/hftype.html).

These are old pocket fixtures for counting/readiness only. The next embedded preparation must recount after including second-shell groups/links and record neighboring metal occupancy. A single-Dy QM region does not inherit the whole-protein two-Dy multiplicity11. Outside-region spectator-metal spin and exchange are absent from a fixed-charge representation and must be recorded as a limitation.

## Analytic gradients and remaining qualification

Native `EnGrad` is the intended route, retaining the complete native dispersion/gCP and ECP gradient components. The existing La embedded outputs demonstrate this route for La; no Dy-specific embedded EnGrad has been executed in this audit. The documented native composite and UKS routes reveal no required numerical-gradient fallback, but **do not mark Dy analytic gradients qualified until its actual output contains the SCF gradient, SHARK ECP gradient, dispersion gradient, gCP completion, matching engrad energy/coordinates, and point-charge derivative file.** Numerical gradients must fail explicitly, not run silently. [Run types](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basics.html).

This is scalar-relativistic, spin-free DFT with an explicit open f shell, not a spin-orbit-coupled lanthanide free-energy model. Retain SCF residuals, total S2 (spin-pure sextet expectation8.75), atom-resolved spin populations, local f occupancy evidence and restart sensitivity. A converged sextet does not prove the spin density stayed on Dy, establish a unique f-orbital occupation, or resolve near-degenerate electronic roots. Do not change multiplicity or select a favorable solution after inspecting the response sign.

If method dependence becomes decisive, a separately versioned scalar-relativistic all-electron treatment can be considered on the smallest subset; it is not needed merely to represent Dy. In particular do not casually append X2C with external-potential transformation: the [relativistic manual](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/relativistic.html) says the optional inclusion of external charges in the X2C model potential is not available for gradients/Hessians.

## Research-only adapter changes required before execution

Current `metal_environment_reference.py` assumes Ca/La only, six tasks, singlet inputs and a Ca/La response formula. `mace_omol_vacuum.parse_endpoint` also hard-codes the element table, 46-electron La ECP, even parity, multiplicity1, and La-only ECP-gradient proof. Do not modify these production assumptions globally just to pass Dy.

Create a separately named research adapter once source preparation is fixed:

1. Enumerate metal/state specifications: atomic number, oxidation hypothesis, input charge, physical multiplicity, native ECP identity/core count and f treatment. Input generation uses the actual multiplicity, not a singlet-only helper.
2. Validate `Nexplicit=sum(Z)-Q-sum(Ncore)`, nonnegative integer alpha/beta from `(Nexplicit±(multiplicity−1))/2`; require output charge, multiplicity, NEL and exact ECP inventory to match. Support Dy atomic number66 explicitly.
3. Require ECP-gradient evidence for both metals, preserve ordinary engrad/pcgrad coordinate/energy/unit checks, and add UKS/S2/local-spin output records. Do not turn unavailable spin evidence into zeros.
4. Use La/Dy response contrasts with explicitly named order; do not import Ca/La bands or solution offsets. Same-element A/B energy differences cancel the element-specific ECP total-energy zero; raw cross-element totals are not an affinity.
5. Keep unsupported spectator-metal parameters, mapping, boundary terms and full-hybrid cross interactions visible. Existing point-charge electronic response can proceed only under its explicitly narrower ledger, not as a claimed complete hybrid force.

## Next bounded step and audit scope

After the Ca/La gates and actual EF3 embedding preparation, run one source's La/Dy embedded origin pair using this exact native recipe and physical states; parse states/analytic-gradient evidence before expanding to the predeclared perturbations/remaining sources. Prefer this direct qualification over changing Hamiltonian speculatively. Freeze geometry, occupancy and spin policy before chemistry.

No Slurm job, energy evaluation, optimization, old vacuum retry, new model install or production edit occurred. One successful installed basis export and six real-file state counts were performed.

Reproduce the library export from repository root:

```bash
/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca_exportbasis \
  -b def2-mTZVPP -a La Dy \
  -o /tmp/nikasha_native_mTZVPP_La_Dy.bas
```
