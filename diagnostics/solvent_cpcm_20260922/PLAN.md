# Matched native GFN2 CPCM-water challenger

Jacob approved the matched-solvent track; parent coordination and authorization
are recorded in `diagnostics/nikasha_next_phase_20260922/PLAN.md`. Preserve the
released native MACE + native GFN2 ALPB scorer, its reference and all old outputs.

## Fixed question and initial scope

Does replacing only the low-level water-continuum transfer provide useful class
information at fixed real PQQ context geometry and chemical state? Initial sources
are consumed 1H4I, 4MAE, canonical Q9Z4J7 and A0A3F2YLY8, each Ca/La: **eight new
native-GFN2 CPCM singlepoints**. Reuse actual matching native MACE vacuum energies
and GFN2 vacuum/ALPB endpoints. No nuclear optimization, folding, DFT or label fit.

Use ORCA 6.1.1 `Native-GFN2-xTB`, default parameters, `NoAutostart`, native xTB
mixer, electronic temperature 300 K, default convergence criteria, and the
previously qualified MaxIter500 ceiling. Only `CPCM(Water)` replaces ALPB.
No functional, basis, dispersion, parameter, dielectric, radius or proton change.
The actual parameter exports and effective solver settings must agree before
an archived vacuum endpoint can be paired with CPCM.

Freeze default Gaussian-vdW water cavity, scale 1.2 and charge-surface density
5.0/Angstrom^2, without DRACO/SMD/CDS. Installed historical native-DFT outputs
confirm water epsilon 80.1510, refractive index 1.3328 and printed radii Ca 2.7720,
La 2.4000, C 2.0400, N 1.8600, O 1.8240 and H 1.3200 Angstrom. These are the
installed defaults, not fitted radii; the actual native-GFN2 CPCM output must
confirm its own boundary. Other element defaults, if present, must be printed and
retained. The metal-specific cavity difference is explicit and immutable.

## Energy and acceptance

`E_new,M = E_MACE,vac,M + (E_GFN2,CPCM,M - E_GFN2,vac,M)`.

`R_new = E_new,Ca - E_new,La`. Positive changes are more La-like on this model's
raw scale. Convert eV/hartree to kcal/mol once. This replaces the ALPB transfer;
it never adds a second solvation term. CPCM and ALPB have different physical and
empirical content. This is an electronic descriptor, not a full binding free
energy. No aquo reference or absolute adaptive decision is available initially;
the old-band display is a transfer check only. A new reference requires all 25
canonical members, with exactly the existing extrema/minimum-gap rule.

Require real converged receipts; exact coordinates, charge/spin, atoms, parameters
and 300 K state; charge closure <=5e-4 e and max|atomic charge|<=4 e. Inspect parsed
CPCM setup, surface radii, reaction-field energy and actual mixer behavior. A
keyword alone or silent zero reaction field does not establish CPCM operation.
Missing/unsupported/failed endpoints remain unavailable with all artifacts.

### Predeclared solver contingency

If the installed CPCM path changes/disables the native mixer or fails to include
the field, do not mix it with old native vacuum energies. Record that result and
prepare a separate technical stage with the **same eight physical states** and
matched ordinary-SCF CPCM/vacuum endpoints (16 calls; CPCM reruns included if
needed). Use genuine explicit-block MORead from the pinned native orbitals,
unchanged parameters/temperature/convergence, and MaxIter500. This is a numerical
compatibility check, not a search for favorable classifications. No stage-two
execution is included in the initial eight-task submission; notify the parent
with observed behavior before using this contingency. Earlier ordinary-SCF
failures remain relevant and are not evidence against the working native method.

## Execution and continuation

Existing task runner, eight concurrent eight-rank endpoints, 64 CPU /128 GiB on
an available explicitly named GPU-partition CPU host, with no GPU requested.
No arbitrary project time/CPU budget; finite manifests and scheduler limits apply.
Record wall/core allocation, failures and preparation/reuse. No shared executor,
dispatch, scoring default, reference, queue priority or running job is changed.

After a coherent pilot, prepare the original 25 canonical +3 crystal +2 unknown
PLM source manifests and the declared 225-fold transfer population. Root coordinates
their execution; complete canonical calibration before evaluating fold decisions.
All cases are already consumed development data; no fresh blind claim.
