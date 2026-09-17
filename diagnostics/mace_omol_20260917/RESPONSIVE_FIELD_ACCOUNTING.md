# Native embedded energy accounting, verified before launch

ORCA6.1 manual, section2.2.6:
https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html#inclusion-of-point-charges
The external point-charge module adds electronic and nuclear core/field
interactions. Interactions among external charges are excluded by default;
explicit `%method DoEQ false end` freezes that convention. The original
native method/basis/ECP/D4/gCP stays unchanged. No second direct term is added.
The input file specifies environmental charge in elementary-charge units and
coordinates in Angstrom. Original coordinates are used at full precision.

Installed executable is the pinnedORCA6.1.1 shared OpenMPI4.1.8 build; its
package directory has license PDFs but no scientific manual. The official
6.1 documentation was consulted online rather than claiming a local manual.
Actual old native outputs underworkspaces/density_embedding_20260916/prepared_v1/
report matching point-charge counts (9094/9087), successful embedding, and
converged states. They are evidence of syntax/accounting capability, not a
substitute for the eight new states. New collection requires native version,
charge count, ECP/electron count, unchanged coordinates, analytic gradients,
D4/gCP, successful SCF and matching immutable execution receipts.

ORCA utility documentation, sections9.2.2 and9.2.10:
https://www.faccts.de/docs/orca/6.1/manual/contents/utilitiesvisualization/utilities.html#orca-vpot
Native vpot reads the geometry/basis and requested density, evaluates an ESP
at user coordinates in Bohr, and writes ordered potentials. CHELPG fits the
molecular ESP under the molecular total-charge constraint. We will independently
check finite potentials, source atom/probe ordering, endpoint/paired fit error,
and direct coupling from the new wavefunction. This documentation alone does
not qualify a new fit, boundary or solvent model. Any unexpected external-field
contamination or potential singularity fails the new representation gate.

The electronic response diagnostic is embedded energy minus vacuum energy
minus the already measured exact vacuum-density coupling. At fixed geometry
and electronic state it should be nonpositive, with the predeclared+.05kcal
numerical allowance. No clamp or favorable branch selection is allowed.

The final candidate adds full GB from the embedded projected QM charges plus
unchanged environmental charges and the reused short MACE context. No old GB
term satisfies that task. Environment-only vacuum energy cancels between
endpoints; the full GB environment-only term cancels but is retained for audit.
The density is not self-consistent with GB. Nuclear geometry and classical
charges do not relax. Reference/class/combined-gradient/relaxation fields stay
null. The plan contains the full expression and frozen acceptance criteria.
