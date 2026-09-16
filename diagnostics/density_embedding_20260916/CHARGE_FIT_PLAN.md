# Cheap charge-extraction check, declared before execution

Under Jacob's [autonomous-research authorization](AGREEMENT.md), test one
alternative charge representation on the same four saved vacuum endpoint
wavefunctions: native ORCA 6.1.1 CHELPG. This adds four standalone utility calls,
no new SCF/DFT endpoint, changed physical state or biological comparison.

Motivation: MBIS population analysis took about 575–748 seconds per endpoint,
while exact-density potentials at all environment sites took a four-CPU job
27 seconds. A global model still needs a charge representation for the current
solvent solver. A cheaper extraction is useful only if its potential quality is
adequate; runtime alone does not validate it.

Use `orca_chelpg saved.gbw`, with native defaults: grid 0.3 Å, outer extent
2.8 Å, COSMO exclusion radii, no dipole constraint. Record actual native echoes
and executable hash. These sampling radii do not change the protein's physical
cavity. Run one configuration uniformly for all four states; no per-protein
choice or tuning against class/partition results.

Evaluate against the already computed quantum potentials at (a) the fixed
exterior ESP probes and (b) every actual environment charge coordinate. Report
RMS/relative RMS, total charge, direct-interaction error and runtime versus MBIS.
The historical exterior criterion remains relative RMS <=0.10 OR absolute RMS
<=0.005 au, but passing it alone does not establish reliable interaction energy.
Report all discrepancies continuously; this is a development comparison, not
calibration or a new global score. No automatically promoted charge scheme.

Retain failed utility attempts and unsupported ECP behavior explicitly. Use
isolated copies of GBW, densities and densitiesinfo; one thread per utility,
four concurrent tasks in an allocation. Expected scale seconds to minutes;
actual cost is recorded, with no CPU/time stopping budget.

Primary references: [standalone utility](https://www.faccts.de/docs/orca/6.1/manual/contents/utilitiesvisualization/utilities.html#orca-chelpg),
[charge model and defaults](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/population.html#chelpg-charges).
