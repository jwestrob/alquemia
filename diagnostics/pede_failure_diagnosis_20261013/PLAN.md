# PedE failure diagnosis from completed artifacts

## Question and scope

Jacob asks whether the two wrong recovered PedE models reflect poor folding,
fine coordination geometry, or a flaw in the scoring/accommodation model.
Compare all five recovered two-Ca PedE sources, including the three correct
ones, with compatible existing PedE references and archived geometry evidence.
These are inspected development cases; no blind validation is claimed.

Three parallel tasks: geometry and actual folding confidence; saved OMOL/solvent
energy decomposition; and review of the implemented physical model plus prior
DFT evidence. Root adds the frozen-calibration/decision comparison across all
125 recovered sources and the 25 canonical calibration sources.

Use existing source, origin and selected-proposal coordinates and exact paired
states. Geometry metrics are descriptive continuous measurements, not new
acceptance cutoffs or classification rules. Record absent confidence fields and
incompatible or missing reference scores explicitly. Gross source admission
does not establish experimental correctness of fine donor positions.

Energy analysis uses the existing expression and unit constants, with no new
energy evaluation. Each relaxation work is separated into OMOL vacuum and
GFN2 ALPB-minus-vacuum terms. Preserve every candidate and both metal rows.
Compare margins under each policy's own frozen bands: raw score movement alone
does not measure changed confidence because static/accommodated bands differ.

## Execution

All structural calculations, scientific algebra and analysis code run through
Slurm. One-CPU saved-artifact jobs on a healthy standard node are sufficient;
the scripts enforce allocation and write distinct receipts. Metadata/code review
and reporting may run on login. No new MACE, xTB, DFT, folding, optimization,
protonation or water-state calculation is included. No labels, bands, bounds,
production defaults or original outputs are changed.

Root owns submissions and integration. Agents own disjoint geometry/,
energetics/ and model_review/ subdirectories. Save actual results and costs;
retain partial failures. A numerical/backend acceptance check does not establish
biological accuracy, and a geometric correlation does not prove causality.

The final opinion must separate observed differences, supported explanations
and unresolved alternatives, then identify the next discriminating test if
needed. Do not invent a rescue classifier or choose favorable samples.
