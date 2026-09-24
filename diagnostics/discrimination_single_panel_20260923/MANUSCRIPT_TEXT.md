# Text for the single-panel Nikasha discrimination figure

Prepared from completed results for the integrated PLM manuscript. This is a
writing deliverable; no new molecular calculations or analyses were performed.

## Figure caption

**Figure X. Nikasha discrimination across alternative predicted structures of
PQQ-enzyme references.** Each row represents one of 25 reference protein groups;
blue and purple indicate Ca-class and La-class references, respectively. Circles
denote structures predicted with Ca and triangles denote structures predicted
with La. Every starting structure is evaluated with both metals using the
consistent-pocket accommodation protocol. The declared comparison comprises
five Ca-conditioned and four La-conditioned structures per protein. Horizontal
positions show the Ca-minus-La composite-energy contrast after subtracting the
midpoint of the protocol's frozen decision band. Values to the left and right
of the grey band support Ca and La classifications, respectively; orange rings
mark inconclusive results. Horizontal lines span the available source scores.
Open diamonds show the equal-weighted mean of the Ca-conditioned and
La-conditioned median contrasts, displayed only when all nine required sources
are available. Right-hand labels give available/declared sources for each
protein. Of 225 declared structures, 207 received the expected classification,
one was inconclusive and 17 were unavailable; no scorable structure received
the opposing class. These are previously examined structural replicas of 25
proteins used during method development. The contrasts are classification
descriptors expressed in model kcal mol−1, rather than measured affinities.

## Results paragraph

We assessed Nikasha's sensitivity to input geometry using alternative
Ca-conditioned and La-conditioned predictions of 25 known-class PQQ proteins.
With consistent pocket membership and bounded local accommodation, 207 of 208
scorable structures received the expected metal class, one was inconclusive
and none received the opposing class (Fig. X); 17 additional structures were
unavailable because of preparation exclusions. Holding the same supported
pockets at their original coordinates and using the same numerical solver
gave 199 correct, one opposing and eight inconclusive calls. Accommodation
therefore resolved seven inconclusive calls and one opposing call across four
protein groups, without losing a correct call in this matched comparison.
These improvements required geometry-dependent score changes and were not
explained solely by the separately calibrated decision bands. The results
support improved discrimination across alternative structures within the
examined PQQ reference population; they do not constitute independent
prospective validation beyond those protein groups.

## Connection to the PLM application

This figure reports the extensively tested ten-fold context policy. The
practical PLM SOP instead defines its context from three predeclared
La-conditioned structures per protein and reports their median classification.
Both methods evaluate Ca and La on a shared pool of candidate geometries within
each source. Keep the practical protocol's separately reported validation and
cohort results distinct from the 225-structure figure; do not describe this panel
as a completed PLM cohort scan or as an identical three-source execution.

The existing [methods/application draft](../plm_pqq_delivery_20260923/MANUSCRIPT.md)
and [SOP](../../docs/PLM_PQQ_SOP.md) provide that material. Update the cohort
paragraph from its actual completed outputs when available.

## Source records

- [Full reference transfer](../strict_native_transfer_20260923/REPORT.md)
- [Matched static comparison](../strict_static_ablation_20260923/REPORT.md)
- [Figure delivery and exact-value check](REPORT.md)
- Figure/data directory: `workspaces/discrimination_transfer_figures_20260923/accommodation_only_v1/`

The original three-panel comparison remains available for a supplement. The
single-panel figure keeps the strongest method, its complete legend, original
axis limits and all original source values.
