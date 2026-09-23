# LanM: last actual test and prepared within-series test

Historical preparation checkpoint below. Subsequent execution is authorized by
EXECUTION_AGREEMENT.md; actual submissions/results are recorded separately.

**Prepared, not executed.** All six real Hans/Mex EF1–3 sites are ready. The
scientific choice and physical/effective Dy state policy remain with Jacob/root;
zero new MACE, GFN2 or DFT endpoints and zero Slurm jobs have run in this track.

## Last LanM calculation

The latest completed LanM scoring found is job1199508, completed
**September15,2026 at20:51:57 PDT** (September16 03:51:57 UTC). Native r2SCAN-3c /
CPCM water on repaired peptide-amide-v3 Hans sites gave the ordered La/Ca scores
**46.2984,45.2519,50.9235kcal/mol**. GGR's comparison score was3.0575. These are
the archived CN8-gauge electronic contrasts; generic-v3 has no canonical-PQQ
absolute classification bands. This tested La/Ca, not lanthanide-series affinity.

Source: `diagnostics/benchmark_set_20260915/SCORING_RESULTS_1199508.md`, its JSON
and scheduler receipt; six Hans endpoints were part of a10endpoint/321s/64CPU
batch (20544allocatedcore-s, no GPU). Later within-LanM molecular tests were not
located; newer mentions of `lanm_qmmm` often name the environment, not a LanM case.

## What prior series work established

| Actual earlier test | Finding |
|---|---|
|August4 Amber12-6-4 La/Dy capability,20trajectories|Failed to produce the four required La-bidentate/Dy-monodentate contrasts; weak EF4 controls also remained positive.|
|August4 Dy outer-water metadynamics,4×15ns|Water stabilizes monodentate Dy by about2.07kcal/mol relative to dry monodentate, but remains about2.22 above bidentate; formal resultNO_CALL.|
|August5 Hans Protenix,16predictions|High confidence did not guarantee metal placement; correctly placed La sites were all incorrectly monodentate.|
|August5–6 Dy QM/MM correction chain|Final production transfer qualification missing; chain stopped with cancelled job1165540. No qualified rescue result.|
|August2 LBT3 lanthanide-series DFT|Different protein: severe geometry/aquo-reference dependence and open-shell failures; no valid series-affinity method to reuse.|

These are preserved failures/limits, not fresh ideas to rerun. Exact source files
are pinned in ARTIFACTS.json.

## Actual prepared inputs

| Protein/source | Ordered sites | Atoms/site | Total charge, La and Dy | Source waters |
|---|---|---:|---:|---|
|Hans,8DQ2 A,La crystal|La201/202/203 →EF1/2/3|50|−1|0/0/0|
|Mex,8FNS A,Nd crystal|Nd201/202/203 →EF1/2/3|44|−1|A326+A345 /A322+A328 /A320+A331|

Hans reuses the actual September15 pH5 prepared states. Mex uses the same pH5
standard-residue preparation and3.3Å source selection, followed by the existing
amide repair. Every deposited heavy atom is preserved exactly; no missing heavy
atoms, terminal atoms or residues were added. Actual deposited water O positions
are retained and added water H coordinates are frozen identically for La and Dy;
their orientational ensemble is not sampled.

The first preparation attempt was rejected because PDBFixer selected the first
Ser52 alternate conformer, while the existing occupancy-aware carver selected
the other deposited conformer. Only Ser52 CA/CB/OG differed (maximum1.594Å
coordinate component atOG); no donor coordinates were repaired. The versioned
second attempt preselects the existing occupancy-aware conformer before adding
H, and its full source-heavy difference is0.0Å. Both attempts are preserved.

## Model, observable and coverage

The native OMOL checkpoint includes both elements and a multiplicity embedding.
Use physical total-cluster multiplicity1 for La and6 for Dy, with exact charge−1.
This is a declared physical-state input; even an actual batch6 check does not
verify localization of f electrons, a ground state or spin–orbit physics.

Installed native GFN2 exports5d/6s/6p shells and three valence electrons for both
La and Dy; no explicit4f shell. Its approximate solvent transfer therefore uses
effective multiplicity1 for both. Hans has150 native valence electrons; Mex136.
All-electron counts are Hans254La/263Dy and Mex234La/243Dy. Physical and effective
parity checks pass separately. No shared singlet-only production adapter changed.

The prepared scope is12MACE plus48GFN2 calls (24cells, each initial followed by
one verified same-cell native continuation). The first HansEF1Dy MACE endpoint
and vacuum cell are capability subsets of those totals. The actual batch6 and
actual Dy effective-electron output checks remain **unrun**. Both matchingGBW
andxtbw are required to activate continuation; output must positively report
`INITIAL GUESS: XTBRESTART`. A failed/unconfirmed continuation has no substituted
score. One continuation does not establish full electronic-branch convergence.

For each site report the balanced electronic exchange
`D=[E_Hans,Dy−E_Hans,La]−[E_Mex,Dy−E_Mex,La]`, with each E the native OMOL vacuum
energy plus native GFN2(ALPB−vacuum). PositiveD means stronger relativeLa
preference in Hans. Aqueous-ion references and element offsets cancel in this
balanced cycle. Each host's composition/water inventory appears unchanged on
both sides. Metal-dependent relaxation, protonation and population effects do
**not** cancel merely because the reference does. EF1–3 remain an ordered vector;
no fitted weighting, summed affinity, site labels, Kd or populations are produced.

The [primary Hans/Mex study](https://www.nature.com/articles/s41586-023-05945-5)
supports the expected direction atpH5: Hans exhibits much stronger light/heavy
discrimination than Mex. Hans main apparent affinities are68pM La versus2.6nM Dy,
with heterogeneous/partially folding Dy response. Exact Mex La/Dy numerical
uncertainties were not recovered here: the supplementary download failed DNS and
the web viewer could not retrieve it. The qualitative direction is supported by
the preserved primary text; no approximated quantitative target is fitted.
Both proteins are consumed development controls, not blind validation.

## Cost and implementation status

Seven focused real-artifact tests pass in2.093s, zero skips. They test actual
site/water/state maps, native archived restart parsing/recipe, path containment,
unchanged shared selector, rejected charge corruption, source-conformer mapping,
balanced algebra from real archived energies and null unrun outcomes. No molecular
success was fabricated. Native/GPU scientific integration is still unrun.

Reference execution receipts suggest a seconds-to-minutes pilot after scheduling:
warm native MACE job1203745 did64larger-context evaluations in5.696s process time,
13s allocated H200 time; native continuation job1210185 did8larger-context cells
in27s on64CPUs (1728core-s). Six comparable GFN waves would be about162s, but Dy
SCF behavior and initial guesses are unmeasured, so this is a scale estimate,
not a guaranteed duration or finite project budget. No GPU-hours regime is
anticipated for44–50atom fixed sites. Actual costs/failures must replace estimates.

Candidate implementation, inputs and finite60call ledger are sealed under
`workspaces/lanm_series_followup_20260923/prepared_v2/`. Native manifests are
contained; continuation manifests are materialized only from real successful
same-cell seeds. Dedicated warm worker, native runner adapter and collector are
runnable after the pending scientific choice; see COMMANDS.md. Baseline/default
and old LanM calculations are unchanged. This readiness package is not a completed
within-lanthanide discrimination result.
