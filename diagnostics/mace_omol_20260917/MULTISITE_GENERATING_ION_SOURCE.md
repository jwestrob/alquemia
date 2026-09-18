# Multisite generating-field ion-template clarification

Before any quantum or metal-score calculation, first bridge preparation failed
on actual Ca template matching. Archived standalone tip3p.xml contains neutral
water only; no ions. Preserve failed multisite_density_prepared_v1 and oldconfig.

Use the installed explicit amber19/tip3p.xml as the generating-field water/ion
parameter source, alongside the unchanged protein.ff19SB.xml. Record its full
path/hash in each new physical bridge. Verify the charges of every real retained
water against its original standalone tip3p.xml by actually typing the same
water subtopology; require exact per-atom charge equality. BackgroundCa must
match the actual CA template and charge+2. No Lennard-Jones energy is calculated
or enters the fixed point-charge field. Positions, waters, protonation, total
charge and ff19SB protein parameters remain unchanged. No ion/La parameters
are invented, no dependency installed, and AMOEBA/GK remains a separate backend.

Same five paired input preparations, no DFT/MACE/solvent/force calls. New outputs
multisite_density_config_v2 and multisite_density_prepared_v2. This is a fixed
explicit parameter-source correction before scores, not score-based tuning.
