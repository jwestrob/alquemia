# Paired bound/aqueous geometry audit and preparation

Declared 21 October 2026. Read-only source/coordinate arithmetic and topology preparation only: no native executable, energy, force, minimization, dynamics or submission.

Audit the pinned `explicit_model_v1/protein_preparation_v2` amidated I10 LBT3 box, source atom joins, water/ion identities, charge, topology, matched La/Dy coordinates and periodic cell. Preserve it unchanged.

Root's explicit new convention is a neutral box for each leg: the same 49.324 Å starting cube and seven background NaCl pairs; bound peptide−4/metal+3 plus one extra Na (Na8/Cl7), aqueous metal+3 plus three extra Cl (Na7/Cl10). This matches background salt-pair concentration, not total ionic strength, ion activities or chemical potentials. Both metal substitutions retain +3; finite-size and ion-atmosphere sensitivity remain later requirements.

Rebuild aqueous water from the exact pinned 500-water, 24.662 Å source tiled twice in each dimension, using the existing box-builder's order and whole-molecule minimum-image reconstruction. Keep the bound metal's translated coordinate. Apply the same O–heavy 2.4 Å and water-H–heavy 1.5 Å overlap cutoffs to that metal; no peptide or peptide cavity is retained. Replace seventeen complete waters, ranked farthest from the metal with stable index tie-break, by seven sodium and ten chloride ions at their oxygen positions; require minimum-image ion–ion separation at least 5 Å. Both aqueous metal endpoints have identical coordinates, topology, water/ion order, and atom counts except metal label/type.

Retain native protein water types402/403, classes73/74, and Na405/class76, Cl414/class85 from the exact `amoebabio09.prm`. Reserved metal types900/901 remain unparameterized placeholders. Save complete atom/molecule provenance mapping, source pins, contact checks and both endpoint XYZ files. Do not write an executable metal/periodic key or claim a qualified Hamiltonian.

Reject duplicate atom coordinates, asymmetric/cross-molecule solvent bonds, ion–water O distances below 2.0 Å or ion–water H distances below 1.2 Å. These gross-overlap gates are preparation checks, not preferred coordination distances or force-field validation. Report actual minima.

Use one local numerical thread. Write a short report with concrete missing model/settings decisions and a finite first stability pilot conditional on a qualified candidate. No validation affinity or reserved molecular outcomes are read.
