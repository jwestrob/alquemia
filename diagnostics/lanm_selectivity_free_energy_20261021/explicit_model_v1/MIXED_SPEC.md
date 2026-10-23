# Fixed hydrated donor follow-up

Question: does the same unqualified La/Dy common-water AMOEBA candidate retain its donor-dependent exchange contrast when the ion has an explicit eight-donor first shell?

Use exactly the chemistry agent's prepared eighteen fixtures: monodentate acetate+7 waters, bidentate acetate+6 waters, or acetamide+7 waters; both metals and ligand radial positions −0.15/0/+0.15 Å. Reuse the original completed aquo8 La/Dy controls. The source panel's contact-screened geometry is frozen; there is no geometry search, optimization, new parameter, solvent model or fitting.

Copy the original explicit-model key/parameter bytes, changing no settings. Map only the metal placeholder to type500 La or501 Dy. Preserve all ligand/water types, bonds, coordinates and inventories. The common-water Dy transfer remains unqualified and La-specific POLPAIR corrections remain La-specific.

One zero-energy preparation binds the twenty-slot manifest. Root then owns one serial CPU Slurm run: eighteen new `analyze E` plus eighteen new `testgrad Y N` invocations, at most36 native program invocations. The two aquo8 energies/gradients are admitted from exact pinned native outputs. No new finite-difference or isolated-fragment calls occur. First native/parser/nonfinite failure stops new scheduling and preserves every slot.

Report all total/component energies and gradients, ligand radial derivatives, and the nine matched contrasts:

`[E(Dy mixed) − E(Dy aquo8)] − [E(La mixed) − E(La aquo8)]`.

This is the balanced metal-differential of a fixed ligand-replacement reaction: free ligand and released-water terms cancel between the metals. **Individual absolute ligand-replacement energies are unavailable**, since no isolated water in this aquo geometry has been evaluated. Never substitute the old, geometrically different isolated water. Fixed-coordinate energies and gradients are not hydration, binding or selectivity free energies. Radial gradients are descriptive and receive no new numerical/accuracy gate.

Keep the exact source, inputs, raw outputs, subprocess CPU/wall cost and unavailable denominator. No ALPB, GPU, dynamics, protein scoring or parameter promotion. Compare with the separately declared reference/MACE matrices before choosing any coupled protein Hamiltonian.
