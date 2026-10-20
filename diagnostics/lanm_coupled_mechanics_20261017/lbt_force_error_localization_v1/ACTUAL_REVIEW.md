# Independent actual localization review

**Job 15797 completes the declared saved-only localization. No binding, coverage or recorded-bookkeeping failure was found.** The common force discrepancy is distributed across peptide heavy atoms. The smaller metal-dependent discrepancy has a more concentrated side-chain pattern, especially in 7CCN.

## Evidence and coverage

| Artifact | SHA256 |
| --- | --- |
| RESULT.json | `c789abb0071a10b1b2dc06db2e24b1c087c44a15cd3fc2682add532cf24ca7d1` |
| LOCALIZE.py | `415027352268cbabb8efdda6e3ea414c420f58e1ff17ba2140b58a9b966f6ee7` |
| CONFIG.json | `7297b790f09b36d64fa6f480856cae003aa23b048dcb25d73e5f25f11b328514` |
| Parent 15720 RESULT | `f42cffcf5717a988ed569db925b840f02dfb8bd3b4c04750b77e33bc231c4790` |
| Inventory | `e2c196d1b5961cc74d50647e0565780928c3c1c103b282d563833f4b7e003957` |

The result is at `workspaces/lanm_coupled_mechanics_20261017/lbt_force_error_localization_v1/run/execution/RESULT.json`. I verified all 47 recorded artifact pins. All eight endpoint metrics, all four matched metrics, inherited caveats and each endpoint's saved vacuum error vectors are exact copies of their parent authorities.

Every pair retains 256 unique physical IDs in order: one metal, 119 hydrogens, 69 backbone heavy atoms and 67 remaining protein heavy atoms. Every pair has all 17 residue records plus the metal. H2/H3 remain hydrogens; OXT remains backbone heavy. The full atom tables and fixed top-ten summaries are present. No endpoint or pair is unavailable.

All recorded checks pass: 777 for each 7CCO pair and 776 for each 7CCN pair. These cover saved difference vectors, endpoint/matched norms, original projections and the squared-norm identity under the unchanged absolute `1e-10` / relative `1e-12` bookkeeping rule. I inspected those saved outcomes without rerunning their arithmetic.

The result records zero molecular calls, model loads, geometry evaluations, optimization and fitting. Application time is 0.286792 seconds. The inspected completion event records COMPLETED, exit 0:0, one second on one CPU; its SHA is `3fbeab08d3a7622018304d603c0d274821eb04f18d01050c1b0330fc2109b458`. That snapshot still marks notification delivery queued/unconfirmed; root owns acknowledgement.

## Reported localization

Here “error” means the saved **MACE-minus-PBE0 vacuum force difference**, not error against an exact physical oracle. Common means `(e_La + e_Dy)/2`; difference means `e_Dy − e_La`. RMS units are kcal/mol/Å. Fractions below are the producer's fractions of the corresponding full squared norm.

| Pair | Common RMS | Difference RMS | Common: backbone fraction | Difference: remaining-heavy fraction |
| --- | ---: | ---: | ---: | ---: |
| 7CCO before | 3.511769 | 0.501500 | 0.525631 | 0.463829 |
| 7CCO after | 3.516169 | 0.476265 | 0.528318 | 0.473215 |
| 7CCN before | 3.513145 | 0.726763 | 0.520255 | 0.860066 |
| 7CCN after | 3.581379 | 0.563471 | 0.525580 | 0.708271 |

**Common part:** Remaining protein heavy atoms separately account for fractions 0.437920–0.446539; hydrogens account for 0.031919–0.035168 and the metal for 0.000313–0.001281. Trp9 is the largest residue contribution in every pair, with fractions 0.126525–0.132165; its CG is also the largest common-error atom. Other large entries span backbone and side chains. Terminal atoms appear among the large entries, but the result does not support a terminal-only explanation or transplant the earlier hybrid's MM-terminal diagnosis onto whole MACE.

**Metal difference:** In 7CCO, Asp7 is the largest residue contribution, with fractions 0.160001 before and 0.139866 after; backbone, other side chains and the metal also contribute. In 7CCN before, Asp7/Glu14/Asn5 lead with fractions 0.315698/0.222661/0.181460. Asp7 OD1 is the largest atom at 6.242802 kcal/mol/Å. In 7CCN after, Asp3/Asp7/Glu11/Glu14 lead with fractions 0.230807/0.141072/0.126079/0.125537; Asp3 OD1 is largest at 4.148825 kcal/mol/Å. This locates a residual discrepancy at chemically relevant side-chain atoms, without assigning its cause.

## Limits and decision

Squared-force fractions are not energies or fractions of deformation work. Group size affects these totals; the saved group RMS values remain available separately. Atom location does not uniquely identify a faulty interaction, parameter, electronic method or boundary. The common mean is a two-metal decomposition at each fixed geometry, not a fitted transferable correction.

Shared ALPB is excluded from these vacuum differences and receives no validation. Physical Dy multiplicity 6 MACE versus frozen-f/effective-singlet DFT, native SCF/gradient limitations, conditional termini and vacuum-dianion assumptions remain. Force accuracy, physical-state qualification and preference claims remain false.

This closes the declared descriptive localization without another molecular calculation. It does not justify fitting or reopening closed model branches. The separate 15545 reference addresses the later model-converged La geometry; these 7CCN “after” records remain the original step-80 endpoint.

Review used hashes, metadata counts/equality and reported values only. No scientific vectors were recomputed, no scientific package was imported, and no job or production/source change was made.
