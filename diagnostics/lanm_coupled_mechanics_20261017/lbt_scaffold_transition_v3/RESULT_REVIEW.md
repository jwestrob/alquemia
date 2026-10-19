# Independent actual-result review: fixed scaffold transition

Reviewed 2026-10-19T07:58:14Z by `/root/mex_health_review`.

**The complete four-cell MM attribution passes its declared implementation gates.**
The saved comparison continues to challenge the hybrid force balance on this
fixed deformation. It identifies both endpoint force offsets and a discrepancy
in their change; it does not identify unique stiffness or equilibrium parameters.

This review reads stored outputs, metadata and hashes only. All numbers below
are copied from the actual result/report; no scientific values were recalculated.

## Actual evidence and admission

- RESULT.json SHA256: `1b32076679677f814319618e3793477deb4c7e1d4df3d8cd7522f203fbdc4da2`.
- REPORT.md SHA256: `e3664c93a563beda704fbc54c4b414ee813478b6379c9530f81988a562fd34ed`.
- FORCE_TRANSITIONS.json SHA256: `8a7007c755c71935209a2ccb21a6344acce6436adade82a6a685ef7f266aa021`.
- Source `d736a589…` and config `bf88807d…` match the reviewed bundle. Their
  hashes, the report/force-transition pins and all four per-cell result pins
  were independently verified.
- Scheduler receipt records job 14787 COMPLETED, exit 0:0, 20 seconds, 1 CPU.
  The application reports OpenMM 8.5.1 Reference, 36 contexts and 816 energy/full-force
  queries: 9 contexts/204 queries and 187 retained interactions for each of the
  four available cells. New DFT/MACE/GFN/ALPB calls, geometries and optimizer
  steps are all zero.
- All four input/class/term stages pass. The three archived MM cells pass replay;
  after-Dy passes the declared new-component common-group control. Every stored
  group/class/full-system batch and all 33-terminal-atom class-sum checks pass,
  with unchanged energy 1e-8 kcal/mol and Cartesian-force 1e-5 kcal/mol/Angstrom
  thresholds. No failed check is hidden by an available row.
- The exact before/after coordinates, both metal reference cells, original
  frozen directions and full 256-atom vectors remain represented. After-Dy MM is
  explicitly component-only. Before-Dy whole-MACE force and after-Dy full-hybrid
  force/work remain unavailable; neither transition is invented.
- V1 and V2 zero-query failures remain preserved through the bound configuration.

## Substantive result

The stored La hybrid model-minus-reference force error has full 256-atom vector RMS
10.141916 before and 8.850041 after, and error-change RMS 4.463637 kcal/mol/Angstrom.
Thus a constant offset alone does not describe the discrepancy. Whole-MACE has
endpoint errors 3.511987/3.516438 and error-change 0.530091 on this same La pair:
its force change is descriptively closer to this reference, while endpoint
differences remain. This comparison supplies no universal MACE accuracy result.

Along the frozen Phe1-CA direction, La reference force changes
3.341027 to −40.851355, whereas hybrid changes 51.063006 to −1.895426.
The MM component changes 52.001900 to −0.969517. The initial MM bonded projection
is 51.505923, with bond/angle contributions 28.246713/22.531745. These stored
values localize the large initial projected force within the current MM model;
they do not assign a unique physical parameter defect.

The finite MM energy decrease is broader than that initial force projection:

| Stored La work, kcal/mol | Value |
|---|---:|
| Retained bonded | −10.682960 |
| MM LJ | −6.180039 |
| MM Coulomb | −3.189036 |
| QM–MM LJ | −3.229524 |
| Total MM | −23.281560 |
| Terminal-touching bonded union | −9.532816 |

All four bonded classes contribute favorable work. The terminal-touching union
contains 187 interactions, with Phe1-associated work −7.378008 and Ala17-associated
work −2.154808; the saved overlap count is zero. The interactions can include
neighboring residues, including CMAP terms. These are whole-interaction
contributions, not energies owned exclusively by terminal atoms or residues.
Bonded and nonbonded contributions both matter; the largest listed interaction
does not establish a single-term explanation.

The retained parent work comparison remains PBE0+ALPB La +12.279593/Dy +11.722656,
whole-MACE+ALPB +11.290076/+11.082401, and full hybrid La −28.244859 kcal/mol.
The new Dy MM work −23.280382 does not supply missing full-hybrid Dy work. Shared
groups 0–2 are physically identical at matched coordinates; this common scaffold
contribution is not a La/Dy preference measurement.

## Limits and decision

No factual or interpretive blocker was found in the generated report. Retain
the reference's failed aggregate printed-SCF-residual flags and unqualified
physical states, conditional free termini, vacuum dianion/post-added solvent,
spin-free frozen-f representation without SOC, nonstationary endpoints and
absence of an independent whole-reference numerical force qualification.

This is a closed, implementation-consistent diagnostic at two fixed endpoints.
It provides no Hessian, force-integrated path work, unique parameter correction,
affinity/preference prediction or production qualification. The original
reporting-resolution and missing-comparator limits remain. No continuation,
refit or model promotion follows automatically from these results.
