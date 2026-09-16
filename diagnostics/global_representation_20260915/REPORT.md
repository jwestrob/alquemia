> **Archived:** see [ARCHIVE.md](ARCHIVE.md) for final failure, cancellation and cost records. The queued status below is historical.

# Whole-system representation: investigation and queued technical retry

## Finding

A whole-system electronic model is worth testing. The existing APBS challenger is not suitable for scoring; its failure does not establish that protein-environment information is useless.

Exact decomposition of saved Ca-minus-La direct interactions (kcal/mol):

| Term | qm33 | qm36 |
|---|---:|---:|
| Asp303 represented in the permanent environment | 61.890346 | -0.342428 |
| All other environment residues | -3.270289 | -0.102747 |
| Total direct interaction | 58.620058 | -0.445175 |
| Total-system reaction-field contrast | -109.255766 | -111.287647 |
| Isolated-core reaction-field contrast, subtracted | -67.999236 | -137.126374 |
| Net environmental correction | 17.363527 | 25.393552 |

Moving Asp303 changes the direct term by -59.065232 and the subtracted reference term by +69.127138; the target reaction-field change is -2.031880. They leave +8.030025 kcal/mol in the correction. With the +2.204093 baseline shift, the final partition jump is +10.234117. The other environment residues account for +3.167543 of the direct-term change. This decomposition is exact accounting, not a controlled causal experiment: charge distributions, atom ownership and reference cavities change together. Asp303 did not physically disappear; most of its side chain moved into QM, while its retained backbone charges were adjusted by the declared boundary scheme.

The refined-grid partition jump remains about 10.03 kcal/mol. Simply refining the grid is unlikely to repair this inconsistency. The all-zero APBS identity failure is a separate unavailable result.

## Concrete alternative

ORCA 6.1 documents native GFN2-xTB and native ALPB. This avoids needing an external otool_xtb installation. Jacob approved the concrete pilot by asking to continue. Job 1198999 was cancelled during startup as its per-rank memory footprint approached node capacity; no energies were produced. Slurm failed to clear all processes and drained the node. The same pair is now queued as job 1199004 with memory-aware MPI sizing. The installed ORCA binary is pinned in the manifest; actual parameter export supports both metals and gives 23,259 orbitals and 25,764 active electrons for each endpoint. Convergence and energies are still unavailable. [ORCA native xTB manual](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/semiempirical.html#native-gfn-xtb-and-gfn2-xtb)

The proposed first model treats every physical atom with the same Hamiltonian. It has no capped QM region, no core/MM charge transfer and no isolated-core solvation subtraction:

`R_global = E_native-GFN2/ALPB(Ca, full system) - E_native-GFN2/ALPB(La, full system)`

Use the complete reported endpoint energy, with all model components and electronic smearing convention retained. There is no inherited r2SCAN aquo offset or classification band. A raw contrast on one development protein cannot demonstrate discrimination. The global charge distribution must be inspected: convergence alone would not establish realistic metal-localized chemistry or exclude spurious charge transfer through the protein.

The [GFN2-xTB primary method paper](https://doi.org/10.1021/acs.jctc.8b01176) describes broad element parameterization and a minimal electronic model. Availability of element parameters is not La/Ca selectivity validation. Actual native parameter export is requested in both inputs; no basis/ECP substitutions or extra dispersion corrections are specified.

ORCA QM/QM2 is an alternative for retaining local high-level accuracy, but it still has a high/low partition and requires a new boundary test. It is not a partition-free global discriminator. It has not been selected or implemented here. [ORCA multiscale manual](https://www.faccts.de/docs/orca/6.1/manual/contents/multiscalesimulations/qmmm-molecules.html)

## Prepared artifact

New preparation protocol: `whole_chain_native_gfn2_alpb_water_vertical_proposed_v1`.

- Real 1H4I deposited catalytic chain A plus the existing oxidized PQQ3- and selected metal; not the biological oligomer. Other deposited chains are explicitly excluded, matching the completed APBS physical system.
- 9,141 atoms; 9,274 source/schema covalent bonds; no synthetic caps or cut bonds. Protein/PQQ/metal are separate connected components as chemically expected. No metal coordination bonds were invented.
- Paired XYZ coordinates identical; source protein coordinates preserved within 1e-9 A. Both old partitions have the same physical-atom hash.
- Protein formal charge -8, PQQ -3, La +3 / Ca +2: total -8 / -9. Both all-electron counts have even parity. The backend's valence-electron convention remains to be checked in actual output.
- Full PQQ naming and bonds mapped to the existing source and versioned PQQ schema, including the existing three hydrogens. Existing source waters remain absent.
- FF19SB is used only to audit original topology/protonation and integer charge. Its partial charges are not inputs to the global electronic model.
- Coordinates, source mappings, settings, implementation/dependency hashes, executable hash and existing-runner policy are recorded. Native method parameters will be exported if run.

Seven real-artifact preparation/algebra/parser/hardware-layout tests pass. These include reconstruction of saved direct terms, paired coordinates/charges, restored Asp303 CA-CB and peptide connectivity, PQQ completeness, rejection of an explicitly corrupted real pair and acceptance by the existing ORCA runner. Native parameter export and MPI startup are now observed; no converged xTB energy or successful scientific result is claimed.

## Decision

The global feasibility pilot is approved; a memory-aware technical retry is queued. Baseline remains unchanged. Numerical credibility, predictive usefulness and affordability of this new global route are all unestablished; runtime at this size must be measured. Whole-protein diagonalization may be the practical bottleneck. No compute or wall-time budget is imposed.

See OPERATIONS.md for exact preparation, test and execution commands. RESULT.json retains artifact hashes and preparation timings. Initial development-only preparation files without runner task IDs are superseded by `prepared_v1/`; no scientific outputs were overwritten.

## Execution failure and recovery

See TECHNICAL_RETRY.md. The initial MPI choice was an execution error: active workers did not establish useful scaling, and native startup replicated large per-rank arrays. Peak batch RSS was 8,053,902,056 KiB. The top-level allocation reports 866 seconds and 297,904 allocated core-seconds; orphan/cleanup activity after cancellation is not fully captured by that figure. No global energies or charge distribution were obtained. This has not tested the scientific discrimination of GFN2-xTB. Current job: 1199004, using retry_memory_v1. It sizes MPI from physical RAM, with a 16-rank-per-endpoint maximum and a 64-GiB-per-rank planning allowance based on measured startup memory. Native parameter export confirms the intended elements/valence accounting; full numerical feasibility and affordability remain unresolved.
