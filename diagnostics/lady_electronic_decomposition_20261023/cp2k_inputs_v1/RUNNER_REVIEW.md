# CP2K quartet runner source review

Reviewed 2026-10-23T21:12:25.539122+00:00.
**Verdict: no remaining source blocker for the declared contained run. Actual
configuration, scheduler wrapper and release binding remain pending root review.**
This review ran no CP2K, MPI, parser, molecular calculation or scheduler command.

## Exact source

- RUN.py: `03d9338b26e8bd507a9fd9bb4523ca6b8102e2597cebd7f582a67fc5f90bd0d7`.
- DECK_PROVENANCE.json: `916c45531a5384a60dd450d14831e27a6d2efc693eb8d43202a2942e6ada6687`.
- Basis: `2ab7190ac0bf52f389fed017672c49f747d75d6d0e92f61dbf4e2f5981b0eb66`.
- Potential: `ba87181563f5ca2b8e6e4ff30f43d274752b9c663325905ee6c1dc41c7beccd5`.
- Conversion provenance: `ff923ac1a23676e0cd799a989e17271bf30e10aefb74f7a503a50899b82553b4`.

## Runtime and preservation findings

1. Exactly one molecular parent launch appears on the executable path. There
   is no retry, fallback, optimizer or resubmission loop. Native syntax checking
   precedes it and failure prevents the molecular call. Installed
   `src/start/cp2k.F` lines 188–191 set `run_it` false for `--check`; this is
   syntax checking, not initialization or molecular qualification.
2. Runtime requires 24 Slurm CPUs, exactly 12 image-local MPI ranks and two
   OpenMP threads, with OpenBLAS/MKL each restricted to one. MPI launches locally
   through the image's MPICH fork/user path; no additional allocation, scheduler
   shell or remote-host list is invoked. Read-only source/build/assets mounts
   and a separate writable `/task` preserve the installed source and inputs.
3. The declared 7200-second molecular parent cap is exact. A timeout sends
   TERM, then KILL if needed, only to the new child process group. There is no
   `scancel`, unrelated PID targeting or queue mutation. Scheduler cgroup cleanup
   remains the outer containment if the Python controller itself is killed.
   Actual process-tree cleanup and performance are not established by source review.
4. Configured geometry-to-task identity is checked. Existing execution or run
   receipt prevents restart; an existing diagnostic geometry destination prevents
   overwrite. Exclusive log creation gives another collision guard. Parent-call
   intent is saved before launch; this counter is an attempted-call denominator,
   not proof that a native SCF started or converged.
5. Native stdout/stderr and `endpoint.out` are retained even on failure. Atomic
   run receipts preserve exceptions, return codes, timeout and wall cost. Original
   files remain in the execution directory; selected `.inp/.out/.log` files and
   JSON receipts are mirrored. Other native files are not deleted. The elapsed
   wall time times 24 is clearly labelled allocated capacity, with actual Slurm
   accounting deferred. Pre-configuration/path/Slurm-environment failures occur
   before the managed receipt block and therefore rely on the batch log; they
   perform no native call. A hard controller kill can leave the last running
   receipt, so actual terminal Slurm state must still be consulted.
6. Non-pilot entries require an admitted pilot record joined to this exact
   configuration hash. A successful runner status explicitly still requires
   scientific review; the runner does not itself admit physical accuracy.

## Deck/assets check

All four deck and geometry hashes match DECK_PROVENANCE. I independently
compared every coordinate token with the archived XYZ: atom order and values
are unchanged, atom 1 is the metal fragment, and every other atom remains in
one complete shell fragment. All four decks bind `/assets/BASIS_LADY` and
`/assets/POTENTIAL_LADY`, the correct La46/Dy55 entry, and exactly one metal
`ECP_SEMI_LOCAL TRUE`. Light atoms use ALL.

The single-token potential names are valid: installed
`src/qs_kind_types.F` lines 1868–1888 derive potential type ECP from the prefix
before the first hyphen while retaining the entire name for library lookup.
The two-token form is optional. Original effective-singlet/electron-count and
H/C/O documented-library-equivalence limitations remain as recorded in
CONVERSION.md. This review does not duplicate the separate ALMO numerical,
energy-ledger or fragment-convergence review.

## Required pending binding

Root must freeze and inspect a configuration whose pin list covers this exact
runner, parser, four decks, both assets/provenance, the installed binary and
runtime image. Its fixed task map, pilot admission path and fresh diagnostic
root must agree with the reviewed wrapper. The wrapper must explicitly request
`--qos=standard_expand`, one node, 24 CPUs and no requeue, retaining scheduler
memory/time defaults. These are pending configuration/release joins, not
claims that this source review submitted or admitted a run.
