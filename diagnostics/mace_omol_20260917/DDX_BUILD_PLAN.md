# Isolated ddX0.9.0 build

Engineering prerequisite only: build/import the maintained continuum solver
without altering either production or MACE environments. No molecule, energy,
force, density integration or biological comparison is executed in this job.
The active autonomous goal authorizes this contained dependency preparation.

Use the verified0.9.0 PyPI source archive and a fresh workspace venv based on
the existing Python3.11.15 executable. Pin NumPy1.26.4, SciPy1.17.1,
pybind11 3.0.1, setuptools82.0.1, wheel0.45.1 and packaging25.0 from downloaded
binary wheels. Installation is offline inside the allocation. Preserve hashes
of every wheel, source archive, executable, plan and build script; retain build
logs, resource receipts, output wheel, linked-library listing and pip freeze.

Request64CPU/32GiB on standard/memory, use up to64parallel build processes,
single-threaded BLAS, native GNU compilers and Release compilation. Do not add
fast-math or native-architecture flags. Disable upstream example/test builds;
no synthetic molecular calculation is presented as scientific evidence.
Smoke validation is import/version/banner only. Record failures as failures;
any corrected build gets a separate attempt directory. No project time limit
or compute budget. A solver pilot requires its own declared physical setup.
