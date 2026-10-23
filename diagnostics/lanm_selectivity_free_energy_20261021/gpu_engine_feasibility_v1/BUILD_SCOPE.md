# Isolated GPU build

Source inspection confirms POLPAIR in CUDA field/energy kernels and existing Ewald/non-Ewald tests. Build pinned Tinker-GPU44bcd7c with its required CPU submodule33dcaf08 in this workspace only, CUDA13.0, explicit H200 architecture90, double precision, four compile processes. No installation or login GPU calls. Actual GPU parity and periodic tests will be separately run on Slurm before using this engine for molecular sampling.
