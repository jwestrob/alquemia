# Parallelize the existing dense source evaluation without changing its arithmetic

The first real multipole bridge evaluations agree with the independent Cartesian
potential to about 1e-15 au, but each takes about 29 seconds. Source inspection
shows that stock ddX's `build_phi_dense` loop over independent cavity points is
serial, even when its model is configured for 64 threads. This is preparation
cost that must count in any future full score.

Build a separate ddX 0.9.0 environment from the pinned original archive and
dependency wheels. Add only an OpenMP parallel-do directive around that point
loop, with the point index, source index, local vector and scalar private.
Each point retains its original source summation order and native multipole
routine. No solver, cavity, source coefficient, physical parameter or arithmetic
expression changes. Preserve the exact patch and verify all other native
source files byte-for-byte. Do not replace the stock or PCM-adapter environments.

After the current bridge finishes, qualify the build on the same real 2FW0
primary Ca/La average/difference distributions: four distributions at 1, 8 and
64 configured threads, twelve native potential calls and twelve source-integral
calls. Reuse the four original stock outputs as references. No continuum solve,
energy, response iteration, DFT, MACE, force or biological score.

Require identical cavity arrays, source moments, coefficients and software
dependency versions. Maximum full-array potential difference <=1e-12 au and
psi difference <=1e-12. Report whether equality is also bitwise; do not require
a speedup to declare numerical agreement. Record per-call wall/CPU time and
memory, allocation costs, and costs of the separate build. The performance test
uses one 64-CPU allocation with sequential thread-count groups. A routine use
decision requires measured speedup and passing representation checks.

This is contained performance engineering under the active goal, not a new
physical model. It preserves the already running stock bridge and resolution
jobs and their complete costs. No existing working environment is modified.
