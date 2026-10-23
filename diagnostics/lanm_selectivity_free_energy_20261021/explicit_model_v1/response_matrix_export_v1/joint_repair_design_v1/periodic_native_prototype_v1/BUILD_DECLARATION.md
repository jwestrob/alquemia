# One isolated compilation

Root authorized one compilation of manifest `c4feea15` after source review `41c11f68`. Use the previously successful `/usr/bin/f95`, RELWITHDEBINFO `-O2 -g -fbacktrace`, static FFTW libraries, and analyze/testgrad targets, with one build worker and one numerical thread per process. Exact argv, pins and environment are in BUILD_CONTEXT.json. Only the build-driver worker count and redundant CMake-manifest lookup change; all frozen source bytes and both fatal periodic guards stay unchanged.

This changes the preparation scope from source-only to compilation-only. Preserve logs, failures, wall/CPU cost and output hashes. No automatic retry, native molecular execution, candidate composition, box evaluation or sampling is released. The existing 57-program vacuum panel may be prepared while compilation runs; its execution awaits root review/release.
