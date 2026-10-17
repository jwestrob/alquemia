# Saved preparation environment recovery

The frozen response adapter is unchanged at bd29908d. Attempt13059 used the native-reference Python environment and failed an upstream exact terminal physical/force replay before molecular calls. Attempt13060 restored the original pool environment and failed a missing Gemmi import before replay. Both failures and allocated costs are preserved separately.

The interpreter binary is identical in both environments, but the installed NumPy versions differ: native1.26.4 versus original pool2.4.6. Static comparison found five saved guard-distance leaves differing in the last printed digits. The precise terminal assertion difference is not yet measured; no scientific tolerance is changed.

Preparation13063 uses the original pool interpreter/dependencies on the original qualified host/12CPU configuration, without requesting a GPU. A workspace-local mini-prefix contains exact copies of the existing Gemmi0.7.5 package and its non-system libz-ng/libstdc++/libgcc dependencies, preserving the extension RPATH layout. No original environment was modified. The per-file overlay hashes and actual successful import/version probe are recorded in the workspace. All molecular calls remain forbidden; future DFT execution still needs full derivative and measured-cost gates.
