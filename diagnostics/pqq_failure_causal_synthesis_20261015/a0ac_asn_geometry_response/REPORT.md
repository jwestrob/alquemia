# Actual donor geometry in the fixed Asn flip

All10 saved sources pass paired Ca/La coordinate, mapped H-parent and input-identity checks. Slurm4626 performs no molecular evaluation.

For the ambiguous Ca-conditioned sample4, the Asn290 oxygen–metal distance changes4.595274311725036→2.562017935679298Å; nitrogen–metal changes2.4460701134677234→4.718439745320599Å. Thus the favorable La vacuum candidate actually presents the oxygen toward the metal. This is consistent with a donor-orientation explanation for part of its local response, without isolating the energetic contribution of one atom.

The correctly called La-conditioned4 undergoes the same qualitative change: oxygen4.557746811748103→2.5321207859434307Å, nitrogen2.413891671140194→4.683172457552463Å. Its original adaptive candidate remains lower in the completed vacuum pool. The other eight sources begin with nearby oxygen; the fixed flip moves it away and their old geometries remain preferred.

The target also introduces an included, nonbonded Asn290HD22–Gly292backboneH contact of1.3821908645320085Å from4.943529710641982Å. La-conditioned4 has the analogous1.3744771216817757Å contact. Both are explicitly represented in the energy model. The nearest-contact field used earlier included covalent bonds; this dedicated analysis excludes directly mapped covalent bonds and distinguishes the actual nonbonded contacts. It does not change any geometry, gate or score. A new low-energy candidate is not thereby a relaxed protein conformation.

[RESULT.json](RESULT.json) gives all ten source distances and named nonbonded contacts; the workspace also contains a joinable TABLE.csv and full pinned result. Actual job4626:1schedulersecond/1CPU; processwall0.2833237440790981s,processCPU0.199004331s,peakRSS41984KiB. Zero new energy calls.
