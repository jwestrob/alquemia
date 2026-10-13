# Local saved-coordinate reader recovery

Job3340 completed all ten raw geometry/confidence/alignment records. Its local reader asserted equality between the 192 scored atoms and 189 physical-map atoms; the map also includes boundary references and represents cap atoms separately through `core_links`. All nine mapped sources stopped at that assertion before local metrics. Their partial statuses remain immutable.

The local-only successor reuses every completed raw/confidence/alignment field. For the same nineteen saved local geometries it replays each archived full-q with the existing `Kinematics` helper, then takes each linked physical atom directly from the actual saved XYZ via `core_links`. Unrepresented physical boundary anchors retain their exact kinematic positions. It records the complete saved/replayed XYZ residual without introducing a new scientific quality threshold. Origin uses the map's exact zero vector. No optimization, energy call, structure repair, new candidate, confidence recomputation or new fit occurs. The previously fixed measurement definitions and all ten source rows remain unchanged; excluded La1 stays raw-only.

One real Slurm CPU. Existing helper, raw result, maps, coordinates and all ten archived proposal receipts are pinned. Root submits; no scientific imports or structural calculations on login.
