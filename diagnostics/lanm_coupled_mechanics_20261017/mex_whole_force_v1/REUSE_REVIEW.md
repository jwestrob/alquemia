# Whole-Mex directional-check reuse review

19 October 2026. Source/hash/metadata review only; no scientific imports, numerical derivative calculation, model initialization or submission. **No metadata reuse blocker found.** Actual arithmetic/native replay and displaced-coordinate reconstruction remain mandatory in Slurm preflight.

## Confirmed reuse

- Four saved whole2049 BEFORE MACE centers are complete in collector14794. The two spin hypotheses for each metal share exactly one solvent receipt: four MACE centers and two ALPB centers.
- Current La and Dy BEFORE physical-coordinate lists and all2049 IDs are exactly equal. The target is physical index2045; spectator Nd IDs occupy2046–2048. This order is the current PANEL order, not the older classical-origin order with target last.
- All16 old La displacement receipts exist at the required coordinates. Their POINT/request/physical/map/XYZ/charge/native-execution links agree. Their fixed-charge file SHA is identical to the current La and Dy center charge file; binary command and fixed-charge rule match the pinned PANEL solvent parent.
- The old12952 Dy center has different coordinates. Its16 displacement receipts are therefore not reusable for this current paired BEFORE geometry, irrespective of sharing atom IDs.
- No previous whole2049 directional qualification is established by14794. Old12952 tests the regional MACE/GFN/MM Hamiltonian; its solvent and geometry artifacts can be reused separately.

## Parent pins

- **collection14794**: `workspaces/lanm_coupled_mechanics_20261017/mex_whole_collection_v1/execution/RESULT.json`; SHA256 `4c0c8b7379a2210f64267f0df94b1df72b54902d7233ec4b9614d4e941349d27`.
- **collection configuration**: `workspaces/lanm_coupled_mechanics_20261017/mex_whole_collection_v1/CONFIG.json`; SHA256 `1536e02656d4887c8828c11c72f31d696b8fd72b4302337195968f0d01c4d6d7`.
- **collection manifest**: `workspaces/lanm_coupled_mechanics_20261017/mex_whole_collection_v1/execution/MANIFEST.json`; SHA256 `acdaa4623a85a6e36fa7edc23d43a94fd94bc51ff9465aea69834a41f1b8d8d0`.
- **PANEL implementation**: `workspaces/lanm_coupled_mechanics_20261017/mex_whole_collection_v1/PANEL.py`; SHA256 `dd611ca03982e4f2548a88dabf9477e33aba7592fc2b41d832e60d13edd2add7`.
- **qualification12952**: `workspaces/lanm_coupled_mechanics_20261017/mex_qualification_v1/execution/RESULT.json`; SHA256 `15fb2a83b3f57e84ed7fcda0f61094dbad241ec265edd42fc6cce88dc4e55987`.
- **FD directions**: `workspaces/lanm_coupled_mechanics_20261017/mex_qualification_v1/execution/FD_DIRECTIONS.json`; SHA256 `8972af45b94e6e1eefdc916030b3fa08fc1fc79a1ae581340a95da2f01c9ee47`.
- **qualification source**: `workspaces/lanm_coupled_mechanics_20261017/mex_preflight_v1/MEX_CARTESIAN.py`; SHA256 `5153440c685dc2c1c335c35bdae42b19ce93f15de38d4706bd5b616365a32a8c`.
- **qualification configuration**: `workspaces/lanm_coupled_mechanics_20261017/mex_preflight_v1/CONFIG.json`; SHA256 `7abb3d5267c7ef260ecf4fc6bd86ffbbdd95de3c9b750f6999d86e7fc8049653`.
- **original source INPUTS**: `workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/Mex8FNS/INPUTS.json`; SHA256 `4b3b98983af0cb5fcdbf7b8f52e35063308201ce9df777a30484acf5146b7930`.
- **native solvent parent**: `workspaces/lanm_coupled_response_20261016/mex_static_v2/execution/RESULT.json`; SHA256 `ff7b27d490de058d6fbc494ff524a144f533e63e754f6dd3282d58c8fc904d0e`.
- **native solvent binary**: `workspaces/lanm_coupled_response_20261016/mex_static_v2/execution/mex_fixed_solvent`; SHA256 `b2bd6f0add6312836d19d02d77ff7051115facbdb9b0732bc7d347763cd98fb6`.

## Central native receipts

| State | Multiplicity | Raw MACE SHA256 | Shared ALPB SHA256 |
| --- | ---: | --- | --- |
| before_La_aligned | 10 | `617e61035b2de1a475bacb4ca136883415ca3ff3e6f35520bd3165ac2608f8c7` | `7c812b2b438aead93124dcfe0d181f8b493c2972e541b5d32262132de7b65a6b` |
| before_Dy_aligned | 15 | `ec5c472039dc88aaa05ed4cdc01c99faf76fc10e4b78ed82db6997f56f5a239b` | `bcb67b34dfa06b9d1102352ffbc27aa0e580264255ae0ed61127f5e94726830b` |
| before_La_one_Nd_reversed | 4 | `d2e6825d8d61e455f0f955693526690a4f11e517593520f67a82b79f9deffb7f` | `7c812b2b438aead93124dcfe0d181f8b493c2972e541b5d32262132de7b65a6b` |
| before_Dy_one_Nd_reversed | 9 | `34de09432a3c376b052a9122635784e0f336171313465270d98860e6b37386ed` | `bcb67b34dfa06b9d1102352ffbc27aa0e580264255ae0ed61127f5e94726830b` |

Raw paths are linked by the exact collection rows. Current before-La POINT is `f6380efc7961da114f62b0fe3073df6884af42eb6592f1e5cb9d1ede90bc8cef`; old Dy start is `944c1db18ba6c371a2b2f9c3ce34e20cce47addeecde7cb39bb231e538c2ab77`. Current before-Dy uses the paired La-coordinate POINT `b93554bf01ad050b800b894194b84b05582c46eadd0cb540eed0a92afb1cf142`.

## Displacement-to-request ordering

Base path: `workspaces/lanm_coupled_mechanics_20261017/mex_qualification_v1/execution/`. The frozen qualification source loops directions in saved insertion order, then h=0.0002/0.0001, then minus/plus. `CARTESIAN_REQUEST_0000/0001` replay the cached La center;0002 is fresh center `point_001_La`. Displacements use requests0003–0018 and POINTs002–017. Do not equate request and point indices.

| La POINT index | Request index | Direction | h (Å) | Sign | POINT SHA256 | ALPB SHA256 |
| ---: | ---: | --- | ---: | --- | --- | --- |
| 002 | 0003 | core_backbone | 0.0002 | minus | `4996e3a6278f96063aeb88da0cd3717c9584bf74c9dfe4433f6f1a2bc33d2ca5` | `eab2bbb4737e369221e49dc0fd93c4eb3c8239eb319a8d1fe796747d678b2f65` |
| 003 | 0004 | core_backbone | 0.0002 | plus | `c0ce9908b9a6d817bdc77f106e78682692d286d7ab896933e9e07ce108f2bc23` | `84cf413f25e28eaf2b6493dbfe32a385b96df70243e20ba6e029e8278d3d66d7` |
| 004 | 0005 | core_backbone | 0.0001 | minus | `dcf4de954b4a55cf33e9e98289384b2f1fb23a9657c8863860c9b29e5b60a5d2` | `f749be4204d8755479d0c19867f74ee7df883e735ee5372d4e7f5c5964c33a02` |
| 005 | 0006 | core_backbone | 0.0001 | plus | `999765c914b0f1570865905e770a3a8764a0c1f9a705dc50d9a0c0a9e116aa2d` | `a31fd64e24cec5ec76fa72d63f8ed7b7c6924ad5846b863f51166c9b1f6fae27` |
| 006 | 0007 | boundary_anchor | 0.0002 | minus | `d574ce419d965daad24ebedc244e038721bfd127d24f9e93ed0270e3b87443a4` | `a33a9aad3c755a4d034772e4208eab11bbc77241203356c5f60cac3e294f5c01` |
| 007 | 0008 | boundary_anchor | 0.0002 | plus | `12f32ca1fa4ab09f7d19589b0b51695415697282a6bc0a8a8ab074c3a6406d12` | `d5045d79c0fd654c2d8046034ce5fa3027fb4ceb3812af999c23370697a97e77` |
| 008 | 0009 | boundary_anchor | 0.0001 | minus | `d5ddc95e885ded4899d35f3b168fbd08d20ed0835243e81509d78361ccd8502c` | `0e9c31c2a58d4f2ca588d389ec6fbf30185a15fcdf44cfca168697057ab9f32e` |
| 009 | 0010 | boundary_anchor | 0.0001 | plus | `1cc24097db5aca16acd862778c0475ef17b832007c16148bfd02d78b1996521d` | `4e59f6c2e110f8580e660ca973b4a2bd8c940d5929ebeb6650d5a0b31bafc141` |
| 010 | 0011 | external_water | 0.0002 | minus | `a3b031fc06e220d72846beec89ad5d8d36babecdd4504521959e3d83c65d0656` | `20e4a34155d2ef6c2fee29a15c5726ffca49ff625a5776b7f4fbed0158ee0e3e` |
| 011 | 0012 | external_water | 0.0002 | plus | `5bebd98eaa27492e3487b88ef2210541292e4c0bf03a83944b1ba5cfa7558fd9` | `678bd71efdd5e17bf53d435feda2585ce5dea8011a989098412c047c80e863a7` |
| 012 | 0013 | external_water | 0.0001 | minus | `68e956e9a533e4d72e16eb5c4a6f16f004aa956536f06919babb6ea737a46812` | `503f9eb811f2352c030b49dda844233073a6947f3d8189cd0337a67ba122c838` |
| 013 | 0014 | external_water | 0.0001 | plus | `6ff3e42db60c95232f5a916214277b8ad5caa1830c33961a3e28f8981a7d2358` | `c66635b8baf2636359e1bacfde7b95684e96f7af079fcceb2a9a62a25aa2d38d` |
| 014 | 0015 | spectator | 0.0002 | minus | `371a447edbd8a36e4cb5d9e57bfa9a26d5d1f4c2e4d89c6ebadac7715ac1f1f2` | `1827cd1a8cbe72dc59a079b08921dc344e76f63ecaf585b0f1b9b20bfe7086dd` |
| 015 | 0016 | spectator | 0.0002 | plus | `27f0026d99e453b530cdee464fbe601026957acbba43767a1bd04dc5e6f039ef` | `05d738a62dca6d6ca2e2a1a777c303d967e441c8069ff3794c9847f2cdc3a23e` |
| 016 | 0017 | spectator | 0.0001 | minus | `7b609b551ebf8467db2658f7b20dd5b019fe3fd5139b57023d5e74062a1c785e` | `92163a8b2e8447a482d7ddc0c2fa7eaba9c8341c31bfc95cecfe362945823f1a` |
| 017 | 0018 | spectator | 0.0001 | plus | `bc0ad0f7842d8f16eb9b82abc4cd52d8ff55538251c82b6f776e5864dbe92008` | `e55ef4f453fc1c51a10d7087136c9ff2a0e402908063bf8130302ee8edfdf15d` |

Each POINT links `cells/La_proposal_La/solvent/RESULT.json`; each request is `CARTESIAN_REQUEST_NNNN.json`. The old Dy requests begin with0019/0020 cache replays,0021 fresh center, then0022–0037 displacements; none supplies the new Dy solvent panel.

## Implementation traps and remaining runtime gates

- Reuse the **La** vectors from pinned FD_DIRECTIONS, not regenerated Dy-specific directions. Original moving IDs are A/84//C, A/83//CA, the complete A/301 water (O/H1/H2), and A/201//ND. The original water vector is normalized in full Cartesian space; h is its full path length.
- Match each saved displacement to original center ±h×vector and verify direction normalization in Slurm. This review establishes archived request/physical equality and frozen loop ordering; it does not compute that geometry formula on login.
- PANEL.admit_solvent requires exact START.physical_mapping and point.mapping identity, raw XYZ symbols/coordinates, fixed charges, native component/gradient replay and assembled-gradient linkage. A new Dy solvent cell needs its own exact Dy physical/mapping/XYZ evidence at the shared La displacement coordinates; no La-symbol solvent reuse for Dy.
- Preserve producer conversion order `raw_gradient * HA / BOHR`. MACE native output is force in eV/Å; its analytic gradient uses the original negative-force conversion. Reassociation of conversions has caused earlier exact-replay failures and is unnecessary.
- The saved request state describes the old local charge0, La1/Dy6 model. It proves provenance of geometry/solvent only. The new native whole model must use charge+2 and aggregate multiplicities10/15 or4/9, with no local-state substitution.
- Keep64 new displaced MACE and16 new Dy ALPB slots, with the16La ALPB reused and no central reruns. Share solvent only across spin hypotheses of the same metal and exact coordinates. Retain failed/unattempted slots; no repeat or fallback.

Metadata checks verified 187 explicitly touched unique pins. Scientific full-gradient/component closure and directional pass/fail remain unexecuted by this review.
