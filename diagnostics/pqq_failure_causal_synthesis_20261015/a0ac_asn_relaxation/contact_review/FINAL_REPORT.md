# Relaxed alternative amide: named contact relief in both Ca4 and La4

All **10 sources, 20 searches and complete common6 pools** are available. Saved-contact job4910 reports all40 origin/fixed-flip/relaxed-Ca/relaxed-La geometries available, with exact paired Ca/La coordinates. No new energy calls occurred in that analysis. The Ca4 target and La4 control both retain the oxygen-facing alternative amide while relieving the specific fixed-flip HD22–Gly292 backbone-H contact. This physical change is shared by both sources; it is not uniquely associated with the problematic prediction.

## Exact named geometry

Distances in Å; Ca4/La4 denote source conditioning and sample4. “Relaxed-Ca/La” denotes the generating optimizer metal, not different distances between exchanged-metal copies of one geometry.

|Source / geometry|Metal–OD1|Metal–ND2|HD21–Gly292 H|HD22–Gly292 H|
|---|---:|---:|---:|---:|
|Ca4 origin|4.59527|2.44607|4.26956|4.94353|
|Ca4 fixed flip|2.56202|4.71844|2.24263|1.38219|
|Ca4 relaxed-Ca|2.32132|4.06994|2.81492|2.16938|
|Ca4 relaxed-La|2.39032|4.15387|2.68477|2.20915|
|La4 origin|4.55775|2.41389|4.07358|4.95671|
|La4 fixed flip|2.53212|4.68317|2.54703|1.37448|
|La4 relaxed-Ca|2.34949|4.06567|3.08587|2.26229|
|La4 relaxed-La|2.41159|4.15112|2.96583|2.27702|

The intact amide retains both physical H atoms attached to ND2 and its inherited bond lengths (approximately1.183–1.188Å here). No independent H/cap relaxation or atom relabeling occurred. Distances establish orientation and contact changes; they do not by themselves prove coordination or allocate an energy contribution.

**Contact relief remains incomplete.** The shortest included nonbonded amide H–H pair becomes HD21–Asp203HB3 at1.59074/1.61641Å for Ca4 relaxed-Ca/La. La4 retains1.78268Å HD21–Asp203HB2 and1.80620Å HD21–Asp203HB3. These are represented atoms in the fixed204-atom context. The inherited severe-contact guard (new heavy<1Å or H-involving<0.55Å) passes; it is not a general steric-quality or equilibrium criterion. Nearest-contact descriptors exclude direct source bonds and caps and query the four amide atoms, not every protein pair. Identity metadata coordinates remain the original source coordinates; reported `distance_A` values describe each saved frame.

## Availability, boundaries and solvent decision

All20 optimizers report success and physical admission. Ca4 displacements are0.70172/0.64598Å and La4 are0.68485/0.62980Å from their flipped starts; none has a boundary flag. Four other searches reach the inherited0.8Å bound: Ca2 both metals, La0's La search and La2's Ca search. Preserve those finite-domain limitations; no bound expansion is justified here.

Saved common4→6 vacuum changes are +16.527363234kcal/mol for Ca4 and +11.463790042 for La4, identical under mathematical/operational selection; the other eight are zero. The frozen sign-blind solvent-review gate passes. No new physical-gate failure was found, but residual contacts, shared source/state assumptions and constrained MACE-only relaxation remain material. Root's separate preflight/review and complete matched solvent stage are required before any composite interpretation; no class improvement, global minimum or causal H-contact energy claim follows.

## Evidence and actual cost

- [Search collection](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_remaining_v1/execution/FINAL_COLLECTION.json), SHA256 `2ed05f5a241f4f1eccbb2a680428a4b032a51449df56ad323e138b83cf32fd82`.
- [Contact result](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_contact_review_v1/execution/RESULT.json), SHA256 `33a74f84eb2c76bc56a85b2f06da6c57d59c060b6269aa8e9cd8e01b33d0657f`; [table](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_contact_review_v1/execution/TABLE.csv), SHA256 `0b115d235413722a958a197f8fe5712f7af520ec1ca268a474124c93d673fa27`.

Actual `sacct` at2026-10-15 04:12:07UTC: scout4831 used570s×7CPUs=3990 allocated CPU-s (30 model forwards); remaining4844 used4197s×7CPUs=29379 allocated CPU-s (236 forwards); contact4910 used2s×1CPU=2 allocated CPU-s, zero model/SCC/DFT calls. All completed0:0. Contact analysis measured0.69345s wall,0.38479s process CPU and44544KiB peak RSS. This report only reads those saved results/receipts and accounting.
