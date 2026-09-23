# What physically changed in the successful combined pilot?

This is a read-only interpretation of completed preparation/score records, not
an additional molecular experiment. It identifies the actual representation
difference behind the hypothesis without presuming that larger contexts are better.

For the remaining difficult A0A3F2YLY8 Ca-conditioned sample1, the old surrounding
context contained142 atoms, while the consistent fragment union contains193.
The added selected units are the Asn285 and Asp327 peptide units and the Cys129
and Trp289 sidechains (chainA); graph closure/caps are handled by the existing
source graph. All added fragment formal charges are zero. The endpoint total
charges remain Ca−3/La−2. The original core and water inventory are retained.
This therefore is not simply moving the score by adding a formal unit charge.

The other prescribed A0A3 fold, Ca-sample3, originally contained179 atoms. Its
only newly selected unit is the Asn285 peptide; the union also has193 atoms.
Thus the fixed-membership rule corrects a concrete inconsistency between these
source representations. It uses the unlabeled union over supported saved folds,
not a fragment selected after seeing this score.

For the A0ACD6B9F2 Ca4 and La4 sources, both old and union contexts contain180
atoms and there are no additional selected fragments or total-charge changes.
That is consistent with the combined method retaining their existing adaptive
response. These cases do not prove that enlarging all pockets helps.

The static A0A3 Ca1 union score changed by+8.740205 model kcal/mol in the earlier
full-union report: native OMOL contributes−5.731701, while the paired solvent
transfer contributes+14.471905. The latest adaptive pool additionally allows
physical donor motions in the fixed context. These coupled changes do not
identify a unique hydrogen bond, cavity term, missing force or mechanical cause.
Removing/adding fragments changes polar interactions, shape, boundaries and
available mapped motions together. Neutral fragments can still change the
metal-dependent solvent response.

The concrete hypothesis is that consistent local membership avoids fold-dependent
omissions, while accommodation addresses strained donor arrangements. The eight
source pilot supports pursuing that combination; full canonical and structural
transfer tests are still needed to establish practical utility. A future three-fold
scanner is not automatically equivalent to this ten-fold union policy.

## Actual source records

- `workspaces/consistent_context_20260922/prepared_v1/PREPARATION.json`
  SHA256 `8f92f9752bf9e978c1b96da041c261deb82d128f531f6f6f4e058e63a19fb8a2`;
  per-source old/new context preparation paths are recorded there.
- `workspaces/union_adaptive_20260923/proposals_v1/manifest.json`;
  contains actual union selection and physical mapping for the combined pilot.
- `diagnostics/consistent_context_20260922/REPORT.md` for completed static-union
  component shifts and the negative broader result of that component alone.
- `diagnostics/union_adaptive_20260923/REPORT.md` for the completed eight-source
  combined experiment, works, ordering, transferred bands and costs.

No atom deletion/addition trial was run for this note; it does not isolate a
causal contribution of any one newly included residue.
