"""Verify only the two predeclared development pairs in the native compilation."""
import os,ast,re
from common import *
assert os.environ.get('SLURM_JOB_ID')
curator=ROOT/'diagnostics/ln_structure_increment_20261030/curate.py'
tree=ast.parse(curator.read_text());node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='raw_records')
namespace={'OLD':ROOT/'workspaces/ln_selectivity_calibration_20261027','re':re}
exec(compile(ast.Module(body=[node],type_ignores=[]),str(curator),'exec'),namespace)
native={r['byte_offset']:r for r in namespace['raw_records']()}
audit=read(W/'EXPERIMENTAL_PAIR_AUDIT.json');verified=[];summaries=[]
for pair in audit['pairs']:
 rows=[]
 for record in pair['records']:
  for metal in ('La','Dy'):
   raw=native[int(record[metal+'_offset'])]
   assert raw['lig_id']==record['lig_id'] and raw['metal']==metal
   assert float(raw['value_raw'])==float(record['logK_'+metal])
   assert raw['I']==float(record['I_M']) and raw['T_C']==25
   assert raw['equilibrium']=='ML=M.L' and raw['stoichiometry']=='1:1'
   assert str(raw['annotation_code'])==record['annotation_code']
   rows.append(raw)
 verified.append(dict(name=pair['name'],records=rows,all_native_checks_pass=True))
 summaries.append(dict(name=pair['name'],ligand_ids=[r['lig_id'] for r in pair['records']],
   temperature_K=298.15,ionic_strength_M=float(pair['records'][0]['I_M']),
   conditional_development_contrast_kcal_mol=pair['D_kcal_conditional'],
   native_K1_1_to_1_LaIII_DyIII_conditions_verified=True,
   ligand_species_label='H4L; constants refer to L in ML=M.L',
   annotations=[r['native_annotation'] for r in pair['records']],
   status=pair['status'],reason=pair['reason'],
   original_study_provenance='Not decoded beyond NIST SRD46 v8 native compilation offsets.',
   missing_evidence=['exact electrolyte anion and shared medium','original studies/adjustments provenance']))
dat=namespace['OLD']/'nist46_dat'
sources=[pin(dat/f) for f in ['F8558_VALUES.DAT','F8548_CODE2.DAT','F8549_CODE3.DAT','F8550_CODE4.DAT','F8552_CODE6.DAT']]
verification_path=(W if RESULTS_DIR==D else RESULTS_DIR)/'EXPERIMENTAL_NATIVE_VERIFICATION.json'
save(verification_path,dict(pairs=verified,curator=pin(curator),source_files=sources,job=os.environ['SLURM_JOB_ID']))
save(RESULTS_DIR/'EXPERIMENTAL_ADMISSION.json',dict(pairs=summaries,admitted_pair=None,
 no_additional_pair_search=True,full_native_verification=pin(verification_path),
 source_audit=pin(W/'EXPERIMENTAL_PAIR_AUDIT.json'),native_records_remain_ignored_per_repository_redistribution_policy=True))
print('Eight native records verified; neither proposed common medium admitted.')
