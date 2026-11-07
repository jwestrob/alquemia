"""Recover bibliography identifiers for only the two predeclared ligand pairs.

NIST native records stay in the ignored workspace. This does not infer a
supporting electrolyte or assign a particular compilation value to a paper.
"""
import os, re, struct
from util import *

assert os.environ.get('SLURM_JOB_ID')
native=ROOT/'workspaces/ln_selectivity_calibration_20261027/nist46_dat'
code=(native/'F8547_CODE1.DAT').read_bytes()
values=(native/'F8558_VALUES.DAT').read_bytes()
citations={}
for line in (native/'F8544_BIB.DAT').read_bytes().decode('latin1').splitlines():
 if '\t' in line:
  key,text=line.split('\t',1);citations[key]=text
links=(native/'F8555_LMBIB.DAT').read_bytes().decode('latin1').splitlines()
rows=[]
for offset in range(0,len(code),256):
 row=code[offset:offset+256]
 if len(row)!=256:continue
 match=re.match(rb'(\d+)\t',row[4:])
 if not match:continue
 ptr=struct.unpack('<f',row[:4])[0]
 assert ptr.is_integer() and ptr>=2
 value_offset=(int(ptr)-1)*15
 raw_id=int.from_bytes(values[value_offset:value_offset+2],'little')
 if raw_id not in (576,597,599,600):continue
 logical_id=int(match[1])
 found=[]
 for line in links:
  fields=line.split()
  if len(fields)!=4 or fields[0]!=str(logical_id) or fields[1] not in ('45','75','59'):continue
  key=fields[2].rstrip('-')
  assert key in citations
  found.append(dict(metal_code=int(fields[1]),bibliography_code=key,
      stored_marker='-' if fields[2].endswith('-') else '',
      bibliography_index=int(fields[3]),citation=citations[key]))
 assert found
 rows.append(dict(native_ligand_id=raw_id,bibliography_ligand_id=logical_id,
    code1_byte_offset=offset,first_value_byte_offset=value_offset,
    identity=row[4:].decode('latin1').strip(),bibliography=found))
assert sorted(r['native_ligand_id'] for r in rows)==[576,597,599,600]
cis=next(r for r in rows if r['native_ligand_id']==599)
assert cis['bibliography_ligand_id']==440
assert {b['bibliography_code'] for b in cis['bibliography'] if b['metal_code'] in (45,75)}=={'87CMF'}
trans=next(r for r in rows if r['native_ligand_id']==600)
assert trans['bibliography_ligand_id']==441
report=Path(os.environ.get('TAURI_REPORT_DIR',str(D)))
out=report if report!=D else W/'provenance'
out.mkdir(exist_ok=True,parents=True)
save(out/'NATIVE_BIBLIOGRAPHY.json',dict(rows=rows,job=os.environ['SLURM_JOB_ID'],
 sources=[pin(native/f) for f in ['F8547_CODE1.DAT','F8558_VALUES.DAT','F8544_BIB.DAT','F8555_LMBIB.DAT']],
 note='Native value identifiers differ from the logical identifiers used by the bibliography. All mappings above are checked through the first-record pointer. No new values or sources are guessed.'))
save(report/'PROVENANCE_STATUS.json',dict(
 native_cis_id=599,bibliography_cis_id=440,native_trans_id=600,bibliography_trans_id=441,
 cis_La_Dy_source=dict(code='87CMF',citation='J. Charlier, E. Merciny and J. Fuger, Analytica Chimica Acta 192 (1987) 95–101',doi='10.1016/S0003-2670(00)85691-3'),
 trans_source_candidates=sorted({b['bibliography_code'] for b in trans['bibliography'] if b['metal_code'] in (45,75)}),
 trans_1984_source=dict(code='84MF',doi='10.1016/S0003-2670(00)84511-0'),
 primary_publisher_metadata_retrieved=True,primary_methods_tables_retrieved=False,
 exact_common_electrolyte_verified=False,experimental_pair_admitted=False,
 remaining_evidence='Original methods and tables for the cis study, and the exact source of the trans 1.0 M compilation row. A K+ background annotation alone does not identify its anion.',
 private_audit=pin(out/'NATIVE_BIBLIOGRAPHY.json')))
print('Bibliography mapping verified for the four predeclared ligands; cis source recovered. Experimental medium remains unverified.')
