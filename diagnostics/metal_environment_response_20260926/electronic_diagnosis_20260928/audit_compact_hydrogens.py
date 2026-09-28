"""Read-only real/source versus synthetic-cap hydrogen geometry audit."""
from pathlib import Path
import json,sys,xml.etree.ElementTree as ET,itertools
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,record,xyz,write_new
m=read_json(ROOT/'workspaces/lanm_series_followup_20260923/prepared_v1/manifest.json')
refs={('Hans8DQ2' if c['protein']=='Hans' else 'Mex8FNS'):c['repair_manifest'] for c in m['cases'] if c['site']=='EF3'}
refs['Hans8FNR']=record(ROOT/'workspaces/lanm_series_followup_20260923/dy_transfer_v1/sources/EF3/amide_v3/repair_manifest.json')
design=read_json(Path(__file__).parent/'COMPACT_FORCE_DESIGN.json')
def key(a):return a['chain'],int(a['resnum']),a.get('insertion_code',''),a['atom']
def angle(a,b,c):
 x=a-b;y=c-b;return float(np.degrees(np.arccos(np.clip(x@y/np.linalg.norm(x)/np.linalg.norm(y),-1,1))))
rows={}
for s,ref in refs.items():
 r=read_json(verify(ref));atoms=xyz(verify(design['sources'][s]['xyz']['La']));q=np.array([a[1:] for a in atoms]);mp=[a for a in r['atom_graph']['source_to_qm'] if 'source' in a];lookup={key(a['source']):a for a in mp};src={}
 for l in verify(r['source_structure']).read_text().splitlines():
  if l.startswith(('ATOM  ','HETATM')):src[(l[21],int(l[22:26]),l[26].strip(),l[12:16].strip())]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
 maxerr=max(float(np.max(abs(q[a['qm_index']]-np.array(a['xyz_A'])))) for a in mp)
 if maxerr>1e-8:raise ValueError('compact mapping coordinates differ')
 templates={a.attrib['name']:a for a in ET.parse(verify(r['topology_definition'])).getroot().find('Residues')}
 caps=r['atom_graph']['cut_bonds_and_caps'];byparent={key(c['retained']):c for c in caps};angles=[];lengths=[];capangles=[]
 for k,a in lookup.items():
  if a['source']['element']!='H':continue
  rn=a['source']['canonical_resname']
  if rn not in templates:continue
  bs=[(b.attrib['atomName1'],b.attrib['atomName2']) for b in templates[rn].findall('Bond')]
  pn=[y if x==k[3] else x for x,y in bs if k[3] in (x,y)]
  if len(pn)!=1:raise ValueError('ambiguous hydrogen parent')
  pk=(*k[:3],pn[0]);hk=k
  if pk not in lookup:continue
  pi=lookup[pk]['qm_index'];hi=a['qm_index'];lengths.append(dict(hydrogen=hk,parent=pk,length_A=float(np.linalg.norm(q[hi]-q[pi])),kind='real_source_H'))
  for x,y in bs:
   if pn[0] not in (x,y):continue
   other=y if x==pn[0] else x;ok=(*k[:3],other)
   if ok==hk or ok not in src or other.startswith('H'):continue
   # Source heavy vector even when corresponding exterior atom replaced by cap.
   av=angle(src[ok],src[pk],src[hk]);angles.append(dict(hydrogen=hk,parent=pk,heavy_neighbor=ok,angle_degrees=av,neighbor_in_QM=ok in lookup,kind='real_source_H'))
   if pk in byparent and key(byparent[pk]['omitted'])==ok:
    c=byparent[pk];ci=c['qm_index'];capangles.append(dict(hydrogen=hk,parent=pk,omitted_heavy=ok,cap_index=ci,cap_parent_real_H_angle_degrees=angle(q[ci],q[pi],q[hi]),source_heavy_parent_real_H_angle_degrees=av,kind='synthetic_cap_substitutes_same_heavy_direction'))
 rows[s]=dict(repair_manifest=ref,source=r['source_structure'],topology=r['topology_definition'],xyz=design['sources'][s]['xyz']['La'],mapped_coordinate_max_error_A=maxerr,real_H_lengths=lengths,real_H_heavy_angles=angles,cap_angles=capangles)
write_new(Path(__file__).parent/'COMPACT_HYDROGEN_AUDIT.json',dict(rows=rows,implementation=record(__file__),molecular_calls=0))
for s,r in rows.items():
 print(s,'realH',len(r['real_H_lengths']),'length range',min(x['length_A'] for x in r['real_H_lengths']),max(x['length_A'] for x in r['real_H_lengths']))
 print('largest angles',sorted(r['real_H_heavy_angles'],key=lambda a:a['angle_degrees'],reverse=True)[:3])
