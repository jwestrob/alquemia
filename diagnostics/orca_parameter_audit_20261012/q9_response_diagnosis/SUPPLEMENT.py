"""Bounded saved-file policy and DFT inventory; no molecular execution."""
import csv, hashlib, json, math, re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
PINS={}
def pin(p):
    p=Path(p).resolve();r={'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()};PINS[str(p)]=r;return r
def load(p):
    if isinstance(p,dict):
        assert pin(p['path'])['sha256']==p['sha256'];p=p['path']
    else:pin(p)
    return json.loads(Path(p).read_text())
def write(name,d):
    with (OUT/name).open('x') as f:json.dump(d,f,indent=2);f.write('\n')
def pdb(p):
    if isinstance(p,dict):assert pin(p['path'])['sha256']==p['sha256'];p=p['path']
    d={}
    for l in Path(p).read_text().splitlines():
        if l[:6].strip() in ['ATOM','HETATM']:
            k=(l[21].strip(),int(l[22:26]),l[26].strip(),l[12:16].strip())
            assert k not in d;d[k]=tuple(float(l[a:b]) for a,b in [(30,38),(38,46),(46,54)])
    return d
def xyz(p):
    if isinstance(p,dict):assert pin(p['path'])['sha256']==p['sha256'];p=p['path']
    t=Path(p).read_text().splitlines();a=[x.split() for x in t[2:] if x.strip()];assert len(a)==int(t[0])
    return [(x[0],tuple(map(float,x[1:4]))) for x in a]

result=load(OUT/'RESULT_V2.json')
discovery=load(ROOT/'workspaces/motion_envelope_20260923/run_v1/DISCOVERY.json')
union=load(ROOT/'workspaces/motion_envelope_20260923/run_v1/selections/q9z4j7-pqq-la_model__envelope_5eba0c9b39ead1ee/UNION.json')
for script in ['second_shell_context.py','motion_envelope.py','consistent_context.py','affordable_peptide.py']:pin(ROOT/'scripts'/script)
targets=[('A',213,'','O'),('A',214,'','N'),('A',349,'','O'),('A',350,'','N'),('A',442,'','O'),('A',591,'','NE1')]
cases=[]
for c in discovery['cases']:
    if c['case_id'] not in union['members']:continue
    pp=load(c['original_preparation']);coords=pdb(pp['source'])
    anchors=[(a['chain'],a['resnum'],a['insertion_code'],a['atom']) for a in c['anchors']]
    rows=[]
    for t in targets:
        nearest=min(anchors,key=lambda a:math.dist(coords[t],coords[a]));d=math.dist(coords[t],coords[nearest])
        recorded=[x for x in c['contacts'] if (x['atom']['chain'],x['atom']['resnum'],x['atom']['insertion_code'],x['atom']['atom'])==t]
        assert bool(recorded)==(d<=4.3)
        if recorded:assert abs(recorded[0]['distance_A']-d)<1e-10
        rows.append({'target':list(t),'nearest_anchor':list(nearest),'distance_A':d,'recorded_selected_contact':bool(recorded)})
    cases.append({'case_id':c['case_id'],'source':pp['source'],'preparation':c['original_preparation'],'anchor_distance_replay':rows})
assert len(cases)==4
policy={'implementation':pin(__file__),'discovery':pin(ROOT/'workspaces/motion_envelope_20260923/run_v1/DISCOVERY.json'),
 'union':pin(ROOT/'workspaces/motion_envelope_20260923/run_v1/selections/q9z4j7-pqq-la_model__envelope_5eba0c9b39ead1ee/UNION.json'),
 'original_radius_A':3.5,'expanded_radius_A':4.3,'selector_elements':['N','O','S'],
 'anchors':'Fixed original-core nonwater N/O within 3.3 A of metal, plus same functional-group O/N partners.',
 'selection_geometry':'Four saved q9 La-conditioned samples 1-4; their fragment union transferred to canonical sample0. Not a swept motion envelope.',
 'omission_mechanisms':[
  {'groups':['Leu485 sidechain','Ala483 sidechain','Phe214 sidechain','Phe442 sidechain'],
   'reason':'No sidechain N/O/S selector. Nonpolar carbon/hydrogen proximity does not independently select these groups.'},
  {'groups':['Glu213 C/O and Phe214 N/H peptide','Trp349 C/O and Asp350 N/H peptide'],
   'reason':'Both polar selector ends are outside 4.3 A of every fixed anchor in all four discovery members. The adjacent selected donor sidechain or peptide does not force the whole residue into the union.'},
  {'groups':['Trp591 sidechain'],
   'reason':'Indole NE1 selector is outside 4.3 A in all four discovery members despite some ring atoms approaching moved Glu hydrogens.'},
  {'groups':['Phe442 peptide versus sidechain'],
   'reason':'Phe442 O selects its peptide in expanded context, but peptide closure/capping does not include its nonpolar sidechain or CA/HA.'}],
 'source_member_replay':cases,'relevant_fragments':union['fragments'],
 'interpretation':'No claim that every omitted atom needs addition. Source parent graft is a coordinate diagnostic, not an evaluated full-protein structure.',
 'pins':list(PINS.values()),'new_molecular_calls':0}
write('OMISSION_POLICY.json',policy)

collection_path=ROOT/'workspaces/accommodation_fold_DFT_20260920/run_v1/final_collection.json'
collection=load(collection_path)
row=next(r for r in collection['rows'] if r['case_id']=='q9z4j7-pqq-la_model__conditioned_La__seed-1_sample-0')
endpoints={}
for metal,e in row['endpoints'].items():
    for k in ['input','output','receipt','manifest','xyz','current_core_xyz']:
        assert pin(e[k]['path'])['sha256']==e[k]['sha256']
    a,b=xyz(e['xyz']),xyz(e['current_core_xyz']);assert [x[0] for x in a]==[x[0] for x in b]
    maximum=max(math.dist(x[1],y[1]) for x,y in zip(a,b));assert maximum==0
    inputtext=Path(e['input']['path']).read_text();outputtext=Path(e['output']['path']).read_text()
    charge,mult=map(int,re.search(r'\*\s*xyzfile\s+(-?\d+)\s+(\d+)',inputtext,re.I).groups())
    energy=float(re.findall(r'FINAL SINGLE POINT ENERGY\s+([-+\d.]+)',outputtext)[-1]);assert energy==e['energy_hartree']
    assert 'ORCA TERMINATED NORMALLY' in outputtext
    endpoints[metal]={'artifacts':{k:e[k] for k in ['input','output','receipt','manifest','xyz','current_core_xyz']},
       'method_headers':[l for l in inputtext.splitlines() if l.strip().startswith('!')],
       'atoms':len(a),'charge':charge,'multiplicity':mult,'energy_hartree':energy,
       'normal_termination':True,'maximum_current_core_coordinate_difference_A':maximum,
       'status':'Complete actual unchanged-core origin; incompatible with 145/192-atom context or either moved proposal.'}
scope=[
 ('diagnostics/accommodation_torsion_20260920/FINAL_DFT_REPORT.md','Completed native response is only 4MAE and two PLM extra-Asp ±0.2 rad paths; explicitly no Glu/coupled/relaxed-candidate validation. No q9.'),
 ('diagnostics/strict_donor_response_20260923/REPORT.md','Reuses the same 18 native DFT points; only GFN2 stopping repeated. No q9.'),
 ('diagnostics/accommodation_nonlinear_20260920/REPORT.md','Native candidate gate failed; zero candidate DFT launched. Four-case pilot excludes q9.'),
 ('diagnostics/accommodation_nonlinear_20260920/NATIVE_VALIDATION_PLAN.md','Proposed-only native candidate validation, not execution evidence.'),
 ('diagnostics/nikasha_recovery_20260922/DELIVERY.md','Separate candidate DFT remains dry-run only.')]
inventory={'implementation':pin(__file__),'collection':pin(collection_path),'case_id':row['case_id'],
 'canonical_coordinate_match_in_collection':row['canonical_coordinate_match'],'actual_origin_endpoints':endpoints,
 'current_context_states':{'original':{'atoms':145,'Ca_charge':-3,'La_charge':-2},'expanded':{'atoms':192,'Ca_charge':-3,'La_charge':-2}},
 'bounded_related_inventory':[{'source':pin(ROOT/p),'finding':finding} for p,finding in scope],
 'exact_compatible_q9_motion_DFT':'Unavailable in this bounded inventory. No actual saved origin-to-original_La or origin-to-expanded_La DFT pair identified for either current context.',
 'not_established':'An unchanged smaller-core origin cannot validate the candidate work, omitted-environment resistance, or the MACE/solvent decomposition. This is not a claim to have exhausted every archive.',
 'new_DFT_calls':0,'new_molecular_calls':0,'pins':list(PINS.values())}
write('DFT_INVENTORY.json',inventory)

with (OUT/'OMITTED_CONTACTS.csv').open('x',newline='') as f:
    fields=['moving','omitted','both_heavy','source_bond_separation_up_to3','excluded_bond_or_angle_pair','one_four_pair','origin_A','original_La_A','expanded_La_A','origin_bin','original_La_bin','expanded_La_bin']
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
    for r in result['omitted_pairs']:
        w.writerow(dict(moving=r['moving']['label'],omitted=r['omitted']['label'],**{k:r[k] for k in fields[2:6]},**dict(zip(fields[6:9],r['distances_A'].values())),**dict(zip(fields[9:],r['bins'].values()))))
print(json.dumps({'member_replay':len(cases),'DFT_origin_atoms':{k:v['atoms'] for k,v in endpoints.items()},'exact_context_motion_DFT':'unavailable','new_molecular_calls':0}))
