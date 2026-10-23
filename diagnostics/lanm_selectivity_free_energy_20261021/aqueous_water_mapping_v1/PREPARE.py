#!/usr/bin/env python3
"""Exact water-only BIO mapping; no numerical engine entry point."""
import argparse
import collections
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import shlex
import shutil
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
WORK=ROOT/'workspaces/lanm_selectivity_free_energy_20261021';EXPLICIT=WORK/'explicit_model_v1'
VALUES=[3.906136716492464,.28782327521600637,3.5865909437359376,.2837564269122712,.20765190997435956,.3173248056679709]
NAMESPACE={'La':900,'Dy':901};NONMETALS=(402,403,405,414)

def require(ok,why):
    if not ok:raise ValueError(why)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def default_config():
    return {'implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),
        'mapping':pin(WORK/'periodic_mapping_v1/prepared_v1/MAPPING.json'),
        'paired_preparation':pin(WORK/'paired_legs_v1/prepared_v1/RESULT.json'),
        'candidate':pin(EXPLICIT/'paired_water_six_parameter_v1/run_v1/execution/RESULT.json'),
        'candidate_review':pin(EXPLICIT/'paired_water_six_parameter_v1/result_review_v1/REVIEW.json'),
        'transfer':pin(EXPLICIT/'paired_water_six_transfer_v1/run_v1/execution/RESULT.json'),
        'transfer_review':pin(EXPLICIT/'paired_water_six_transfer_v1/actual_transfer_review_v1/REVIEW.json'),
        'inventory_helper':pin(HERE.parent/'periodic_mapping_v1/INVENTORY.py'),
        'resolver':pin(HERE.parent/'alchemical_path_v1/PATH.py'),
        'values':VALUES,'namespace':NAMESPACE,'common_metal_mass_Da':138.91,
        'native_calls_authorized':0,'donor_parameter_binding':False}

def load(cfg):
    require(cfg==default_config(),'Configuration differs from frozen inputs')
    expected={'mapping':'edaa4cc6cbeb8b82670e0314423432667102640eb4cf0df84e2eb0bc300f34e7',
        'candidate':'6027e319a261eb85171fb06e30aa3e519411f6371a7905eeca0e8a79812c48fc',
        'candidate_review':'a3147adcd7448012437e3ac275872fa58fe65e98a114619a95ec01e17a0df187',
        'transfer':'e3eeb870aac22f408c8b38a5f050404cfee5a8f4c0468633c07122fbc5e07254',
        'transfer_review':'3f7b3a9b375cf3cdad0580a75852ce9488a664232aa1da8b1254b547ac2b4ce4'}
    for k,v in expected.items():require(cfg[k]['sha256']==v,'Wrong actual authority: '+k)
    for p in cfg.values():
        if isinstance(p,dict) and 'sha256' in p:verify(p)
    mapping=read(cfg['mapping']);paired=read(cfg['paired_preparation']);candidate=read(cfg['candidate']);transfer=read(cfg['transfer'])
    cr=read(cfg['candidate_review']);tr=read(cfg['transfer_review']);bound=read(transfer['bound'])
    require(candidate['job_id']=='17605' and candidate['candidate_frozen'] and candidate['development_checks_passed'] and candidate['error'] is None,'Passing water candidate required')
    require(candidate['selected']['parameters']==VALUES and candidate['selected']['score']['development_checks_passed'],'Selected water values changed')
    require(cr['result']==cfg['candidate'] and cr['all_selected_raw_replays_exact'] and cr['exact_evaluator_replay'],'Candidate raw review differs')
    require(transfer['job_id']=='17612' and transfer['stage_A_passed'] and transfer['error'] is None and not transfer['global_stop'],'Completed numerical transfer required')
    require(tr['actual_result']==cfg['transfer'] and tr['status']=='passed_actual_raw_numerical_and_comparison_review','Transfer review differs')
    require(bound['candidate_result']==cfg['candidate'] and bound['parameters']==VALUES,'Transfer candidate differs')
    require(mapping['paired_preparation']==cfg['paired_preparation'] and paired['status']=='geometry_prepared_hamiltonian_unassigned','Wrong paired preparation')
    require(mapping['protein_parameters']==paired['protein_parameters'] and mapping['protein_parameters']['sha256']=='eda92c3683846e340c9b07991ed15824c00023d3f3b64effee08d1b53583b6ca','BIO parameter file differs')
    verify(mapping['protein_parameters']);inv=module(cfg['inventory_helper'],'aqueous_inventory');resolver=module(cfg['resolver'],'aqueous_resolver')
    bio=inv.parse(mapping['protein_parameters']);check_namespace(bio)
    sources={};intrinsic={}
    for metal,oldtype in (('La',500),('Dy',501)):
        gid='aquo8_central_'+metal;row=transfer['rows'][gid];value=candidate['selected']['model'][gid];files=bound['inputs'][gid]['files']
        require(row['status']=='complete' and row['cartesian_check']['fd_passed'],'Aquo numerical check missing')
        require(all(row['model'][k]==value[k] for k in ('energy','forces','native')),'Transfer is not selected candidate E/F')
        for k,p in files.items():verify(p);require(p['sha256']==value['provenance']['effective_input'][k]['sha256'],'Effective candidate bytes differ')
        for k in ('analyze_output','gradient_output'):verify(value['native'][k])
        require(value['native']['finite_difference'] is None,'Expected saved selected analytic record')
        pp=resolver.parse_parameters(verify(files['parameters']).read_text(),verify(files['key']).read_text());i=0 if metal=='La' else 2
        pair=resolver.resolve_pair(pp,oldtype,36)
        require([pair['Rmin_pair_A'],pair['thole'],pair['epsilon_kcal_mol']]==[VALUES[i],VALUES[i+1],VALUES[4 if metal=='La' else 5]],'Resolved selected water pair differs')
        intrinsic[metal]={k:copy.deepcopy(pp[k][oldtype]) for k in ('atom','vdw','polarize','metal_multipoles')}
        sources[metal]={'gid':gid,'files':files,'model':value,'resolved_water_pair':pair,'parameters':pp}
    return mapping,paired,bio,intrinsic,sources,inv,resolver

def check_namespace(bio):
    for n in NAMESPACE.values():
        require(n not in bio['atom'] and n not in bio['vdw'] and n not in bio['polarize'] and n not in bio['multipole'],'Metal namespace occupied')
        require(all(v['class']!=n for v in bio['atom'].values()),'Metal class occupied')
        require(all(n not in map(int,v['tokens'][1:3]) for v in bio['pairs']),'Metal pair namespace occupied')

def render_key(intrinsic):
    lines=['parameters amoebabio09.prm','polar-eps 1.0e-8','digits 10','openmp-threads 1',
        '# Frozen17605 aqueous mapping; BIO nonmetal records remain unchanged.']
    for metal,n in NAMESPACE.items():
        d=intrinsic[metal];a=d['atom'];v=d['vdw'];p=d['polarize'];i=0 if metal=='La' else 2
        require(d['metal_multipoles']==[3.]+[0.]*9 and not p['groups'] and a['valence']==0,'Intrinsic isolated-ion state differs')
        lines += [f'atom {n} {n} {metal} "{metal}III frozen water mapping; common kinetic mass" {a["atomic_number"]} 138.91 0',
            f'multipole {n} 0 0 3.0',' 0.0 0.0 0.0',' 0.0',' 0.0 0.0',' 0.0 0.0 0.0',
            f'polarize {n} {p["alpha_A3"]:.17g} {p["thole"]:.17g}',
            f'vdw {n} {v["diameter_A"]:.17g} {v["epsilon_kcal_mol"]:.17g}',
            f'vdwpr 73 {n} {VALUES[i]:.17g} {VALUES[4 if metal=="La" else 5]:.17g}',
            f'polpair 402 {n} {VALUES[i+1]:.17g}']
    return '\n'.join(lines)+'\n'

def validate_key(text,intrinsic,bio):
    check_namespace(bio)
    require(text==render_key(intrinsic),'Unselected key override, recipient or parameter change')
    # Read actual emitted definitions, not a relabeled copy of the input objects.
    p={'settings':{},'atom':{},'vdw':{},'polarize':{},'vdwpr':{},'polpair':{},'metal_multipoles':{}}
    for typ in NONMETALS:
        a=bio['atom'][typ];v,=bio['vdw'][a['class']];pol,=bio['polarize'][typ]
        p['atom'][typ]={'class':a['class'],'atomic_number':a['atomic_number'],'mass':float(a['mass']),'valence':a['valence']}
        t=v['tokens'];p['vdw'][a['class']]={'diameter_A':float(t[2]),'epsilon_kcal_mol':float(t[3]),'reduction':float(t[4]) if len(t)>4 else 0.}
        t=pol['tokens'];p['polarize'][typ]={'alpha_A3':float(t[2]),'thole':float(t[3]),'groups':list(map(int,t[4:]))}
    lines=text.splitlines()
    for i,line in enumerate(lines):
        t=shlex.split(line,comments=True)
        if not t:continue
        if t[0]=='atom':p['atom'][int(t[1])]={'class':int(t[2]),'atomic_number':int(t[5]),'mass':float(t[6]),'valence':int(t[7])}
        elif t[0]=='vdw':p['vdw'][int(t[1])]={'diameter_A':float(t[2]),'epsilon_kcal_mol':float(t[3]),'reduction':0.}
        elif t[0]=='polarize':p['polarize'][int(t[1])]={'alpha_A3':float(t[2]),'thole':float(t[3]),'groups':list(map(int,t[4:]))}
        elif t[0] in ('vdwpr','polpair'):p[t[0]][tuple(sorted(map(int,t[1:3])))]=tuple(map(float,t[3:]))
        elif t[0]=='multipole':p['metal_multipoles'][int(t[1])]=[float(t[-1])]+[float(v) for row in lines[i+1:i+5] for v in row.split()]
    require(set(p['vdwpr'])=={(73,900),(73,901)} and set(p['polpair'])=={(402,900),(402,901)},'Unexpected pair coverage')
    return p

def map_aquo(text,metal):
    lines=text.splitlines();require(int(lines[0].split()[0])==25 and len(lines)==26,'Expected isolated aquo8')
    old=500 if metal=='La' else 501;mapping={old:NAMESPACE[metal],36:402,37:403};out=[lines[0]+'; BIO atom-type mapping']
    types=[];bonds=[];positions=[]
    for i,line in enumerate(lines[1:],1):
        t=line.split();typ=int(t[5]);require(int(t[0])==i and typ in mapping,'Unexpected source type/order')
        require((typ==old)==(i==1) and (i!=1 or len(t)==6),'Metal is not isolated first atom')
        t[5]=str(mapping[typ]);out.append(' '.join(t));types.append(int(t[5]));bonds.append([int(x)-1 for x in t[6:]]);positions.append([float(x) for x in t[2:5]])
    require(collections.Counter(types)=={NAMESPACE[metal]:1,402:8,403:16},'Wrong aquo composition')
    return '\n'.join(out)+'\n',types,bonds,positions

def global_settings(text):
    block=text.split('Literature References',1)[0];out={}
    for line in block.splitlines():
        t=line.split()
        if len(t)==2 and t[0]!='forcefield' and not t[0].startswith('#'):out[t[0]]=t[1]
    return out

def prepare(config_path,out):
    wall=time.monotonic();cpu=time.process_time();cfg=json.loads(Path(config_path).read_text());mapping,paired,bio,intrinsic,sources,inv,resolver=load(cfg)
    key=render_key(intrinsic);resolved=validate_key(key,intrinsic,bio);biofile=verify(mapping['protein_parameters']);smallfile=verify(mapping['small_parameters'])
    require(global_settings(biofile.read_text())==global_settings(smallfile.read_text()),'Global nonbonded/bonded conventions differ')
    require(mapping['water_comparison']['different_terms']==['ureybrad_74_73_74'],'Unrecorded water difference')
    source_ub=mapping['water_comparison']['small_normalized']['ureybrad_74_73_74'];bio_ub=mapping['water_comparison']['protein']['ureybrad_74_73_74']
    require(source_ub['normalized_tokens']==['ureybrad','74','73','74','-7.60','1.5337'] and bio_ub['normalized_tokens']==['ureybrad','74','73','74','-7.60','1.5537'],'Urey difference changed')
    pairs={};cells=[];boxes=[];out.mkdir(parents=True,exist_ok=False)
    for metal,n in NAMESPACE.items():
        pairs[metal]={str(t):resolver.resolve_pair(resolved,n,t) for t in NONMETALS}
        i=0 if metal=='La' else 2;w=pairs[metal]['402'];require([w['Rmin_pair_A'],w['thole'],w['epsilon_kcal_mol']]==[VALUES[i],VALUES[i+1],VALUES[4 if metal=='La' else 5]],'Mapped water pair differs')
        si=sources[metal];oldtype=500 if metal=='La' else 501;siH=resolver.resolve_pair(si['parameters'],oldtype,37);bh=pairs[metal]['403']
        require(all(siH[k]==bh[k] for k in ('Rmin_pair_A','epsilon_kcal_mol','thole')),'WaterH native mixing changed')
        for t in (403,405,414):require(pairs[metal][str(t)]['VDW_origin']=='native_CUBIC_MEAN_HHG' and pairs[metal][str(t)]['POL_origin']=='native_minimum_nonzero','Unselected H/salt pair override')
        for kind in ('aquo8','aqueous'):
            cid='aquo8_central_'+metal+'_BIO' if kind=='aquo8' else metal+'_aqueous_BIO';directory=out/cid;directory.mkdir()
            if kind=='aquo8':
                xyz,types,bonds,positions=map_aquo(verify(si['files']['xyz']).read_text(),metal);(directory/'system.xyz').write_text(xyz)
            else:
                source=paired['aqueous_boxes'][metal];shutil.copyfile(verify(source),directory/'system.xyz');lines=(directory/'system.xyz').read_text().splitlines();require(int(lines[0].split()[0])==11958 and len(lines)==11960,'Aqueous box dimensions changed')
                types=[int(line.split()[5]) for line in lines[2:]];require(collections.Counter(types)=={n:1,402:3980,403:7960,405:7,414:10},'Aqueous composition changed');require(pin(directory/'system.xyz')['sha256']==source['sha256'],'Aqueous coordinates/topology changed')
            (directory/'system.key').write_text(key);shutil.copyfile(biofile,directory/'amoebabio09.prm')
            files={k:pin(directory/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoebabio09.prm'))};require(files['parameters']['sha256']==mapping['protein_parameters']['sha256'],'BIO bytes changed')
            row={'id':cid,'metal':metal,'atoms':25 if kind=='aquo8' else 11958,'files':files,**files,'parameter_filename':'amoebabio09.prm','types':types,'native_calls':0,'numerical_execution_authorized':False}
            if kind=='aquo8':
                row.update(reference_native={k:si['model']['native'][k] for k in ('analyze_output','gradient_output')},reference_xyz=si['files']['xyz'],reference_key=si['files']['key'],reference_parameters=si['files']['parameters'],si_model=si['model'],si_result=cfg['transfer'],candidate_result=cfg['candidate'],bonds_zero_based=bonds,positions_A=positions,atom_map=[{'source_index1':j+1,'mapped_index1':j+1,'source_type':500 if metal=='La' else 501,'mapped_type':n} if j==0 else {'source_index1':j+1,'mapped_index1':j+1,'source_type':36 if t==402 else 37,'mapped_type':t} for j,t in enumerate(types)])
                row['reference_native']['derivative_test']=False;cells.append(row)
            else:row.update(source=paired['aqueous_boxes'][metal],prepared_only=True,box_length_A=paired['box_length_A'],formal_charge=0);boxes.append(row)
    result={'status':'aqueous_water_mapping_prepared_no_calls','implementation':pin(__file__),'configuration':pin(config_path),'sources':{k:v for k,v in cfg.items() if isinstance(v,dict) and 'sha256'in v},'values':VALUES,'intrinsic_source_definitions':intrinsic,'mapped_namespace':NAMESPACE,'common_metal_mass_Da':138.91,'Dy_original_mass_Da':intrinsic['Dy']['atom']['mass'],'BIO_parameters_preserved_byte_for_byte':mapping['protein_parameters'],'cells':cells,'aqueous_endpoints':boxes,'resolved_metal_pairs':pairs,'water_model_difference':{'only_mapped_water_record_difference':'ureybrad_74_73_74','SI':source_ub,'BIO':bio_ub,'global_settings_equal':True,'analytic_energy_force_correction_computed':False,'physical_difference_erased':False},'donor_parameters_bound':False,'periodic_numerical_policy_bound':False,'native_calls':0,'energy_calls':0,'force_calls':0,'dynamics_steps':0,'physical_accuracy_qualified':False,'local_threads':1,'elapsed_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu}
    save(out/'MANIFEST.json',result);print(json.dumps({'manifest':pin(out/'MANIFEST.json'),'cells':len(cells),'aqueous_endpoints':len(boxes),'seconds':result['elapsed_seconds']}));return result

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True);a=s.add_parser('config');a.add_argument('--out',required=True);a=s.add_parser('prepare');a.add_argument('--config',required=True);a.add_argument('--out',required=True);a=p.parse_args()
    if a.action=='config':save(a.out,default_config())
    else:prepare(a.config,Path(a.out).resolve())
