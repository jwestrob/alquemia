#!/usr/bin/env python3
"""Metadata freeze or Slurm-only two-bond radial geometry sensitivity. No energies."""
import argparse, copy, hashlib, importlib.util, json, os, resource, shutil, socket, sys, time, traceback
from pathlib import Path
import xml.etree.ElementTree as ET
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/donor_H_sensitivity'
WORK=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015'
PROTOCOL='A0AC_Asn290_two_H_ff19SB_radial_same_path_sensitivity_v1'
METALS=('Ca','La'); CANDIDATES=('origin','amide_flip','relaxed_flip_La'); HYDROGENS=('HD21','HD22')
SOURCES=tuple('a0acd6b9f2-pqq-la_model__conditioned_'+z+'__seed-1_sample-4' for z in METALS)
def require(v,why):
    if not v: raise ValueError(why)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def check(p):
    require(pin(p['path'])==p,'Changed input: '+p['path']);return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    spec=importlib.util.spec_from_file_location(name,check(p));m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def source_key(a):return tuple(a.get(k,'') for k in ('chain','resnum','insertion_code','resname','atom','element'))
def walk_pins(v):
    if isinstance(v,dict):
        if set(v)=={'path','sha256'}:check(v)
        else:
            for x in v.values():walk_pins(x)
    elif isinstance(v,list):
        for x in v:walk_pins(x)

def freeze(output):
    old=read(WORK/'asn_relaxation_v2/CONFIG.json');pp=pin(WORK/'asn_relaxation_v2/preparation/PREPARATION.json');prep=read(pp)
    cp=pin(WORK/'asn_relaxation_remaining_v1/execution/FINAL_COLLECTION.json');collection=read(cp)
    require(collection['all20_searches_and_common6_complete'] and len(collection['rows'])==10,'Actual complete common6 required')
    sources=[]
    for sid in SOURCES:
        p,=[r for r in prep['rows'] if r['source_case_id']==sid];r,=[r for r in collection['rows'] if r['source_case_id']==sid]
        require(p['status']=='prepared' and r['status']=='complete','Source unavailable')
        sources.append({'source_case_id':sid,'prepared':p,'baseline':{z:{g:r['matrix'][z][g] for g in CANDIDATES} for z in METALS}})
    archived=ROOT/'workspaces/mace_hydrogen_20260916/pilot_v2/medium/implementation'
    ff=Path('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/lib/python3.11/site-packages/openmm/app/data/amber19/protein.ff19SB.xml')
    forcefield=pin(ff);require(forcefield['sha256']=='086bfdb1c05e7d5cb1d330e6c39fab5056a95e89492cdbd3044c746db13e9d09','Wrong documented ff19SB')
    kernel={n:pin(archived/n) for n in ('mace_hydrogen.py','mace_response_trace.py','affordable_common.py')}
    require(kernel['mace_hydrogen.py']['sha256']=='e187838f57d54d7d9e7d52e1b250e5410a5c76f48fff10ee6d14ea07111f5ed1','Wrong archived radial kernel')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);(out/'implementation').mkdir()
    for name in ('PREPARE.py','RUN.sbatch','PLAN.md'):shutil.copy2(HERE/name,out/'implementation'/name)
    c={'protocol':PROTOCOL,'implementation':{n:pin(out/'implementation'/n) for n in ('PREPARE.py','RUN.sbatch')},'plan':pin(out/'implementation/PLAN.md'),
       'python':old['python'],'preparation':pp,'collection':cp,'sources':sources,'source_ids':list(SOURCES),'candidate_order':list(CANDIDATES),
       'forcefield':forcefield,'radial_kernel':kernel,'source_helper':old['source_preparation_helper'],'source_parser':old['source_parser'],
       'settings':{'replay_tolerance_A':1e-12,'subset_atom_names':['ND2',*HYDROGENS],'new_overlap_heavy_A':1.0,'new_overlap_H_A':.55,
                   'serialization':'.17g; all other coordinate lines byte-identical','geometry_only':True},
       'declared_sources':2,'declared_endpoints':12,'output':str(out/'execution'),'new_energy_calls':0}
    walk_pins({'sources':sources,'kernel':kernel});write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def template_bonds(c,context,named,entries,state):
    xml=ET.parse(check(c['forcefield'])).getroot();template,=[x for x in xml.findall('./Residues/Residue') if x.get('name')=='ASN']
    types={x.get('name'):x.attrib for x in xml.findall('./AtomTypes/Type')};atoms={x.get('name'):x.get('type') for x in template.findall('Atom')}
    meta=context['source_atom_metadata'];phys={a['atom']:i for i,a in enumerate(meta) if a and source_key(a)[:4]==('A',290,'','ASN')}
    require(set(phys)==set(named)==set(atoms),'Source/local/ASN-template atom inventory differs')
    require(all(meta[phys[n]]['element']==types[atoms[n]]['element']==entries[named[n]]['source']['element'] for n in atoms),'ASN template element differs')
    bonds=[tuple(b) for b in context['bonds']];rev={v:k for k,v in phys.items()}
    internal={frozenset((rev[i],rev[j])) for i,j in bonds if i in rev and j in rev}
    expected={frozenset((b.get('atomName1'),b.get('atomName2'))) for b in template.findall('Bond')}
    require(internal==expected,'Source ASN internal graph differs from exact template')
    external=[(rev[i],j) if i in rev else (rev[j],i) for i,j in bonds if (i in rev)!=(j in rev)]
    require(sorted(n for n,j in external)==sorted(x.get('atomName') for x in template.findall('ExternalBond')),'ASN is not the supported interior template')
    for n,j in external:
        require(meta[j] and meta[j]['chain']=='A' and meta[j]['atom']==('C' if n=='N' else 'N'),'External source bond is not a peptide attachment')
    parents={tuple(h):tuple(p) for h,p in state['H_covalent_parents']};rows=[]
    for h in HYDROGENS:
        require(parents[source_key(entries[named[h]]['source'])]==source_key(entries[named['ND2']]['source']),'Frozen H parent differs')
        require([j if i==phys[h] else i for i,j in bonds if phys[h] in (i,j)]==[phys['ND2']],'H has unsupported covalent neighbors')
        classes={types[atoms[n]]['class'] for n in ('ND2',h)}
        term,=[x for x in xml.findall('./HarmonicBondForce/Bond') if {x.get('class1'),x.get('class2')}==classes]
        rows.append({'atom_a_id':'A/290//ND2','atom_b_id':'A/290//'+h,'forcefield_equilibrium_A':10*float(term.get('length')),
            'template':'ASN','atom_names':['ND2',h],'atom_types':[atoms['ND2'],atoms[h]],'bond_parameter':dict(term.attrib),
            'units':'XML length nm converted once to angstrom; no force-field energy evaluated'})
    return rows,{'atom_types':atoms,'internal_bonds':[sorted(x) for x in sorted(expected,key=lambda x:sorted(x))],
                 'external_source_bonds':[{'atom':n,'neighbor':meta[j]} for n,j in external],'neutral_state_verified':True}

def prepare(s,d,c,np,helper,base,kernel):
    p=s['prepared'];r=p['original_row'];sid=s['source_case_id'];state=read(r['state_signature']);cp=read(r['source_preparation'])
    require(read(p['admission'])['pass'] and cp['source']==r['protonated_source']==p['contact']['source'],'Source/admission differs')
    require(cp['new_atom_count']==204 and state['endpoint_charges']=={'Ca':-4,'La':-3} and state['multiplicities']=={'Ca':1,'La':1} and not state['water_inventory'],'Frozen dry paired state differs')
    require(['A:ASN290',0,None] in state['core_microstates'],'Neutral Asn state absent')
    maps={z:read(r['maps'][z]) for z in METALS};require(maps['Ca']==maps['La'],'Paired source maps differ')
    entries={x['qm_index']:x for x in cp['mapping']['source_to_qm']};require(set(entries)==set(range(1,204)),'Incomplete204 source map')
    named={x['source']['atom']:i for i,x in entries.items() if x.get('source') and source_key(x['source'])[:4]==('A',290,'','ASN')}
    require(named==r['named_Asn_indices_zero_based'],'Named source/local Asn mapping differs')
    bonds,template=template_bonds(c,maps['Ca']['context'],named,entries,state);moved=[named[h] for h in HYDROGENS];fixed=[i for i in range(204) if i not in moved]
    caps=set(p['contact']['caps']);require(len(caps)==14 and not caps.intersection(moved),'Physical H/cap membership differs')
    ids=read(p['contact']['identities'])['centers'];saved=np.load(check(p['contact']['origin']),allow_pickle=False)
    union=saved['positions'].copy();symbols=saved['symbols'].tolist();saved.close()
    require(len(ids)==len(union)==9546 and len(symbols)==len(union),'Frozen source union differs')
    physical=helper.physical_source(check(r['protonated_source']),base);lookup={source_key(a['source']):a for a in physical}
    source_ids=[source_key(a['source']) for a in ids if 'source' in a]
    require(len(source_ids)==len(set(source_ids))==len(physical)==9529 and set(source_ids)==set(lookup),'Source identity union duplicates/omits atoms')
    for i,a in enumerate(ids):
        if i>=204:require(np.array_equal(union[i],np.asarray(lookup[source_key(a['source'])]['xyz'])),'Exterior source coordinate changed')
        if i in moved:require(a['kind']=='physical_source' and a['source']==entries[i]['source'],'Moved H is not the exact source site')
    d.mkdir();write(d/'FF19SB_TARGETS.json',{'forcefield':c['forcefield'],'template_match':template,'bonds':bonds,'subset_policy':'Only actual ND2,HD21,HD22 are passed to the unchanged archived repair kernel'})
    matrix={z:{} for z in METALS};admissions=[]
    for candidate in CANDIDATES:
        original={z:s['baseline'][z][candidate] for z in METALS};pair={z:base.xyz_atoms(check(original[z]['xyz'])) for z in METALS}
        require(all(len(v)==204 and v[0]['element']==z for z,v in pair.items()),'Wrong paired204 endpoint')
        before=np.asarray([a['xyz'] for a in pair['Ca']]);require(np.array_equal(before,np.asarray([a['xyz'] for a in pair['La']])) and all(pair['Ca'][i]['element']==pair['La'][i]['element']==symbols[i] for i in range(1,204)),'Paired coordinate/element mismatch')
        view=[{'id':'A/290//'+n,'source_index':named[n],'element':entries[named[n]]['source']['element'],'xyz_A':before[named[n]].tolist()} for n in ('ND2',*HYDROGENS)]
        repaired,moves=kernel.repair(view,bonds);after=before.copy()
        for a in repaired:
            if a['element']=='H':after[a['source_index']]=a['xyz_A']
        require(len(moves)==2 and np.array_equal(after[fixed],before[fixed]),'Projection changed the wrong coordinates')
        radial=[]
        for h,b in zip(HYDROGENS,bonds):
            i=named[h];j=named['ND2'];oldvec=before[i]-before[j];newvec=after[i]-after[j];target=b['forcefield_equilibrium_A'];tol=c['settings']['replay_tolerance_A']
            length=float(np.linalg.norm(newvec));expected=before[j]+oldvec*(target/np.linalg.norm(oldvec))
            require(abs(length-target)<=tol and np.max(np.abs(after[i]-expected))<=tol,'Radial target/direction replay failed')
            radial.append({'hydrogen':h,'parent':'ND2','local_index':i,'old_length_A':float(np.linalg.norm(oldvec)),'new_length_A':length,'target_A':target})
        before_union=union.copy();before_union[:204]=before;after_union=before_union.copy();after_union[:204]=after
        contacts=helper.contact_check(np,before_union,after_union,symbols,moved,ids,caps)
        admission={'status':'prepared' if contacts['pass'] else 'unavailable','pass':contacts['pass'],'source_case_id':sid,'candidate':candidate,
            'radial_checks':radial,'all_other_coordinates_fixed':True,'cap_coordinates_and_anchors_fixed':True,'all_protons_retained':True,
            'paired_coordinates_equal':True,'moved_indices_zero_based':moved,'contacts':contacts,'ff19SB_targets':pin(d/'FF19SB_TARGETS.json'),
            'old_geometry_only_reference':True,'no_optimizer_domain_claim':True,'chemical_state':r['state_signature']}
        ap=d/(candidate+'_ADMISSION.json');write(ap,admission);admissions.append(pin(ap))
        for z in METALS:
            old=original[z];require(old['status']=='complete','Old MACE cell unavailable');walk_pins(old)
            result=read(old['receipt']);require(result['energy_eV']==old['energy_eV'],'Old MACE energy join differs')
            if 'full_source_admission' in old:require(read(old['full_source_admission'])['pass'],'Saved relaxed geometry was not admitted')
            lines=check(old['xyz']).read_text().splitlines(keepends=True)
            for i in moved:lines[i+2]='H '+' '.join(format(v,'.17g') for v in after[i])+'\n'
            xp=d/(z+'_'+candidate+'_radial_H.xyz')
            with xp.open('x') as f:f.write(''.join(lines))
            actual=base.xyz_atoms(xp);require(all(actual[i]['line']==pair[z][i]['line'] for i in fixed) and np.array_equal(np.asarray([a['xyz'] for a in actual]),after),'Emitted fixed-coordinate/replay check failed')
            matrix[z][candidate]={'status':admission['status'],'xyz':pin(xp),'charge':state['endpoint_charges'][z],'multiplicity':1,'atom_count':204,
                'geometry_admission':pin(ap),'original':old,'original_xyz':old['xyz'],'original_MACE':old['receipt'],'state_signature':r['state_signature']}
    passed=all(read(a)['pass'] for a in admissions)
    return {'status':'prepared' if passed else 'unavailable','source_case_id':sid,'case_key':r['case_key'],'record_id':r['record_id'],
        'conditioning':r['conditioning'],'sample':4,'maps':r['maps'],'source_preparation':r['source_preparation'],'protonated_source':r['protonated_source'],
        'state_signature':r['state_signature'],'moved_indices_zero_based':moved,'matrix':matrix,'admissions':admissions,'original_contact':p['contact'],
        'maps_role':'Unchanged source identity/topology reference; new XYZ carries the two-H coordinate derivation, not an old-map q0 replay',
        'coordinate_derivation':{'kernel':c['radial_kernel']['mace_hydrogen.py'],'targets':pin(d/'FF19SB_TARGETS.json'),'changed_atom_names':list(HYDROGENS)}}

def execute(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before numerical imports')
    require(int(os.environ.get('SLURM_JOB_NUM_NODES','0'))==1 and int(os.environ.get('SLURM_CPUS_ON_NODE','0').split('(')[0])>=1,'One-node CPU allocation required')
    config={'path':str(Path(path).resolve()),'sha256':sha};c=read(config);require(pin(__file__)==c['implementation']['PREPARE.py'] and pin(sys.executable)==c['python'],'Frozen driver/interpreter differs')
    out=Path(c['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time();rows=[]
    receipt={'status':'failed','config':config,'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'new_energy_calls':0}
    try:
        walk_pins(c);require(c['source_ids']==list(SOURCES) and c['candidate_order']==list(CANDIDATES),'Finite scope changed')
        sys.path.insert(0,str(Path(c['radial_kernel']['mace_hydrogen.py']['path']).parent))
        kernel=load('_asn_radial_kernel',c['radial_kernel']['mace_hydrogen.py']);helper=load('_asn_contact_helper',c['source_helper']);base=load('_asn_source_parser',c['source_parser'])
        import numpy as np
        require(kernel.RECIPE=='ff19sb_protein_H_radial_bond_projection_v1','Archived repair recipe differs')
        for s in c['sources']:
            try:rows.append(prepare(s,out/s['source_case_id'],c,np,helper,base,kernel))
            except Exception as exc:
                rows.append({'status':'unavailable','source_case_id':s['source_case_id'],'reason':str(exc),'traceback':traceback.format_exc(),
                             'matrix':{z:{g:{'status':'unavailable','reason':str(exc),'original':s['baseline'][z][g]} for g in CANDIDATES} for z in METALS}})
        complete=len(rows)==2 and all(r['status']=='prepared' for r in rows);status='complete' if complete else 'partial'
        write(out/'PREPARATION.json',{'protocol':PROTOCOL,'config':config,'status':status,'declared_sources':2,'declared_endpoints':12,'rows':rows,'new_energy_calls':0,
            'interpretation':'Same-path two-H geometry sensitivity only; no candidate-pool or class decision'})
        receipt.update(status=status,preparation=pin(out/'PREPARATION.json'))
    except BaseException:receipt['traceback']=traceback.format_exc();raise
    finally:
        receipt.update(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu,max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();sp=p.add_subparsers(dest='op',required=True);f=sp.add_parser('freeze');f.add_argument('--output',required=True)
    e=sp.add_parser('execute');e.add_argument('--config',required=True);e.add_argument('--config-sha',required=True);a=p.parse_args()
    if a.op=='freeze':freeze(a.output)
    else:sys.exit(execute(a.config,a.config_sha))
