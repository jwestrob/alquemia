#!/usr/bin/env python3
"""Freeze byte pins on login; prepare three exact two-Ca inputs and map1FLG on Slurm."""
from pathlib import Path
import argparse, copy, hashlib, json, os, resource, socket, subprocess, sys, time, traceback, urllib.request
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = Path('/groups/banfield/projects/multienv/corkscrew/supplementary_structures')
ADH9 = Path('/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue/adh9/af3_comparison')
THREADS = ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
PROTOCOL = 'Nikasha_Ca_two_site_occupancy_AF3_v1'

def pin(path):
    p=Path(path).resolve(); h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}

def verify(rec):
    if pin(rec['path']) != rec:raise ValueError('pin changed: '+rec['path'])
    return Path(rec['path'])

def read(path):return json.loads(Path(path).read_text())

def write(path,data):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n')

def require_slurm(threads=1):
    if not os.environ.get('SLURM_JOB_ID') or not os.environ.get('SLURM_JOB_NODELIST'):
        raise RuntimeError('molecular preparation requires real Slurm allocation')
    for key in THREADS:os.environ[key]=str(threads)

def command_prefix(m,gpu=False):
    return [m['singularity']['path'],'exec',*(['--nv'] if gpu else []),'--bind','/groups,/shared','--pwd','/app/alphafold',m['container']['path']]

def freeze(workspace):
    out=Path(workspace).resolve();out.mkdir(parents=True,exist_ok=False)
    prior=read(ADH9/'manifest.json'); targets=[]
    for acc,kind in [('q9z4j7','Ca-binding-controls'),('c5axv8','La-binding-controls'),('q88jh5','Ca-binding-controls')]:
        p=SOURCE/kind/'Ca-alone'/(acc+'-pqq-ca')/(acc+'-pqq-ca_data.json')
        targets.append({'accession':acc,'source_input':pin(p),'new_state_id':acc+'__occupancy2Ca_v1',
                        'new_source_ids':[acc+'-pqq-la_model__conditioned_2Ca_v1__seed-1_sample-'+str(i) for i in range(5)],
                        'original_source_ids':[acc+'-pqq-la_model__conditioned_Ca__seed-1_sample-'+str(i) for i in range(5)]})
    canonical=ROOT/'workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3/CANONICAL_COLLECTION.json'
    m={'protocol_id':PROTOCOL,'targets':targets,'model_seeds':[1],'num_diffusion_samples':5,'num_recycles':10,
       'new_MSA_searches':0,'prepared_output':str(out/'prepared'),'fold_output':str(out/'folds'),
       'original_reference':pin(canonical),'source_placement_audit':pin(ROOT/'workspaces/ca_coverage_recovery_20261012/pqq15/run_v1/execution_v1/AUDIT.json'),
       'original_AF3_manifest':pin(ADH9/'manifest.json'),'original_AF3_runner':pin(ADH9/'run_fold.py'),
       'container':{'path':prior['container'],'sha256':prior['container_sha256']},
       'weights':{'path':prior['weights'],'sha256':prior['weights_sha256']},
       'singularity':pin('/usr/local/bin/singularity'),'python':pin('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python'),
       'implementation':{name:pin(HERE/name) for name in ('PREPARE.py','RUN_AF3.py','PREPARE.sbatch','FOLD.sbatch','PLAN.md')},
       'geometry_code':{name:pin(ROOT/'scripts'/name) for name in ('carve_with_pqq.py','coordination_policy.py','detect_pqq.py','pqq_microstates.py')},
       'reference_url':'https://files.rcsb.org/download/1FLG.cif','scout_accession':'q9z4j7',
       'declared_models':15,'new_scores':0,'new_fit':False,'production_changed':False}
    write(out/'manifest.json',m);print(json.dumps({'manifest':pin(out/'manifest.json'),'executed':False}))

def reference_audit(m,out):
    import gemmi
    sys.path.insert(0,str(ROOT/'scripts'))
    import carve_with_pqq as base
    p=out/'1FLG.cif'
    with urllib.request.urlopen(m['reference_url'],timeout=60) as response:p.write_bytes(response.read())
    write(out/'1FLG_DOWNLOAD.json',{'url':m['reference_url'],'file':pin(p),'retrieved_epoch':time.time()})
    st=gemmi.read_structure(str(p));model=st[0]
    if len(st)!=1:raise ValueError('expected single1FLG model')
    q9=read(verify(m['targets'][0]['source_input']))
    query=q9['sequences'][0]['protein']['sequence'];mappings=[];by_chain={}
    scoring=gemmi.AlignmentScoring();scoring.match=2;scoring.mismatch=-1;scoring.gapo=-6;scoring.gape=-1
    for chain in model:
        residues=[r for r in chain if r.name in base.STANDARD_AA]
        if not residues:continue
        seq=''.join(gemmi.find_tabulated_residue(r.name).one_letter_code for r in residues)
        al=gemmi.align_string_sequences(list(query),list(seq),[],scoring)
        a,b=al.add_gaps(query,1),al.add_gaps(seq,2);qi=ri=0;mapping=[]
        for qa,ra in zip(a,b):
            if qa!='-':qi+=1
            if ra!='-':ri+=1
            if qa!='-' and ra!='-':
                r=residues[ri-1]
                mapping.append({'source_position':qi,'source_residue':qa,'crystal_chain':chain.name,
                                'crystal_resnum':r.seqid.num,'crystal_icode':r.seqid.icode.strip(),
                                'crystal_resname':r.name,'identity':qa==ra})
        by_chain[chain.name]={(r['crystal_resnum'],r['crystal_icode']):r for r in mapping}
        mappings.append({'chain':chain.name,'alignment_score':al.score,'identity_percent':al.calculate_identity(),
                         'alignment':al.formatted(query,seq),'mapping':mapping,
                         'six_donor_region':[r for r in mapping if r['source_position'] in (45,48,51)]})
    ions=[]
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.element.name!='Ca':continue
                site=base.select_metal_sites(model,site_chain=chain.name,site_resnum=residue.seqid.num,
                    site_icode=residue.seqid.icode.strip(),site_resname=residue.name,site_atom=atom.name)[0]
                nearby=base._prepare_nearby_pqq(model,atom.pos,base.DEFAULT_PQQ_MICROSTATE)
                contacts=base._classify_contacts(model,site,nearby)[0]
                rows=[]
                for contact in contacts:
                    row=contact.metadata();row['Q9Z4J7_mapping']=by_chain.get(row['chain'],{}).get((row['resnum'],row['icode']))
                    rows.append(row)
                pqq=[(c.name,r) for c in model for r in c if r.name=='PQQ']
                ions.append({'site':site.label(),'typed_CN':len(contacts),'typed_contacts':rows,
                             'PQQ_associations':[k.label() for k in nearby],
                             'nearest_PQQ_atom_A':min(atom.pos.dist(a.pos) for _,r in pqq for a in r)})
    if len(ions)!=4 or len(mappings)!=2:raise ValueError('1FLG fourCa/two-protein-chain inventory differs')
    result={'reference':pin(p),'URL':m['reference_url'],'expected_entry':'1FLG','Q9_input':m['targets'][0]['source_input'],
            'alignment_policy':{'algorithm':'Gemmi global affine','match':2,'mismatch':-1,'gap_open':-6,'gap_extend':-1,'free_gaps':[]},
            'mappings':mappings,'Ca_sites':ions,'typed_cutoff_A':3.1,'PQQ_association_cutoff_A':4.0,
            'scope':'descriptive reference mapping, not a prediction or scoring admission'}
    write(out/'1FLG_AUDIT.json',result);return pin(out/'1FLG_AUDIT.json')

def execute(path):
    require_slurm();m=read(path)
    for rec in m['implementation'].values():verify(rec)
    for rec in m['geometry_code'].values():verify(rec)
    verify(m['singularity']);verify(m['container']);verify(m['python']);verify(m['original_reference'])
    out=Path(m['prepared_output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    receipt={'manifest':pin(path),'job_id':os.environ['SLURM_JOB_ID'],'hostname':socket.gethostname(),
             'status':'running','new_energy_calls':0,'new_fold_calls':0,'numerical_threads':1}
    write(out/'START.json',receipt)
    try:
        targets=[]
        for target in m['targets']:
            old=read(verify(target['source_input']));new=copy.deepcopy(old)
            assert old['modelSeeds']==[1] and old['dialect']=='alphafold3' and old['version']==2
            assert len(old['sequences'])==3 and old['sequences'][0]['protein']['id']=='A'
            assert old['sequences'][1]=={'ligand':{'id':'B','smiles':'[Ca+2]'}}
            assert old['sequences'][2]['ligand']['id']=='C'
            extra=copy.deepcopy(old['sequences'][1]);extra['ligand']['id']='D';new['sequences'].append(extra)
            restored=copy.deepcopy(new);restored['sequences'].pop();assert restored==old
            ip=out/'inputs'/(target['accession']+'.json');write(ip,new)
            protein=old['sequences'][0]['protein']
            semantic={key:hashlib.sha256(json.dumps(protein[key],sort_keys=True,separators=(',',':')).encode()).hexdigest()
                      for key in ('sequence','modifications','unpairedMsa','pairedMsa','templates')}
            validation='import json,sys;from pathlib import Path;from alphafold3.common import folding_input;p=Path(sys.argv[1]);x=folding_input.Input.from_json(p.read_text(),json_path=p);print(json.dumps({"name":x.name,"seeds":x.rng_seeds,"chains":[c.id for c in x.chains]}))'
            command=command_prefix(m)+['python','-c',validation,str(ip)]
            proc=subprocess.run(command,capture_output=True,text=True)
            (out/(target['accession']+'_schema.log')).write_text(proc.stdout+proc.stderr)
            if proc.returncode:raise RuntimeError('upstream input schema rejected '+target['accession'])
            targets.append({**target,'input':pin(ip),'original_name_retained':old['name'],'protein_fields_sha256':semantic,
                            'only_change':'append separateD ligand identical toB exceptid','all_original_fields_identical':True,
                            'schema_command':command,'schema_log':pin(out/(target['accession']+'_schema.log'))})
            del old,new,restored,protein
        reference=reference_audit(m,out)
        prepared={'protocol_id':PROTOCOL,'manifest':pin(path),'targets':targets,'reference_audit':reference,
                  'status':'prepared','new_scores':0,'new_folds':0,'original_bands':m['original_reference']}
        write(out/'PREPARED.json',prepared);receipt.update(status='complete',prepared=pin(out/'PREPARED.json'))
    except Exception as exc:
        receipt.update(status='failed',reason=str(exc),traceback=traceback.format_exc());raise
    finally:
        usage=resource.getrusage(resource.RUSAGE_SELF);child=resource.getrusage(resource.RUSAGE_CHILDREN)
        receipt.update(wall_seconds=time.monotonic()-start,CPU_seconds=usage.ru_utime+usage.ru_stime+child.ru_utime+child.ru_stime,
                       peak_RSS_KiB=usage.ru_maxrss,child_peak_RSS_KiB=child.ru_maxrss)
        write(out/'EXECUTION.json',receipt)

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('freeze');a.add_argument('--workspace',required=True)
    a=s.add_parser('execute');a.add_argument('--manifest',required=True)
    args=p.parse_args();freeze(args.workspace) if args.op=='freeze' else execute(args.manifest)
