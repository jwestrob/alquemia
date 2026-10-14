#!/usr/bin/env python3
"""Finite fixed-coordinate packing-context preparation; science requires Slurm."""
import argparse, copy, hashlib, json, os, resource, socket, sys, time, traceback, types
from collections import Counter
from pathlib import Path
import xml.etree.ElementTree as ET

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ORDER=('origin','adaptive_Ca','adaptive_La')
METALS=('Ca','La')
VALENCE={'H':1,'C':4,'N':5,'O':6,'S':6,'Ca':2,'La':3}
NEUTRAL=('ALA','VAL','LEU','ILE','PHE','MET')

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p, ('changed input',p)
    return Path(p['path'])
def read(p): return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')

def freeze(review,output):
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    review=Path(review).resolve();scope=read(review)
    ids=['q9z4j7-pqq-la_model__conditioned_2Ca_v1__seed-1_sample-'+str(i) for i in (0,3)]
    requests=[];inputs={}
    def add(p):
        if p:verify(p);inputs[p['path']]=p
    for cid in ids:
        item=next(c for c in scope['cases'] if c['record']['source_case_id']==cid)
        r=item['record'];bm=read(r['baseline_manifest']);case=bm['cases'][r['baseline_manifest_case_index']]
        n=case['native_source_record'];p=n['source']['origin_row']['source'];core=read(p['original_core']['parent'])
        assert p['status']=='prepared' and p['mapping_eligible'] and p['paired_maps_equal']
        old=read(r['source_preparation'])
        for z in METALS:
            add(r['maps'][z])
            for q in ORDER:add(r['matrix'][z][q]['xyz'])
            for key in ('xyz','input'):add(p['original_core']['endpoints'][z][key])
        for x in (r['baseline_manifest'],r['source_preparation'],r['state_signature'],p['original_core']['parent'],p['spectator_admission'],p['source']['source_geometry_audit'],core['source_structure'],*n['proposals'].values()):add(x)
        if core.get('protonation_manifest'):add(core['protonation_manifest']);add(read(core['protonation_manifest'])['output'])
        requests.append({'case_id':cid,'role':'control' if cid.endswith('sample-0') else 'target',
            'original_review_identity':{'case_key':item['case_key'],'record_id':item['record_id']},
            'record':r,'prepared_source':p,'native_proposals':n['proposals'],'old_context_audit':r['source_preparation']})
    sel=ROOT/'workspaces/motion_envelope_20260923/run_v2/SELECTION.json';cfg=read(sel)['config'];add(cfg['topology'])
    modules={p.name:pin(p) for p in sorted((ROOT/'scripts').glob('*.py'))}
    m={'protocol':'Nikasha_fixed_pool_allheavy_packing4p3_v1','driver':pin(__file__),'plan':pin(HERE/'PLAN.md'),
       'review':pin(review),'config_source':pin(sel),'topology':cfg['topology'],'requests':requests,'inputs':list(inputs.values()),
       'script_pins':modules,'python':pin('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python'),
       'settings':{'contact_A':4.3,'moved_A':1e-8,'old_replay_A':1e-12,'kinematic_replay_A':1e-8,
          'bond_preservation_A':1e-8,'minimum_distance_A':.45,'spectator_boundary_A':4.3},
       'output':str(out/'execution'),'declared_sources':2,'declared_new_cells':12,
       'new_energy_calls':0,'new_optimizations':0,'numerical_threads':1}
    write(out/'manifest.json',m);print(json.dumps(pin(out/'manifest.json')))

def source_id(a):return (a['chain_index'],a['residue_index'],a['atom'])
def atom_id(a):
    if a['kind']=='source':return ('source',*source_id(a['source']))
    if a['kind']=='opaque_cofactor':return ('cofactor',a['fragment'],a['name'])
    if a['kind']=='sigma_link_H':return ('cap',source_id(a['retained']),source_id(a['omitted']))
    raise ValueError('unsupported mapped atom kind')

def execute(path):
    assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_JOB_NODELIST'),'actual Slurm required'
    cpus=int(os.environ['SLURM_CPUS_ON_NODE']);assert 1<=cpus<=12 and len(os.sched_getaffinity(0))>=1
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','OPENMM_CPU_THREADS'):os.environ[k]='1'
    m=read(path)
    for p in [m['driver'],m['plan'],m['review'],m['config_source'],m['python'],*m['inputs'],*m['script_pins'].values()]:verify(p)
    out=Path(m['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    receipt={'manifest':pin(path),'job_id':os.environ['SLURM_JOB_ID'],'hostname':socket.gethostname(),
       'allocated_CPUs':cpus,'status':'running','new_energy_calls':0,'new_optimizations':0}
    write(out/'START.json',receipt)
    try:
        import numpy as np
        sys.path.insert(0,str(ROOT/'scripts'))
        import second_shell_context as shell, consistent_context as context
        from mace_site_kinematics import Kinematics
        from affordable_common import xyz,cache_key
        from mace_hybrid import check_atoms,write_xyz
        settings=m['settings'];templates={r.attrib['name']:r for r in ET.parse(verify(m['topology'])).findall('.//Residues/Residue')}
        neutral_checks=[]
        def sidechain(g,rkey):
            r=g.residues[rkey];rn=r.canonical_resname
            if rn not in NEUTRAL:return shell.sidechain(g,rkey)
            assert rkey in g.supported_topology_residues and rn in templates
            template=templates[rn];tnames={a.attrib['name'] for a in template.findall('Atom')}
            names=tnames-{n for n in tnames if n.startswith('H')}-{'N','CA','C','O','OXT'}
            assert 'CB' in names
            observed={n for n in g.by_residue[rkey] if g.atoms[g.by_residue[rkey][n]].element.name not in ('H','D')}-{'N','CA','C','O','OXT'}
            assert observed==names, ('nonpolar heavy membership',r.source_dict(),observed,names)
            expected=Counter()
            for bond in template.findall('Bond'):
                a,b=bond.attrib['atomName1'],bond.attrib['atomName2']
                if a.startswith('H') and b in names:expected[b]+=1
                if b.startswith('H') and a in names:expected[a]+=1
            actual=Counter(p[2] for h,p in g.hparents.items() if p[:2]==rkey and p[2] in names)
            assert actual==expected, ('nonpolar hydrogen state',r.source_dict(),actual,expected)
            nodes={g.key(rkey,n) for n in names}
            neutral_checks.append({'source':r.source_dict(),'heavy_atoms':sorted(names),'template_H_parents':dict(expected),'source_H_parents':dict(actual),'formal_charge':0,'template':m['topology']})
            return nodes,0
        globals_=dict(context.forced_expansion.__globals__);globals_['sidechain']=sidechain
        expand=types.FunctionType(context.forced_expansion.__code__,globals_,'packing_expansion')
        prepared=[];union=set();discovery=[]
        for req in m['requests']:
            p=req['prepared_source'];r=req['record'];s=shell.parent_state(p['original_core'],m['topology'],require_endpoint_receipts=False)
            old=read(req['old_context_audit']);graph=s['graph'];km=read(r['maps']['Ca'])['context'];kin=Kinematics(km)
            original,check=context.forced_expansion(s,old['added_fragments'])
            old_error={}
            for z in METALS:
                saved=xyz(verify(r['matrix'][z]['origin']['xyz']))
                assert [x[0] for x in saved]==[x[0] for x in original[z]]
                old_error[z]=float(np.max(np.abs(np.asarray([a[1:] for a in saved])-np.asarray([a[1:] for a in original[z]]))))
                assert old_error[z]<=settings['old_replay_A']
            coords={};full={};replays={}
            for q in ORDER:
                vec=np.zeros(len(km['modes'])) if q=='origin' else np.asarray(read(req['native_proposals'][q.removeprefix('adaptive_')])['proposal']['full_q'])
                physical,rendered,_,_=kin.evaluate(vec);saved=xyz(verify(r['matrix']['Ca'][q]['xyz']))
                replays[q]=float(np.max(np.abs(rendered-np.asarray([a[1:] for a in saved]))));assert replays[q]<=settings['kinematic_replay_A']
                for i,kind,*rest in km['core_links']:
                    if kind=='source':physical[rest[0]]=saved[i][1:]
                coords[q]=physical
                full[q]={k:np.array([a.pos.x,a.pos.y,a.pos.z]) for k,a in graph.atoms.items()}
                for i,meta in enumerate(km['source_atom_metadata']):
                    if meta:
                        k=source_id(meta);assert k in full[q]
                        full[q][k]=physical[i]
            moved=[i for i,meta in enumerate(km['source_atom_metadata']) if meta and meta['element'] not in ('H','D') and
                   any(float(np.linalg.norm(coords[q][i]-coords['origin'][i]))>settings['moved_A'] for q in ORDER)]
            assert moved,'no real saved movement'
            oldids={source_id(a['source']) for a in old['mapping']['source_to_qm'] if 'source' in a}
            active_selector=p['source']['metal'];active_residue=graph.locate(active_selector)
            active_key=graph.key(active_residue.key,active_selector['atom'])
            assert graph.atoms[active_key].element.name=='Ca'
            fragments={context.fragment_id(f) for f in old['added_fragments']};contacts=[];excluded=[]
            def fragment(kind,rkey):
                f={'kind':kind,'source':graph.residues[rkey].source_dict()};fragments.add(context.fragment_id(f));return f
            def previous(rkey):
                prev=graph.amide_prev.get(rkey);assert prev is not None,('missing preceding peptide',rkey)
                return fragment('peptide',prev[:2])
            for k,atom in graph.atoms.items():
                if k in oldids or k==active_key or atom.element.name in ('H','D'):continue
                distances=[(float(np.linalg.norm(full[q][k]-coords[q][i])),q,i) for q in ORDER for i in moved]
                d,q,i=min(distances)
                if d>settings['contact_A']:continue
                row={'source':graph.meta[k],'nearest_moved':km['source_atom_metadata'][i],'distance_A':d,'candidate':q}
                rn=graph.meta[k]['canonical_resname']
                if rn=='HOH':excluded.append({**row,'reason':'fixed_water_inventory'});continue
                assert k[:2] in graph.supported_topology_residues,('unsupported nearby packing species',row)
                selected=[]
                if k[2] in ('C','O','OXT'):selected.append(fragment('peptide',k[:2]))
                elif k[2]=='N':selected.append(previous(k[:2]))
                elif k[2]=='CA':
                    selected.extend([previous(k[:2]),fragment('peptide',k[:2])])
                    if rn not in ('GLY','PRO'):selected.append(fragment('sidechain',k[:2]))
                elif rn=='PRO':selected.append(previous(k[:2]))
                else:selected.append(fragment('sidechain',k[:2]))
                contacts.append({**row,'selects':selected})
            union|=fragments
            discovery.append({'case_id':req['case_id'],'moved_heavy_atoms':[km['source_atom_metadata'][i] for i in moved],
                'old_replay_max_A':old_error,'kinematic_replay_max_A':replays,'contacts':contacts,'excluded':excluded,
                'source_fragment_ids':sorted(fragments)})
            prepared.append((req,s,old,full,old_error,replays))
        fragments=[context.fragment_record(x) for x in sorted(union)]
        union_record={'rule':m['settings'],'fragments':fragments,'source_ids':[r['case_id'] for r in m['requests']],
            'discovery':discovery,'selection_uses_energies_or_labels':False}
        write(out/'UNION.json',union_record);rows=[]
        for req,s,old,full,old_error,replays in prepared:
            d=out/req['case_id'];d.mkdir();row={'case_id':req['case_id'],'role':req['role'],'status':'unavailable',
                'original_review_identity':req['original_review_identity'],'original_record':req['record'],
                'selection_id':req['record']['selection_id'],'source':req['prepared_source']['source'],
                'original_context':req['old_context_audit'],'old_matrix':req['record']['matrix'],'union':pin(out/'UNION.json')}
            try:
                origin,audit=expand(s,fragments);g=s['graph'];oldmap={atom_id(a):a for a in old['mapping']['source_to_qm']}
                charges={z:req['prepared_source']['original_core']['endpoints'][z]['charge']+audit['added_formal_charge'] for z in METALS}
                assert charges['La']==charges['Ca']+1
                signature=context.state_signature(s,audit,charges);write(d/'state_signature.json',signature)
                write(d/'context_preparation.json',audit);matrix={z:{} for z in METALS};maps={};checks={}
                for q in ORDER:
                    paired={};mapping=copy.deepcopy(audit['mapping']);oldcells={z:xyz(verify(req['record']['matrix'][z][q]['xyz'])) for z in METALS}
                    assert oldcells['Ca'][1:]==oldcells['La'][1:]
                    for z in METALS:
                        new=[oldcells[z][0]]+list(origin[z][1:]);shared=[];newcaps=[]
                        for a in mapping['source_to_qm']:
                            i=a['qm_index'];prior=oldmap.get(atom_id(a))
                            if prior is not None:
                                new[i]=oldcells[z][prior['qm_index']];shared.append(i)
                            elif a['kind']=='source':
                                k=source_id(a['source']);new[i]=(a['source']['element'],*full[q][k].tolist())
                            elif a['kind']=='sigma_link_H':
                                x=full[q][source_id(a['retained'])];y=full[q][source_id(a['omitted'])]
                                pos=x+a['length_A']*(y-x)/np.linalg.norm(y-x);new[i]=('H',*pos.tolist());newcaps.append(i)
                            else:raise ValueError('new unmatched cofactor atom')
                        assert all(atom_id(a) in {atom_id(v) for v in mapping['source_to_qm']} for a in old['mapping']['source_to_qm'] if a['kind']!='sigma_link_H')
                        check=check_atoms(new,charges[z]);assert set(a[0] for a in new)<=set(VALENCE)
                        electrons=sum(VALENCE[a[0]] for a in new)-charges[z];assert electrons%2==0
                        pos=np.array([a[1:] for a in new]);dist=np.linalg.norm(pos[:,None]-pos[None,:],axis=2);np.fill_diagonal(dist,np.inf)
                        assert float(dist.min())>=settings['minimum_distance_A']
                        # Neighbor occupancy remains in the full source, outside this fixed local model.
                        heavypos=pos[[i for i,a in enumerate(new) if a[0]!='H']]
                        # The exact spectator identity is recorded by the original admission.
                        audit_rows=read(req['prepared_source']['source']['source_geometry_audit'])['rows']
                        sel=next(v['spectator_Ca'] for v in audit_rows if v['new_source_id']==req['case_id'])
                        hits=[k for k,meta in g.meta.items() if (meta['chain'],meta['resnum'],meta['atom'])==(sel['chain'],sel['resnum'],sel['atom'])]
                        assert len(hits)==1;minimum=float(np.min(np.linalg.norm(heavypos-full[q][hits[0]],axis=1)))
                        assert minimum>settings['spectator_boundary_A']
                        xp=d/(z+'_'+q+'.xyz');write_xyz(xp,new)
                        matrix[z][q]={'xyz':pin(xp),'charge':charges[z],'multiplicity':1,'uhf':0,'atom_count':len(new),
                            'valence_electrons':electrons,'all_electron_count':check['all_electron_count'],'status':'prepared_uncomputed'}
                        paired[z]=new
                        checks[z+'_'+q]={'shared_real_and_surviving_caps_exact':True,'shared_atom_count':len(shared),
                            'new_caps':newcaps,'spectator_minimum_A':minimum,'minimum_allatom_distance_A':float(dist.min())}
                    assert paired['Ca'][1:]==paired['La'][1:]
                    for a in mapping['source_to_qm']:a['xyz_A']=list(paired['La'][a['qm_index']][1:])
                    mapping['candidate']=q;write(d/(q+'_mapping.json'),mapping);maps[q]=pin(d/(q+'_mapping.json'))
                    bond_error=max((abs(float(np.linalg.norm(full[q][source_id(b['a'])]-full[q][source_id(b['b'])]))-
                                        float(np.linalg.norm(full['origin'][source_id(b['a'])]-full['origin'][source_id(b['b'])]))) for b in mapping['retained_bonds']),default=0.)
                    assert bond_error<=settings['bond_preservation_A'];checks[q+'_bond_preservation_A']=bond_error
                composition=[[a[0] for a in xyz(verify(matrix['Ca'][q]['xyz']))] for q in ORDER];assert composition[0]==composition[1]==composition[2]
                row.update(status='prepared',matrix=matrix,maps=maps,state_signature=pin(d/'state_signature.json'),state_key=cache_key(signature),
                    preparation=pin(d/'context_preparation.json'),source_preparation=pin(d/'context_preparation.json'),charge_ledger=audit['source_formal_charge_ledger'],checks=checks,
                    original_replay_max_A=old_error,kinematic_replay_max_A=replays)
            except Exception as e:row.update(reason=str(e),traceback=traceback.format_exc())
            rows.append(row);write(d/'PREPARATION.json',row)
        admitted=[r for r in rows if r['status']=='prepared']
        if len(admitted)==2:assert admitted[0]['state_key']==admitted[1]['state_key'],'shared fragment/state identity mismatch'
        result={'protocol':m['protocol'],'manifest':pin(path),'status':'prepared' if len(admitted)==2 else 'partial_or_unavailable',
          'rows':rows,'union':pin(out/'UNION.json'),'neutral_sidechain_template_checks':neutral_checks,
          'declared_sources':2,'declared_new_cells':12,'new_energy_calls':0,'new_optimizations':0,'same_state_across_sources':len(admitted)==2}
        write(out/'PREPARATION.json',result);receipt.update(status=result['status'],result=pin(out/'PREPARATION.json'))
    except Exception as e:receipt.update(status='failed',reason=str(e),traceback=traceback.format_exc());raise
    finally:
        u=resource.getrusage(resource.RUSAGE_SELF);receipt.update(wall_seconds=time.monotonic()-start,CPU_seconds=u.ru_utime+u.ru_stime,peak_RSS_KiB=u.ru_maxrss)
        write(out/'EXECUTION.json',receipt)

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='command',required=True)
    f=s.add_parser('freeze');f.add_argument('--review',required=True);f.add_argument('--output',required=True)
    f=s.add_parser('execute');f.add_argument('--manifest',required=True);a=p.parse_args()
    freeze(a.review,a.output) if a.command=='freeze' else execute(a.manifest)
