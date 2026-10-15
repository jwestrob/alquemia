"""Slurm-only checks on the five actual stationary-axis cap fixtures; no energies."""
import hashlib, itertools, json, os, sys
from pathlib import Path

PARENT_SHA='be3356652accc7ce0d413180cd05eb465e5a5b716106b1adcaac107728ade0ec'
REBASE_SHA='d896d60f2c4a3e0ba685c83785672c3d9c61f621f5e8a9269450b6b0821d2c4d'
FIXTURES={
 'p15279-pqq-la_model__conditioned_Ca__seed-1_sample-4':('context_8e39d5d61eeabb1222d871f1','1cfb8cae57dc8c355200892b42852255d0ce6fb4d9cf0f086e16baef680ad5db'),
 'p15279-pqq-la_model__conditioned_Ca__seed-1_sample-3':('context_922a5cee4869dd437a14d56f','100ca35a9dc17a66c169b5f40d2137add424fe31be1c9bcff10303c7610f98d5'),
 'atq70401.1-pqq-la_model__conditioned_Ca__seed-1_sample-2':('context_a2a0096c29ff3a1c6f246e82','227afab5b50b1bd339e9ecda96e754c7774be4850c4910249de6bc85f201369b'),
 'p38539-pqq-la_model__conditioned_Ca__seed-1_sample-0':('context_f18dd62cf49978f15c75f7f1','6daa19766108c86dd4a11d3d86264c4087143434c0e7c178a3724b859decfaf4'),
 'p16027-pqq-la_model__conditioned_Ca__seed-1_sample-0':('context_fabc02c52734b6bff84b9728','6c7840ca07e2ddfe87a1ae78c63362e0a3f3e1f208c4e3f6c18c4bdb3fd73f68')}

def require(value,reason):
    if not value:raise ValueError(reason)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def check(p):require(pin(p['path'])==p,'Changed fixture: '+p['path']);return Path(p['path'])
def read(p):return json.loads(check(p).read_text())

def check_row(row,np,Kinematics,angular,generic):
    """Call after prepare_row in the same Slurm preparation; caller saves this receipt."""
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before numerical fixture checks')
    sid=row['source_case_id']
    if sid not in FIXTURES:return {'status':'not_applicable','source_case_id':sid,'reason':'Not one of the five actual stationary-axis failures','new_energy_calls':0}
    cid,sha=FIXTURES[sid];require(row['context_id']==cid and row['status']=='prepared','Wrong actual fixture/context')
    original=row['original_row']['maps'];require(all(original[z]['sha256']==sha for z in ('Ca','La')),'Original fixture mapping changed')
    old=read(original['Ca']);require(read(original['La'])==old,'Original paired mapping differs')
    kin=Kinematics(read(row['mapping'])['context']);tasks=row['task_templates'];require(len(tasks)==2 and {t['metal'] for t in tasks}=={'Ca','La'},'Actual task pair missing')
    task=next(t for t in tasks if t['metal']=='Ca');other=next(t for t in tasks if t['metal']=='La');active=task['active_indices']
    require(len(active)==4 and all(task[k]==other[k] for k in ('active_indices','active_mode_ids','mode_count','mapping','contact')),'Actual shared four-mode pair differs')
    require(task['mapping']==row['mapping'] and task['contact']==row['contact'],'Task mapping/contact differs')
    require(kin.modes==old['context']['modes'] and kin.links==old['context']['core_links'],'Physical modes/cap links changed')
    require([kin.modes[i]['id'] for i in active]==task['active_mode_ids'],'Selected mode identity differs')
    cap_links=[v for v in kin.links if v[1]=='cap'];moving={i for k in active for i in kin.modes[k]['moving_indices']}
    affected=[v for v in cap_links if v[2] in moving or v[3] in moving]
    require(len(affected)==1 and affected[0][:5]==[140,'cap',117,116,1.09],'Actual offending cap fixture differs')
    proof=[]
    for link in cap_links:
        for anchor in link[2:4]:
            for k in active:
                mode=kin.modes[k];listed=anchor in mode['moving_indices'];axis=mode['unit']=='radian' and anchor in mode['axis_indices']
                require(not listed or axis,'Genuinely moving cap anchor is still unsupported')
                if listed:proof.append({'cap_qm_index':link[0],'physical_anchor_index':anchor,'source':kin.data['source_atom_metadata'][anchor],
                                       'mode_id':mode['id'],'axis_indices':mode['axis_indices'],'stationary_axis_endpoint':True})
    zero=np.zeros(len(kin.modes));p0,c0,_,_=kin.evaluate(zero);caps=sorted(v[0] for v in cap_links);anchors=sorted({i for v in cap_links for i in v[2:4]})
    contact=row['contact'];require(caps==contact['caps'],'Declared cap inventory changed')
    with np.load(check(contact['origin']),allow_pickle=False) as data:union=data['positions'].copy();symbols=data['symbols'].tolist()
    ids=read(contact['identities'])['centers'];n=row['atom_count']
    require(np.max(np.abs(c0-union[:n]))<=1e-12,'Actual rebased q0 differs from source-contact input')
    bound=angular.SETTINGS['angular_bound_radian'];require(bound==.8,'Inherited angular bound differs')
    points=[('zero',[0.]*4)]
    for i in range(4):
        for sign in (-1,1):
            q=[0.]*4;q[i]=sign*bound;points.append(('individual_'+str(i)+'_'+str(sign),q))
    points.extend(('coupled_'+''.join('p' if s>0 else 'm' for s in signs),[s*bound for s in signs]) for signs in itertools.product((-1,1),repeat=4))
    trials=[]
    for name,q in points:
        full=angular.full_q(task,np.asarray(q));p,core,_,_=kin.evaluate(full)
        anchor_error=float(np.max(np.abs(p[anchors]-p0[anchors])));cap_error=float(np.max(np.abs(core[caps]-c0[caps])))
        require(anchor_error<=1e-12 and cap_error<=1e-12,'Stationary axis/cap invariance failed at '+name)
        local=angular.geometry_check(kin,full,symbols[:n]);physical=angular.physical_status(kin,task,np.asarray(q))
        require(local['jacobian_max_error_A_per_unit']<=1e-7 and local['bond_max_error_A']<=1e-9 and local['fixed_atoms_unchanged'],'Original Jacobian/bond/fixed-atom checks failed at '+name)
        changed=union.copy();changed[:n]=core;exterior=generic.contact_check(np,union,changed,symbols,contact['moving_indices'],ids,set(caps),n)
        admitted=bool(local['pass'] and physical['physical_feasible'] and exterior['pass']);final=None
        if admitted:final=angular.final_geometry(kin,task,np.asarray(q),symbols[:n]);require(final['physical_feasible'],'Original final admission failed at '+name)
        if name=='zero':require(admitted,'Actual rebased q0 is not physically admissible')
        trials.append({'name':name,'active_q_radian':q,'anchor_max_error_A':anchor_error,'cap_max_error_A':cap_error,
           'original_local_geometry':local,'physical_status':physical,'original_final_geometry':final,
           'full_source_contact':{k:exterior[k] for k in ('pass','new_severe_contacts','policy','nearest_omitted_source_contacts')},
           'physical_admission_pass':admitted,'energy_evaluation_performed':False})
    return {'status':'passed','source_case_id':sid,'context_id':cid,'implementation':pin(__file__),'original_maps':original,'rebased_mapping':row['mapping'],
       'parent_rebase_sha256':PARENT_SHA,'axis_aware_rebase_sha256':REBASE_SHA,'kinematics':pin(sys.modules[Kinematics.__module__].__file__),
       'angular_geometry':pin(angular.__file__),'full_source_contact_helper':pin(generic.__file__),'source_contact':contact,
       'cap_anchor_proof':proof,'trials':trials,'declared_trials':25,'cap_and_anchor_tolerance_A':1e-12,
       'jacobian_tolerance_A_per_unit':1e-7,'bond_tolerance_A':1e-9,'unchanged_angular_bound_radian':bound,
       'interpretation':'Cap invariance qualification only. Trial points failing original physical/contact gates remain inadmissible; none is an energy candidate.',
       'new_energy_calls':0,'new_optimizations':0}
