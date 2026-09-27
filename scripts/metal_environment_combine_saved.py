"""Combine strictly matched saved native/classical components, no molecular calls."""
from pathlib import Path
import argparse,json
import numpy as np
from affordable_common import read_json,record,verify,write_new,HA_TO_KCAL,BOHR_TO_A,InvalidArtifact
from affordable_response import cap_jacobians
from metal_environment_reference import read_pcgrad
from metal_environment_force_checks import xyz_data,pc_data,task_id,rotation
from mace_omol_vacuum import parse_endpoint
GROUPS=('retained_bonded','MM_LJ','MM_Coulomb','QM_MM_LJ')


def matches(config_name,x,task,particles,caps,env_atoms):
    _,symbols,q=xyz_data(verify(task['xyz']));_,charges,pc=pc_data(verify(task['pointcharges']));ids={p['id']:i for i,p in enumerate(particles)};errors=[];seen=set()
    for i,p in enumerate(particles):
        qi=p.get('qm_index')
        if qi is None:continue
        if symbols[qi]!=(task['metal'] if p['id']=='metal' else p['element']):raise InvalidArtifact('real QM element mismatch')
        errors.append(float(abs(q[qi]-x[i]).max()));seen.add(qi)
    for cap in caps:
        ia,ib=[ids[cap[k]] for k in ('retained_source_id','omitted_source_id')];ra,rb=x[ia],x[ib];l=cap['length_A'];orig_a=np.array(cap['retained_xyz_A']);orig_b=np.array(cap['omitted_xyz_A']);offset=np.array(cap['xyz_A'])-(orig_a+l*(orig_b-orig_a)/np.linalg.norm(orig_b-orig_a))
        if config_name=='rigid':offset=offset@rotation().T
        expected=ra+l*(rb-ra)/np.linalg.norm(rb-ra)+offset;errors.append(float(abs(expected-q[cap['qm_index']]).max()));seen.add(cap['qm_index'])
    if seen!=set(range(len(q))):raise InvalidArtifact('unmapped or duplicate electronic atom')
    if len(env_atoms)!=len(pc):raise InvalidArtifact('external atom count mismatch')
    for i,p in enumerate(env_atoms):
        pi=ids[p['id']]
        if charges[i]!=p['charge_e'] or particles[pi]['charge_e']!=charges[i]:raise InvalidArtifact('permanent charge mismatch')
        errors.append(float(abs(pc[i]-x[pi]).max()))
    error=max(errors)
    if error>1e-11:raise InvalidArtifact('unmatched molecular coordinates: '+str(error))
    return error


def physical_gradient(x,qg,pg,particles,caps,env_atoms):
    ids={p['id']:i for i,p in enumerate(particles)};g=np.zeros_like(x)
    for i,p in enumerate(particles):
        if p.get('qm_index') is not None:g[i]+=qg[p['qm_index']]
    for c in caps:
        ia,ib=[ids[c[k]] for k in ('retained_source_id','omitted_source_id')];ja,jb=cap_jacobians(x[ia],x[ib],c['length_A']);g[ia]+=ja.T@qg[c['qm_index']];g[ib]+=jb.T@qg[c['qm_index']]
    for i,p in enumerate(env_atoms):g[ids[p['id']]]+=pg[i]
    return g


def combine(classical_path,scout_path,force_path):
    from ggr_sensitivity import executed
    classical=read_json(classical_path);cm=read_json(verify(classical['manifest']));ledger=read_json(verify(cm['ledger']));config=read_json(verify(ledger['inputs']['source_manifest']));particles=read_json(verify(ledger['artifacts']['particles.json']));caps=read_json(verify(ledger['artifacts']['caps.json']));env=read_json(verify(config['environments']['A']['atoms']));ct={t['task_id']:t for t in cm['tasks']}
    if classical['completed']!=32 or not all(c['passed'] for c in classical['checks']):raise InvalidArtifact('complete passed classical matrix required')
    scouts=read_json(scout_path);forces=read_json(force_path);sm,sactual=executed(verify(scouts['manifest']));fm,factual=executed(verify(forces['manifest']))
    if sm['inputs']!=fm['inputs'] or sm['inputs']!=ledger['inputs']['source_manifest']:raise InvalidArtifact('scientific source preparation mismatch')
    if sm['method']!=fm['method']:raise InvalidArtifact('electronic methods differ')
    st={t['task_id']:t for t in sm['tasks']};ft={t['task_id']:t for t in fm['tasks']};rows={};qgrads={};cforces={};maxerror=0.
    pairs=[]
    for metal in ('Ca','La'):
        for name in ('A','B'):pairs.append((metal,name,scouts['rows'][metal+'_'+name],sactual[metal+'_'+name],st[metal+'_'+name]))
        for mode in ('repeat','rigid'):tid=task_id(metal,mode,0);pairs.append((metal,mode,forces['rows'][tid],factual[tid],ft[tid]))
        for mode in ('hydroxyl','boundary'):
            for h in (-.001,.001,-.0005,.0005):
                tid=task_id(metal,'mm' if mode=='hydroxyl' else mode,h);pairs.append((metal,f'{mode}_{h:+.4f}',forces['rows'][tid],factual[tid],ft[tid]))
    for metal,name,qr,observed,qt in pairs:
        key=metal+'_'+name;cr=classical['rows'][key];receipt=Path(verify(classical['manifest'])).parent/'results'/key/'receipt.json'
        if read_json(receipt)!=cr or cr['manifest']!=classical['manifest'] or cr['model']!=cm['models'][metal]:raise InvalidArtifact('classical receipt/model mismatch')
        verify(cr['model'])
        if qr['status']!='complete' or observed['status']!='complete' or cr['status']!='complete':raise InvalidArtifact('required shared cell unavailable')
        if any(qr[k]!=observed[k] for k in ('energy_hartree','output','receipt')):raise InvalidArtifact('electronic receipt mismatch')
        parsed=parse_endpoint(qt,verify(qr['output']),verify(qr['engrad']),permanent_field=True)
        if parsed['energy_hartree']!=qr['energy_hartree']:raise InvalidArtifact('electronic energy mismatch')
        x=np.load(verify(ct[key]['coordinates']));err=matches(name,x,qt,particles,caps,env);maxerror=max(maxerror,err)
        qg=np.array(parsed['gradient_kcal_mol_per_A']);pg=read_pcgrad(verify(qr['pointcharge_gradient']),len(env))*HA_TO_KCAL/BOHR_TO_A
        qgrads[key]=physical_gradient(x,qg,pg,particles,caps,env);ff=np.load(verify(cr['forces']))['forces_kcal_mol_A']
        if ff.shape!=(4,len(particles),3) or not np.isfinite(ff).all():raise InvalidArtifact('classical force dimensions/values mismatch')
        cforces[key]=ff.sum(axis=0)
        rows[key]={'quantum_energy_Hartree':qr['energy_hartree'],'classical_energies_kcal_mol':cr['energies_kcal_mol'],'classical_sum_kcal_mol':sum(cr['energies_kcal_mol'][k] for k in GROUPS),'quantum_output':qr['output'],'quantum_receipt':qr['receipt'],'quantum_engrad':qr['engrad'],'quantum_pcgrad':qr['pointcharge_gradient'],'classical_receipt':record(receipt),'classical_forces':cr['forces'],'coordinates':ct[key]['coordinates'],'coordinate_max_residual_A':err}
    def work(a,b):
        return {'electronic':(rows[b]['quantum_energy_Hartree']-rows[a]['quantum_energy_Hartree'])*HA_TO_KCAL,**{g:rows[b]['classical_energies_kcal_mol'][g]-rows[a]['classical_energies_kcal_mol'][g] for g in GROUPS}}
    responses={metal:work(metal+'_A',metal+'_B') for metal in ('Ca','La')}
    for values in responses.values():values['combined']=sum(values.values())
    delta={k:responses['La'][k]-responses['Ca'][k] for k in responses['Ca']};fd=[];repeat_rigid=[];missing=[]
    for metal in ('Ca','La'):
        base=metal+'_A';g=qgrads[base]-cforces[base]
        for mode in ('hydroxyl','boundary'):
            spec=cm['modes'][mode];i=spec['physical_index'];v=np.array(spec['direction']);anq=float(qgrads[base][i]@v);anc=float(-cforces[base][i]@v);analytic=anq+anc
            for h in (.001,.0005):
                changes=work(f'{metal}_{mode}_{-h:+.4f}',f'{metal}_{mode}_{h:+.4f}');numeric=sum(changes.values())/(2*h)
                fd.append({'metal':metal,'mode':mode,'step_A':h,'analytic_electronic_kcal_mol_A':anq,'analytic_classical_kcal_mol_A':anc,'analytic_combined_kcal_mol_A':analytic,'FD_electronic_kcal_mol_A':changes['electronic']/(2*h),'FD_classical_kcal_mol_A':sum(changes[k] for k in GROUPS)/(2*h),'FD_combined_kcal_mol_A':numeric,'signed_residual_kcal_mol_A':numeric-analytic,'physical_index':i,'direction':v.tolist(),'qualification':'residual_only_original_native_gates_not_replaced'})
        for mode in ('repeat','rigid'):
            key=metal+'_'+mode;changes=work(base,key);expected=g@rotation().T if mode=='rigid' else g;residual=(qgrads[key]-cforces[key])-expected
            repeat_rigid.append({'metal':metal,'mode':mode,'energy_components_shift_kcal_mol':changes,'combined_energy_shift_kcal_mol':sum(changes.values()),'maximum_physical_combined_gradient_residual_kcal_mol_A':float(abs(residual).max()),'original_electronic_gate':next(c for c in forces['checks'] if c['metal']==metal and c['mode']==mode)})
        missing.append({'metal':metal,'mode':'metal','status':'unavailable_quantum_displacements_not_evaluated','classical_cells_available':4,'combined_finite_difference':None})
    return {'protocol_id':'nikasha_saved_dry_additive_components_comparison_v1','inputs':{'classical':record(classical_path),'scout':record(scout_path),'force_checks':record(force_path),'ledger':cm['ledger']},'implementation':record(__file__),'common_source':sm['inputs'],'matched_configurations':len(rows),'physical_atom_count':len(particles),'maximum_coordinate_correspondence_error_A':maxerror,'rows':rows,'per_metal_A_B_works_kcal_mol':responses,'La_minus_Ca_response_kcal_mol':delta,'directional_FD_residuals':fd,'repeat_rigid_residuals':repeat_rigid,'unavailable_comparisons':missing,'native_qualification_status':forces['status'],'full_hybrid_qualified':False,'solvent_present':False,'new_molecular_evaluations':0,'classification':None,'limitations':['Saved dry additive candidate only; native rigid qualification remains failed.','Protein MM strain and the capped-covalent boundary approximation remain unresolved.','Hydration-fit metal12-6 cross terms and GAFF2 PQQ LJ are declared candidates, not validated specificity models.','No quantum metal-displacement pairs; no full force validation, solvent, relaxation, affinity or classifier claim.']}


def main():
    p=argparse.ArgumentParser();p.add_argument('--classical',required=True);p.add_argument('--scout',required=True);p.add_argument('--force-checks',required=True);p.add_argument('--output',required=True);a=p.parse_args();r=combine(a.classical,a.scout,a.force_checks);write_new(a.output,r);print(json.dumps({k:v for k,v in r.items() if k not in ('rows','inputs')},indent=2))
if __name__=='__main__':main()
