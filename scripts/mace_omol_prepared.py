"""Opt-in charge-ablated descriptor for explicit, source-backed whole-chain inputs."""
from __future__ import annotations
import argparse
import copy
import importlib.metadata
import json
from pathlib import Path
import sys
import numpy as np
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import check_atoms,write_xyz,EV_TO_KCAL,accepted_attempt
from mace_global_prepare import protein,cofactor_and_waters,FF,WATER_FF,POLICY
from mace_omol_panel_prepare import ACTIVE_POLICY,POLICY as LEGACY_CANONICAL_POLICY
from mace_omol_ablation import ADAPTER,COMPONENT
from mace_omol_ablation_run import PROTOCOL,SEMANTICS,descriptor_model,validate as development_validate,collect as development_collect
from mace_omol_intact import EVALUATION,energy_task,geometry,collect as collect_endpoints
from mace_omol_panel import ADAPTER_TASK
from mace_omol_source_prepare import POLICY as SOURCE_POLICY, validate_bridge

ASSEMBLY='deposited_chain_A_with_declared_cofactor_and_site_waters'


@cached_file_checks
def audit_preparation(path):
    """Replay only the established preparation arithmetic; no energy calculation."""
    p=read_json(path)
    from mace_omol_multisite import POLICY as MULTISITE_POLICY, audit as audit_multisite
    if p.get('policy_id')==MULTISITE_POLICY:
        return audit_multisite(path)
    if p.get('status')!='prepared' or p.get('policy_id') not in (POLICY,ACTIVE_POLICY,LEGACY_CANONICAL_POLICY,SOURCE_POLICY):
        raise InvalidArtifact('unsupported whole-chain preparation protocol')
    if p['policy_id']==SOURCE_POLICY:validate_bridge(p)
    if p.get('assembly')!=ASSEMBLY or p.get('forcefield')!=record(FF) or p.get('water_forcefield')!=record(WATER_FF):
        raise InvalidArtifact('unsupported assembly or changed preparation forcefield')
    for key in ('source','source_preparation','forcefield','water_forcefield'):verify(p[key])
    row=p['source_audit_row']
    if row['source_structure']!=p['source'] or row['preparation_manifest']!=p['source_preparation']:
        raise InvalidArtifact('preparation provenance does not match source row')
    atoms,q,details=protein(row,allow_terminal_completion=True,check_peptide_connectivity=True)
    core=read_json(verify(p['source_preparation']))
    if core.get('pqq') and not row['case_id'].startswith('PQQ'):
        raise InvalidArtifact('recorded PQQ must not be omitted by a non-PQQ preparation identifier')
    extra,cq,waters,moves=cofactor_and_waters(row,core,expected_water_count=len(p['explicit_waters']))
    metal=xyz(verify(row['endpoints']['La']['xyz']))[0]
    if metal[0]!='La' or not np.allclose(metal[1:],row['selected_site']['xyz_A'],rtol=0,atol=1e-6):
        raise InvalidArtifact('selected metal identity/coordinate mapping differs')
    expected=[*atoms,*extra,{'id':'metal','element':'M','xyz_A':list(metal[1:]),'radius_A':1.8,'kind':'selected_metal'}]
    if (p['physical_atoms']!=expected or p['preparation_details']!=details
            or p['protein_charge_e']!=q or p['cofactor_charge_e']!=cq
            or p['explicit_waters']!=waters or p['water_H_moves']!=moves):
        raise InvalidArtifact('whole-chain preparation no longer matches exact source replay')
    if len({a['id'] for a in expected})!=len(expected) or any('cap' in a['kind'] for a in expected):
        raise InvalidArtifact('duplicate/synthetic fragment atom in whole-chain input')
    # Existing explicit ligand/water exclusions remain visible; do not reduce a
    # multisite source to the singleton model merely because a preparation did.
    metals=[r for r in details['source_residue_inventory'] if r['chain']=='A' and
            r['resname'] in ('CA','LA','CE','PR','ND','SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','Y')]
    if len(metals)!=1:raise InvalidArtifact('single metal-bearing chain required; multisite input unsupported')
    for m in ('La','Ca'):
        e=p['endpoints'][m];coords=xyz(verify(e['xyz']));wanted=[(m if a['element']=='M' else a['element'],*a['xyz_A']) for a in expected]
        charge=q+cq+(3 if m=='La' else 2)
        if (coords!=wanted or e['charge']!=charge or e['spin_multiplicity']!=1
                or e['state']!=check_atoms(coords,charge) or not -100<=charge<=100):
            raise InvalidArtifact('paired coordinates, physical charge or electron state differ')
    return {'status':'pass','preparation':record(path),'atoms':len(expected),'policy_id':p['policy_id'],
            'assembly':ASSEMBLY,'protein_charge_e':q,'cofactor_charge_e':cq,
            'explicit_water_count':len(waters),'recorded_source_exclusions':[
                r for r in details['source_residue_inventory'] if not r['selected_protein']],
            'scope':'exact replay of recorded preparation, including its declared ligand/water omissions',
            'new_model_forwards':0}


def qualification(path):
    c=read_json(path);mp=verify(c['manifest']);m=read_json(mp)
    if m.get('stage')!='ablation_development' or m.get('protocol_id')!=PROTOCOL:
        raise InvalidArtifact('qualified charge-feature descriptor required')
    development_validate(mp)
    if c!=development_collect(mp) or not c['canonical_extension_permitted']:
        raise InvalidArtifact('completed passing descriptor qualification required')
    return c,m


def tasks(preparation):
    p=read_json(preparation);result=[]
    for position in ('bound','detached'):
        for metal in ('La','Ca'):
            e=p['endpoints'][metal]
            result.append(energy_task({'task_id':f'{metal}_{position}_primary','case_id':p['case_id'],
                'metal':metal,'kind':'full','position':position,'variant':'primary',
                'metal_index':len(p['physical_atoms'])-1,'source_xyz':e['xyz'],'preparation':record(preparation),
                'assembly':p['assembly'],'microstate':p['microstate'],'explicit_waters':p['explicit_waters'],
                'evidence':p['evidence'],'energy_component':COMPONENT,'charge_feature_adapter':ADAPTER,
                'output_semantics':SEMANTICS,'edge_adapter':ADAPTER_TASK,
                'require_isolated_metal':position=='detached',**{k:v for k,v in e.items() if k!='xyz'}}))
    return result


def reuse(collection,wanted,parent):
    if collection is None:return {},{}
    saved=read_json(collection);mp=verify(saved['manifest']);source=read_json(mp)
    if (source.get('protocol_id')!=PROTOCOL or source.get('model')!=parent['model']
            or source.get('software')!=parent['software']):
        raise InvalidArtifact('reuse requires the identical qualified descriptor model/software')
    for ref in source['implementation'].values():verify(ref)
    refs={};values={}
    for t in wanted:
        matches=[s for s in source['tasks'] if s.get('preparation')==t['preparation'] and s['metal']==t['metal']
                 and s.get('position')==t['position'] and s.get('variant')=='primary' and s['kind']=='full']
        if len(matches)!=1:raise InvalidArtifact('explicit reuse lacks a unique compatible endpoint')
        old=matches[0]
        for key in ('source_xyz','preparation','charge','spin_multiplicity','state','metal_index','assembly','microstate',
                    'explicit_waters','energy_component','charge_feature_adapter','edge_adapter','energy_only','capture_native_readout'):
            if old[key]!=t[key]:raise InvalidArtifact('reuse changed descriptor physical/method field: '+key)
        if xyz(verify(old['xyz']))!=geometry(t):raise InvalidArtifact('reuse coordinates differ')
        r=saved['rows'][old['task_id']]
        if not any(accepted_attempt(a,old,mp)==r for a in sorted((mp.parent/'execution'/old['task_id']).glob('attempt_*'))):
            raise InvalidArtifact('explicit reuse lacks an actual matching execution receipt')
        refs[t['task_id']]={'source_collection':record(collection),'source_task_id':old['task_id']};values[t['task_id']]=r
    return refs,values


def scoring_tasks(preparation,factorization,parent):
    wanted=tasks(preparation)
    if factorization is not None:
        from mace_omol_factorization import verified
        ref=verified(factorization)
        if ref['model']!=parent['model'] or ref['software']!=parent['software']:
            raise InvalidArtifact('factorization reference belongs to another descriptor model/software')
        wanted=[t for t in wanted if t['position']=='bound']
    return wanted


@cached_file_checks
def prepare(preparation,development_collection,agreement,output,reuse_collection=None,factorization=None):
    from mace_omol import common,seal
    audit=audit_preparation(preparation);_,parent=qualification(development_collection)
    wanted=scoring_tasks(preparation,factorization,parent);refs,_=reuse(reuse_collection,wanted,parent)
    _,out,m=common(verify(parent['inventory']),verify(parent['software']),agreement,output,'ablation_prepared')
    for name in ('mace_omol_ablation.py','mace_omol_products.py','mace_omol_edges.py','mace_omol_readout.py'):
        if m['implementation'][name]['sha256']!=parent['implementation'][name]['sha256']:
            raise InvalidArtifact('descriptor adapter changed after qualification')
    m.update(protocol_id=PROTOCOL,model=descriptor_model(verify(parent['software'])),preparation=record(preparation),
             preparation_audit=audit,development_collection=record(development_collection),reused=refs,
             reuse_collection=record(reuse_collection) if reuse_collection else None,
             factorization=record(factorization) if factorization else None,
             energy_evaluation=EVALUATION,output_semantics=SEMANTICS,
             evidence_use=read_json(preparation).get('evidence_use','unrecorded'),
             driver_python=record(sys.executable),
             driver_versions={k:importlib.metadata.version(k) for k in ('openmm','numpy')},
             tasks=[t for t in wanted if t['task_id'] not in refs])
    for t in m['tasks']:
        if t['position']=='bound':t['xyz']=t['source_xyz']
        else:
            path=out/(t['task_id']+'.xyz');write_xyz(path,geometry(t));t['xyz']=record(path)
    return seal(out,m)


@cached_file_checks
def validate(manifest):
    from mace_omol import SCHEMA,TOL
    m=read_json(manifest);audit=audit_preparation(verify(m['preparation']));_,parent=qualification(verify(m['development_collection']))
    if (m['driver_python']!=record(sys.executable) or
            m['driver_versions']!={k:importlib.metadata.version(k) for k in ('openmm','numpy')}):
        raise InvalidArtifact('use the recorded preparation driver environment; MACE workers use their separate pinned environment')
    wanted=scoring_tasks(verify(m['preparation']),verify(m['factorization']) if m.get('factorization') else None,parent)
    refs,_=reuse(verify(m['reuse_collection']) if m['reuse_collection'] else None,wanted,parent)
    if (m['schema_version']!=SCHEMA or m['protocol_id']!=PROTOCOL or m['stage']!='ablation_prepared'
            or m['model']!=descriptor_model(verify(m['software'])) or m['software']!=parent['software']
            or m['inventory']!=parent['inventory'] or m['tolerances']!=TOL or m['energy_evaluation']!=EVALUATION
            or m['output_semantics']!=SEMANTICS or m['reused']!=refs or m['preparation_audit']!=audit):
        raise InvalidArtifact('prepared descriptor manifest differs from qualified method/source')
    for ref in [m['agreement'],*m['implementation'].values()]:verify(ref)
    for name in ('mace_omol_ablation.py','mace_omol_products.py','mace_omol_edges.py','mace_omol_readout.py'):
        if m['implementation'][name]['sha256']!=parent['implementation'][name]['sha256']:
            raise InvalidArtifact('prepared descriptor adapter changed')
    selected=[t for t in wanted if t['task_id'] not in refs]
    if len(selected)!=len(m['tasks']) or len(selected)+len(refs)!=(2 if m.get('factorization') else 4):
        raise InvalidArtifact('declared scoring state inventory changed')
    for t,e in zip(m['tasks'],selected):
        payload={k:v for k,v in t.items() if k not in ('xyz','cache_key')}
        if payload!=e or xyz(verify(t['xyz']))!=geometry(t):raise InvalidArtifact('prepared scoring state differs')
        if t['cache_key']!=cache_key({'task':{k:v for k,v in t.items() if k!='cache_key'},'model':m['model'],
                                    'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('prepared scoring cache differs')
    return {'status':'pass','tasks':len(selected),'reused_endpoints':len(refs),'manifest':record(manifest)}


@cached_file_checks
def collect(manifest):
    m=read_json(manifest);result=collect_endpoints(manifest)
    _,parent=qualification(verify(m['development_collection']))
    wanted=scoring_tasks(verify(m['preparation']),verify(m['factorization']) if m.get('factorization') else None,parent)
    refs,values=reuse(verify(m['reuse_collection']) if m['reuse_collection'] else None,wanted,parent)
    result['reused_rows']={k:{**ref,'result':values[k]} for k,ref in refs.items()}
    rows={**result['rows'],**values};terms={}
    if m.get('factorization'):
        from mace_omol_factorization import verified,FACTORIZATION
        ref=verified(verify(m['factorization']));bound={metal:rows[metal+'_bound_primary'] for metal in ('La','Ca')}
        complete=all(r['status']=='computed' for r in bound.values())
        value=(bound['Ca']['energy_eV']-bound['La']['energy_eV']-ref['Ca_minus_La_disconnected_atom_model_eV'])*EV_TO_KCAL if complete else None
        terms=None
        result.update(score_evaluation=FACTORIZATION,factorization_reference=m['factorization'],
                      bound_model_eV={metal:r['energy_eV'] for metal,r in bound.items()},
                      disconnected_atom_model_eV=ref['disconnected_atom_model_eV'])
    else:
        for metal in ('La','Ca'):
            pair={p:rows[f'{metal}_{p}_primary'] for p in ('bound','detached')}
            terms[metal]=pair['bound']['energy_eV']-pair['detached']['energy_eV'] if all(r['status']=='computed' for r in pair.values()) else None
        complete=all(v is not None for v in terms.values())
        value=(terms['Ca']-terms['La'])*EV_TO_KCAL if complete else None
        result['score_evaluation']='four_full_bound_detached_states_v1'
    checks=[{'name':key+'_accounting','error_model_kcal':r['native_readout']['component_sum_error_kcal_mol'],
             'pass':abs(r['native_readout']['component_sum_error_kcal_mol'])<=.01}
            for key,r in rows.items() if r['status']=='computed']
    result.update(status='complete' if complete else 'incomplete',protocol_id=PROTOCOL,output_semantics=SEMANTICS,
                  score_unit='kcal_equivalent_model_units',bound_minus_detached_model_eV=terms,
                  R_mask_model_kcal=value,
                  numerical_gate_pass=complete and all(c['pass'] for c in checks),checks=checks,
                  decision_status='compatible_calibration_not_supplied',calibrated_class=None,
                  baseline_changed=False,production_promotion=False)
    return result


@cached_file_checks
def report(manifest,output,calibration=None):
    validate(manifest);result=collect(manifest);m=read_json(manifest)
    if calibration is not None:
        from mace_omol_mask_calibration import decision
        result.update(decision(result,m,calibration))
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    result.update(preparation=m['preparation'],preparation_audit=m['preparation_audit'],
                  report_implementation=record(__file__))
    write_new(out/'result.json',result)
    (out/'REPORT.md').write_text('# Prepared whole-chain MACE descriptor\n\n'
        f"Status: {result['status']}. Score: {result['R_mask_model_kcal']} kcal-equivalent model units.\n\n"
        'Modified learned descriptor, not a quantum energy or binding free energy.\n'
        f"Decision: {result['calibrated_class']}; status: {result['decision_status']}. Baseline unchanged.\n")
    return {'status':result['status'],'result':record(out/'result.json')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('audit');a.add_argument('--preparation',required=True)
    a=sub.add_parser('prepare')
    for key in ('preparation','development-collection','agreement','output'):a.add_argument('--'+key,required=True)
    a.add_argument('--reuse-collection');a.add_argument('--factorization')
    a=sub.add_parser('report')
    for key in ('manifest','output'):a.add_argument('--'+key,required=True)
    a.add_argument('--calibration')
    args=vars(p.parse_args());command=args.pop('command')
    result=audit_preparation(args['preparation']) if command=='audit' else {'prepare':prepare,'report':report}[command](**args)
    print(json.dumps(result,indent=2))
