"""Frozen author-domain multisite challenge; reuses the prepared MACE interface."""
from __future__ import annotations
import argparse
import importlib.metadata
import json
from pathlib import Path
import random
import resource
import time
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from mace_file_checks import cached_file_checks
from mace_omol_multisite import AUTHOR_POLICY, prepare as prepare_site
from mace_omol_prepared import prepare as prepare_score, collect as collect_score, validate

CASES = {
 'A0A7': {'file':'s001','sites':6,'La_Kd_uM':17.,'La_CI95_uM':2.,'Ca_CD_proxy_uM':2000.,
          'accession':'A0A7L4YJY0',
          'sequence':'STSEDQYYDTNYDGQVDTVVTDTDGNGVYDAAVYDTDGNGVADTVAYDSDENGVVDTVGFDYNEDGVVDEVVTDYNEDGYADSSSSS'},
 'HEW5': {'file':'s003','sites':8,'La_Kd_uM':5.2,'La_CI95_uM':1.7,'Ca_CD_proxy_uM':750.,
          'accession':'A0A6P0HEW5',
          'sequence':'MGSGPSSTEYDADGDGYVDTRESDTDGDGYVDTIETDTDGDGWVDTVATDTDGDGYIDTVATDTDGDGYADVVETDTDGDGYTDEVAYDADGDGYIDTVEADTDGDGYTDTVVHDGS'},
 'RTX': {'file':'s008','sites':8,'La_Kd_uM':40.,'La_CI95_uM':4.,'Ca_CD_proxy_uM':250.,
         'accession':None,
         'sequence':'GSARDDVLIGDAGANVLNGLAGNDVLSGGAGDDVLLGDEGSDLLSGDAGNDDLFGGQGDDTYLFGVGYGHDTIYESGGGHDTIRINAGADQLWFARQGNDLEIRILGTDDALTVHDWYRDADHRVEIIHAANQAVDQAGIEKLVEAMAQYPD'}}


def structure_inventory(path):
    import gemmi
    s = gemmi.read_structure(str(path))
    if len(s) != 1: raise InvalidArtifact('exactly one author model required')
    seq=''.join(gemmi.find_tabulated_residue(r.name).one_letter_code for c in s[0] if c.name=='A' for r in c)
    return {'source':record(path),'chain_A_sequence':seq,
            'residues':[{'chain':c.name,'resid':str(r.seqid),'name':r.name,'atoms':len(r)} for c in s[0] for r in c],
            'metals':[{'chain':c.name,'resid':str(r.seqid),'atom':a.name,'xyz_A':list(a.pos)}
                      for c in s[0] for r in c if r.name=='CA' for a in r],
            'geometry_origin':'author_AF3_Ca_conditioned_prediction'}


def normalize(source, cfg, out):
    inv=structure_inventory(source)
    if inv['chain_A_sequence'] != cfg['sequence'] or len(inv['metals']) != cfg['sites']:
        raise InvalidArtifact('author sequence or metal count differs from declared domain')
    for r in inv['residues']:
        if r['chain']!='A' and (r['name']!='CA' or r['atoms']!=1):
            raise InvalidArtifact('unsupported author assembly/cofactor')
    metal_map={(r['chain'],r['resid'],r['atom']):len(cfg['sequence'])+i+1 for i,r in enumerate(inv['metals'])}
    mapping={}; lines=[]; sites=[]
    for line in Path(source).read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            chain=line[21]; resid=int(line[22:26]); ic=line[26].strip(); atom=line[12:16].strip()
            if line[16].strip() or ic: raise InvalidArtifact('unsupported author alternate/insertion')
            old=f'{chain}/{resid}//{atom}'
            if chain!='A':
                newnum=metal_map[(chain,str(resid),atom)]
                line=line[:21]+'A'+f'{newnum:4d}'+line[26:]
            else: newnum=resid
            mapping[old]=f'A/{newnum}//{atom}'
        lines.append(line)
    dest=out/'normalized.pdb'; dest.write_text('\n'.join(lines)+'\n')
    receipt={'policy_id':AUTHOR_POLICY,'original':record(source),'normalized':record(dest),
             'source_to_normalized':mapping,'source_heavy_coordinates_changed':False,
             'assembly_changed':False,'chain_A_sequence':cfg['sequence'],
             'experimental_construct_sequence':'MPVP'+cfg['sequence'],
             'construct_difference':'experimental_N_terminal_MPVP_scar_absent_from_author_model'}
    write_new(out/'normalization.json',receipt)
    for row in inv['metals']:
        sites.append({'site_key':row['chain'],'metal_selector':{'chain':'A','resnum':metal_map[(row['chain'],row['resid'],row['atom'])],
                      'icode':'','resname':'CA','atom':row['atom']},
                      'source_metal_xyz_A':row['xyz_A'],'first_shell_waters':[],
                      'original_metal_selector':row})
    return dest, sites


def protonate(source, output):
    """One fixed native hydrogen preparation; never rebuild internal heavy atoms."""
    import openmm as mm
    from openmm import app,unit
    from mace_global_prepare import FF, STANDARD, atom_id, terminal_position
    pdb=app.PDBFile(str(source)); model=app.Modeller(pdb.topology,pdb.positions)
    metal_atoms=[a for a in model.topology.atoms() if a.residue.name=='CA']
    original={atom_id(a):(a.element.symbol,list(pdb.positions[a.index].value_in_unit(unit.angstrom)))
              for a in pdb.topology.atoms()}
    if any(a.element.symbol=='H' for a in pdb.topology.atoms()): raise InvalidArtifact('author source unexpectedly contains H')
    model.delete(metal_atoms)
    residues=list(model.topology.residues())
    if any(r.name not in STANDARD for r in residues): raise InvalidArtifact('nonstandard author protein')
    positions=np.array(model.positions.value_in_unit(unit.angstrom)); additions=[]
    for a,b in model.topology.bonds():
        if a.residue!=b.residue and {a.name,b.name}=={'C','N'}:
            distance=float(np.linalg.norm(positions[a.index]-positions[b.index]))
            if not 1. <= distance <= 1.8: raise InvalidArtifact('broken author peptide geometry')
    last=residues[-1]; named={a.name:a for a in last.atoms()}
    if 'OXT' not in named:
        if not {'C','CA','O'}<=named.keys(): raise InvalidArtifact('missing terminal heavy atoms')
        v=terminal_position(*(positions[named[k].index] for k in ('C','CA','O')))
        oxygen=model.topology.addAtom('OXT',app.element.oxygen,last); model.topology.addBond(named['C'],oxygen)
        model.positions=np.concatenate((positions,v[None,:]))*.1*unit.nanometer
        additions.append({'id':atom_id(oxygen),'xyz_A':v.tolist(),'recipe':'reflect_carbonyl_about_C_CA_v1'})
    random.seed(20260918); np.random.seed(20260918)
    ff=app.ForceField(str(FF))
    variants=model.addHydrogens(ff,pH=6.0,platform=mm.Platform.getPlatformByName('Reference'))
    # The original source has OXT in this pilot. If completion is needed, the
    # multisite protein builder repeats the declared exact terminal construction.
    if additions: model.delete([a for a in model.topology.atoms() if atom_id(a) in {v['id'] for v in additions}])
    with output.open('x') as h: app.PDBFile.writeFile(model.topology,model.positions,h,keepIds=True)
    text=output.read_text().splitlines(); text=[l for l in text if not l.startswith('END')]
    text.extend(l for l in Path(source).read_text().splitlines() if l.startswith('HETATM'))
    text.append('END'); output.write_text('\n'.join(text)+'\n')
    final=app.PDBFile(str(output)); actual={atom_id(a):(a.element.symbol,list(final.positions[a.index].value_in_unit(unit.angstrom)))
             for a in final.topology.atoms() if a.element.symbol!='H'}
    if set(original)!=set(actual) or any(v[0]!=actual[k][0] or not np.allclose(v[1],actual[k][1],rtol=0,atol=1e-9) for k,v in original.items()):
        raise InvalidArtifact('hydrogen preparation changed original heavy atoms')
    receipt={'protocol_id':'openmm_author_domain_pH6_fixed_seed_v1','source':record(source),'output':record(output),
             'ph':6.0,'add_missing_residues':False,'detected_missing_residues':{},'repaired_missing_atom_count':0,
             'terminal_completion':additions,'hydrogen_variants':variants,'random_seed':20260918,
             'platform':'Reference','forcefield':record(FF),'openmm':importlib.metadata.version('openmm'),
             'source_heavy_coordinates_preserved':True,'water_inventory':[]}
    dest=output.with_suffix('.json');write_new(dest,receipt);return dest


@cached_file_checks
def prepare(sources,development_collection,factorization,agreement,output):
    start=time.monotonic(); out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    source=Path(sources).resolve()
    write_new(out/'inventory.json',{'structures':[structure_inventory(p) for p in sorted(source.glob('*.pdb'))],
               'supplement':record(source/'SC-016-D5SC02315G-s007.pdf'),
               'supplementary_download':record(source/'supplementary_download.json'),
               'selection':'A0A7_HEW5_RTX_all_sites;other_four_inventory_only',
               'exposure':'identifiers absent from prior scored result/diagnostic records; labels and models inspected before new scoring'})
    proteins=[]; receipts={}; failures=[]
    for name,cfg in CASES.items():
        d=out/name;d.mkdir()
        try:
            normalized,sites=normalize(source/('SC-016-D5SC02315G-'+cfg['file']+'.pdb'),cfg,d)
            pm=protonate(normalized,d/'protonated.pdb')
            proteins.append({'panel_key':name,'source':record(normalized),'sites':sites})
            receipts[name]=pm
        except Exception as exc:
            failures.append({'domain':name,'stage':'source_preparation','reason':str(exc),'exception':type(exc).__name__})
    write_new(out/'panel.json',{'proteins':proteins,'agreement':record(agreement)})
    rows=[]
    for domain in proteins:
        name=domain['panel_key'];cfg=CASES[name];d=out/name
        evidence={'evidence_stratum':'qualified_protein_level_La_ITC_vs_Ca_CD_folding_proxy',
             'doi':'10.1039/D5SC02315G','table':'S6','group':name,
             'La_ITC_Kd_uM':cfg['La_Kd_uM'],'La_ITC_CI95_halfwidth_uM':cfg['La_CI95_uM'],
             'Ca_CD_folding_concentration_uM':cfg['Ca_CD_proxy_uM'],'Ca_fitted_Kd':None,
             'direction':'qualified_La_like_protein_level','site_specific_affinity_label':None,
             'assay_pH':6.,'geometry':'author_Ca_conditioned_AF3_domain_prediction_missing_MPVP_scar'}
        recipe={'policy_id':AUTHOR_POLICY,'panel':record(out/'panel.json'),'panel_key':name,'case_prefix':'KHOURY_'+name,
             'protonation':record(receipts[name]),'site_order':[s['site_key'] for s in domain['sites']],
             'retain_source_acetyl':False,'retained_water_oxygen_ids':[],
             'author_normalization':record(d/'normalization.json'),'evidence':evidence}
        write_new(d/'recipe.json',recipe)
        for site in recipe['site_order']:
            row={'domain':name,'site':site}
            try:
                p=prepare_site(d/'recipe.json',site,agreement,d/site/'prepared')
                m=prepare_score(verify(p['preparation']),development_collection,agreement,d/site/'score',factorization=factorization)
                row.update(status='prepared',preparation=p['preparation'],manifest=record(d/site/'score'/'manifest.json'),dry_run=m)
            except Exception as exc:
                row.update(status='unsupported',reason=str(exc),exception=type(exc).__name__)
            rows.append(row)
            print(json.dumps({k:row.get(k) for k in ('domain','site','status','reason')}),flush=True)
    result={'protocol_id':'khoury_author_domains_masked_MACE_all_sites_v1','agreement':record(agreement),
            'inventory':record(out/'inventory.json'),'rows':rows,'source_failures':failures,
            'declared_domains':list(CASES),'declared_sites':22,'declared_endpoint_forwards':44,
            'prepared_endpoint_forwards':sum(2 for r in rows if r['status']=='prepared'),
            'primary_summary':'unweighted_mean_all_declared_sites;unavailable_if_any_site_unavailable',
            'direction_screen_model_kcal':.02,'factorization':record(factorization),
            'development_collection':record(development_collection),'implementation':record(__file__),
            'wall_seconds':time.monotonic()-start,'process_CPU_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,
            'new_model_forwards':0,'baseline_changed':False}
    write_new(out/'manifest.json',result)
    (out/'manifests.txt').write_text(''.join(r['manifest']['path']+'\n' for r in rows if r['status']=='prepared'))
    return {'manifest':record(out/'manifest.json'),'prepared_endpoints':result['prepared_endpoint_forwards'],'failures':failures}


@cached_file_checks
def recover(manifest,output):
    """Retry failed preparation only, preserving failed artifacts and all inputs."""
    p=read_json(manifest);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    rows=[]
    for old in p['rows']:
        row=dict(old)
        if old['status']!='prepared':
            existing=Path(manifest).parent/old['domain']/old['site']/'score/manifest.json'
            if existing.exists():
                try:
                    checked=validate(existing)
                    row={'domain':old['domain'],'site':old['site'],'status':'prepared',
                         'preparation':read_json(existing)['preparation'],'manifest':record(existing),
                         'dry_run':checked,'prior_failure':old,'recovery':'validated_existing_preparation_no_new_preparation'}
                    rows.append(row);continue
                except InvalidArtifact:
                    pass
            recipe=Path(manifest).parent/old['domain']/'recipe.json'
            d=out/old['domain']/old['site']
            prepared=prepare_site(recipe,old['site'],verify(p['agreement']),d/'prepared')
            checked=prepare_score(verify(prepared['preparation']),verify(p['development_collection']),
                verify(p['agreement']),d/'score',factorization=verify(p['factorization']))
            row={'domain':old['domain'],'site':old['site'],'status':'prepared',
                 'preparation':prepared['preparation'],'manifest':record(d/'score/manifest.json'),
                 'dry_run':checked,'prior_failure':old}
        rows.append(row)
    p.update(rows=rows,parent_preparation=record(manifest),recovery_implementation=record(__file__),
             prepared_endpoint_forwards=sum(2 for r in rows if r['status']=='prepared'))
    write_new(out/'manifest.json',p)
    (out/'manifests.txt').write_text(''.join(r['manifest']['path']+'\n' for r in rows if r['status']=='prepared'))
    return {'manifest':record(out/'manifest.json'),'prepared_endpoints':p['prepared_endpoint_forwards']}


@cached_file_checks
def collect(manifest,ggr_development,ggr_2fw0,ggr_2fvy,output):
    p=read_json(manifest);verify(p['agreement']);rows=[];vectors={k:[] for k in CASES}
    for row in p['rows']:
        r={**row,'R_mask_model_kcal':None,'numerical_gate_pass':False}
        if row['status']=='prepared':
            mp=verify(row['manifest']);validate(mp);c=collect_score(mp)
            r.update(result=c,R_mask_model_kcal=c['R_mask_model_kcal'],numerical_gate_pass=c['numerical_gate_pass'],status=c['status'])
        rows.append(r);vectors[row['domain']].append({'site':row['site'],'score':r['R_mask_model_kcal'],'status':r['status']})
    means={k:float(np.mean([r['score'] for r in v])) if len(v)==CASES[k]['sites'] and all(r['score'] is not None for r in v) else None for k,v in vectors.items()}
    from mace_hybrid import EV_TO_KCAL
    calcium_checks=[]
    for domain,cfg in CASES.items():
        values=[r['result']['bound_model_eV']['Ca'] for r in rows if r['domain']==domain
                and r.get('result',{}).get('bound_model_eV',{}).get('Ca') is not None]
        spread=(max(values)-min(values))*EV_TO_KCAL if len(values)==cfg['sites'] else None
        calcium_checks.append({'domain':domain,'all_Ca_energy_spread_model_kcal':spread,
                               'tolerance_model_kcal':.01,'pass':spread is not None and spread<=.01})
    dev=read_json(ggr_development); qualified=read_json(verify(p['development_collection']))
    if dev['manifest']!=qualified['manifest'] or not dev['numerical_gate_pass']:
        raise InvalidArtifact('GGR development comparator belongs to another qualification')
    ggr={'GGR_1GLG':dev['scores']['GGR_1GLG']['R_mask_model_kcal']}
    for case,path in (('GGR_2FW0',ggr_2fw0),('GGR_2FVY',ggr_2fvy)):
        old=read_json(path)
        if (old['protocol_id']!=dev['protocol_id'] or old['status']!='complete' or not old['numerical_gate_pass']
                or old['factorization_reference']!=p['factorization']
                or read_json(verify(old['preparation']))['case_id']!=case):
            raise InvalidArtifact('GGR comparator method/reference/preparation mismatch')
        ggr[case]=old['R_mask_model_kcal']
    comparisons=[{'domain':k,'comparator':g,'margin_model_kcal':None if value is None else value-gv,
                  'direction_pass':None if value is None else value-gv>.02} for k,value in means.items() for g,gv in ggr.items()]
    result={'manifest':record(manifest),'rows':rows,'ordered_site_vectors':vectors,'primary_domain_means':means,
          'ggr_scores':ggr,'ggr_sources':[record(x) for x in (ggr_development,ggr_2fw0,ggr_2fvy)],
          'comparisons':comparisons,'passed_comparisons':sum(r['direction_pass'] is True for r in comparisons),
          'comparison_denominator':9,'independent_new_domain_groups':3,'new_site_denominator':22,
          'all_Ca_permutation_checks':calcium_checks,
          'numerical_gate_pass':len(rows)==22 and all(r['numerical_gate_pass'] for r in rows) and all(r['pass'] for r in calcium_checks),
          'status':'complete' if all(v is not None for v in means.values()) else 'incomplete',
          'evidence_stratum':'qualified_protein_level_La_ITC_vs_Ca_CD_proxy;author_domain_construct_mismatch',
          'baseline_changed':False,'calibration':None,'broad_affinity_validation':False,'report_implementation':record(__file__)}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);write_new(out/'result.json',result)
    return {'status':result['status'],'means':means,'passed_comparisons':result['passed_comparisons'],'result':record(out/'result.json')}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);sp=ap.add_subparsers(dest='operation',required=True)
    p=sp.add_parser('prepare')
    for name in ('sources','development-collection','factorization','agreement','output'):p.add_argument('--'+name,required=True)
    p=sp.add_parser('collect')
    for name in ('manifest','ggr-development','ggr-2fw0','ggr-2fvy','output'):p.add_argument('--'+name,required=True)
    p=sp.add_parser('recover')
    for name in ('manifest','output'):p.add_argument('--'+name,required=True)
    args=vars(ap.parse_args());operation=args.pop('operation');print(json.dumps(globals()[operation](**args),indent=2))
