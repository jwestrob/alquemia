#!/usr/bin/env python3
"""Source/metadata inventory only. No numerical engine or parameter binding."""
import argparse,collections,hashlib,json,shlex,time
from pathlib import Path
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):assert pin(p['path'])==p;return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def parse(p):
    lines=verify(p).read_text().splitlines();out={k:{} for k in ('atom','vdw','polarize','multipole')};out['pairs']=[];out['water_terms']=[]
    for n,l in enumerate(lines):
        t=l.split()
        if not t:continue
        key=t[0].lower()
        if key=='atom':
            a=shlex.split(l);out['atom'][int(a[1])]={'class':int(a[2]),'label':a[3],'description':a[4],
                'atomic_number':int(a[5]),'mass':a[6],'valence':int(a[7]),'line':n+1,'text':l}
        elif key in ('vdw','polarize','multipole'):
            entry={'line':n+1,'tokens':t,'text':l}
            if key=='multipole':entry['text']='\n'.join(lines[n:n+5]);entry['components']=[s.split() for s in lines[n+1:n+5]]
            out[key].setdefault(int(t[1]),[]).append(entry)
        elif key in ('vdwpr','polpair'):out['pairs'].append({'line':n+1,'tokens':t,'text':l})
        elif key in ('bond','angle','strbnd','ureybrad'):
            out['water_terms'].append({'line':n+1,'tokens':t,'text':l})
    return out
def run(config,out):
    start=time.monotonic();c=json.loads(Path(config).read_text());assert c['implementation']==pin(__file__);verify(c['specification'])
    r=read(c['paired_preparation']);assert r['status']=='complete' and r['metal_parameters_selected'] is False
    m=read(r['atom_mapping']);bio=parse(r['protein_parameters']);small=parse(c['small_parameters']);key=parse(c['small_key'])
    for p in (r['parent_bound_result'],r['source_water']):verify(p)
    for group in ('bound_boxes','aqueous_boxes'):
        for p in r[group].values():verify(p)
    assert all(t not in bio['atom'] for t in (900,901))
    assert all(t not in {a['class'] for a in bio['atom'].values()} for t in (900,901))
    inventories={};classes={}
    for leg in ('bound','aqueous'):
        counts=collections.Counter(a['type'] for a in m[leg] if 'type'in a);inventories[leg]={}
        for typ,count in sorted(counts.items()):
            a=bio['atom'][typ];cl=a['class'];assert typ in bio['polarize'] and typ in bio['multipole'] and cl in bio['vdw']
            contexts=[x['source_join'] for x in m[leg] if x.get('type')==typ and x['role']=='peptide']
            inventories[leg][str(typ)]={'count':count,'atom':a,'vdw_class_records':bio['vdw'][cl],
                 'polarize_type_records':bio['polarize'][typ],'multipole_type_records':bio['multipole'][typ],
                 'peptide_contexts':contexts,'pair_policy':'native BIO mixing unless explicitly named in future qualified map'}
            classes.setdefault(str(cl),{'types':set(),'bound_count':0,'aqueous_count':0})['types'].add(typ)
            classes[str(cl)][leg+'_count']+=count
    for cl in classes:classes[cl]['types']=sorted(classes[cl]['types'])
    active=set(map(int,inventories['bound']))
    chemical={
        'carboxylate_O':{'small_type':403,'small_class':403,'protein_types':[143,161],'protein_class':31},
        'amide_O':{'small_type':105,'small_class':55,'protein_types':[5,11,153],'protein_class':5},
        'amide_C':{'small_type':104,'small_class':53,'protein_types':[3,9,152],'protein_class':3},
        'water_O':{'small_type':36,'small_class':34,'protein_types':[402],'protein_class':73},
        'water_H':{'small_type':37,'small_class':35,'protein_types':[403],'protein_class':74}}
    for name,v in chemical.items():
        assert sorted(t for t in active if bio['atom'][t]['class']==v['protein_class'])==v['protein_types']
        source=key if v['small_type']in key['atom'] else small
        v['source_atom']=source['atom'][v['small_type']]
        pol=key['polarize'].get(v['small_type'],small['polarize'].get(v['small_type']))
        mp=key['multipole'].get(v['small_type'],small['multipole'].get(v['small_type']))
        vd=key['vdw'].get(v['small_class'],small['vdw'].get(v['small_class']))
        v.update(source_polarize=pol,source_multipoles=mp,source_vdw=vd,
                 bound_atoms=sum(int(inventories['bound'][str(t)]['count']) for t in v['protein_types']),
                 transfer_qualified=False)
    rules=[]
    for p in key['pairs']:
        t=p['tokens'];a,b=map(int,t[1:3]);assert b==500
        if t[0]=='vdwpr':
            cls={55:5,53:3,403:31}[a];rules.append({'operator':'VDWPR','source_class_pair':[a,500],
                'destination_class_pair':[cls,900],'source_values':t[3:],'source_line':p,
                'applies_to_types':classes[str(cls)]['types'],'selected':False})
        else:
            types={105:[5,11,153],403:[143,161]}[a]
            for typ in types:rules.append({'operator':'POLPAIR','source_type_pair':[a,500],
                'destination_type_pair':[typ,900],'source_values':t[3:],'source_line':p,'selected':False})
    optional=[]
    for metal,dest in [('La',900),('Dy',901)]:
        for name in ('water_O','carboxylate_O','amide_O'):
            v=chemical[name];optional.append({'metal':metal,'family':name,'candidate_source':None,
                'VDWPR_destination_class_pair':[v['protein_class'],dest],
                'POLPAIR_destination_type_pairs':[[t,dest] for t in v['protein_types']],
                'numerical_values':None,'status':'unselected future qualified candidate only'})
    def water_records(table,tc,cc):
        result={}
        for typ,new in tc.items():
            for kind in ('polarize','multipole'):
                vals=table[kind][typ];assert len(vals)==1;v=vals[0];t=v['tokens'];tokens=t[:]
                end=len(t) if kind=='polarize' else len(t)-1
                for i in ([1]+list(range(4,end)) if kind=='polarize' else range(1,end)):
                    z=int(tokens[i]);tokens[i]=str((1 if z>=0 else -1)*tc.get(abs(z),abs(z)))
                result[kind+'_'+str(new)]={'normalized_tokens':tokens,'source':v}
        for cls,new in cc.items():
            v=table['vdw'][cls];assert len(v)==1;t=v[0]['tokens'][:];t[1]=str(new);result['vdw_'+str(new)]={'normalized_tokens':t,'source':v[0]}
        for v in table['water_terms']:
            t=v['tokens'];n=2 if t[0]=='bond' else 3;indices=list(map(int,t[1:n+1]))
            if all(i in cc for i in indices):
                norm=[t[0]]+[str(cc[i]) for i in indices]+t[n+1:];k=t[0]+'_'+'_'.join(norm[1:n+1]);result[k]={'normalized_tokens':norm,'source':v}
        return result
    sw=water_records(small,{36:402,37:403},{34:73,35:74});bw=water_records(bio,{402:402,403:403},{73:73,74:74})
    assert set(sw)==set(bw)
    differences=[]
    for k in sw:
        match=sw[k]['normalized_tokens']==bw[k]['normalized_tokens']
        if k.startswith('multipole'):match=match and sw[k]['source']['components']==bw[k]['source']['components']
        if not match:differences.append(k)
    assert differences==['ureybrad_74_73_74'],differences
    result={'status':'mapping_inventory_complete','implementation':pin(__file__),'config':pin(config),
        'specification':c['specification'],'paired_preparation':c['paired_preparation'],'protein_parameters':r['protein_parameters'],
        'small_parameters':c['small_parameters'],'small_key':c['small_key'],'atom_mapping':r['atom_mapping'],
        'bound_boxes':r['bound_boxes'],'aqueous_boxes':r['aqueous_boxes'],'all_represented_types':inventories,
        'all_represented_classes':classes,'chemical_family_map':chemical,'published_La_pair_remapping':rules,
        'future_candidate_pair_destinations':optional,'water_comparison':{'small_normalized':sw,'protein':bw,'different_terms':differences},
        'reserved_metal_types_and_classes':{'La':900,'Dy':901,'unoccupied_in_entire_protein_namespace':True},
        'collision_examples':{'small_acetate_403':bio['atom'][403],'small_amide_O105':bio['atom'][105],'small_amide_C104':bio['atom'][104],
                              'small_carboxyl_C402':bio['atom'][402],'small_water_O36':bio['atom'][36]},
        'native_calls':0,'energy_calls':0,'force_calls':0,'dynamics_steps':0,'candidate_frozen':False,
        'executable_key_written':False,'transfer_qualified':False,'reserved_validation_outcomes_opened':False,
        'local_threads':1,'elapsed_seconds':time.monotonic()-start}
    out.mkdir(parents=True,exist_ok=False);save(out/'MAPPING.json',result)
    print(json.dumps({'mapping':pin(out/'MAPPING.json'),'seconds':result['elapsed_seconds'],'types':len(active),'classes':len(classes)}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.config,a.out)
