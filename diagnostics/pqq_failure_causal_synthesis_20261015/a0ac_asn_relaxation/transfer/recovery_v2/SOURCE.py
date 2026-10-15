"""Scoped source parser: Ce is allowed only at the declared exchanged ion."""
from COMMON import *

def matches(src,sel):
    return all(src.get(k,'')==sel.get(k,'') for k in ('chain','resnum','resname','atom')) and src.get('insertion_code','')==sel.get('insertion_code',sel.get('icode','')) and ('element' not in sel or src['element']==sel['element'])

class DeclaredSource:
    def __init__(self,old,scope,base):
        self.old=old;self.scope=scope;self.base=base
        path=check(scope['protonated_source']);atoms=[];seen=set();models=0
        for line in path.read_text().splitlines():
            if line.startswith('MODEL '):models+=1
            if not line.startswith(('ATOM  ','HETATM')):continue
            a={'chain':line[21:22].strip(),'resnum':int(line[22:26]),'insertion_code':line[26:27].strip(),
               'resname':line[17:20].strip(),'atom':line[12:16].strip(),'altloc':line[16:17].strip(),
               'element':line[76:78].strip().capitalize()}
            declared_ce=a['element']=='Ce' and matches(a,scope['active_selector'])
            require((a['element'] in ('H','C','N','O','S','Ca','La') or declared_ce) and not a['altloc'],'Unsupported source element/alternate location')
            k=base.key(a);require(k not in seen,'Duplicate source identity');seen.add(k)
            atoms.append({'source':a,'element':a['element'],'xyz':[float(line[i:j]) for i,j in ((30,38),(38,46),(46,54))]})
        require(models<=1 and atoms,'Missing or multiple source models')
        active,=[v for v in atoms if matches(v['source'],scope['active_selector'])]
        require(active['element'] in ('Ca','La','Ce'),'Declared source ion is not an allowed exchanged element')
        cp=read(scope['source_preparation']);state=read(scope['state_signature']);maps={z:read(scope['maps'][z]) for z in METALS}
        require(cp['source']==scope['protonated_source'] and maps['Ca']==maps['La'],'Source/preparation/paired map differs')
        context=maps['Ca']['context']
        require(context['metal_index']==0 and context['source_atom_metadata'][0] is None and context['physical_ids'][0]=='metal','Source ion is not the declared endpoint index0')
        pair={z:base.xyz_atoms(check(scope['original_xyz'][z])) for z in METALS}
        require(all(len(v)==cp['new_atom_count'] and v[0]['element']==z and base.close(v[0]['xyz'],active['xyz'],1e-12) for z,v in pair.items()),'Declared source ion does not match paired endpoint index0')
        require(pair['Ca'][0]['xyz']==pair['La'][0]['xyz'] and [(v['element'],v['xyz']) for v in pair['Ca'][1:]]==[(v['element'],v['xyz']) for v in pair['La'][1:]],'Paired original endpoint coordinates differ')
        require(all(base.key(v['source'])!=base.key(active['source']) for v in cp['mapping']['source_to_qm'] if 'source' in v),'Declared ion also appears in nonmetal source map')
        require(all(v['element']!='Ce' or v is active for v in atoms),'Unselected source Ce is unsupported')
        self.atoms=atoms;self.active=active
        self.evidence={'source':scope['protonated_source'],'selector':scope['active_selector'],'selector_evidence':scope['active_selector_evidence'],
          'native_source_identity':active['source'],'native_source_element':active['element'],'native_source_xyz_A':active['xyz'],
          'endpoint_index':0,'endpoint_symbols':['Ca','La'],'original_xyz':scope['original_xyz'],'maps':scope['maps'],
          'state_signature':scope['state_signature'],'endpoint_charges':state['endpoint_charges'],'multiplicities':state['multiplicities'],
          'only_declared_ion_extended':True,'source_bytes_unchanged':True,'electronic_Ce_evaluation':False,'coordinate_tolerance_A':1e-12}

    def __getattr__(self,name):return getattr(self.old,name)

    def physical_source(self,path,base):
        require(pin(path)==self.scope['protonated_source'] and base is self.base,'Scoped source parser called with another source/parser')
        return copy.deepcopy(self.atoms)

    def accepts_declared_source_ion(self,atom):
        return atom==self.active and atom['element']=='Ce'
