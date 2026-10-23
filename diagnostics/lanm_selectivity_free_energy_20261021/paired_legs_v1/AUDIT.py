#!/usr/bin/env python3
"""Independent saved-coordinate, charge and source-water audit; no native calls."""
from collections import Counter
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import sys
import time


def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def readxyz(p):
    lines=Path(p).read_text().splitlines(); n=int(lines[0].split()[0])
    data=[line.split() for line in lines[-n:]]
    return [{'i':int(t[0]),'label':t[1],'pos':tuple(map(float,t[2:5])),
             'type':int(t[5]),'bond':tuple(map(int,t[6:]))}for t in data]


def main():
    start=time.monotonic(); p=Path(sys.argv[1]).resolve(); out=Path(sys.argv[2]).resolve()
    result=json.loads(p.read_text()); checks=0; pins=[]
    def check(ok,what):
        nonlocal checks
        checks+=1
        if not ok: raise AssertionError(what)
    def verify(x):
        check(pin(x['path'])==x,x['path']); pins.append(x)
        return Path(x['path'])
    for name in ('implementation','config','parent_bound_result','source_water','protein_parameters','atom_mapping'):
        verify(result[name])
    cfg=json.loads(Path(result['config']['path']).read_text()); verify(cfg['specification']); verify(cfg['bound_charge_review'])
    review=json.loads(Path(cfg['bound_charge_review']['path']).read_text()); verify(review['native_parameter_output'])
    check(Decimal(review['sum_printed_monopoles_e'])==Decimal('-4'),'peptide charge')
    maps=json.loads(Path(result['atom_mapping']['path']).read_text())
    source=readxyz(result['source_water']['path']); L=result['box_length_A']; L0=L/2
    charges={}; prm=Path(result['protein_parameters']['path']).read_text().splitlines()
    for line in prm:
        t=line.split()
        if t[:1]==['multipole'] and int(t[1])in(402,403,405,414): charges[int(t[1])]=Decimal(t[-1])
    check(charges[402]+2*charges[403]==0,'water exact neutral')
    check(charges[405]==1 and charges[414]==-1,'ion charges')
    stats={}; max_source_error=0.
    for leg,key,metal_i in [('bound','bound_boxes',258),('aqueous','aqueous_boxes',1)]:
        endpoints={metal:readxyz(verify(v))for metal,v in result[key].items()}
        a,b=endpoints['La'],endpoints['Dy']; check(len(a)==len(b),'endpoint count')
        for x,y in zip(a,b):
            check(x['i']==y['i'] and x['pos']==y['pos'] and x['bond']==y['bond'],'identity mapping')
            check((x['type'],y['type'])==(900,901) if x['i']==metal_i else x==y,'nonmetal endpoint identity')
        mapping=maps[leg]; check([x['index_1based']for x in mapping]==list(range(1,len(a)+1)),'mapping coverage')
        count=Counter(x['type']for x in a); check(count[402]*2==count[403],'water stoichiometry')
        formal=Decimal(-4 if leg=='bound' else 0)+3+count[405]*charges[405]+count[414]*charges[414]
        check(formal==0,'neutral leg')
        oh=[]; angles=[]
        for r,m in zip(a,mapping):
            for neighbor in r['bond']:check(r['i']in a[neighbor-1]['bond'],'reciprocal connectivity')
            if r['type']not in (402,403,405,414):continue
            j=m['source_water_index_zero_based']*3
            oxygen=source[j]['pos']; shift=[(v-.5)*L0 for v in m['tile']]
            if r['type']in(402,405,414):expected=[oxygen[k]+shift[k]for k in range(3)]
            else:
                oi=r['bond'][0]; oxygenrow=a[oi-1]
                # Bond order, separately parsed from written XYZ, identifies H1/H2.
                hk=oxygenrow['bond'].index(r['i'])+1; h=source[j+hk]['pos']
                expected=[oxygen[k]+shift[k]+h[k]-oxygen[k]-L0*round((h[k]-oxygen[k])/L0)for k in range(3)]
            error=max(abs(r['pos'][k]-expected[k]-L*round((r['pos'][k]-expected[k])/L))for k in range(3))
            max_source_error=max(max_source_error,error);check(error<1e-8,'saved water/ion source geometry')
            if r['type']==402:
                check(len(r['bond'])==2,'water connectivity')
                vectors=[[a[h-1]['pos'][k]-r['pos'][k]-L*round((a[h-1]['pos'][k]-r['pos'][k])/L)for k in range(3)]for h in r['bond']]
                norms=[math.sqrt(sum(v*v for v in vec))for vec in vectors];oh.extend(norms)
                cos=sum(x*y for x,y in zip(*vectors))/(norms[0]*norms[1]);angles.append(math.degrees(math.acos(max(-1,min(1,cos)))))
        stats[leg]={'atoms':len(a),'waters':count[402],'Na':count[405],'Cl':count[414],
                    'exact_formal_charge_e':str(formal),'metal_index_1based':metal_i,
                    'water_OH_min_max_A':[min(oh),max(oh)],'water_HOH_min_max_deg':[min(angles),max(angles)]}
    data={'status':'passed','implementation':pin(__file__),'prepared_result':pin(p),
          'checks':checks,'pins':pins,'local_threads':1,'native_calls':0,
          'max_source_coordinate_error_A':max_source_error,'legs':stats,
          'scope':'saved topology, source-water coordinates, endpoint identity and neutrality; no Hamiltonian qualification',
          'elapsed_seconds':time.monotonic()-start}
    with out.open('x')as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:data[k]for k in ('status','checks','max_source_coordinate_error_A','legs','elapsed_seconds')}))


if __name__=='__main__':main()
