#!/usr/bin/env python3
"""Restricted, source-pinned Tinker parameter renderer; no native execution."""
import math
import shlex

TYPE=502
LAMBDAS=(0.,.25,.5,.75,1.,.4999,.5001,.49995,.50005)
SETTINGS={'vdwtype':'BUFFERED-14-7','radiusrule':'CUBIC-MEAN',
          'radiustype':'R-MIN','radiussize':'DIAMETER','epsilonrule':'HHG','vdwindex':'CLASS'}


def require(condition,message):
    if not condition:raise ValueError(message)


def parse_parameters(prm_text,key_text):
    """Read relevant values in native parameter-file then key override order."""
    result={'settings':{'vdwindex':'CLASS'},'atom':{},'vdw':{},'polarize':{},'vdwpr':{},'polpair':{},'metal_multipoles':{}}
    for text in (prm_text,key_text):
        lines=text.splitlines()
        for i,line in enumerate(lines):
            tokens=line.split()
            if not tokens or tokens[0].lower() not in set(SETTINGS)|{'atom','vdw','polarize','vdwpr','polpair','multipole'}:continue
            a=shlex.split(line,comments=True)
            if not a:continue
            keyword=a[0].lower()
            if keyword in SETTINGS:result['settings'][keyword]=a[1].upper()
            elif keyword=='atom':
                result['atom'][int(a[1])]={'class':int(a[2]),'symbol':a[3],'atomic_number':int(a[5]),
                                         'mass':float(a[6]),'valence':int(a[7])}
            elif keyword=='vdw':
                result['vdw'][int(a[1])]={'diameter_A':float(a[2]),'epsilon_kcal_mol':float(a[3]),
                                        'reduction':float(a[4]) if len(a)>4 else 0.}
            elif keyword=='polarize':
                # Current candidate uses explicit real-valued Thole, integer group types.
                require(len(a)>=4 and any(x in a[3].lower()for x in ('.','e')),'Unsupported polarize record')
                require(all(v.lstrip('+-').isdigit() for v in a[4:]),'Unsupported direct Thole damping')
                result['polarize'][int(a[1])]={'alpha_A3':float(a[2]),'thole':float(a[3]),
                                             'groups':[int(v)for v in a[4:]]}
            elif keyword in ('vdwpr','polpair'):
                pair=tuple(sorted((int(a[1]),int(a[2]))))
                require(len(a)==(5 if keyword=='vdwpr' else 4),'Unsupported pair record')
                result[keyword][pair]=tuple(float(x)for x in a[3:])
            elif keyword=='multipole' and int(a[1]) in (500,501,502):
                require(len(a)==5 and a[2:4]==['0','0'],'Unsupported metal multipole axes')
                tails=[shlex.split(v,comments=True)for v in lines[i+1:i+5]]
                require([len(v)for v in tails]==[3,1,2,3],'Incomplete metal multipole')
                result['metal_multipoles'][int(a[1])]=[float(a[4])]+[float(v)for row in tails for v in row]
    require(result['settings']==SETTINGS,'Unsupported VDW conventions')
    require(TYPE not in result['atom'] and TYPE not in result['vdw'] and TYPE not in result['polarize'],
            'Artificial type/class already exists')
    require(all(atom['class']!=TYPE for atom in result['atom'].values()),'Artificial class already used by another type')
    require(all(TYPE not in pair for name in ('vdwpr','polpair')for pair in result[name]),'Artificial pair already exists')
    for t,z in ((500,57),(501,66)):
        require(result['atom'][t]['class']==t and result['atom'][t]['atomic_number']==z and
                result['atom'][t]['valence']==0,'Endpoint metal identity')
        require(result['metal_multipoles'][t]==[3.]+[0.]*9,'Endpoint metal charge/multipoles')
        require(result['polarize'][t]['alpha_A3']>0 and result['polarize'][t]['thole']>0 and
                not result['polarize'][t]['groups'],'Metal polarization state')
    return result


def resolve_pair(parameters,metal_type,other_type):
    atom=parameters['atom'];mc,oc=atom[metal_type]['class'],atom[other_type]['class']
    m,o=parameters['vdw'][mc],parameters['vdw'][oc]
    require(m['diameter_A']>0 and o['diameter_A']>0 and m['epsilon_kcal_mol']>0 and o['epsilon_kcal_mol']>0,
            'Unsupported nonpositive VDW endpoint')
    radii=[m['diameter_A']/2,o['diameter_A']/2]
    pair=tuple(sorted((mc,oc)))
    if pair in parameters['vdwpr']:
        radius,epsilon=parameters['vdwpr'][pair];vdw_origin='explicit_VDWPR'
    else:
        radius=2*sum(x**3 for x in radii)/sum(x*x for x in radii)
        em,eo=m['epsilon_kcal_mol'],o['epsilon_kcal_mol']
        epsilon=4*em*eo/(math.sqrt(em)+math.sqrt(eo))**2
        vdw_origin='native_CUBIC_MEAN_HHG'
    pair=tuple(sorted((metal_type,other_type)))
    if pair in parameters['polpair']:
        thole,=parameters['polpair'][pair];pol_origin='explicit_POLPAIR'
    else:
        values=[parameters['polarize'][metal_type]['thole'],parameters['polarize'][other_type]['thole']]
        thole=min(values) or max(values);pol_origin='native_minimum_nonzero'
    require(radius>0 and epsilon>0 and thole>0,'Unsupported nonpositive resolved pair')
    return {'nonmetal_type':other_type,'nonmetal_class':oc,'Rmin_pair_A':radius,
            'epsilon_kcal_mol':epsilon,'thole':thole,'VDW_origin':vdw_origin,'POL_origin':pol_origin}


def interpolate(a,b,lam):
    require(math.isfinite(lam)and 0<=lam<=1,'Lambda outside [0,1]')
    return (1-lam)*a+lam*b


def endpoint_resolution(parameters,nonmetal_types):
    return {str(t):{'La':resolve_pair(parameters,500,t),'Dy':resolve_pair(parameters,501,t)}
            for t in sorted(set(nonmetal_types))}


def render_key(prm_text,key_text,nonmetal_types,lam):
    p=parse_parameters(prm_text,key_text);resolved=endpoint_resolution(p,nonmetal_types)
    alpha=interpolate(p['polarize'][500]['alpha_A3'],p['polarize'][501]['alpha_A3'],lam)
    thole=interpolate(p['polarize'][500]['thole'],p['polarize'][501]['thole'],lam)
    # One isolated metal: atomic VDW values allocate its class. Every actual
    # metal/nonmetal pair is explicitly overridden below; there is no self pair.
    diam=interpolate(p['vdw'][500]['diameter_A'],p['vdw'][501]['diameter_A'],lam)
    eps=interpolate(p['vdw'][500]['epsilon_kcal_mol'],p['vdw'][501]['epsilon_kcal_mol'],lam)
    extra=[f'# Numerical substitution path lambda={lam:.8f}; candidate unqualified',
           'atom 502 502 La "Artificial La-to-Dy configurational path" 57 138.91 0',
           'multipole 502 0 0 3.000000000000000',
           ' 0.0 0.0 0.0',' 0.0',' 0.0 0.0',' 0.0 0.0 0.0',
           f'polarize 502 {alpha:.16g} {thole:.16g}',f'vdw 502 {diam:.16g} {eps:.16g}']
    classes={};current={}
    for t,ends in resolved.items():
        a,b=ends['La'],ends['Dy'];cl=a['nonmetal_class']
        pair={name:interpolate(a[name],b[name],lam)for name in ('Rmin_pair_A','epsilon_kcal_mol','thole')}
        if cl in classes:require(classes[cl]==(pair['Rmin_pair_A'],pair['epsilon_kcal_mol']),'Conflicting class resolution')
        else:
            extra.append(f"vdwpr {cl} 502 {pair['Rmin_pair_A']:.16g} {pair['epsilon_kcal_mol']:.16g}")
            classes[cl]=(pair['Rmin_pair_A'],pair['epsilon_kcal_mol'])
        extra.append(f"polpair {t} 502 {pair['thole']:.16g}")
        current[t]=dict(pair,nonmetal_class=cl)
    return key_text.rstrip()+'\n\n'+'\n'.join(extra)+'\n',{'lambda':lam,'alpha_A3':alpha,
            'charge':3,'fixed_mass_u':138.91,'metal_type':502,'metal_class':502,
            'actual_pair_parameters':current,'resolved_endpoints':resolved}


def render_xyz(text):
    lines=text.splitlines();n=int(lines[0].split()[0]);rows=[v.split()for v in lines[1:]if v.strip()]
    require(len(rows)==n and [int(v[0])for v in rows]==list(range(1,n+1)),'XYZ dimensions/index order')
    require(int(rows[0][5])==500 and len(rows[0])==6,'Expected isolated La atom first')
    require(all(int(v[5])not in (500,501,502)for v in rows[1:]),'Multiple metal atoms')
    require(all('1'not in v[6:]for v in rows[1:]),'Nonmetal declares bond to isolated metal')
    # Preserve every nonmetal line and all coordinate strings exactly.
    first=rows[0][:];first[1]='La';first[5]='502'
    return lines[0]+'\n'+' '.join(first)+'\n'+'\n'.join(lines[2:])+'\n',rows
