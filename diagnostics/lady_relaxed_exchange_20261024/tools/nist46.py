import sys,re
D=sys.argv[1]
def lines(fn):
    raw=open(f'{D}/{fn}','rb').read().decode('latin1')
    return raw.split('\r\n') if '\r\n' in raw else raw.split('\r')
def table(fn, width=None):
    raw=open(f'{D}/{fn}','rb').read().decode('latin1')
    # entries begin with index number
    out={}
    for m in re.finditer(r'(?:^|\s)(\d{1,4})\s{1,2}(\S.*?)(?=\s+\d{1,4}\s{1,2}\S|\r|\Z)',raw[12:]):
        pass
    return raw
code2={}
for L in lines('F8548_CODE2.DAT')[1:]:
    m=re.match(r'\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)',L)
    if m: code2[int(m.group(1))]=m.group(4)
code2[1]='?'
code3={}
for L in lines('F8549_CODE3.DAT')[1:]:
    m=re.match(r'\s*(\d+)\s+(\S+)\s+(\S+)',L)
    if m: code3[int(m.group(1))]=m.group(2)+' '+m.group(3)
raw4=open(f'{D}/F8550_CODE4.DAT','rb').read().decode('latin1')
code4={int(i):v.strip() for i,v in re.findall(r'(\d+) ([KHS\*] [\dtv\*]+ [\d\.\*]+)',raw4)}
raw5=open(f'{D}/F8551_CODE5.DAT','rb').read().decode('latin1')
code5={int(i):v for i,v in re.findall(r'(\d+) ([#+\-]?\d+)',raw5[12:])}
code6={}
for L in lines('F8552_CODE6.DAT'):
    m=re.match(r'\s*(\d+)\s+(\S+)\s+(.*)',L)
    if m: code6[int(m.group(1))]=m.group(2)+': '+m.group(3).strip()
c1=lines('F8547_CODE1.DAT')
lig={}
for i,L in enumerate(c1):
    m=re.search(r'([A-Z]\d{3,4})\t(.*?)\t',L)
    if m: lig[i+1]=(m.group(1),m.group(2))
v=open(f'{D}/F8558_VALUES.DAT','rb').read()
n=int(v[:6]); recs=[]
p=15
while p+15<=len(v):
    r=v[p:p+15]; p+=15
    recs.append((int.from_bytes(r[0:2],'little'),r[2],int.from_bytes(r[3:5],'little'),r[5],r[6:13].decode('latin1').strip(),r[13],r[14]))
mode=sys.argv[2]
if mode=='find':
    for k,(c,nm) in lig.items():
        if re.search(sys.argv[3],nm): print(k,c,nm)
elif mode=='show':
    ligs=set(int(x) for x in sys.argv[3].split(','))
    mets=set(sys.argv[4].split(','))
    for (L,M,E,C,val,c5,c6) in recs:
        if L in ligs and code2.get(M,'?').split('/')[0] in mets:
            print(L, code2.get(M), code3.get(E,E), code4.get(C,C), val, code5.get(c5,c5), code6.get(c6,c6) if c6>1 else '')
print('nrec',len(recs),'hdr',n, file=sys.stderr)
if mode=='probe':
    from collections import defaultdict
    t=defaultdict(dict)
    for (L,M,E,C,val,c5,c6) in recs:
        key=(code2.get(M),code3.get(E,'').split()[0] if E in code3 else E,code4.get(C))
        try: t[L][key]=float(val.strip('()'))
        except: pass
    for L,d in t.items():
        ca=d.get(('Ca/2+','ML','K 25 0.1')); la=d.get(('La/3+','ML','K 25 0.1')); h=d.get(('H/+','HL','K 25 0.1'))
        if ca and la and 10.3<ca<11 and 14.8<la<16: print(L,ca,la,h)
if mode=='probeac':
    from collections import defaultdict
    t=defaultdict(dict)
    for (L,M,E,C,val,c5,c6) in recs:
        key=(code2.get(M),code3.get(E,'').split()[0] if E in code3 else E,code4.get(C))
        try: t[L][key]=float(val.strip('()'))
        except: pass
    for L,d in t.items():
        h=d.get(('H/+','HL','K 25 0.0')); la=[k for k in d if k[0]=='La/3+']
        if h and 4.70<h<4.80 and la: print(L,h,len(d),d.get(('Cu/2+','ML','K 25 0.1')),d.get(('H/+','HL','H 25 0.0')))
