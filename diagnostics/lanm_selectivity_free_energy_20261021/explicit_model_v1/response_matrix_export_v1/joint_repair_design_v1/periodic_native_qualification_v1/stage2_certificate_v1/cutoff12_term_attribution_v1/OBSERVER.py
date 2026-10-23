"""Insertion-only source generation and strict small term-trace reader."""
from decimal import Decimal, localcontext
from pathlib import Path
import re
TERMS=[('bond','ebond1','eb'),('angle','eangle1','ea'),('strbnd','estrbnd1','eba'),
 ('urey','eurey1','eub'),('angang','eangang1','eaa'),('opbend','eopbend1','eopb'),
 ('opdist','eopdist1','eopd'),('improp','eimprop1','eid'),('imptor','eimptor1','eit'),
 ('tors','etors1','et'),('pitors','epitors1','ept'),('strtor','estrtor1','ebt'),
 ('angtor','eangtor1','eat'),('tortor','etortor1','ett'),('charge','echarge1','ec'),
 ('chgdpl','echgdpl1','ecd'),('dipole','edipole1','ed'),('mpole','empole1','em'),
 ('polar','epolar1','ep'),('chgtrn','echgtrn1','ect'),('rxnfld','erxnfld1','erxf'),
 ('vdw',None,'ev'),('repel','erepel1','er'),('disp','edisp1','edsp'),
 ('solv','esolv1','es'),('metal','emetal1','elf'),('geom','egeom1','eg'),('extra','extra1','ex')]
BEGIN='c     TERM_OBSERVER_INSERT_BEGIN\n';END='c     TERM_OBSERVER_INSERT_END\n'
def need(ok,msg):
 if not ok:raise ValueError(msg)
def marked(s):return BEGIN+s+END
def inverse(s):return re.sub(re.escape(BEGIN)+r'.*?'+re.escape(END),'',s,flags=re.S)
def render(old):
 s=old
 def add(anchor,code):
  nonlocal s
  need(s.count(anchor)==1,'Unique source anchor '+anchor);s=s.replace(anchor,anchor+marked(code))
 add('      einter = 0.0d0\n',"      call term_observe (0,'INITIAL',0.0d0,.true.)\n")
 for i,(name,call,energy) in enumerate(TERMS,1):
  anchor=('      if (use_'+name+')  call '+call+'\n') if call else "         if (vdwtyp .eq. 'GAUSSIAN')  call egauss1\n      end if\n"
  add(anchor,f"      call term_observe ({i},'{name}',{energy},use_{name})\n")
 add('      call watfour2 (derivs)\n',"      call term_observe (29,'FINAL',energy,.true.)\n")
 s+=marked(HELPER)
 need(inverse(s)==old,'Exact insertion inverse')
 need(all(len(l)<=72 for l in HELPER.splitlines() if not l.startswith('c')),'Fixed form width')
 return s
HELPER="""      subroutine term_observe (idx,label,value,enabled)
      use virial, only: vir
      use energi
      use vdwpot, only: use_vcorr,vdwtyp
      implicit none
      integer idx,i,j,u,last,active
      character*(*) label
      real*8 value
      logical enabled
      save u,last
      data last /-1/
      if (idx.ne.last+1) stop 'TERM OBSERVER ORDER'
      if (idx.eq.0) then
         open(newunit=u,file='term_trace.dat',status='new')
         write(u,'(a)') 'TERM_TRACE_V1'
         write(u,'(a)') 'UNITS energy=kcal/mol virial=kcal/mol'
         write(u,'(a,1x,l1,1x,a)') 'VDW',use_vcorr,trim(vdwtyp)
      end if
      active=0
      if (enabled) active=1
      write(u,10) idx,trim(label),active,value,
     &            ((vir(i,j),j=1,3),i=1,3)
   10 format(i3,1x,a,1x,i1,10(1x,es26.17e3))
      last=idx
      if (idx.eq.29) then
         write(u,20) eb,ea,eba,eub,eaa,eopb,eopd,eid,eit,et,
     &       ept,ebt,eat,ett,ec,ecd,ed,em,ep,ect,erxf,ev,er,edsp,
     &       es,elf,eg,ex
   20    format('FINAL_TERMS',28(1x,es26.17e3))
         write(u,'(a)') 'END snapshots=30 gradient_calls=1'
         close(u)
      else
         flush(u)
      end if
      end
"""
def parse(path):
 lines=Path(path).read_text().splitlines()
 need(len(lines)==35,'Complete single trace')
 need(lines[:2]==['TERM_TRACE_V1','UNITS energy=kcal/mol virial=kcal/mol'],'Trace header')
 vdw=lines[2].split();need(len(vdw)==3 and vdw[0]=='VDW' and vdw[1] in ('T','F'),'VDW flags')
 labels=['INITIAL']+[t[0] for t in TERMS]+['FINAL'];rows=[]
 for i,(label,line) in enumerate(zip(labels,lines[3:33])):
  v=line.split();need(len(v)==13 and v[:2]==[str(i),label] and v[2] in ('0','1'),'Trace ordering/width')
  need(all(Decimal(x).is_finite() for x in v[3:]),'Finite trace')
  rows.append(dict(label=label,enabled=v[2]=='1',energy_token=v[3],virial_tokens=v[4:]))
 final=lines[33].split();need(len(final)==29 and final[0]=='FINAL_TERMS','Final energy components')
 need(final[1:]==[r['energy_token'] for r in rows[1:-1]],'No component changed after its snapshot')
 need(lines[34]=='END snapshots=30 gradient_calls=1','One gradient termination')
 need(rows[0]['enabled'] and rows[-1]['enabled'] and all(Decimal(x)==0 for x in [rows[0]['energy_token'],*rows[0]['virial_tokens']]),'Initial zero accumulators')
 return dict(rows=rows,vdw_long_range_correction=vdw[1]=='T',vdw_type=vdw[2])
def checks(trace,gradient,old,helper):
 need(gradient['E_G_virial_token_sha256']==old['E_G_virial_token_sha256'],'Exact saved total E/G/V tokens')
 for key in ('polarization_token','electric_token','dielec_token','role','state_id'):
  need(gradient[key]==old[key],'Exact saved '+key)
 rows=trace['rows'];need(rows[-1]['energy_token']==gradient['energy_token'] and rows[-1]['virial_tokens']==gradient['virial_tokens'],'Final observer equals native report tokens')
 need(rows[19]['energy_token']==gradient['polarization_token'],'Polarization energy tokens')
 need(not trace['vdw_long_range_correction'] and trace['vdw_type']=='BUFFERED-14-7','Fixed baseline VDW mode')
 with localcontext() as c:
  c.prec=70;energies=[Decimal(r['energy_token']) for r in rows[1:-1]]
  error=abs(sum(energies)-Decimal(gradient['energy_token']))
  bound=helper.arithmetic(len(energies),sum(map(abs,energies)),Decimal(gradient['energy_token']))+sum(helper.rounding(r['energy_token']) for r in rows[1:])
  need(error<=bound,'Term energy sum arithmetic closure')
  increments=[[Decimal(b)-Decimal(a) for a,b in zip(rows[i-1]['virial_tokens'],r['virial_tokens'])] for i,r in enumerate(rows) if i]
  need(all(x==0 for x in increments[-1]),'No unassigned post-term virial increment')
  need(all(sum(r[j] for r in increments)==Decimal(rows[-1]['virial_tokens'][j]) for j in range(9)),'Virial telescoping closure')
  return dict(exact_total_token_parity=True,energy_sum_error=float(error),energy_sum_bound=float(bound),virial_closure=True,post_term_virial_zero=True,
              term_virial_increments_tokens=[[str(x) for x in r] for r in increments[:-1]])
