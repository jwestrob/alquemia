"""Report-only saved parent cycles; no additional fit rows or numerical calls."""
import math
import numpy as np

PARENTS={
 'aquo8_central_La':('ion_La','aquo8_central_La__shell',8),
 'aquo8_central_Dy':('ion_Dy','aquo8_central_Dy__shell',8),
 'aquo9_central_La':('ion_La','aquo9_central_La__shell',9),
 'aquo9_central_Dy':('ion_Dy','aquo9_central_Dy__shell',9)}
CONTRAST={'aquo9_central_Dy':1,'aquo9_central_La':-1,
          'aquo8_central_Dy':-1,'aquo8_central_La':1}

def total(terms,pred,key):
 values=[c*np.asarray(pred[id][key]) for id,c in terms.items()]
 if values[0].ndim==0:return math.fsum(float(x) for x in values)
 return np.array([math.fsum(float(x.flat[j]) for x in values)
                  for j in range(values[0].size)]).reshape(values[0].shape)

def comparison(value,target):
 v=np.asarray(value);t=np.asarray(target)
 return {'prediction':v.tolist(),'target':t.tolist(),'error':(v-t).tolist()}

def report(pred,bundle):
 parents={};hartree=bundle['units']['hartree_kcal_mol']
 energy_fields={'C':'cluster_interaction_hartree','P':'pair_sum_hartree',
                'B':'cooperative_remainder_hartree'}
 for parent,(ion,shell,n) in PARENTS.items():
  C={parent:1,ion:-1,shell:-1};P={}
  pairs=[r for r in bundle['pair_cycles'] if r['id'].startswith(parent+'__pair_')]
  if len(pairs)!=n:raise ValueError('exact parent pair-cycle denominator '+parent)
  for row in pairs:
   for id,c in row['model_terms'].items():P[id]=P.get(id,0)+c
  B=C.copy()
  for id,c in P.items():B[id]=B.get(id,0)-c
  rows={}
  for name,terms in [('C',C),('P',P),('B',B)]:
   reference=bundle['water_common_basis_components']['per_parent_electronic'][parent]
   moment=bundle['water_parent_components_and_cycles'][parent][name]
   rows[name]={'energy_kcal_mol':comparison(total(terms,pred,'energy'),reference[energy_fields[name]]*hartree),
               'mu_Debye':comparison(total(terms,pred,'mu'),moment['reference_mu_Debye']),
               'model_terms':terms,'reference_mu_print_bound_Debye_per_component':moment['reference_print_bound_Debye_per_component']}
  parents[parent]=rows
 def contrast(terms):
  answer={}
  for component in ('C','P','B'):
   answer[component]={}
   for quantity in ('energy_kcal_mol','mu_Debye'):
    values={id:parents[id][component][quantity] for id in terms}
    answer[component][quantity]=comparison(total(terms,values,'prediction'),total(terms,values,'target'))
  return {'parent_terms':terms,'components':answer}
 return {'parents':parents,'four_parent_contrast':contrast(CONTRAST),
         'Dy_minus_La_by_coordination':{str(n):contrast({f'aquo{n}_central_Dy':1,f'aquo{n}_central_La':-1}) for n in (8,9)},
         'method':'Water electronic common-parent-basis reference; D4 excluded. Signed conjugate moments at the declared common origin.',
         'fit_rows_added':0,'report_only':True,
         'interpretation':'C=P+B; component errors, metal differences and their cancellation remain distinct. Consumed data, no exchange free-energy validation.'}
