"""Declared anchor chain rule; inherited slope algorithm at two fixed new steps."""
from decimal import Decimal,localcontext
from types import FunctionType
STEPS=(Decimal('0.0000125'),Decimal('0.00000625'))
def corrected_projection(helper,central,direction):
 with localcontext() as ctx:
  ctx.prec=70
  raw,raw_print=helper.projection(central,direction,None)
  sums=[sum((helper.number(t) for t in central['gradient_tokens'][j::3]),Decimal(0)) for j in range(3)]
  printing=[sum((helper.rounding(t) for t in central['gradient_tokens'][j::3]),Decimal(0)) for j in range(3)]
  S=[[Decimal.from_float(float(x)) for x in row] for row in direction['strain_matrix']]
  a=[Decimal.from_float(float(x)) for x in direction['anchor_A']];L=Decimal.from_float(float(direction['L0_A']))
  helper.require(L>0 and len(central['gradient_tokens'])%3==0,'Finite anchor projection dimensions')
  coeff=[sum((S[i][j]*a[j] for j in range(3)),Decimal(0))/L for i in range(3)]
  adjustment=-sum((c*g for c,g in zip(coeff,sums)),Decimal(0))
  allowance=sum((abs(c)*e for c,e in zip(coeff,printing)),Decimal(0))
  return dict(raw=raw,raw_print=raw_print,corrected=raw+adjustment,corrected_print=raw_print+allowance,adjustment=adjustment,adjustment_print=allowance,gradient_sums=sums,coefficients=coeff)

def checks(helper,name,analytic,printing,energies):
 # Code object and original gate/rounding logic remain unchanged; only prospectively
 # declared STEPS changes in this private globals view. No module/source mutation.
 fn=FunctionType(helper.slope_checks.__code__,dict(helper.slope_checks.__globals__,STEPS=STEPS),helper.slope_checks.__name__)
 return fn(name,analytic,printing,energies)
