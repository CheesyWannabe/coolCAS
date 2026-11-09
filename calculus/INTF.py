from INTE import *
if gc:
  gc.collect()
def lin(g):
  if free(g):
    return None
  d=D(g)
  if not free(d) or isn(d,0):
    return None
  b=S(Sb(g,Mu(d,X)))
  if not free(b):
    return None
  return d,b
def mulx(a,b):
  P1=a[1] if a[0]=='+' else [a]
  P2=b[1] if b[0]=='+' else [b]
  return S(('+',[Mu(x,y) for x in P1 for y in P2]))
def expand(e):
  t=e[0]
  if t=='+':
    return S(('+',[expand(a) for a in e[1]]))
  if t=='^' and e[1][0]=='+' and isn(e[2]) and e[2][1].q==1 and e[2][1]>1 and e[2][1]<9:
    b=expand(e[1])
    r=b
    for i in range(e[2][1].p-1):
      r=mulx(r,b)
    return r
  if t=='*':
    r=ONE
    for a in e[1]:
      r=mulx(r,expand(a))
    return r
  return e
RULES=[]

def I(e,dp=0):
  if dp>8:
    return None
  e=S(e)
  if free(e):
    return S(Mu(e,X))
  q=toR(e)
  if q:
    v=ratint(q[0],q[1])
    if v is not None:
      return v
  t=e[0]
  if t=='+':
    out=[]
    for a in e[1]:
      v=I(a,dp+1)
      if v is None:
        return None
      out.append(v)
    return S(('+',out))
  fs=e[1] if t=='*' else [e]
  cs=[a for a in fs if free(a)]
  if cs:
    v=I(prod([a for a in fs if not free(a)]),dp+1)
    if v is None:
      return None
    return S(Mu(prod(cs),v))
  for rule in RULES:
    v=rule(e,fs,dp)
    if v is not None:
      return v
  x=expand(e)
  if P(x)!=P(e):
    return I(x,dp+1)
  return None
