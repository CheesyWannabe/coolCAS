from INTM import *
if gc:
  gc.collect()

Y=('k','y')
Z=('k','z')
TA=('k','a')
TB=('k','b')
TC=('k','c')
TS=('k','s')
SYS=[1]
KV['a']=0.3
KV['b']=0.7
KV['c']=0.9
KV['s']=0.6

def tidy(e):
  try:
    r=expand(e)
    if len(P(r))<len(P(e)):
      return r
  except Exception:
    pass
  return e

def setsys(k):
  SYS[0]=k
  MAP.clear()
  ALIAS.clear()
  MAP['x']='r'
  if k==3:
    MAP['y']='p'
    MAP['z']='t'
  else:
    MAP['y']='t'

def aliases():
  ALIAS.clear()
  ALIAS['r']='x'
  if SYS[0]==3:
    ALIAS['p']='y'
    ALIAS['t']='z'
  else:
    ALIAS['t']='y'

def angles():
  if SYS[0]==3:
    return [Z,Y]
  return [Y]

SMP=((0.7,1.1,0.4),(1.9,0.5,2.2),(1.3,2.6,-0.9),(2.4,0.3,1.6))

def indep(e,a):
  vs=[]
  for s in SMP:
    for w in (0.37,1.9,4.1):
      KV[a[1]]=w
      try:
        KV['y'],KV['z']=s[1],s[2]
        KV[a[1]]=w
        vs.append((s,ev(e,s[0])))
      except Exception:
        return False
    a0=vs[-3][1]
    for q in vs[-3:]:
      if abs(q[1]-a0)>1e-9*(1+abs(a0)):
        return False
  return True

def has(e,a):
  t=e[0]
  if t=='k':
    return e==a
  if t=='+' or t=='*':
    for b in e[1]:
      if has(b,a):
        return True
    return False
  if t=='^':
    return has(e[1],a) or has(e[2],a)
  if t=='f':
    return has(e[2],a)
  return False

def kill(e,a):
  if has(e,a) and indep(e,a):
    return S(rep(e,a,ZERO))
  return e

def simp(e):
  t=e[0]
  if t=='+' or t=='*':
    e=S((t,[simp(b) for b in e[1]]))
  elif t=='^':
    e=S(('^',simp(e[1]),simp(e[2])))
  elif t=='f':
    e=S(('f',e[1],simp(e[2])))
  if e[0]!='+':
    return e
  for a in angles():
    e=kill(e,a)
    if e[0]!='+':
      return e
    gs={}
    order=[]
    for b in e[1]:
      if not has(b,a):
        continue
      fs=b[1] if b[0]=='*' else [b]
      k=P(S(prod([f for f in fs if not has(f,a)])))
      if k not in gs:
        gs[k]=[]
        order.append(k)
      gs[k].append(b)
    for k in order:
      L=gs[k]
      if len(L)<2:
        continue
      sm=S(('+',L))
      if indep(sm,a):
        rest=[b for b in e[1] if not any([b is q for q in L])]
        e=S(('+',rest+[S(rep(sm,a,ZERO))]))
        if e[0]!='+':
          return e
  return e

def nonneg(f):
  for s in SMP:
    try:
      KV['y'],KV['z']=s[1],s[2]
      if SYS[0]!=3:
        KV['y']=s[1]*2
      if ev(f,s[0])<-1e-12:
        return False
    except Exception:
      return False
  return True

def posroot(e):
  t=e[0]
  if t=='+' or t=='*':
    return S((t,[posroot(b) for b in e[1]]))
  if t=='f':
    return S(('f',e[1],posroot(e[2])))
  if t!='^':
    return e
  b=posroot(e[1])
  k=e[2]
  if not (isn(k) and k[1].q!=1):
    return S(('^',b,k))
  fs=b[1] if b[0]=='*' else [b]
  out=[]
  keep=[]
  for f in fs:
    if f[0]=='^' and isn(f[2]) and nonneg(f[1]):
      out.append(Pw(f[1],N(f[2][1]*k[1])))
    elif nonneg(f):
      out.append(Pw(f,k))
    else:
      keep.append(f)
  if keep:
    out.append(Pw(prod(keep),k))
  return S(prod(out))

def conv(f):
  g=rep(rep(rep(f,X,TA),Y,TB),Z,TC)
  if SYS[0]==3:
    xr=Mu(X,F('sin',Y),F('cos',Z))
    yr=Mu(X,F('sin',Y),F('sin',Z))
    zr=Mu(X,F('cos',Y))
  else:
    xr=Mu(X,F('cos',Y))
    yr=Mu(X,F('sin',Y))
    zr=Z
  g=S(rep(rep(rep(g,TA,xr),TB,yr),TC,zr))
  g=simp(g)
  g=posroot(g)
  return tidy(simp(g))

def jac():
  if SYS[0]==3:
    return Mu(Pw(X,N(2)),F('sin',Y))
  return X
