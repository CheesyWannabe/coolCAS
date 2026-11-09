from INTM import *
if gc:
  gc.collect()

def mono(vals):
  up=True
  dn=True
  for i in range(len(vals)-1):
    if not vals[i+1]>vals[i]:
      up=False
    if not vals[i+1]<vals[i]:
      dn=False
  if up:
    return 1
  if dn:
    return -1
  return 0

def bnd(g,qs,vals):
  if numc(g):
    return 0,None
  m=mono(vals)
  if m==0:
    return None
  h=inv(g,W)
  if h is None:
    return None
  h=S(h)
  for i in range(1,len(qs)-1):
    KV['w']=vals[i]
    if abs(ev(h)-qs[i])>1e-7*max(1,abs(qs[i])):
      return None
  return m,h

def pick(cands,ps,big):
  for c in cands:
    good=True
    for p in ps:
      KV['w']=p
      v=ev(c)
      for d in cands:
        w=ev(d)
        if (big and w>v+1e-9) or (not big and w<v-1e-9):
          good=False
    if good:
      return c
  return None

def swapreg(od,lm):
  p,q=od
  a,b=lm[1]
  lo,hi=lm[0]
  if not (numc(a) and numc(b)):
    return None
  if numc(lo) and numc(hi):
    return [q,p],[(a,b),(lo,hi)]
  sv=KV['w']
  try:
    gl=S(sw(lo,q))
    gh=S(sw(hi,q))
    av=ev(a)
    bv=ev(b)
    qs=[av+(bv-av)*i/16 for i in range(17)]
    lv=[ev(gl,s) for s in qs]
    hv=[ev(gh,s) for s in qs]
    for i in range(17):
      if lv[i]>hv[i]+1e-12:
        return None
    L=bnd(gl,qs,lv)
    H=bnd(gh,qs,hv)
    if L is None or H is None:
      return None
    lows=[a]
    ups=[b]
    if L[0]==1:
      ups.append(L[1])
    elif L[0]==-1:
      lows.append(L[1])
    if H[0]==1:
      lows.append(H[1])
    elif H[0]==-1:
      ups.append(H[1])
    if lv[0]<=lv[-1]:
      P0=S(subst(gl,a))
    else:
      P0=S(subst(gl,b))
    if hv[0]>=hv[-1]:
      P1=S(subst(gh,a))
    else:
      P1=S(subst(gh,b))
    p0=ev(P0)
    p1=ev(P1)
    ps=[p0+(p1-p0)*(i+0.5)/12 for i in range(12)]
    ql=pick(lows,ps,True)
    qh=pick(ups,ps,False)
    if ql is None or qh is None:
      return None
    nl=[(S(rep(ql,W,nd(p))),S(rep(qh,W,nd(p)))),(P0,P1)]
  except Exception:
    return None
  finally:
    KV['w']=sv
  return [q,p],nl

def same(f,od,lm,od2,lm2):
  for g in (ONE,A(ONE,Mu(N(R(3,7)),X),Mu(N(R(5,9)),('k','y')),Mu(N(R(2,11)),('k','z')))):
    v1=numint(g,od,lm,4)
    v2=numint(g,od2,lm2,4)
    if v1 is None or v2 is None or abs(v1-v2)>1e-6*max(1,abs(v1)):
      return False
  return True
