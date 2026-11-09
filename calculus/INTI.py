from INTH import *
if gc:
  gc.collect()
def wrep(e,K):
  t=e[0]
  if t=='x':
    raise ValueError("x")
  if t=='n' or t=='k':
    return e
  if t=='+' or t=='*':
    return (t,[wrep(a,K) for a in e[1]])
  if t=='^':
    return ('^',wrep(e[1],K),wrep(e[2],K))
  if e[1] in ('sin','cos','tan') and P(e[2])==K:
    d=Iv(A(ONE,Pw(X,N(2))))
    if e[1]=='sin':
      return Mu(N(2),X,d)
    if e[1]=='cos':
      return Mu(Sb(ONE,Pw(X,N(2))),d)
    return Mu(N(2),X,Iv(Sb(ONE,Pw(X,N(2)))))
  raise ValueError("fn")
def findtrig(e):
  t=e[0]
  if t=='f':
    if e[1] in ('sin','cos','tan'):
      return e[2]
    return findtrig(e[2])
  if t=='+' or t=='*':
    for a in e[1]:
      g=findtrig(a)
      if g:
        return g
  if t=='^':
    return findtrig(e[1]) or findtrig(e[2])
  return None
def weier(e,dp):
  g=findtrig(e)
  if g is None:
    return None
  L=lin(g)
  if not L:
    return None
  try:
    h=S(Mu(N(2),wrep(e,P(g)),Iv(A(ONE,Pw(X,N(2))))))
  except Exception:
    return None
  q=toR(h)
  if q is None:
    return None
  H=ratint(q[0],q[1])
  if H is None:
    return None
  return S(Mu(Iv(L[0]),subst(H,F('tan',S(Mu(HALF,g))))))
def roots_of(e,K,out):
  t=e[0]
  if t=='^' and e[2][0]=='n' and e[2][1].q>1 and P(e[1])==K:
    out.append(e[2][1].q)
    return
  if t=='+' or t=='*':
    for a in e[1]:
      roots_of(a,K,out)
  elif t=='^':
    roots_of(e[1],K,out)
    roots_of(e[2],K,out)
  elif t=='f':
    roots_of(e[2],K,out)
def findroot(e):
  t=e[0]
  if t=='^' and e[2][0]=='n' and e[2][1].q>1 and lin(e[1]):
    return e[1]
  if t=='+' or t=='*':
    for a in e[1]:
      g=findroot(a)
      if g:
        return g
  if t=='^':
    return findroot(e[1]) or findroot(e[2])
  if t=='f':
    return findroot(e[2])
  return None
def rrep(e,K,n,xe):
  t=e[0]
  if t=='^' and e[2][0]=='n' and P(e[1])==K:
    return Pw(X,N(e[2][1]*n))
  if t=='x':
    return xe
  if t=='n' or t=='k':
    return e
  if t=='+' or t=='*':
    return (t,[rrep(a,K,n,xe) for a in e[1]])
  if t=='^':
    return ('^',rrep(e[1],K,n,xe),rrep(e[2],K,n,xe))
  return ('f',e[1],rrep(e[2],K,n,xe))
def rootsub(e,dp):
  L=findroot(e)
  if L is None:
    return None
  K=P(L)
  ds=[]
  roots_of(e,K,ds)
  n=1
  for d in ds:
    n=lcm(n,d)
  a,b=lin(L)
  xe=Mu(Sb(Pw(X,N(n)),b),Iv(a))
  h=S(Mu(rrep(e,K,n,xe),N(n),Pw(X,N(n-1)),Iv(a)))
  H=I(h,dp+1)
  if H is None:
    return None
  return S(subst(H,Pw(L,N(R(1,n)))))
def erep(e,r):
  t=e[0]
  if t=='x':
    raise ValueError("x")
  if t=='n' or t=='k':
    return e
  if t=='+' or t=='*':
    return (t,[erep(a,r) for a in e[1]])
  if t=='^':
    return ('^',erep(e[1],r),erep(e[2],r))
  if e[1]=='exp':
    L=lin(e[2])
    if not L or not isn(L[0]):
      raise ValueError("arg")
    return Mu(F('exp',L[1]),Pw(X,N(L[0][1]/r)))
  return ('f',e[1],erep(e[2],r))
def exprates(e,out):
  t=e[0]
  if t=='f':
    if e[1]=='exp':
      L=lin(e[2])
      if L and isn(L[0]):
        out.append(L[0][1])
    exprates(e[2],out)
  elif t=='+' or t=='*':
    for a in e[1]:
      exprates(a,out)
  elif t=='^':
    exprates(e[1],out)
    exprates(e[2],out)
def expsub(e,dp):
  rs=[]
  exprates(e,rs)
  if not rs:
    return None
  g=0
  l=1
  for r in rs:
    g=gcd(g,r.p)
    l=lcm(l,r.q)
  r=R(g,l)
  try:
    h=S(Mu(erep(e,r),Iv(Mu(N(r),X))))
  except Exception:
    return None
  H=I(h,dp+1)
  if H is None:
    return None
  return S(subst(H,F('exp',Mu(N(r),X))))
RULES.append(lambda e,fs,dp:rootsub(e,dp))
RULES.append(lambda e,fs,dp:expsub(e,dp))
RULES.append(lambda e,fs,dp:weier(e,dp))
