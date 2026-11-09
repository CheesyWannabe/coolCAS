from COCV import *

def solvefor(E,v):
  if v=='x':
    r=inv(E,ZERO)
    return None if r is None else posroot(S(r))
  T=Y if v=='y' else Z
  E2=rep(rep(E,X,TS),T,X)
  r=inv(S(E2),ZERO)
  if r is None:
    return None
  return posroot(S(rep(rep(S(r),X,T),TS,X)))

def angsol(h,a,lo,hi):
  out=[]
  n=360
  def fv(t):
    KV[a[1]]=t
    return ev(h,1.0)
  prev=None
  for i in range(n+1):
    t=lo+(hi-lo)*i/n
    try:
      v=fv(t)
    except Exception:
      prev=None
      continue
    if v==0:
      out.append(t)
    elif prev is not None and prev[1]*v<0:
      a_,b_=prev[0],t
      for k in range(70):
        m=(a_+b_)/2
        if fv(a_)*fv(m)<=0:
          b_=m
        else:
          a_=m
      out.append((a_+b_)/2)
    prev=(t,v)
  res=[]
  for t in out:
    if not any([abs(t-q)<1e-7 for q in res]):
      res.append(t)
  return res

def surface():
  s=input("equation (eg z=2):")
  if '=' not in s:
    print("need an = sign")
    return
  l,r=s.split('=')
  E=S(Sb(parse(l),parse(r)))
  E=conv(E)
  print("becomes: "+mp(E)+" = 0")
  pref=['x','y'] if SYS[0]==1 else ['z','x','y'] if SYS[0]==2 else ['x','y','z']
  Es=[E]
  for k in (1,2):
    try:
      Es.append(tidy(expand(S(Mu(E,Pw(X,N(-k)))))))
    except Exception:
      Es.append(None)
  for v in pref:
    T=X if v=='x' else Y if v=='y' else Z
    for k in range(3):
      Ek=Es[k]
      if Ek is None:
        continue
      if v=='x' and free(Ek):
        continue
      if v!='x' and not has(Ek,T):
        continue
      try:
        w=solvefor(Ek,v)
      except Exception:
        w=None
      if w is not None:
        if k>0:
          print("(divide by "+mv('x')+("^"+str(k) if k>1 else "")+")")
        print(mv(v)+" = "+mp(w))
        return
  for k in (1,2):
    Ek=Es[k]
    if Ek is None or not free(Ek):
      continue
    for a in angles():
      if has(Ek,a):
        if a==Y and SYS[0]==3:
          rs=angsol(Ek,a,0.0,pi)
        else:
          rs=angsol(Ek,a,0.0,2*pi)
        nm=mv('y' if a==Y else 'z')
        print("(r cancels) "+nm+" = "+", ".join([ident(t) or str(t) for t in rs]))
        return
  print("can't solve; use form above")

def run(s):
  setsys(s)
  try:
    surface()
  finally:
    MAP.clear()
    ALIAS.clear()
