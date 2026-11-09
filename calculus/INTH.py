from INTG import *
if gc:
  gc.collect()
def nconst(a,b):
  if hassym(a) or hassym(b):
    return None
  vs=[]
  for x in (0.37,0.83,1.41,2.29,-0.71,3.1):
    try:
      w=ev(b,x)
      if abs(w)>1e-12:
        vs.append(ev(a,x)/w)
    except Exception:
      pass
  if len(vs)<3:
    return None
  for v in vs:
    if abs(v-vs[0])>1e-9*max(1,abs(vs[0])):
      return None
  r=rat(vs[0],1e-9*max(1,abs(vs[0])),1000)
  if r is None:
    return None
  return N(r)
def usub(fs,dp):
  for i in range(len(fs)):
    f=fs[i]
    rest=prod(fs[:i]+fs[i+1:])
    cs=[]
    if f[0]=='f' and f[1]!='abs':
      if not lin(f[2]):
        cs.append((f[2],('f',f[1],X)))
      cs.append((f,X))
    if f[0]=='^':
      b,k=f[1],f[2]
      if free(k) and not free(b) and not lin(b):
        cs.append((b,Pw(X,k)))
      if free(b) and not free(k) and not lin(k):
        cs.append((k,Pw(b,X)))
    for g,o in cs:
      dg=D(g)
      if isn(dg,0):
        continue
      q=S(Mu(rest,Iv(dg)))
      if not free(q):
        q=nconst(rest,dg)
        if q is None:
          continue
      H=I(o,dp+1)
      if H is None:
        continue
      return S(Mu(q,subst(H,g)))
  return None
def parts(fs,dp):
  for i in range(len(fs)):
    f=fs[i]
    ok=f[0]=='f' and f[1] in ('ln','atan','asin','acos')
    if f[0]=='^' and f[1][0]=='f' and f[1][1]=='ln' and isn(f[2]) and f[2][1].q==1 and f[2][1]>0:
      ok=True
    if not ok:
      continue
    V=I(prod(fs[:i]+fs[i+1:]),dp+1)
    if V is None:
      continue
    W=I(S(Mu(V,D(f))),dp+1)
    if W is None:
      continue
    return S(Sb(Mu(V,f),W))
  if len(fs)==2:
    a,b=fs
    if b[0]=='f' and b[1]=='exp':
      a,b=b,a
    if a[0]=='f' and a[1]=='exp' and b[0]=='f' and b[1] in ('sin','cos'):
      L1=lin(a[2])
      L2=lin(b[2])
      if L1 and L2:
        p,c=L1[0],L2[0]
        dn=Iv(A(Pw(p,N(2)),Pw(c,N(2))))
        if b[1]=='sin':
          w=Sb(Mu(p,b),Mu(c,F('cos',b[2])))
        else:
          w=A(Mu(p,b),Mu(c,F('sin',b[2])))
        return S(Mu(a,dn,w))
  for i in range(len(fs)):
    f=fs[i]
    tab=(f[0]=='f' and f[1] in ('sin','cos','exp') and lin(f[2])) or (f[0]=='^' and free(f[1]) and lin(f[2]))
    if not tab:
      continue
    Pe=prod(fs[:i]+fs[i+1:])
    if toP(Pe) is None:
      continue
    out=[]
    sg=1
    Pk=Pe
    Tk=I(f,dp+1)
    for it in range(15):
      if Tk is None:
        return None
      out.append(Mu(N(sg),Pk,Tk))
      Pk=D(Pk)
      if isn(Pk,0):
        return S(('+',out))
      Tk=I(Tk,dp+1)
      sg=-sg
  return None
RULES.append(lambda e,fs,dp:usub(fs,dp))
RULES.append(lambda e,fs,dp:parts(fs,dp))
