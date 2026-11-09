from INTD import *
if gc:
  gc.collect()
def lnabs(u):
  return F('ln',F('abs',u))
def prim(Fp):
  l=1
  for c in Fp:
    l=lcm(l,c.q)
  return pexpr(pscal(Fp,R(l)))
def term(Fp,j,Nm):
  d=pdeg(Fp)
  Nm=Nm+[R(0)]*(d-len(Nm))
  Fe=pexpr(Fp)
  if d==1:
    c=N(Nm[0])
    if j==1:
      return Mu(c,lnabs(prim(Fp)))
    return Mu(c,N(R(1,1-j)),Pw(Fe,N(1-j)))
  if d==2:
    p,q=Fp[1],Fp[0]
    a,b=Nm[1],Nm[0]
    k=q-p*p/4
    t=A(Mu(N(2),X),N(p))
    out=[]
    if a!=0:
      if j==1:
        if k>0:
          out.append(Mu(N(a/2),F('ln',prim(Fp))))
        else:
          out.append(Mu(N(a/2),lnabs(prim(Fp))))
      else:
        out.append(Mu(N(a/2/(1-j)),Pw(Fe,N(1-j))))
    bb=b-a*p/2
    if bb!=0:
      if k>0:
        J=Mu(Pw(N(k),N(R(-1,2))),F('atan',Mu(t,Iv(Mu(N(2),Sq(N(k)))))))
      else:
        sd=Sq(N(-k*4))
        J=Mu(N(R(1,2)),Pw(N(-k),N(R(-1,2))),Sb(lnabs(Sb(t,sd)),lnabs(A(t,sd))))
      for n in range(2,j+1):
        J=A(Mu(N(R(1)/(k*4*(n-1))),t,Pw(Fe,N(1-n))),Mu(N(R(2*n-3)/(k*2*(n-1))),J))
      out.append(Mu(N(bb),J))
    return ('+',out)
  if j>1:
    return None
  if d==3:
    c=Fp[0]
    s=1
    if c<0:
      s=-1
    a=S(Mu(N(s),Pw(N(R(abs(c.p),c.q)),N(R(1,3)))))
    L1=lnabs(A(X,a))
    L2=F('ln',A(Pw(X,N(2)),Ng(Mu(a,X)),Pw(a,N(2))))
    T=Mu(sq(3),F('atan',Mu(Sb(Mu(N(2),X),a),Iv(Mu(a,sq(3))))))
    out=[]
    if Nm[2]!=0:
      out.append(Mu(N(Nm[2]/3),lnabs(Fe)))
    if Nm[0]!=0:
      out.append(Mu(N(Nm[0]/3),Pw(a,N(-2)),A(L1,Mu(N(R(-1,2)),L2),T)))
    if Nm[1]!=0:
      out.append(Mu(N(Nm[1]/3),Iv(a),A(Ng(L1),Mu(HALF,L2),T)))
    return ('+',out)
  c=Fp[0]
  k=Pw(N(c),N(R(1,4)))
  s2=sq(2)
  u=Mu(s2,X,Iv(k))
  G=Sb(F('ln',A(Pw(X,N(2)),Mu(s2,k,X),Pw(k,N(2)))),F('ln',A(Pw(X,N(2)),Ng(Mu(s2,k,X)),Pw(k,N(2)))))
  H=A(F('atan',A(u,ONE)),F('atan',Sb(u,ONE)))
  out=[]
  if Nm[0]!=0:
    out.append(Mu(N(Nm[0]/2),Iv(Mu(s2,Pw(k,N(3)))),A(Mu(HALF,G),H)))
  if Nm[1]!=0:
    out.append(Mu(N(Nm[1]/2),Pw(N(c),N(R(-1,2))),F('atan',Mu(Pw(X,N(2)),Pw(N(c),N(R(-1,2)))))))
  if Nm[2]!=0:
    out.append(Mu(N(Nm[2]/2),Iv(Mu(s2,k)),A(Mu(N(R(-1,2)),G),H)))
  if Nm[3]!=0:
    out.append(Mu(N(Nm[3]/4),F('ln',Fe)))
  return ('+',out)
def polyint(p):
  return ('+',[Mu(N(p[i]/(i+1)),Pw(X,N(i+1))) for i in range(len(p)) if p[i]!=0])
def ratint(Pn,Qd):
  if pdeg(Pn)<0:
    return ZERO
  lc=Qd[-1]
  Pn=pscal(Pn,R(1)/lc)
  Qd=pscal(Qd,R(1)/lc)
  g=pgcd(Pn,Qd)
  if pdeg(g)>0:
    Pn=pdiv(Pn,g)[0]
    Qd=pdiv(Qd,g)[0]
  q,r=pdiv(Pn,Qd)
  out=[polyint(q)]
  if pdeg(r)<0:
    return S(('+',out))
  fs=[]
  for sf,k in yun(Qd):
    L=splitf(sf)
    if L is None:
      return None
    for f in L:
      if pdeg(f)>2 and k>1:
        return None
      fs.append((f,k))
  cols=[]
  meta=[]
  for f,k in fs:
    for j in range(1,k+1):
      base=pdiv(Qd,ppow(f,j))[0]
      for t in range(pdeg(f)):
        cols.append(pmul(base,[R(0)]*t+[R(1)]))
        meta.append((f,j,t))
  n=len(cols)
  M=[]
  for i in range(n):
    row=[]
    for c in cols:
      row.append(c[i] if i<len(c) else R(0))
    row.append(r[i] if i<len(r) else R(0))
    M.append(row)
  sol=solve(M,n)
  i=0
  for f,k in fs:
    for j in range(1,k+1):
      nm=sol[i:i+pdeg(f)]
      i+=pdeg(f)
      allz=True
      for v in nm:
        if v!=0:
          allz=False
      if allz:
        continue
      tt=term(f,j,nm)
      if tt is None:
        return None
      out.append(tt)
  return S(('+',out))
