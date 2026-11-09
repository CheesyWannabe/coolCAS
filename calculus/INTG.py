from INTF import *
if gc:
  gc.collect()
def table(e,fs,dp):
  t=e[0]
  if t=='f':
    n,g=e[1],e[2]
    if n=='ln' and g[0]=='f' and g[1]=='abs':
      L=lin(g[2])
      if L:
        return S(Mu(Iv(L[0]),Sb(Mu(g[2],e),g[2])))
    L=lin(g)
    if not L:
      return None
    if n=='sin':
      r=Ng(F('cos',g))
    elif n=='cos':
      r=F('sin',g)
    elif n=='tan':
      r=Ng(lnabs(F('cos',g)))
    elif n=='exp':
      r=e
    elif n=='ln':
      r=Sb(Mu(g,e),g)
    elif n=='atan':
      r=Sb(Mu(g,e),Mu(HALF,F('ln',A(ONE,Pw(g,N(2))))))
    elif n=='asin':
      r=A(Mu(g,e),Sq(Sb(ONE,Pw(g,N(2)))))
    elif n=='acos':
      r=Sb(Mu(g,e),Sq(Sb(ONE,Pw(g,N(2)))))
    elif n=='abs':
      r=Mu(HALF,g,e)
    else:
      return None
    return S(Mu(Iv(L[0]),r))
  if t!='^':
    return None
  b,k=e[1],e[2]
  if free(k):
    L=lin(b)
    if L:
      if isn(k,-1):
        return S(Mu(Iv(L[0]),lnabs(b)))
      k1=S(A(k,ONE))
      return S(Mu(Pw(b,k1),Iv(Mu(L[0],k1))))
    if b[0]=='f' and isn(k) and k[1].q==1:
      L=lin(b[2])
      kk=k[1].p
      if L:
        g=b[2]
        if b[1]=='cos' and kk==-2:
          return S(Mu(Iv(L[0]),F('tan',g)))
        if b[1]=='sin' and kk==-2:
          return S(Ng(Mu(Iv(L[0]),F('cos',g),Iv(F('sin',g)))))
        if b[1]=='tan' and kk>=2:
          v=I(Pw(b,N(kk-2)),dp+1)
          if v is None:
            return None
          return S(Sb(Mu(Iv(Mu(L[0],N(kk-1))),Pw(b,N(kk-1))),v))
    if isn(k) and k[1].q==2 and abs(k[1].p)==1:
      Q=toP(b)
      if Q and pdeg(Q)==2:
        al=Q[2]
        h=-Q[1]/(al*2)
        kc=Q[0]-Q[1]*Q[1]/(al*4)
        if kc==0:
          return None
        s=Sq(N(R(abs(al.p),al.q)))
        u=Mu(s,A(X,N(-h)))
        sk=Sq(N(R(abs(kc.p),kc.q)))
        if al<0:
          if not kc>0:
            return None
          W=F('asin',Mu(u,Iv(sk)))
        else:
          W=lnabs(A(u,Sq(b)))
        if k[1].p<0:
          return S(Mu(Iv(s),W))
        return S(Mu(Iv(s),A(Mu(HALF,u,Sq(b)),Mu(N(kc/2),W))))
  if free(b) and not free(k):
    L=lin(k)
    if L:
      try:
        if ev(b)>0 and abs(ev(b)-1)>1e-12:
          return S(Mu(e,Iv(Mu(L[0],F('ln',b)))))
      except Exception:
        pass
  return None
def trinfo(f):
  if f[0]=='f' and f[1] in ('sin','cos','tan'):
    return f[1],f[2],1
  if f[0]=='^' and f[1][0]=='f' and f[1][1] in ('sin','cos','tan') and isn(f[2]) and f[2][1].q==1:
    return f[1][1],f[1][2],f[2][1].p
  return None
def trig(fs,dp):
  info=[]
  for f in fs:
    v=trinfo(f)
    if v is None:
      return None
    info.append(v)
  if len(info)==2 and P(info[0][1])!=P(info[1][1]):
    (n1,a1,k1),(n2,a2,k2)=info
    if k1!=1 or k2!=1 or n1=='tan' or n2=='tan':
      return None
    if n1=='cos' and n2=='sin':
      n1,a1,n2,a2=n2,a2,n1,a1
    sm=A(a1,a2)
    df=Sb(a1,a2)
    if n1=='sin' and n2=='cos':
      g=Mu(HALF,A(F('sin',sm),F('sin',df)))
    elif n1=='sin':
      g=Mu(HALF,Sb(F('cos',df),F('cos',sm)))
    else:
      g=Mu(HALF,A(F('cos',df),F('cos',sm)))
    return I(g,dp+1)
  key=P(info[0][1])
  m=0
  n=0
  for v in info:
    if P(v[1])!=key:
      return None
    if v[0]=='sin':
      m+=v[2]
    elif v[0]=='cos':
      n+=v[2]
    else:
      m+=v[2]
      n-=v[2]
  g=info[0][1]
  L=lin(g)
  if not L:
    return None
  ia=Iv(L[0])
  if n%2:
    h=Mu(Pw(X,N(m)),Pw(Sb(ONE,Pw(X,N(2))),N(R(n-1,2))))
    H=I(h,dp+1)
    if H is None:
      return None
    return S(Mu(ia,subst(H,F('sin',g))))
  if m%2:
    h=Ng(Mu(Pw(X,N(n)),Pw(Sb(ONE,Pw(X,N(2))),N(R(m-1,2)))))
    H=I(h,dp+1)
    if H is None:
      return None
    return S(Mu(ia,subst(H,F('cos',g))))
  if m<0 or n<0:
    return None
  pc=pmul(ppow([R(1,2),R(-1,2)],m//2),ppow([R(1,2),R(1,2)],n//2))
  c2=F('cos',S(Mu(N(2),g)))
  h=('+',[Mu(N(pc[i]),Pw(c2,N(i))) for i in range(len(pc)) if pc[i]!=0])
  return I(h,dp+1)
RULES.append(table)
RULES.append(lambda e,fs,dp:trig(fs,dp))
