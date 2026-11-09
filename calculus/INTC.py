from INTB import *
if gc:
  gc.collect()
def D(e):
  t=e[0]
  if t=='x':
    return ONE
  if t=='n' or t=='k':
    return ZERO
  if t=='+':
    return S(('+',[D(a) for a in e[1]]))
  if t=='*':
    L=e[1]
    out=[]
    for i in range(len(L)):
      d=D(L[i])
      if isn(d,0):
        continue
      out.append(('*',L[:i]+[d]+L[i+1:]))
    return S(('+',out))
  if t=='^':
    b,x=e[1],e[2]
    if free(x):
      return S(Mu(x,Pw(b,Sb(x,ONE)),D(b)))
    if free(b):
      return S(Mu(e,F('ln',b),D(x)))
    return S(Mu(e,A(Mu(D(x),F('ln',b)),Mu(x,D(b),Iv(b)))))
  n,a=e[1],e[2]
  if n=='sin':
    g=F('cos',a)
  elif n=='cos':
    g=Ng(F('sin',a))
  elif n=='tan':
    g=Pw(F('cos',a),N(-2))
  elif n=='exp':
    g=e
  elif n=='ln':
    if a[0]=='f' and a[1]=='abs':
      return S(Mu(D(a[2]),Iv(a[2])))
    g=Iv(a)
  elif n=='abs':
    g=Mu(a,Iv(e))
  elif n=='asin':
    g=Pw(Sb(ONE,Pw(a,N(2))),N(R(-1,2)))
  elif n=='acos':
    g=Ng(Pw(Sb(ONE,Pw(a,N(2))),N(R(-1,2))))
  else:
    g=Iv(A(ONE,Pw(a,N(2))))
  return S(Mu(g,D(a)))
def subst(e,v):
  t=e[0]
  if t=='x':
    return v
  if t=='n' or t=='k':
    return e
  if t=='+' or t=='*':
    return (t,[subst(a,v) for a in e[1]])
  if t=='^':
    return ('^',subst(e[1],v),subst(e[2],v))
  return ('f',e[1],subst(e[2],v))
def Z():
  return [R(0)]
def ptrim(p):
  while len(p)>1 and p[-1]==0:
    p.pop()
  return p
def pdeg(p):
  if len(p)==1 and p[0]==0:
    return -1
  return len(p)-1
def padd(a,b):
  n=max(len(a),len(b))
  r=[]
  for i in range(n):
    x=R(0)
    if i<len(a):
      x=x+a[i]
    if i<len(b):
      x=x+b[i]
    r.append(x)
  return ptrim(r)
def pscal(a,c):
  return ptrim([x*c for x in a])
def psub(a,b):
  return padd(a,pscal(b,R(-1)))
def pmul(a,b):
  r=[R(0)]*(len(a)+len(b)-1)
  for i in range(len(a)):
    for j in range(len(b)):
      r[i+j]=r[i+j]+a[i]*b[j]
  return ptrim(r)
def ppow(a,k):
  r=[R(1)]
  for i in range(k):
    r=pmul(r,a)
  return r
def pdiv(a,b):
  db=pdeg(b)
  q=[R(0)]*max(1,len(a)-len(b)+1)
  r=a[:]
  while pdeg(r)>=db:
    d=pdeg(r)-db
    c=r[-1]/b[-1]
    q[d]=c
    r=psub(r,[R(0)]*d+pscal(b,c))
  return ptrim(q),r
def pgcd(a,b):
  while pdeg(b)>=0:
    a,b=b,pdiv(a,b)[1]
  return pscal(a,R(1)/a[-1])
def pder(a):
  if len(a)==1:
    return Z()
  return ptrim([a[i]*i for i in range(1,len(a))])
def pval(a,x):
  s=R(0)
  for c in reversed(a):
    s=s*x+c
  return s
def pexpr(p,v=X):
  L=[]
  for i in range(len(p)):
    if p[i]!=0:
      L.append(Mu(N(p[i]),Pw(v,N(i))))
  return S(('+',L))
def toP(e):
  t=e[0]
  if t=='n':
    return [e[1]]
  if t=='x':
    return [R(0),R(1)]
  if t=='+' or t=='*':
    s=Z() if t=='+' else [R(1)]
    for a in e[1]:
      p=toP(a)
      if p is None:
        return None
      s=padd(s,p) if t=='+' else pmul(s,p)
    return s
  if t=='^' and e[2][0]=='n' and e[2][1].q==1 and not e[2][1]<0 and e[2][1].p<40:
    p=toP(e[1])
    if p is None:
      return None
    return ppow(p,e[2][1].p)
  return None
def toR(e):
  t=e[0]
  if t=='n' or t=='x':
    return toP(e),[R(1)]
  if t=='+':
    n,d=Z(),[R(1)]
    for a in e[1]:
      q=toR(a)
      if q is None:
        return None
      n,d=padd(pmul(n,q[1]),pmul(q[0],d)),pmul(d,q[1])
    return n,d
  if t=='*':
    n,d=[R(1)],[R(1)]
    for a in e[1]:
      q=toR(a)
      if q is None:
        return None
      n,d=pmul(n,q[0]),pmul(d,q[1])
    return n,d
  if t=='^' and e[2][0]=='n' and e[2][1].q==1 and abs(e[2][1].p)<40:
    q=toR(e[1])
    if q is None:
      return None
    k=e[2][1].p
    if k<0:
      return ppow(q[1],-k),ppow(q[0],-k)
    return ppow(q[0],k),ppow(q[1],k)
  return None
