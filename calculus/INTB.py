from INTA import *
if gc:
  gc.collect()
def S(e):
  t=e[0]
  if t=='n' or t=='x' or t=='k':
    return e
  if t=='f':
    return Sf(e[1],S(e[2]))
  if t=='^':
    return Sp(S(e[1]),S(e[2]))
  if t=='*':
    return Sm([S(a) for a in e[1]])
  return Sa([S(a) for a in e[1]])
def xdeg(e):
  if free(e):
    return -1.0
  if e[0]=='x':
    return 1.0
  if e[0]=='^' and e[1][0]=='x' and e[2][0]=='n':
    return e[2][1].f()
  if e[0]=='*':
    for a in e[1]:
      d=xdeg(a)
      if d>0:
        return d
  return 0.0
def Sa(L):
  fl=[]
  for a in L:
    if a[0]=='+':
      fl.extend(a[1])
    else:
      fl.append(a)
  c=R(0)
  T=[]
  for a in fl:
    if a[0]=='n':
      c=c+a[1]
      continue
    k,r=split(a)
    key=P(r)
    for t in T:
      if t[0]==key:
        t[1]=t[1]+k
        break
    else:
      T.append([key,k,r])
  out=[]
  for t in T:
    k,r=t[1],t[2]
    if k==0:
      continue
    if k==1:
      out.append(r)
    elif r[0]=='*':
      out.append(('*',[N(k)]+r[1]))
    else:
      out.append(('*',[N(k),r]))
  if c!=0:
    out.append(N(c))
  if not out:
    return ZERO
  if len(out)==1:
    return out[0]
  out=sorted(out,key=lambda a:-xdeg(a))
  return ('+',out)
def fac(n):
  d={}
  p=2
  while p*p<=n and p<5000:
    while n%p==0:
      d[p]=d.get(p,0)+1
      n=n//p
    p+=1
  if n>1:
    d[n]=d.get(n,0)+1
  return d
def cls(f):
  if f[0]=='^' and f[1][0]=='n' and f[2][0]=='n':
    return 1
  if free(f):
    return 2
  if f[0]=='x' or (f[0]=='^' and f[1][0]=='x'):
    return 3
  if f[0]=='+' or (f[0]=='^' and f[1][0]=='+'):
    return 4
  return 5
def content(L):
  g=0
  l=1
  for a in L:
    k=split(a)[0]
    if a[0]=='n':
      k=a[1]
    g=gcd(g,k.p)
    l=lcm(l,k.q)
  return R(g,l)
def Sm(L):
  fl=[]
  for a in L:
    if a[0]=='*':
      fl.extend(a[1])
    else:
      fl.append(a)
  c=R(1)
  sd={}
  bs=[]
  ex=[]
  for a in fl:
    t=a[0]
    if t=='n':
      c=c*a[1]
      continue
    if t=='^' and a[1][0]=='n' and a[2][0]=='n':
      b=a[1][1]
      r=a[2][1]
      ok=True
      if b<0:
        if r.q%2==1:
          if r.p%2:
            c=-c
          b=-b
        else:
          ok=False
      if ok and b!=0:
        for p,k in fac(b.p).items():
          sd[p]=sd.get(p,R(0))+r*k
        for p,k in fac(b.q).items():
          sd[p]=sd.get(p,R(0))-r*k
        continue
    if t=='f' and a[1]=='exp':
      ex.append(a[2])
      continue
    if t=='^':
      b,k=a[1],a[2]
    else:
      b,k=a,ONE
    key=P(b)
    for z in bs:
      if z[0]==key:
        z[2].append(k)
        break
    else:
      bs.append([key,b,[k]])
  if c==0:
    return ZERO
  gr={}
  for p in sd:
    r=sd[p]
    k=r.p//r.q
    fr=r-k
    c=c*rpow(R(p),k)
    if fr!=0:
      kk=(fr.p,fr.q)
      gr[kk]=gr.get(kk,1)*p
  sur=[]
  for kk in sorted(gr.keys()):
    sur.append(('^',N(gr[kk]),N(R(kk[0],kk[1]))))
  oth=[]
  def put(v):
    if v[0]=='n':
      return v[1]
    if v[0]=='*':
      cc=R(1)
      for w in v[1]:
        if w[0]=='n':
          cc=cc*w[1]
        else:
          oth.append(w)
      return cc
    oth.append(v)
    return R(1)
  if ex:
    c=c*put(Sf('exp',Sa(ex)))
  for z in bs:
    c=c*put(Sp(z[1],Sa(z[2])))
  if c==0:
    return ZERO
  oth=sorted(oth,key=lambda f:(cls(f),P(f)))
  res=[]
  if c!=1:
    res.append(N(c))
  res=res+sur+oth
  if len(res)==2 and res[0][0]=='n' and res[1][0]=='+':
    return Sa([Sm([res[0],t]) for t in res[1][1]])
  return prod(res)
def Sp(b,x):
  if x[0]=='n':
    r=x[1]
    if r==0:
      return ONE
    if r==1:
      return b
    if b[0]=='n':
      if b[1]==0 and r>0:
        return ZERO
      if r.q==1:
        if b[1]==0 and r<0:
          raise ZeroDivisionError("div 0")
        return N(rpow(b[1],r.p))
      return Sm([('^',b,x)])
    if b[0]=='^' and r.q==1:
      return Sp(b[1],S(Mu(b[2],x)))
    if b[0]=='*' and r.q==1:
      return Sm([Sp(f,x) for f in b[1]])
    if b[0]=='+' and r<0:
      c=content(b[1])
      if c!=1 and (r.q==1 or c>0):
        w=Sa([Sm([N(R(1)/c),t]) for t in b[1]])
        return Sm([Sp(N(c),x),('^',w,x)])
  if b[0]=='f' and b[1]=='exp':
    return Sf('exp',S(Mu(b[2],x)))
  if b==EC:
    return Sf('exp',x)
  if isn(b,1):
    return ONE
  return ('^',b,x)
def sq(n):
  return ('^',N(n),HALF)
def sin15(k):
  k=k%24
  s=1
  if k>=12:
    k-=12
    s=-1
  if k>6:
    k=12-k
  V=[ZERO,Mu(N(R(1,4)),Sb(sq(6),sq(2))),HALF,Mu(HALF,sq(2)),Mu(HALF,sq(3)),Mu(N(R(1,4)),A(sq(6),sq(2))),ONE]
  return S(Mu(N(s),V[k]))
def near(v):
  k=int(floor(v+0.5))
  if abs(v-k)<1e-9:
    return k
  return None
def Sf(n,a):
  if n=='exp':
    if isn(a,0):
      return ONE
    if a[0]=='f' and a[1]=='ln':
      return a[2]
    if a[0]=='*' and len(a[1])==2 and a[1][0][0]=='n' and a[1][0][1].q==1 and a[1][1][0]=='f' and a[1][1][1]=='ln':
      return Sp(a[1][1][2],a[1][0])
    return ('f','exp',a)
  if n=='ln':
    if a[0]=='n':
      r=a[1]
      if not r>0:
        raise ValueError("ln dom")
      if r==1:
        return ZERO
      d=fac(r.p)
      for p,k in fac(r.q).items():
        d[p]=d.get(p,0)-k
      if len(d)==1:
        for p in d:
          if d[p]==1:
            return ('f','ln',N(p))
      return Sa([('*',[N(d[p]),('f','ln',N(p))]) for p in d])
    if a==EC:
      return ONE
    if a[0]=='f' and a[1]=='exp':
      return a[2]
    if a[0]=='^' and a[1][0]=='n' and a[2][0]=='n' and a[1][1]>0:
      return S(Mu(a[2],Sf('ln',a[1])))
    if a[0]=='*' and free(a):
      ok=True
      for f in a[1]:
        if not ((f[0]=='n' and f[1]>0) or (f[0]=='^' and f[1][0]=='n' and f[1][1]>0 and f[2][0]=='n')):
          ok=False
      if ok:
        return Sa([Sf('ln',f) for f in a[1]])
    return ('f','ln',a)
  if n=='abs':
    if numc(a):
      try:
        if ev(a)<0:
          return S(Ng(a))
        return a
      except Exception:
        pass
    if a[0]=='f' and a[1]=='exp':
      return a
    return ('f','abs',a)
  if n in ('sin','cos','tan','asin','atan'):
    k,r=split(a)
    if k<0:
      w=Sf(n,Sm([N(-k),r]))
      if n=='cos':
        return w
      return S(Ng(w))
  if numc(a):
    try:
      v=ev(a)
      if n in ('sin','cos','tan'):
        k=near(v*12/pi)
        if k is not None:
          if n=='sin':
            return sin15(k)
          if n=='cos':
            return sin15(k+6)
          if k%2==0 or k%3==0:
            c=sin15(k+6)
            if not isn(c,0):
              return S(Mu(sin15(k),Iv(c)))
      if n in ('asin','acos','atan'):
        if n=='asin':
          w=asin(v)
        elif n=='acos':
          w=acos(v)
        else:
          w=atan(v)
        k=near(w*12/pi)
        if k is not None:
          return S(Mu(N(R(k,12)),PI))
    except Exception:
      pass
  return ('f',n,a)
