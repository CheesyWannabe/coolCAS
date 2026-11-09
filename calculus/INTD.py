from INTC import *
if gc:
  gc.collect()
def yun(f):
  out=[]
  fd=pder(f)
  a=pgcd(f,fd)
  b=pdiv(f,a)[0]
  c=pdiv(fd,a)[0]
  d=psub(c,pder(b))
  i=1
  while pdeg(b)>0:
    a=pgcd(b,d)
    if pdeg(a)>0:
      out.append((a,i))
    b=pdiv(b,a)[0]
    c=pdiv(d,a)[0]
    d=psub(c,pder(b))
    i+=1
  return out
def divs(n):
  L=[]
  i=1
  while i*i<=n:
    if n%i==0:
      L.append(i)
      if i*i!=n:
        L.append(n//i)
    i+=1
  return L
def rroots(f):
  L=1
  for c in f:
    L=L*c.q//gcd(L,c.q)
  z=[(c*L).p for c in f]
  out=[]
  if z[0]==0:
    out.append(R(0))
    i=0
    while z[i]==0:
      i+=1
    z=z[i:]
  a0,an=abs(z[0]),abs(z[-1])
  if a0>10**7 or an>10**7:
    return out
  for p in divs(a0):
    for q in divs(an):
      for s in (1,-1):
        r=R(s*p,q)
        if pval(f,r)==0:
          new=True
          for w in out:
            if w==r:
              new=False
          if new:
            out.append(r)
  return out
def cm(a,b):
  return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def cdv(a,b):
  d=b[0]*b[0]+b[1]*b[1]
  return ((a[0]*b[0]+a[1]*b[1])/d,(a[1]*b[0]-a[0]*b[1])/d)
def croots(f):
  c=[x.f() for x in f]
  n=len(c)-1
  z=[]
  w=(0.4,0.9)
  p=(1.0,0.0)
  for i in range(n):
    z.append(p)
    p=cm(p,w)
  for it in range(400):
    mx=0.0
    for i in range(n):
      v=(c[n],0.0)
      for k in range(n-1,-1,-1):
        v=cm(v,z[i])
        v=(v[0]+c[k],v[1])
      d=(1.0,0.0)
      for j in range(n):
        if j!=i:
          d=cm(d,(z[i][0]-z[j][0],z[i][1]-z[j][1]))
      s=cdv(v,d)
      z[i]=(z[i][0]-s[0],z[i][1]-s[1])
      mx=max(mx,abs(s[0])+abs(s[1]))
    if mx<1e-14:
      break
  return z
def rat(x,t,Q=10000):
  s=1
  if x<0:
    s,x=-1,-x
  h0,h1,k0,k1=0,1,1,0
  y=x
  for i in range(40):
    a=int(y)
    h0,h1=h1,a*h1+h0
    k0,k1=k1,a*k1+k0
    if k1>Q:
      return None
    if abs(x-h1/k1)<=t:
      return R(s*h1,k1)
    f=y-a
    if f<=0:
      return None
    y=1/f
  return None
def quads(f):
  z=croots(f)
  re=[]
  cx=[]
  for w in z:
    if abs(w[1])<1e-7:
      re.append(w[0])
    elif w[1]>0:
      cx.append(w)
  if len(re)%2:
    return None
  cand=[]
  for w in cx:
    cand.append([(-2*w[0]),(w[0]*w[0]+w[1]*w[1])])
  def match(L):
    if not L:
      return [[]]
    out=[]
    for j in range(1,len(L)):
      rest=L[1:j]+L[j+1:]
      for m in match(rest):
        out.append([[-(L[0]+L[j]),L[0]*L[j]]]+m)
    return out
  for m in match(re):
    qs=[]
    ok=True
    for p in cand+m:
      b=rat(p[0],1e-7)
      c=rat(p[1],1e-7)
      if b is None or c is None:
        ok=False
        break
      qs.append([c,b,R(1)])
    if not ok:
      continue
    g=f
    for q in qs:
      a,r=pdiv(g,q)
      if pdeg(r)>=0:
        ok=False
        break
      g=a
    if ok:
      return qs
  return None
def splitf(f):
  out=[]
  for r in rroots(f):
    l=[-r,R(1)]
    out.append(l)
    f=pdiv(f,l)[0]
  d=pdeg(f)
  if d<=0:
    return out
  if d==2:
    return out+[f]
  if d==3:
    if f[1]==0 and f[2]==0:
      return out+[f]
    return None
  if d>8:
    return None
  qs=quads(f)
  if qs:
    return out+qs
  if d==4 and f[1]==0 and f[2]==0 and f[3]==0 and f[0]>0:
    return out+[f]
  return None
def solve(M,n):
  for c in range(n):
    p=c
    while M[p][c]==0:
      p+=1
    M[c],M[p]=M[p],M[c]
    iv=R(1)/M[c][c]
    M[c]=[v*iv for v in M[c]]
    for r in range(n):
      if r!=c and M[r][c]!=0:
        k=M[r][c]
        M[r]=[M[r][i]-k*M[c][i] for i in range(n+1)]
  return [M[i][n] for i in range(n)]
