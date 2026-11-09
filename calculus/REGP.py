from INTM import *
if gc:
  gc.collect()

Y=('k','y')

def rd(prompt):
  s=input(prompt)
  if s=="":
    return None
  return readf(s)

def num(e):
  return ev(e)

def pt(s):
  a,b=s.split(',')
  return (S(parse(a)),S(parse(b)))

def line(p,q,v):
  x1,y1=p
  x2,y2=q
  if v=='x':
    return S(A(y1,Mu(Sb(y2,y1),Iv(Sb(x2,x1)),Sb(X,x1))))
  return S(A(x1,Mu(Sb(x2,x1),Iv(Sb(y2,y1)),Sb(Y,y1))))

def hull(P0):
  cx=sum([num(p[0]) for p in P0])/len(P0)
  cy=sum([num(p[1]) for p in P0])/len(P0)
  return sorted(P0,key=lambda p:atan2(num(p[1])-cy,num(p[0])-cx))

def slabs(V,v):
  k=0 if v=='x' else 1
  cs=[]
  for p in V:
    new=True
    for c in cs:
      if abs(num(c)-num(p[k]))<1e-12:
        new=False
    if new:
      cs.append(p[k])
  cs=sorted(cs,key=num)
  n=len(V)
  out=[]
  for i in range(len(cs)-1):
    a,b=cs[i],cs[i+1]
    m=(num(a)+num(b))/2
    es=[]
    for j in range(n):
      p,q=V[j],V[(j+1)%n]
      lo=min(num(p[k]),num(q[k]))
      hi=max(num(p[k]),num(q[k]))
      if lo<=m<=hi and hi-lo>1e-12:
        es.append(line(p,q,v))
    if len(es)!=2:
      return None
    Ls=es
    if v=='x':
      vals=[ev(L,m) for L in Ls]
    else:
      vals=[ev(rep(L,Y,X),m) for L in Ls]
    if vals[0]>vals[1]:
      Ls=[Ls[1],Ls[0]]
    out.append((Ls[0],Ls[1],a,b))
  return out

def polyint(f,V):
  V=hull(V)
  sx=slabs(V,'x')
  sy=slabs(V,'y')
  if sx is None and sy is None:
    print("Polygon must be convex")
    return
  use='x'
  pcs=sx
  if sx is None or (sy is not None and len(sy)<len(sx)):
    use='y'
    pcs=sy
  od=['y','x'] if use=='x' else ['x','y']
  print("Setup ("+str(len(pcs))+" piece"+("s" if len(pcs)>1 else "")+"):")
  for lo,hi,a,b in pcs:
    print(" "+use+": "+P(a)+" to "+P(b))
    print("  "+od[0]+": "+P(lo)+" to "+P(hi))
  total(f,od,[[(lo,hi),(a,b)] for lo,hi,a,b in pcs])

def total(f,od,lms):
  G=ZERO
  ok=True
  vt=0.0
  for lm in lms:
    if len(lms)>1:
      print("piece "+od[1]+"="+P(lm[1][0])+".."+P(lm[1][1])+":")
    g,v=report(f,od,lm)
    if v is None:
      print("Can't evaluate")
      return
    vt+=v
    if g is None:
      ok=False
    else:
      G=S(A(G,g))
  if len(lms)>1:
    print("TOTAL:")
    if ok:
      print("= "+P(lncomb(G)))
      print("~ "+str(ev(G)))
    else:
      print("~ "+str(vt))
      s=ident(vt)
      if s:
        print("maybe = "+s)

def lin(eq):
  l,r=eq.split('=')
  E=S(Sb(parse(l),parse(r)))
  c0=put(E,ZERO,ZERO)
  a=S(Sb(put(E,ONE,ZERO),c0))
  b=S(Sb(put(E,ZERO,ONE),c0))
  if abs(ev(put(E,N(2),N(3)))-ev(S(A(Mu(N(2),a),Mu(N(3),b),c0))))>1e-9:
    raise ValueError("not a line")
  return a,b,c0

def put(e,xe,ye):
  return S(subst(rep(e,Y,ye),xe))

def lines2poly(Ls):
  V=[]
  for i in range(len(Ls)):
    for j in range(i+1,len(Ls)):
      a1,b1,c1=Ls[i]
      a2,b2,c2=Ls[j]
      d=S(Sb(Mu(a1,b2),Mu(a2,b1)))
      if abs(ev(d))<1e-12:
        continue
      x=S(Mu(Sb(Mu(b1,c2),Mu(b2,c1)),Iv(d)))
      y=S(Mu(Sb(Mu(a2,c1),Mu(a1,c2)),Iv(d)))
      V.append((x,y))
  return V

def run(m):
  f=rd("integrand f(x,y)=")
  if f is None:
    return
  if m=='1':
    n=int(input("number of vertices="))
    V=[pt(input("vertex "+str(i+1)+" x,y=")) for i in range(n)]
    polyint(f,V)
  else:
    n=int(input("number of lines="))
    Ls=[lin(input("line "+str(i+1)+" (eg y=2x-1): ")) for i in range(n)]
    V=lines2poly(Ls)
    if len(V)!=n:
      print("Use vertices (mode 1)")
      return
    print("vertices: "+", ".join(["("+P(p[0])+","+P(p[1])+")" for p in hull(V)]))
    polyint(f,V)
