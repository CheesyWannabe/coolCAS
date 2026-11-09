from INTK import *
if gc:
  gc.collect()

Y=('k','y')

def readf(s):
  try:
    e=parse(s)
  except Exception:
    print("Syntax error: check")
    print(" brackets, * and names")
    return None
  try:
    return S(e)
  except MemoryError:
    print("Out of memory")
  except Exception as ex:
    print("Error: "+repr(ex))
  return None

def tidy(e):
  try:
    r=expand(e)
    if len(P(r))<len(P(e)):
      return r
  except Exception:
    pass
  return e

def dy(e):
  return S(unsw(D(S(sw(e,'y'))),'y'))

def at(e,x,y):
  KV['y']=y
  return ev(e,x)

PTS=((0.37,0.61),(1.23,0.29),(0.71,1.47),(2.1,1.1),(0.53,2.3),(1.7,0.83),(-0.6,0.45),(0.45,-0.7))

def zero(e):
  n=0
  for x,y in PTS:
    try:
      v=at(e,x,y)
    except Exception:
      continue
    if abs(v)>1e-8:
      return False
    n+=1
  return n>=3

def xfree(e,y0):
  if not hassym(e) or free(e):
    pass
  vs=[]
  for y in (0.61,1.29,0.47):
    w=[]
    for x in (0.37,1.23,2.1):
      try:
        w.append(at(e,x,y))
      except Exception:
        pass
    if len(w)>=2 and max(w)-min(w)>1e-8*(1+abs(w[0])):
      return False
  return True

def onlyy(e,x0):
  if free(e):
    return e
  e2=tidy(e)
  if free(e2):
    return e2
  return S(subst(e,x0))

def integ(e,v):
  KV['y']=0.41
  KV['u']=0.53
  f=S(sw(e,v))
  F_=one(f,v)
  if F_ is None and f[0]=='+':
    out=[]
    for t in f[1]:
      w=one(t,v)
      if w is None:
        return None
      out.append(w)
    F_=S(('+',out))
  if F_ is None:
    return None
  return S(unsw(F_,v))

def one(f,v):
  try:
    F_=I(f)
    if F_ is not None and not chk(f,F_):
      F_=None
  except Exception:
    F_=None
  if F_ is None and hassym(f):
    F_=scale(f,Y if v=='x' else U,v)
  return F_

def scale(f,K,v):
  for s in (K,Iv(K)):
    try:
      h=S(Mu(subst(f,Mu(X,Iv(s))),Iv(s)))
      if hassym(h):
        continue
      H=I(h)
      if H is None:
        continue
      F_=S(subst(H,Mu(X,s)))
      if chk(f,F_):
        print("(sub t="+P(S(unsw(Mu(X,s),v)))+")")
        return F_
    except Exception:
      pass
  return None

def put(e,xe,ye):
  return S(subst(rep(e,Y,ye),xe))

def main():
  print("POTENT: find f(x,y)")
  print("from fx, fy (+ f(x0,y0)=c)")
  while True:
    s=input("fx=")
    if s=="":
      break
    Pe=readf(s)
    if Pe is None:
      continue
    Qe=readf(input("fy="))
    if Qe is None:
      continue
    sa=input("x0 (blank=no cond)=")
    x0=y0=c0=None
    if sa!="":
      try:
        x0=S(parse(sa))
        y0=S(parse(input("y0=")))
        c0=S(parse(input("f(x0,y0)=")))
        ev(x0)
        ev(y0)
        ev(c0)
      except Exception:
        print("Can't read condition")
        continue
    if gc:
      gc.collect()
    try:
      solve(Pe,Qe,x0,y0,c0)
    except MemoryError:
      print("Out of memory")
    except Exception as ex:
      print("Error: "+repr(ex))

def solve(Pe,Qe,x0,y0,c0):
  Py=tidy(dy(Pe))
  Qx=tidy(D(Qe))
  print("Py="+P(Py))
  print("Qx="+P(Qx))
  if not zero(S(Sb(Py,Qx))):
    print("Py != Qx: not exact,")
    print("no such f(x,y) exists")
    return
  print("Py = Qx: f exists")
  F1=integ(Pe,'x')
  if F1 is None:
    print("Can't integrate fx dx")
    return
  F1=tidy(F1)
  print("int fx dx = "+P(F1))
  print("f = "+P(F1)+" + g(y)")
  gp=tidy(S(Sb(Qe,dy(F1))))
  if not xfree(gp,y0):
    print("g'(y) depends on x?")
    print("check fx, fy")
    return
  xs=x0
  if xs is None:
    xs=ONE
  try:
    gp=onlyy(gp,xs)
  except Exception:
    gp=onlyy(gp,N(R(1,2)))
  print("g'(y) = fy - d/dy(...)")
  print("      = "+P(gp))
  if isn(gp,0):
    g=ZERO
  else:
    g=integ(gp,'y')
    if g is None:
      print("Can't integrate g'(y)")
      return
  g=tidy(g)
  print("g(y) = "+P(g)+" + C")
  f=tidy(S(A(F1,g)))
  if not (zero(S(Sb(D(f),Pe))) and zero(S(Sb(dy(f),Qe)))):
    print("Check failed, sorry")
    return
  if x0 is None:
    print("f(x,y) = "+P(lncomb(f))+" + C")
    return
  v=put(f,x0,y0)
  C=S(Sb(c0,v))
  print("f("+P(x0)+","+P(y0)+") = "+P(lncomb(v))+" + C = "+P(c0))
  print("C = "+P(lncomb(C)))
  f=tidy(S(A(f,C)))
  print("f(x,y) = "+P(lncomb(f)))

main()
