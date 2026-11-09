from CORE import *
if gc:
  gc.collect()

Y=('k','y')
Z=('k','z')
TT=('k','#')
KV['#']=0.37
VX=[X,Y,Z]

def tidy(e):
  try:
    r=expand(e)
    if len(P(r))<len(P(e)):
      return r
  except Exception:
    pass
  return e

def ex(s):
  return S(parse(s))

def vec(s):
  return [ex(a) for a in s.split(',')]

def vs(v):
  return "("+", ".join([P(c) for c in v])+")"

def dv(e,i):
  if i==0:
    return tidy(S(D(e)))
  v='y' if i==1 else 'z'
  return tidy(S(unsw(D(S(sw(e,v))),v)))

def dsym(e,s):
  e1=rep(rep(e,X,TT),s,X)
  return S(rep(rep(D(S(e1)),X,s),TT,X))

def put3(e,p):
  e=rep(e,Y,p[1])
  if len(p)>2:
    e=rep(e,Z,p[2])
  return S(subst(e,p[0]))

def dot(a,b):
  return S(('+',[Mu(a[i],b[i]) for i in range(len(a))]))

def eqn(s):
  l,r=s.split('=')
  return S(Sb(parse(l),parse(r)))

def plane(n,p):
  lhs=dot(n,VX[:len(n)])
  return P(lhs)+" = "+P(dot(n,p))

def normal(E):
  n=[dv(E,i) for i in range(3)]
  for c in n:
    if not numc(c):
      raise ValueError("not a plane")
  return n,S(Ng(put3(E,[ZERO,ZERO,ZERO])))

def m_line():
  t=input("line: 1 two points 2 pt+dir =")
  Pp=vec(input("point P (x,y,z)="))
  if t=='2':
    d=vec(input("direction (a,b,c)="))
  else:
    Q=vec(input("point Q (x,y,z)="))
    d=[S(Sb(Q[i],Pp[i])) for i in range(3)]
  n,k=normal(eqn(input("plane (eg 2x+3y=1):")))
  print("r(t) = "+vs(Pp)+" + t*"+vs(d))
  nd=dot(n,d)
  r0=S(Sb(k,dot(n,Pp)))
  print("n.d = "+P(nd))
  if abs(ev(nd))<1e-12:
    if abs(ev(r0))<1e-12:
      print("line lies in the plane")
    else:
      print("parallel, no intersection")
    return
  t0=S(Mu(r0,Iv(nd)))
  print("t = "+P(t0))
  print("point = "+vs([S(A(Pp[i],Mu(t0,d[i]))) for i in range(3)]))

def ev3(e,x,y,z):
  KV['y']=y
  KV['z']=z
  return ev(e,x)

def gauss(M):
  n=len(M)
  for c in range(n):
    p=max(range(c,n),key=lambda r:abs(M[r][c]))
    if abs(M[p][c])<1e-14:
      return None
    M[c],M[p]=M[p],M[c]
    for r in range(n):
      if r!=c:
        f=M[r][c]/M[c][c]
        M[r]=[M[r][i]-f*M[c][i] for i in range(n+1)]
  return [M[i][n]/M[i][i] for i in range(n)]

def snap(v):
  s=ident(v)
  if s is None:
    q=rat(v,1e-9*max(1,abs(v)),1000)
    if q is None:
      return None
    s=str(q)
  try:
    return ex(s)
  except Exception:
    return None

def m_par():
  F=eqn(input("surface (eg x^2+z=4):"))
  n,k=normal(eqn(input("plane (eg 6x-y+2z=5):")))
  g=[dv(F,i) for i in range(3)]
  print("grad F = "+vs(g))
  print("need grad F = L*"+vs(n)+", F=0")
  H=[[dv(g[i],j) for j in range(3)] for i in range(3)]
  nn=[ev(c) for c in n]
  pts=[]
  for a in (-3.1,-0.7,1.3,3.7):
    for b in (-2.9,0.6,2.2,4.1):
      for c in (-2.3,0.9,3.3):
        x,y,z=a,b,c
        try:
          gv=[ev3(q,x,y,z) for q in g]
        except Exception:
          continue
        L=sum([gv[i]*nn[i] for i in range(3)])/sum([q*q for q in nn])
        ok=False
        for it in range(40):
          try:
            gv=[ev3(q,x,y,z) for q in g]
            Fv=ev3(F,x,y,z)
            M=[[ev3(H[i][j],x,y,z) for j in range(3)]+[-nn[i],-(gv[i]-L*nn[i])] for i in range(3)]
            M.append(gv+[0.0,-Fv])
          except Exception:
            break
          s=gauss(M)
          if s is None:
            break
          x+=s[0]
          y+=s[1]
          z+=s[2]
          L+=s[3]
          if abs(s[0])+abs(s[1])+abs(s[2])+abs(s[3])<1e-13*(1+abs(x)+abs(y)+abs(z)):
            ok=True
            break
        if not ok:
          continue
        if any([abs(p[0]-x)+abs(p[1]-y)+abs(p[2]-z)<1e-6 for p in pts]):
          continue
        pts.append((x,y,z,L))
  if not pts:
    print("no such point found")
    return
  for x,y,z,L in pts:
    e=[snap(x),snap(y),snap(z)]
    if None in e:
      print("point ~ (%.10g, %.10g, %.10g)" % (x,y,z))
      continue
    print("point = "+vs(e))
    Lv=snap(L)
    if Lv is not None:
      print("  (L = "+P(Lv)+")")
    print("  tangent plane: "+plane(n,e))

def m_dir():
  f=rd("f=")
  if f is None:
    return
  p=vec(input("point (x,y[,z])="))
  d=vec(input("direction (a,b[,c])="))
  n=len(p)
  g=[dv(f,i) for i in range(n)]
  print("grad f = "+vs(g))
  g0=[put3(q,p) for q in g]
  print("grad at point = "+vs(g0))
  L=S(Sq(dot(d,d)))
  print("|v| = "+P(L))
  print("unit u = "+vs([S(Mu(c,Iv(L))) for c in d]))
  r=S(Mu(dot(g0,d),Iv(L)))
  print("D_u f = grad.u = "+P(r))
  print("~ "+str(ev(r)))
  print("(max rate |grad| = "+P(S(Sq(dot(g0,g0))))+")")

def m_path():
  t=input("1 P(x,y,z)  2 grad parts =")
  if t=='2':
    g=[rd("dP/dx="),rd("dP/dy="),rd("dP/dz= (blank=2D)")]
    if g[2] is None:
      g=g[:2]
  else:
    Pf=rd("P(x,y,z)=")
    g=[dv(Pf,i) for i in range(3)]
    print("grad P = "+vs(g))
  ALIAS['t']='x'
  try:
    r=[ex(input(c+"(t)=")) for c in ('x','y','z')[:len(g)]]
  finally:
    ALIAS.clear()
  t0=ex(input("at t="))
  pt=[S(subst(c,t0)) for c in r]
  vel=[S(subst(D(c),t0)) for c in r]
  print("r(t0) = "+vs(pt))
  print("r'(t0) = "+vs(vel))
  g0=[put3(q,pt) for q in g]
  print("grad at r(t0) = "+vs(g0))
  res=dot(g0,vel)
  print("dP/dt = grad.r' = "+P(res))

def m_chain():
  ALIAS['u']='#u'
  ALIAS['v']='#v'
  ALIAS['w']='#w'
  try:
    z=rd("z(u,v,w)=")
  finally:
    ALIAS.clear()
  if z is None:
    return
  qs=[]
  for nm in ('u','v','w'):
    s=('k',nm)
    if P(z).find(nm)<0:
      continue
    e=rd(nm+"(x,y)= (blank=partials)")
    if e is not None:
      qs.append((s,nm,e,dv(e,0),dv(e,1)))
    else:
      qs.append((s,nm,None,rd(" d"+nm+"/dx="),rd(" d"+nm+"/dy=")))
  for k,var in ((3,'x'),(4,'y')):
    print("dz/d"+var+" = "+" + ".join(["z_"+q[1]+"*"+q[1]+"_"+var for q in qs]))
    out=[]
    for q in qs:
      zq=dsym(z,q[0])
      print(" z_"+q[1]+"="+P(zq)+", "+q[1]+"_"+var+"="+P(q[k]))
      out.append(Mu(zq,q[k]))
    res=S(('+',out))
    for q in qs:
      if q[2] is not None:
        res=rep(res,q[0],q[2])
    print("= "+P(tidy(expand(S(res)))))

def m_tan():
  E=eqn(input("surface (eg z=x^2+y^2):"))
  p=vec(input("point (a,b,c)="))
  n=[put3(dv(E,i),p) for i in range(3)]
  print("normal = grad = "+vs(n))
  print("tangent plane: "+plane(n,p))
  print("normal line: "+vs(p)+" + t*"+vs(n))

def rd(p):
  s=input(p)
  if s=="":
    return None
  return readf(s)

def main():
  print("VECTOR: gradients & geometry")
  print("1 line meets plane")
  print("2 tangent plane // a plane")
  print("3 directional derivative")
  print("4 rate along a path r(t)")
  print("5 chain rule z(u,v,w)")
  print("6 tangent plane at a point")
  M={'1':m_line,'2':m_par,'3':m_dir,'4':m_path,'5':m_chain,'6':m_tan}
  while True:
    m=input("mode (blank=quit)=")
    if m=="":
      break
    if m not in M:
      continue
    if gc:
      gc.collect()
    try:
      M[m]()
    except MemoryError:
      print("Out of memory")
    except Exception as ex_:
      print("Error: "+repr(ex_))
    ALIAS.clear()

main()
