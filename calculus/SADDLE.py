from CORE import *
if gc:
  gc.collect()

Y=('k','y')

def dy(e):
  return S(unsw(D(S(sw(e,'y'))),'y'))

def tidy(e):
  try:
    r=expand(e)
    if len(P(r))<len(P(e)):
      return r
  except Exception:
    pass
  return e

def at(e,x,y):
  KV['y']=y
  return ev(e,x)

def grad(G,x,y):
  return abs(at(G[0],x,y))+abs(at(G[1],x,y))

def newton(G,H,x,y,lim):
  for it in range(60):
    try:
      a=at(G[0],x,y)
      b=at(G[1],x,y)
      A=at(H[0],x,y)
      B=at(H[1],x,y)
      C=at(H[2],x,y)
    except Exception:
      return None
    h=A*A+2*B*B+C*C
    dt=A*C-B*B
    if abs(dt)>1e-12*h:
      dx=(C*a-B*b)/dt
      dv=(A*b-B*a)/dt
    else:
      p=A*A+B*B
      q=B*(A+C)
      r=B*B+C*C
      lm=1e-10*h+1e-300
      p+=lm
      r+=lm
      u=A*a+B*b
      w=B*a+C*b
      dt=p*r-q*q
      if dt==0:
        return None
      dx=(r*u-q*w)/dt
      dv=(p*w-q*u)/dt
    x-=dx
    y-=dv
    if abs(x)>lim or abs(y)>lim:
      return None
    if abs(dx)+abs(dv)<1e-14*(1+abs(x)+abs(y)):
      break
  try:
    g=grad(G,x,y)
  except Exception:
    return None
  if g>1e-7:
    return None
  return x,y

def snap(v):
  s=ident(v)
  if s is None:
    r=rat(v,1e-3*max(1,abs(v)),100)
    if r is None:
      return None
    s=str(r)
  try:
    return S(parse(s))
  except Exception:
    return None

def put(e,xe,ye):
  return S(subst(rep(e,Y,ye),xe))

def short(e,v):
  s=P(lncomb(e))
  if len(s)>24:
    return str(v)
  return s

G0=[]
H0=[]

def kind(f,x,y,Dv,A,sc):
  if Dv<-1e-9*sc:
    return "SADDLE POINT"
  if Dv>1e-9*sc:
    if A>0:
      return "local min"
    return "local max"
  f0=at(f,x,y)
  up=0
  dn=0
  for r in (0.01,0.001):
    for k in range(24):
      t=k*pi/12+0.1
      try:
        d=at(f,x+r*cos(t),y+r*sin(t))-f0
      except Exception:
        continue
      if d>1e-13*max(1,abs(f0)):
        up+=1
      elif d<-1e-13*max(1,abs(f0)):
        dn+=1
  w="D=0, test fails;"
  try:
    B=at(H0[1],x,y)
    C=at(H0[2],x,y)
    for vx,vy in ((-B,A),(C,-B)):
      m=abs(vx)+abs(vy)
      if m==0:
        continue
      for t in (0.01,-0.01,0.003,-0.003,0.001,-0.001):
        d=at(f,x+t*vx/m,y+t*vy/m)-f0
        if d>1e-13*max(1,abs(f0)):
          up+=1
        elif d<-1e-13*max(1,abs(f0)):
          dn+=1
    vs=[(-B,A),(C,-B),(1,0),(0,1)]
    for vx,vy in vs:
      m=abs(vx)+abs(vy)
      if m==0:
        continue
      cu=True
      for h in (0.01,-0.01,0.03):
        if grad(G0,x+h*vx/m,y+h*vy/m)>1e-12:
          cu=False
      if cu:
        w="ON A LINE/CURVE of critical pts;"
        break
  except Exception:
    pass
  if up and dn:
    return w+" looks like SADDLE"
  if up:
    return w+" looks like min"
  if dn:
    return w+" looks like max"
  return w+" flat here"

def main():
  print("SADDLE / critical points")
  print("f(x,y), use x y ^ sqrt e")
  while True:
    s=input("f(x,y)=")
    if s=="":
      break
    f=readf(s)
    if f is None:
      continue
    try:
      fx=tidy(D(f))
      fy=tidy(dy(f))
      fxx=tidy(D(fx))
      fxy=tidy(dy(fx))
      fyy=tidy(dy(fy))
    except MemoryError:
      print("Out of memory")
      continue
    except Exception as ex:
      print("Error: "+repr(ex))
      continue
    print("fx="+P(fx))
    print("fy="+P(fy))
    print("fxx="+P(fxx))
    print("fyy="+P(fyy))
    print("fxy="+P(fxy))
    r=input("search |x|,|y|<= (blank=5)=")
    try:
      R=5.0
      if r!="":
        R=abs(ev(S(parse(r))))
    except Exception:
      print("Bad range")
      continue
    print("searching...")
    G=(fx,fy)
    H=(fxx,fxy,fyy)
    G0[:]=G
    H0[:]=H
    pts=[]
    n=9
    for i in range(n):
      for j in range(n):
        x0=R*(2*(i+0.5)/n-1)+0.0123*R
        y0=R*(2*(j+0.5)/n-1)-0.0071*R
        p=newton(G,H,x0,y0,20*R+20)
        if p is None:
          continue
        x,y=p
        if abs(x)>R*1.0001 or abs(y)>R*1.0001:
          continue
        xe=snap(x)
        ye=snap(y)
        if xe is not None and ye is not None:
          xv=ev(xe)
          yv=ev(ye)
          try:
            if grad(G,xv,yv)<=max(1e-10,grad(G,x,y)):
              x,y=xv,yv
            else:
              xe=ye=None
          except Exception:
            xe=ye=None
        new=True
        for q in pts:
          if abs(q[0]-x)+abs(q[1]-y)<1e-5*(1+abs(x)+abs(y)):
            new=False
            break
        if new:
          pts.append((x,y,xe,ye))
    if not pts:
      print("No critical points found")
      print("in this box")
      continue
    pts=sorted(pts,key=lambda q:(q[0],q[1]))
    ns=0
    shown=0
    cv=0
    for x,y,xe,ye in pts:
      A=at(fxx,x,y)
      B=at(fxy,x,y)
      C=at(fyy,x,y)
      Dv=A*C-B*B
      sc=max(1,abs(A*C),B*B)
      k=kind(f,x,y,Dv,A,sc)
      if "CURVE" in k:
        cv+=1
        if cv>1:
          continue
      if "SADDLE" in k:
        ns+=1
      shown+=1
      if shown>16:
        continue
      fv=at(f,x,y)
      if xe is not None:
        try:
          cs="("+P(xe)+", "+P(ye)+")"
          fs=short(put(f,xe,ye),fv)
          ds=short(S(Sb(Mu(put(fxx,xe,ye),put(fyy,xe,ye)),Pw(put(fxy,xe,ye),N(2)))),Dv)
        except Exception:
          xe=None
      if xe is None:
        cs="(%.10g, %.10g)" % (x,y)
        fs="%.10g" % fv
        ds="%.6g" % Dv
      print(cs+": "+k)
      print("  f="+fs+"  D="+ds)
    if shown>16:
      print("...+"+str(shown-16)+" more (smaller box?)")
    if cv>1:
      print("Infinitely many critical")
      print("points on a curve; "+str(cv)+" found")
    print(str(ns)+" saddle point(s) of "+str(shown))

main()
