from CORE import *
if gc:
  gc.collect()

Y=('k','y')

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

def put(e,xe,ye):
  return S(subst(rep(e,Y,ye),xe))

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

def res(F,x,y):
  return abs(at(F[0],x,y))+abs(at(F[1],x,y))

def newton2(F,J,x,y,lim):
  for it in range(60):
    try:
      a=at(F[0],x,y)
      b=at(F[1],x,y)
      A=at(J[0],x,y)
      B=at(J[1],x,y)
      C=at(J[2],x,y)
      E=at(J[3],x,y)
    except Exception:
      return None
    h=A*A+B*B+C*C+E*E
    dt=A*E-B*C
    if abs(dt)>1e-12*h:
      dx=(E*a-B*b)/dt
      dv=(A*b-C*a)/dt
    else:
      p=A*A+C*C
      q=A*B+C*E
      r=B*B+E*E
      lm=1e-10*h+1e-300
      p+=lm
      r+=lm
      u=A*a+C*b
      w=B*a+E*b
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
    if res(F,x,y)>1e-7:
      return None
  except Exception:
    return None
  return x,y

def jac(F):
  return (D(F[0]),dy(F[0]),D(F[1]),dy(F[1]))

def solve2(F,R,n=9,ok=None):
  J=jac(F)
  pts=[]
  for i in range(n):
    for j in range(n):
      x0=R*(2*(i+0.5)/n-1)+0.0123*R
      y0=R*(2*(j+0.5)/n-1)-0.0071*R
      p=newton2(F,J,x0,y0,20*R+20)
      if p is None:
        continue
      x,y=p
      if abs(x)>R*1.0001 or abs(y)>R*1.0001:
        continue
      xe=snap(x)
      ye=snap(y)
      if xe is not None and ye is not None:
        try:
          xv=ev(xe)
          yv=ev(ye)
          if res(F,xv,yv)<=max(1e-10,res(F,x,y)):
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
  return sorted(pts,key=lambda q:(q[0],q[1]))

def pstr(p):
  x,y,xe,ye=p
  if xe is not None:
    return "("+P(xe)+", "+P(ye)+")"
  return "(%.8g, %.8g)" % (x,y)

def val(e,p):
  x,y,xe,ye=p
  v=at(e,x,y)
  if xe is not None:
    try:
      w=S(lncomb(put(e,xe,ye)))
      if len(P(w))<=26 and abs(ev(w)-v)<1e-8*(1+abs(v)):
        return w,v
    except Exception:
      pass
  return None,v

def vstr(w,v):
  if w is not None:
    return P(w)
  g=ident(v)
  if g:
    return g+" (~%.10g)" % v
  return "%.10g" % v
