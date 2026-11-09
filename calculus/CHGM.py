from INTM import *
if gc:
  gc.collect()

Y=('k','y')
UU=('k','u')
VV=('k','v')

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

def put(e,xe,ye):
  return S(subst(rep(e,Y,ye),xe))

def hasx(e):
  return not free(e)

def hasy(e):
  return P(e).find('y')>=0

def at(e,x,y):
  KV['y']=y
  return ev(e,x)

def inverse(u,v):
  for a,b,A_,B_ in ((u,v,UU,VV),(v,u,VV,UU)):
    try:
      if not hasy(a) and hasx(a):
        xe=inv(a,A_)
        if xe is not None:
          xe=S(xe)
          w=S(rep(subst(b,xe),Y,X))
          ye=inv(w,B_)
          if ye is not None:
            return xe,S(ye)
      if not hasx(a) and hasy(a):
        ye=inv(S(rep(a,Y,X)),A_)
        if ye is not None:
          ye=S(ye)
          w=put(b,X,ye)
          xe=inv(w,B_)
          if xe is not None:
            return S(xe),ye
    except Exception:
      pass
  try:
    p=[S(D(u)),S(dy(u)),S(D(v)),S(dy(v))]
    for q in p:
      if not numc(q):
        return None
    a,b,c,d=p
    u0=put(u,ZERO,ZERO)
    v0=put(v,ZERO,ZERO)
    dt=S(Sb(Mu(a,d),Mu(b,c)))
    U1=Sb(UU,u0)
    V1=Sb(VV,v0)
    xe=S(Mu(Sb(Mu(d,U1),Mu(b,V1)),Iv(dt)))
    ye=S(Mu(Sb(Mu(a,V1),Mu(c,U1)),Iv(dt)))
    return xe,ye
  except Exception:
    return None

def run(f,u,v):
  J=tidy(S(Sb(Mu(D(u),dy(v)),Mu(dy(u),D(v)))))
  J4=(D(u),dy(u),D(v),dy(v))
  print("J=d(u,v)/d(x,y)="+P(J))
  print("dA = du dv / |J|")
  o=input("order (blank=dudv)=").replace(' ','')
  if o=="":
    o="dudv"
  od=[c for c in o if c in 'uv']
  if len(od)!=2 or od[0]==od[1]:
    print("order must be dudv or dvdu")
    return
  ALIAS['u']='x'
  ALIAS['v']='y'
  try:
    lm=[]
    for w in od:
      lo=S(parse(input(w+" from=")))
      hi=S(parse(input(w+" to=")))
      lm.append((lo,hi))
  finally:
    ALIAS.clear()
  hs=input("point in region x,y (blank=auto)=")
  hint=None
  if hs!="":
    a,b=hs.split(',')
    hint=(ev(S(parse(a))),ev(S(parse(b))))
  ur=None
  vr=None
  KV['u']=0.53
  for i in range(2):
    try:
      r=(ev(rep(rep(lm[i][0],X,ONE),Y,ONE),1.0),ev(rep(rep(lm[i][1],X,ONE),Y,ONE),1.0))
    except Exception:
      r=(0.5,1.5)
    if od[i]=='u':
      ur=r
    else:
      vr=r
  uc=(ur[0]+ur[1])/2
  vc=(vr[0]+vr[1])/2
  iv=inverse(u,v)
  p0=None
  if iv is not None:
    try:
      KV['u']=uc
      KV['v']=vc
      xv=ev(iv[0])
      yv=ev(iv[1])
      if abs(at(u,xv,yv)-uc)+abs(at(v,xv,yv)-vc)>1e-7*(1+abs(uc)+abs(vc)):
        iv=None
      elif hint and abs(xv-hint[0])+abs(yv-hint[1])>1e-6+0.5*(abs(xv)+abs(yv)):
        print("(inverse is on another branch)")
        iv=None
      else:
        p0=(xv,yv)
    except Exception:
      iv=None
  F_=None
  if p0 is None:
    F_=__import__('CHGF')
    p0=F_.start(u,v,J4,uc,vc,hint)
  if p0 is None:
    print("Can't find region in x,y")
    return
  sg=1
  try:
    if at(J,p0[0],p0[1])<0:
      sg=-1
  except Exception:
    pass
  Ja=S(Mu(N(sg),J))
  h=tidy(S(Mu(f,Iv(Ja))))
  print("f/|J| = "+P(h))
  g=None
  if iv is not None:
    xe,ye=iv
    print("x = "+P(xe))
    print("y = "+P(ye))
    g=tidy(put(h,xe,ye))
  if g is None:
    print("(no inverse; fitting g)")
    if F_ is None:
      F_=__import__('CHGF')
    g=F_.fit(h,u,v,J4,ur,vr,p0)
  if g is None:
    print("Can't write g(u,v) exactly")
    F_.numeric(f,u,v,J4,Ja,od,lm,p0)
    return
  F_=None
  print("g(u,v) = "+P(g))
  ws=["int["+od[i]+"="+P(lm[i][0]).replace('x','u').replace('y','v')+".."+P(lm[i][1]).replace('x','u').replace('y','v')+"]" for i in (1,0)]
  print(ws[0]+" "+ws[1])
  print(" ("+P(g)+") d"+od[0]+" d"+od[1])
  gE=S(rep(rep(g,UU,X),VV,Y))
  oe=['x' if w=='u' else 'y' for w in od]
  MAP['x']='u'
  MAP['y']='v'
  try:
    report(gE,oe,lm)
  finally:
    MAP.clear()

def main(s):
  f=readf(s)
  if f is None:
    return
  u=readf(input("u="))
  v=readf(input("v="))
  if u is None or v is None:
    return
  if gc:
    gc.collect()
  run(f,u,v)

