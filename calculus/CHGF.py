from CHGM import *

def ninv(u,v,J4,u0,v0,x,y):
  for it in range(50):
    try:
      a=at(u,x,y)-u0
      b=at(v,x,y)-v0
      A_=at(J4[0],x,y)
      B_=at(J4[1],x,y)
      C=at(J4[2],x,y)
      E=at(J4[3],x,y)
    except Exception:
      return None
    dt=A_*E-B_*C
    if dt==0:
      return None
    dx=(E*a-B_*b)/dt
    dv=(A_*b-C*a)/dt
    x-=dx
    y-=dv
    if abs(dx)+abs(dv)<1e-14*(1+abs(x)+abs(y)):
      break
  try:
    if abs(at(u,x,y)-u0)+abs(at(v,x,y)-v0)>1e-9*(1+abs(u0)+abs(v0)):
      return None
  except Exception:
    return None
  return x,y

def start(u,v,J4,u0,v0,hint):
  if hint:
    return ninv(u,v,J4,u0,v0,hint[0],hint[1])
  best=None
  for i in range(7):
    for j in range(7):
      p=ninv(u,v,J4,u0,v0,-4.9+1.63*i,-4.7+1.61*j)
      if p is None:
        continue
      sc=0 if (p[0]>=0 and p[1]>=0) else 1
      if best is None or sc<best[0]:
        best=(sc,p)
  if best:
    return best[1]
  return None

def solve(A_,M):
  n=len(A_)
  for c in range(n):
    p=max(range(c,n),key=lambda r:abs(M[r][c]))
    M[c],M[p]=M[p],M[c]
    if abs(M[c][c])<1e-300:
      return None
    for r in range(n):
      if r!=c:
        k=M[r][c]/M[c][c]
        M[r]=[M[r][i]-k*M[c][i] for i in range(n+1)]
  return [M[i][n]/M[i][i] for i in range(n)]

def fit(h,u,v,J4,ur,vr,hint):
  mons=[(i,j) for i in range(4) for j in range(4) if i+j<=3]
  rows=[]
  pt=hint
  for i in range(5):
    for j in range(5):
      u0=ur[0]+(ur[1]-ur[0])*(i+0.5)/5
      v0=vr[0]+(vr[1]-vr[0])*(j+0.37)/5
      p=start(u,v,J4,u0,v0,pt)
      if p is None:
        p=start(u,v,J4,u0,v0,None)
      if p is None:
        continue
      pt=p
      try:
        rows.append((u0,v0,at(h,p[0],p[1])))
      except Exception:
        pass
  if len(rows)<len(mons)+3:
    return None
  n=len(mons)
  M=[[0.0]*(n+1) for k in range(n)]
  for u0,v0,hv in rows:
    b=[u0**i*v0**j for i,j in mons]
    for r in range(n):
      for c in range(n):
        M[r][c]+=b[r]*b[c]
      M[r][n]+=b[r]*hv
  cf=solve(mons,M)
  if cf is None:
    return None
  terms=[]
  for k in range(n):
    if abs(cf[k])<1e-9:
      continue
    q=rat(cf[k],1e-7*max(1,abs(cf[k])),1000)
    if q is None:
      return None
    terms.append(Mu(N(q),Pw(UU,N(mons[k][0])),Pw(VV,N(mons[k][1]))))
  g=S(('+',terms)) if terms else ZERO
  for u0,v0,hv in rows:
    KV['u']=u0
    KV['v']=v0
    if abs(ev(g)-hv)>1e-7*(1+abs(hv)):
      return None
  return g

def numeric(f,u,v,J4,Ja,od,lm,p0):
  print("numeric, please wait...")
  ND=nodes(3)
  st=[p0]
  def gval(uu,vv):
    q=ninv(u,v,J4,uu,vv,st[0][0],st[0][1])
    if q is None:
      raise ValueError("inv")
    st[0]=q
    return at(f,q[0],q[1])/abs(at(Ja,q[0],q[1]))
  a0,b0=ev(lm[1][0]),ev(lm[1][1])
  s=0.0
  for d1,w1,m1 in ND:
    c1=(a0+b0)/2
    h1=(b0-a0)/2
    for o in ([c1] if m1 else [b0-h1*d1,a0+h1*d1]):
      KV['y']=o
      a1=ev(lm[0][0],o)
      b1=ev(lm[0][1],o)
      c2=(a1+b1)/2
      h2=(b1-a1)/2
      t=0.0
      for d2,w2,m2 in ND:
        for i_ in ([c2] if m2 else [b1-h2*d2,a1+h2*d2]):
          if od[0]=='u':
            t+=w2*gval(i_,o)
          else:
            t+=w2*gval(o,i_)
      s+=w1*t*h2
  s*=h1
  print("~ "+str(s))
  g=ident(s)
  if g:
    print("maybe = "+g)
