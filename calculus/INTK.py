from INTJ import *
if gc:
  gc.collect()

U=('k','u')
MAP={}

def mv(v):
  return MAP.get(v,v)

def mp(e):
  s=P(e)
  if not MAP:
    return s
  return ''.join([MAP.get(ch,ch) for ch in s])

def rep(e,o,n):
  if e==o:
    return n
  t=e[0]
  if t=='+' or t=='*':
    return (t,[rep(a,o,n) for a in e[1]])
  if t=='^':
    return ('^',rep(e[1],o,n),rep(e[2],o,n))
  if t=='f':
    return ('f',e[1],rep(e[2],o,n))
  return e

def sw(e,v):
  if v=='x':
    return e
  return rep(rep(e,X,U),('k',v),X)

def unsw(e,v):
  if v=='x':
    return e
  return rep(rep(e,X,('k',v)),U,X)

def vs(e,out):
  t=e[0]
  if t=='x':
    out.append('x')
  elif t=='k':
    out.append(e[1])
  elif t=='+' or t=='*':
    for a in e[1]:
      vs(a,out)
  elif t=='^':
    vs(e[1],out)
    vs(e[2],out)
  elif t=='f':
    vs(e[2],out)
  return out

def step(G,v,a,b):
  f=S(sw(G,v))
  try:
    Fx=I(f)
  except Exception:
    Fx=None
  if Fx is None or not chk(f,Fx):
    print("d"+mv(v)+": no elementary F")
    return None
  print("d"+mv(v)+": F="+mp(lncomb(S(unsw(Fx,v))))+"+C")
  ex,nq,st=definite(f,Fx,S(sw(a,v)),S(sw(b,v)))
  if ex is None:
    print(" "+mv(v)+"="+mp(a)+".."+mp(b)+": no exact value")
    return None
  r=S(unsw(ex,v))
  try:
    r2=expand(r)
    if len(P(r2))<len(P(r)):
      r=r2
  except Exception:
    pass
  return r

XV=[0.0]
def setv(v,p):
  if v=='x':
    XV[0]=p
  else:
    KV[v]=p

def nodes(m):
  h=1.0/(2**m)
  L=[]
  k=0
  while True:
    t=k*h
    et=exp(t)
    u=pi/4*(et-1/et)
    if u>30:
      break
    eu=exp(u)
    w=pi/4*(et+1/et)/((eu+1/eu)/2)**2*h
    if k>0 and w<1e-40:
      break
    L.append((2/(exp(2*u)+1),w,k==0))
    k+=1
  return L

def nint(f,od,lm,i,ND):
  v=od[i]
  a=ev(lm[i][0],XV[0])
  b=ev(lm[i][1],XV[0])
  c=(a+b)/2
  hw=(b-a)/2
  s=0.0
  for d,w,mid in ND:
    if mid:
      ps=[c]
    else:
      ps=[b-hw*d,a+hw*d]
    for p in ps:
      setv(v,p)
      try:
        if i==0:
          s+=w*ev(f,XV[0])
        else:
          s+=w*nint(f,od,lm,i-1,ND)
      except Exception:
        pass
  return s*hw

W=('k','w')
KV['w']=0.47
MI={'exp':'ln','ln':'exp','sin':'asin','asin':'sin','tan':'atan','atan':'tan','cos':'acos','acos':'cos'}

def nd(v):
  if v=='x':
    return X
  return ('k',v)

def inv(g,w):
  if g==X:
    return w
  t=g[0]
  if t=='+' or t=='*':
    dp=[a for a in g[1] if not free(a)]
    if len(dp)!=1:
      return None
    c=[a for a in g[1] if free(a)]
    if not c:
      return inv(dp[0],w)
    if t=='+':
      return inv(dp[0],Sb(w,('+',c)))
    return inv(dp[0],Mu(w,Iv(prod(c))))
  if t=='^':
    b,k=g[1],g[2]
    if free(k):
      return inv(b,Pw(w,Iv(k)))
    if free(b):
      return inv(k,Mu(F('ln',w),Iv(F('ln',b))))
    return None
  if t=='f' and g[1] in MI:
    return inv(g[2],F(MI[g[1]],w))
  return None




