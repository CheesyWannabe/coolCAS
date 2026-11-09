from INTI import *
if gc:
  gc.collect()
ALIAS={}
FN=['arcsin','arccos','arctan','asin','acos','atan','sqrt','cbrt','sin','cos','tan','sec','csc','cot','exp','abs','log','ln','pi']
def tok(s):
  s=s.replace('**','^').replace(' ','').lower()
  T=[]
  i=0
  while i<len(s):
    ch=s[i]
    if ch.isdigit() or ch=='.':
      j=i
      while j<len(s) and (s[j].isdigit() or s[j]=='.'):
        j+=1
      T.append(('n',s[i:j]))
      i=j
      continue
    if ch.isalpha():
      w=None
      for f in FN:
        if s[i:i+len(f)]==f:
          w=f
          break
      if w is None:
        if ch in ALIAS:
          w=ALIAS[ch]
          T.append(('w',w))
          i+=1
          continue
        if ch not in 'xyze':
          raise ValueError("bad name")
        w=ch
      T.append(('w',w))
      i+=len(w)
      continue
    if ch in '+-*/^()':
      T.append(('o',ch))
      i+=1
      continue
    raise ValueError("bad char")
  out=[]
  for t in T:
    if out:
      p=out[-1]
      ef=p[0]=='n' or (p[0]=='w' and (p[1] in ('x','y','z','e','pi') or p[1][0]=='#')) or p==('o',')')
      sf=t[0]=='w' or t==('o','(')
      if ef and sf:
        out.append(('o','*'))
    out.append(t)
  return out
def parse(s):
  T=tok(s)
  ps=[0]
  def pk():
    if ps[0]<len(T):
      return T[ps[0]]
    return None
  def tk():
    t=T[ps[0]]
    ps[0]+=1
    return t
  def expr():
    a=term_()
    while pk() in (('o','+'),('o','-')):
      o=tk()
      b=term_()
      a=A(a,b) if o[1]=='+' else Sb(a,b)
    return a
  def term_():
    a=unary()
    while pk() in (('o','*'),('o','/')):
      o=tk()
      b=unary()
      a=Mu(a,b) if o[1]=='*' else Mu(a,Iv(b))
    return a
  def unary():
    if pk()==('o','-'):
      tk()
      return Ng(unary())
    if pk()==('o','+'):
      tk()
      return unary()
    return power()
  def power():
    a=atom()
    if pk()==('o','^'):
      tk()
      return Pw(a,unary())
    return a
  def atom():
    t=tk()
    if t[0]=='n':
      s=t[1]
      if '.' in s:
        w,f=s.split('.')
        return N(R(int((w or '0')+f),10**len(f)))
      return N(int(s))
    if t==('o','('):
      a=expr()
      if tk()!=('o',')'):
        raise ValueError("paren")
      return a
    if t[0]=='w':
      w=t[1]
      if w=='x':
        return X
      if w=='e':
        return EC
      if w=='y' or w=='z':
        return ('k',w)
      if w[0]=='#':
        return ('k',w[1:])
      if w=='pi':
        return PI
      if tk()!=('o','('):
        raise ValueError("need (")
      a=expr()
      if tk()!=('o',')'):
        raise ValueError("paren")
      if w=='sqrt':
        return Sq(a)
      if w=='cbrt':
        return Pw(a,N(R(1,3)))
      if w=='sec':
        return Iv(F('cos',a))
      if w=='csc':
        return Iv(F('sin',a))
      if w=='cot':
        return Mu(F('cos',a),Iv(F('sin',a)))
      if w=='log':
        w='ln'
      if w[:3]=='arc':
        w='a'+w[3:]
      return F(w,a)
    raise ValueError("syntax")
  e=expr()
  if ps[0]!=len(T):
    raise ValueError("syntax")
  return e
ep=1.0
while 1.0+ep/2>1.0:
  ep=ep/2
TOL=max(1e-7,ep*1000)
def chk(f,Fx):
  try:
    dF=D(Fx)
  except Exception:
    return False
  n=0
  for x in (0.31,0.77,1.37,2.13,-0.6,3.3,-1.7,0.05):
    try:
      a=ev(dF,x)
      b=ev(f,x)
    except Exception:
      continue
    if abs(a-b)>TOL*max(1,abs(b)):
      return False
    n+=1
  return n>=2
def quad(f,a,b,m):
  c=(a+b)/2
  hw=(b-a)/2
  h=1.0/(2**m)
  s=0.0
  bad=False
  k=0
  while True:
    t=k*h
    et=exp(t)
    u=pi/4*(et-1/et)
    if u>30:
      break
    eu=exp(u)
    w=pi/4*(et+1/et)/((eu+1/eu)/2)**2
    if k>0 and w<1e-18:
      break
    if k==0:
      ps=[c]
    else:
      d=hw*2/(exp(2*u)+1)
      ps=[b-d,a+d]
    for p in ps:
      try:
        s+=w*ev(f,p)
      except Exception:
        if min(abs(p-a),abs(p-b))>1e-6*abs(b-a):
          bad=True
    k+=1
  return s*h*hw,bad
def vat(G,v,sg):
  L=G[1] if G[0]=='+' else [G]
  out=[]
  x0=ev(v)
  for t in L:
    try:
      w=S(subst(t,v))
      ev(w)
      out.append(w)
    except Exception:
      vs=[ev(t,x0+sg*h*max(1,abs(x0))) for h in (1e-5,1e-7,1e-9)]
      if abs(vs[2])<1e-4 and abs(vs[2])<=abs(vs[0]):
        continue
      raise ValueError("limit")
  return S(('+',out))
def ident(v):
  t=1e-9*max(1,abs(v))
  if abs(v)<t:
    return "0"
  r=rat(v,t,1000)
  if r:
    return str(r)
  r=rat(v/pi,t/pi,1000)
  if r:
    return P(S(Mu(N(r),PI)))
  for c in (ONE,PI):
    u=v/ev(c)
    r=rat(u*u,2*abs(u)*t,1000)
    if r and r>0:
      w=S(Mu(Sq(N(r)),c))
      if w[0]!='n':
        if v<0:
          w=S(Ng(w))
        return P(w)
  for d in range(1,13):
    for n in (2,3,5,6,7,10):
      for b in range(1,13):
        for sg in (1,-1):
          r=v*d-sg*b*sqrt(n)
          a=int(floor(r+0.5))
          if a!=0 and abs(r-a)<d*t:
            return P(S(Mu(Iv(N(d)),A(N(a),Mu(N(sg*b),Sq(N(n)))))))
  r=rat(exp(v),exp(v)*t,1000) if abs(v)<30 else None
  if r and r!=1:
    return P(S(F('ln',N(r))))
  return None
def lncomb(e):
  if e[0]!='+':
    return e
  oth=[]
  ls=[]
  for a in e[1]:
    k,r=split(a)
    if r[0]=='f' and r[1]=='ln' and r[2][0]=='n' and r[2][1].q==1:
      ls.append((k,r[2][1].p))
    else:
      oth.append(a)
  if len(ls)<2:
    return e
  Dn=1
  for k,p in ls:
    Dn=Dn*k.q//gcd(Dn,k.q)
  nu=1
  de=1
  for k,p in ls:
    m=k*Dn
    if abs(m.p)>40:
      return e
    if m>0:
      nu*=p**m.p
    else:
      de*=p**(-m.p)
  a=R(nu,de)
  s=1
  if a<1:
    a=R(de,nu)
    s=-1
  if a!=1:
    oth.append(('*',[N(R(s,Dn)),('f','ln',N(a))]))
  if not oth:
    return ZERO
  if len(oth)==1:
    return oth[0]
  return ('+',oth)
def run(s):
  f=S(parse(s))
  try:
    Fx=I(f)
  except Exception:
    Fx=None
  if Fx is not None and not chk(f,Fx):
    Fx=None
  return f,Fx
def definite(f,Fx,a,b):
  av=ev(a)
  bv=ev(b)
  nq=None
  st=0
  tl=0.05
  try:
    nq,b1=quad(f,av,bv,6)
    n2,b2=quad(f,av,bv,5)
    d=abs(n2-nq)
    if b1 or b2 or not abs(nq)<1e15:
      st=0
    elif d<1e-9*max(1,abs(nq)):
      tl=1e-6
      st=2
    elif d<1e-4*max(1,abs(nq)):
      st=1
  except Exception:
    nq=None
  ex=None
  if Fx is not None:
    try:
      G=Fx
      cs=[]
      if Fx[0]=='*':
        cs=[w for w in Fx[1] if free(w)]
        sm=[w for w in Fx[1] if not free(w)]
        if len(sm)==1:
          G=sm[0]
        else:
          cs=[]
      ex=S(Mu(prod(cs),Sb(vat(G,b,-1),vat(G,a,1))))
      xv=ev(ex)
      if nq is not None and abs(xv-nq)>tl*max(1,abs(nq)):
        ex=None
    except Exception:
      ex=None
  return ex,nq,st
