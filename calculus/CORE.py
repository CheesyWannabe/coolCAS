from INTC import *
if gc:
  gc.collect()

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

def mulx(a,b):
  P1=a[1] if a[0]=='+' else [a]
  P2=b[1] if b[0]=='+' else [b]
  return S(('+',[Mu(x,y) for x in P1 for y in P2]))

def expand(e):
  t=e[0]
  if t=='+':
    return S(('+',[expand(a) for a in e[1]]))
  if t=='^' and e[1][0]=='+' and isn(e[2]) and e[2][1].q==1 and e[2][1]>1 and e[2][1]<9:
    b=expand(e[1])
    r=b
    for i in range(e[2][1].p-1):
      r=mulx(r,b)
    return r
  if t=='*':
    r=ONE
    for a in e[1]:
      r=mulx(r,expand(a))
    return r
  return e

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

U=('k','u')

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
