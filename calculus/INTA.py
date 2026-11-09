try:
  import gc
except ImportError:
  gc=None
from math import *
def gcd(a,b):
  a,b=abs(a),abs(b)
  while b:
    a,b=b,a%b
  return a
class R:
  def __init__(s,p,q=1):
    if q==0:
      raise ZeroDivisionError("div 0")
    if q<0:
      p,q=-p,-q
    g=gcd(p,q)
    if g>1:
      p,q=p//g,q//g
    s.p=p
    s.q=q
  def __add__(s,o):
    o=rr(o)
    return R(s.p*o.q+o.p*s.q,s.q*o.q)
  def __sub__(s,o):
    o=rr(o)
    return R(s.p*o.q-o.p*s.q,s.q*o.q)
  def __mul__(s,o):
    o=rr(o)
    return R(s.p*o.p,s.q*o.q)
  def __truediv__(s,o):
    o=rr(o)
    return R(s.p*o.q,s.q*o.p)
  def __neg__(s):
    return R(-s.p,s.q)
  def __eq__(s,o):
    o=rr(o)
    return s.p==o.p and s.q==o.q
  def __ne__(s,o):
    return not s.__eq__(o)
  def __lt__(s,o):
    o=rr(o)
    return s.p*o.q<o.p*s.q
  def __gt__(s,o):
    o=rr(o)
    return s.p*o.q>o.p*s.q
  def f(s):
    return s.p/s.q
  def __str__(s):
    if s.q==1:
      return str(s.p)
    return str(s.p)+"/"+str(s.q)
def rr(o):
  if isinstance(o,R):
    return o
  return R(o)
def rpow(r,k):
  if k<0:
    return R(r.q**(-k),r.p**(-k))
  return R(r.p**k,r.q**k)
X=('x',)
PI=('k','pi')
EC=('k','e')
def N(r):
  return ('n',rr(r))
ZERO=N(0)
ONE=N(1)
HALF=N(R(1,2))
def A(*a):
  return ('+',list(a))
def Mu(*a):
  return ('*',list(a))
def Pw(b,x):
  return ('^',b,x)
def F(n,a):
  return ('f',n,a)
def Ng(a):
  return ('*',[N(-1),a])
def Sb(a,b):
  return ('+',[a,Ng(b)])
def Iv(a):
  return ('^',a,N(-1))
def Sq(a):
  return ('^',a,HALF)
def prod(L):
  if not L:
    return ONE
  if len(L)==1:
    return L[0]
  return ('*',L)
def isn(e,v=None):
  if e[0]!='n':
    return False
  return v is None or e[1]==v
KV={'y':0.41,'z':0.29,'u':0.53}

def hassym(e):
  t=e[0]
  if t=='k':
    return e[1]!='pi' and e[1]!='e'
  if t=='n' or t=='x':
    return False
  if t=='+' or t=='*':
    for a in e[1]:
      if hassym(a):
        return True
    return False
  if t=='^':
    return hassym(e[1]) or hassym(e[2])
  return hassym(e[2])

def numc(e):
  return free(e) and not hassym(e)

def free(e):
  t=e[0]
  if t=='x':
    return False
  if t=='n' or t=='k':
    return True
  if t=='+' or t=='*':
    for a in e[1]:
      if not free(a):
        return False
    return True
  if t=='^':
    return free(e[1]) and free(e[2])
  return free(e[2])
def ev(e,x=0.0):
  t=e[0]
  if t=='n':
    return e[1].f()
  if t=='x':
    return x
  if t=='k':
    if e[1]=='pi':
      return pi
    if e[1]=='e':
      return exp(1)
    return KV[e[1]]
  if t=='+':
    s=0.0
    for a in e[1]:
      s+=ev(a,x)
    return s
  if t=='*':
    s=1.0
    for a in e[1]:
      s*=ev(a,x)
    return s
  if t=='^':
    b=ev(e[1],x)
    if e[2][0]=='n':
      r=e[2][1]
      if r.q==1:
        return b**r.p
      if b<0:
        if r.q%2==1:
          v=(-b)**r.f()
          if r.p%2:
            return -v
          return v
        raise ValueError("dom")
      return b**r.f()
    p=ev(e[2],x)
    if b<0:
      raise ValueError("dom")
    return b**p
  n=e[1]
  a=ev(e[2],x)
  if n=='sin':
    return sin(a)
  if n=='cos':
    return cos(a)
  if n=='tan':
    return tan(a)
  if n=='exp':
    return exp(a)
  if n=='ln':
    if a<=0:
      raise ValueError("dom")
    return log(a)
  if n=='abs':
    return abs(a)
  if n=='asin':
    return asin(a)
  if n=='acos':
    return acos(a)
  if n=='atan':
    return atan(a)
  raise ValueError("fn")
def split(a):
  if a[0]=='*' and a[1][0][0]=='n':
    r=a[1][1:]
    if len(r)==1:
      return a[1][0][1],r[0]
    return a[1][0][1],('*',r)
  return R(1),a
def lcm(a,b):
  return a*b//gcd(a,b)
def negx(a):
  if a[0]=='n':
    return N(-a[1])
  k,r=split(a)
  if k==-1:
    return r
  if r[0]=='*':
    return ('*',[N(-k)]+r[1])
  return ('*',[N(-k),r])
def P(e):
  t=e[0]
  if t=='n':
    return str(e[1])
  if t=='x':
    return 'x'
  if t=='k':
    return e[1]
  if t=='f':
    return Pf(e)
  if t=='^':
    return Pp(e)
  if t=='*':
    return Pm(e[1])
  s=''
  Ls=[P(a) for a in e[1]]
  if Ls[0][0]=='-':
    for i in range(1,len(Ls)):
      if Ls[i][0]!='-':
        Ls=[Ls[i]]+Ls[:i]+Ls[i+1:]
        break
  for q in Ls:
    if s and q[0]!='-':
      s+='+'
    s+=q
  return s
def wrap(e):
  q=P(e)
  if e[0]=='+':
    return '('+q+')'
  return q
def jn(L,start):
  s=start
  for a in L:
    q=wrap(a)
    if s=='':
      s=q
    elif s==start and not q[0].isdigit():
      s+=q
    else:
      s+='*'+q
  return s
CLEAR=[1]

def Pm(L):
  c=R(1)
  num=[]
  den=[]
  for a in L:
    if a[0]=='n':
      c=c*a[1]
    elif a[0]=='f' and a[1]=='exp' and ((a[2][0]=='n' and a[2][1]<0) or split(a[2])[0]<0):
      den.append(('f','exp',negx(a[2])))
    elif a[0]=='^' and a[2][0]=='n' and a[2][1]<0:
      if a[2][1]==-1:
        den.append(a[1])
      else:
        den.append(('^',a[1],N(-a[2][1])))
    else:
      num.append(a)
  sg=''
  if c.p<0:
    sg='-'
  p=abs(c.p)
  st=''
  if p!=1 or not num:
    st=str(p)
  ns=jn(num,st)
  if c.q==1 and not den:
    return sg+ns
  st=''
  if c.q!=1:
    st=str(c.q)
  ds=jn(den,st)
  if len(den)+(c.q!=1)>1:
    ds='('+ds+')'
  elif CLEAR[0] and ('^' in ds or '*' in ds) and not (ds[0]=='(' and ds[-1]==')' and ds.count('(')==1):
    ds='('+ds+')'
  if CLEAR[0] and len(num)+(p!=1)>1:
    ns='('+ns+')'
  return sg+ns+'/'+ds
def Pp(e):
  b,x=e[1],e[2]
  if x[0]=='n':
    r=x[1]
    if r<0:
      return Pm([e])
    if r==R(1,2):
      return 'sqrt('+P(b)+')'
  bs=P(b)
  if b[0] in ('+','*','^') or (b[0]=='n' and (b[1]<0 or b[1].q!=1)) or (b[0]=='f' and b[1]=='exp'):
    bs='('+bs+')'
  xs=P(x)
  if not (x[0]=='x' or x[0]=='k' or (x[0]=='n' and x[1].q==1 and not x[1]<0)):
    xs='('+xs+')'
  return bs+'^'+xs
def Pf(e):
  n,a=e[1],e[2]
  if n=='abs':
    return '|'+P(a)+'|'
  if n=='ln' and a[0]=='f' and a[1]=='abs':
    return 'ln|'+P(a[2])+'|'
  if n=='exp':
    q=P(a)
    if q=='1':
      return 'e'
    if a[0]=='x' or a[0]=='k' or (a[0]=='n' and a[1].q==1 and not a[1]<0):
      return 'e^'+q
    return 'e^('+q+')'
  return n+'('+P(a)+')'
