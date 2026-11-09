from CORE import *
if gc:
  gc.collect()

LOG=[]
LH=[0]
AP=[ZERO,1,False]

def say(s):
  LOG.append(s)

def tidy(e):
  try:
    r=expand(e)
    if len(P(r))<len(P(e)):
      return r
  except Exception:
    pass
  return e

def apr():
  a,s,inf=AP
  if inf:
    return [1e6,1e4,1e3,100.0,30.0]
  av=ev(a)
  return [av+s*h for h in (1e-7,1e-5,1e-3,1e-2)]

def sgn(e):
  for x in apr():
    try:
      v=ev(e,x)
      if v>0:
        return 1
      if v<0:
        return -1
    except Exception:
      pass
  return 0

def deabs(e):
  t=e[0]
  if t=='+' or t=='*':
    return (t,[deabs(a) for a in e[1]])
  if t=='^':
    return ('^',deabs(e[1]),deabs(e[2]))
  if t=='f':
    g=deabs(e[2])
    if e[1]=='abs':
      s=sgn(g)
      if s:
        w=S(Mu(N(s),g))
        say("near the point |"+P(g)+"|="+P(w))
        return w
    return ('f',e[1],g)
  return e

def vpow(e):
  t=e[0]
  if t=='^':
    return not free(e[2]) or vpow(e[1])
  if t=='+' or t=='*':
    for a in e[1]:
      if vpow(a):
        return True
    return False
  if t=='f':
    return vpow(e[2])
  return False

def haslog(e):
  t=e[0]
  if t=='f':
    return e[1] in ('ln','atan','asin','acos') or haslog(e[2])
  if t=='+' or t=='*':
    for a in e[1]:
      if haslog(a):
        return True
    return False
  if t=='^':
    return haslog(e[1]) or haslog(e[2])
  return False

def hashalf(e):
  t=e[0]
  if t=='^':
    return (isn(e[2]) and e[2][1].q==2) or hashalf(e[1])
  if t=='+' or t=='*':
    for a in e[1]:
      if hashalf(a):
        return True
    return False
  if t=='f':
    return hashalf(e[2])
  return False

def isz(r):
  if r is None or r[0]!='v':
    return False
  try:
    return abs(ev(r[1]))<1e-13
  except Exception:
    return False

def sv(r):
  v=ev(r[1])
  if v>0:
    return 1
  return -1

def nd(e):
  fs=e[1] if e[0]=='*' else [e]
  nu=[]
  de=[]
  for f in fs:
    if f[0]=='^' and isn(f[2]) and f[2][1]<0:
      de.append(S(Pw(f[1],N(-f[2][1]))))
    elif f[0]=='f' and f[1]=='exp' and not free(f[2]) and split(f[2])[0]<0:
      de.append(S(F('exp',Ng(f[2]))))
    else:
      nu.append(f)
  return prod(nu),prod(de),len(de)>0

def lim(e,dp=0):
  if dp>30:
    return None
  e=S(e)
  if len(P(e))>400:
    return None
  if free(e):
    return ('v',e)
  a,s,inf=AP
  if e==X:
    if inf:
      return ('i',1)
    return ('v',a)
  if not inf and not vpow(e):
    try:
      v=S(subst(e,a))
      if numc(v) and abs(ev(v))<1e12:
        return ('v',v)
    except Exception:
      pass
  t=e[0]
  if t=='+':
    return lsum(e[1],dp)
  Nm,Dn,hd=nd(e)
  if hd:
    return quot(Nm,Dn,dp)
  if t=='*':
    return lprod(e[1],dp)
  if t=='^':
    return lpow(e,dp)
  if t=='f':
    return lfun(e,dp)
  return None

def lead(e):
  if free(e):
    if numc(e) and abs(ev(e))>1e-14:
      return (e,R(0))
    return None
  t=e[0]
  if t=='x':
    return (ONE,R(1))
  if t=='^' and isn(e[2]):
    L=lead(e[1])
    if not L:
      return None
    k=e[2][1]
    if k.q!=1 and not ev(L[0])>0:
      return None
    return (S(Pw(L[0],e[2])),L[1]*k)
  if t=='*':
    c=ONE
    p=R(0)
    for f in e[1]:
      L=lead(f)
      if not L:
        return None
      c=Mu(c,L[0])
      p=p+L[1]
    return (S(c),p)
  if t=='+':
    Ls=[]
    for a in e[1]:
      L=lead(a)
      if not L:
        return None
      Ls.append(L)
    pm=Ls[0][1]
    for L in Ls:
      if L[1]>pm:
        pm=L[1]
    c=S(('+',[L[0] for L in Ls if L[1]==pm]))
    if not numc(c) or abs(ev(c))<1e-12:
      return None
    return (c,pm)
  return None

def mono(L):
  return S(Mu(L[0],Pw(X,N(L[1]))))

def rlead(Nm,Dn):
  Ln=lead(Nm)
  Ld=lead(Dn)
  if not Ln or not Ld:
    return None
  if Ln[1]>0 and Ld[1]>0:
    say("inf/inf, dominant terms:")
    say(" ("+P(mono(Ln))+")/("+P(mono(Ld))+")")
  p=Ln[1]-Ld[1]
  c=S(Mu(Ln[0],Iv(Ld[0])))
  if p>0:
    return ('i',sv(('v',c)))
  if p<0:
    return ('v',ZERO)
  return ('v',c)

def lhop(Nm,Dn,dp,form):
  if LH[0]>=8:
    return None
  LH[0]+=1
  n2=tidy(D(Nm))
  d2=tidy(D(Dn))
  say(form+" -> L'Hopital:")
  say(" ("+P(n2)+")/("+P(d2)+")")
  if isn(d2,0):
    return None
  mk=len(LOG)
  r=lim(S(Mu(n2,Iv(d2))),dp+1)
  if r is not None and r[0] in ('v','i'):
    return r
  del LOG[mk-2:]
  say(form+": L'Hopital inconclusive")
  say(" (f'/g' has no limit)")
  if Nm[0]=='+':
    ts=[S(Mu(t,Iv(Dn))) for t in Nm[1]]
    say("split: "+" + ".join(["("+P(t)+")" for t in ts]))
    r=lsum(ts,dp+1)
    if r is not None and r[0] in ('v','i'):
      return r
  if Nm[0]=='*':
    fs=Nm[1]
    for i in range(len(fs)):
      q=S(Mu(fs[i],Iv(Dn)))
      rest=fs[:i]+fs[i+1:]
      m2=len(LOG)
      say("split: ("+P(q)+")*"+P(S(prod(rest))))
      r=lprod([q]+rest,dp+1)
      if r is not None and r[0] in ('v','i'):
        return r
      del LOG[m2:]
  return None

def quot(Nm,Dn,dp):
  a,s,inf=AP
  if inf:
    r=rlead(Nm,Dn)
    if r:
      return r
  LN=lim(Nm,dp+1)
  LD=lim(Dn,dp+1)
  if LN is None or LD is None:
    return None
  zn=isz(LN)
  zd=isz(LD)
  if zn and zd:
    if not inf:
      pn=toP(Nm)
      pd=toP(Dn)
      if pn and pd:
        g=pgcd(pn,pd)
        if pdeg(g)>0:
          say("0/0, factor & cancel ("+P(pexpr(g))+")")
          return lim(S(Mu(pexpr(pdiv(pn,g)[0]),Iv(pexpr(pdiv(pd,g)[0])))),dp+1)
    return lhop(Nm,Dn,dp,"0/0")
  if LN[0]=='i' and LD[0]=='i':
    return lhop(Nm,Dn,dp,"inf/inf")
  if zd:
    if LN[0] in ('v','i'):
      g=sgn(S(Mu(Nm,Iv(Dn))))
      if g:
        return ('i',g)
      return None
    return ('d',)
  if LD[0]=='i':
    if LN[0]=='v':
      return ('v',ZERO)
    if LN[0]=='b':
      say("squeeze: bounded/inf -> 0")
      return ('v',ZERO)
    return None
  if LN[0]=='i':
    if LD[0]=='v':
      return ('i',LN[1]*sv(LD))
    return ('d',)
  if LN[0]=='v' and LD[0]=='v':
    return ('v',S(Mu(LN[1],Iv(LD[1]))))
  if LN[0]=='b' and LD[0]=='v':
    return ('b',)
  return None

def lprod(fs,dp):
  Rs=[]
  for f in fs:
    r=lim(f,dp+1)
    if r is None or r[0]=='d':
      return None
    Rs.append(r)
  zs=[]
  ins=[]
  oth=[]
  ov=[]
  nb=0
  sg=1
  for i in range(len(fs)):
    r=Rs[i]
    if isz(r):
      zs.append(fs[i])
    elif r[0]=='i':
      ins.append(fs[i])
      sg*=r[1]
    elif r[0]=='b':
      nb+=1
    else:
      oth.append(fs[i])
      ov.append(r[1])
      sg*=sv(r)
  if zs and ins:
    if nb:
      return None
    Z=S(prod(zs))
    In=S(prod(ins))
    tr=[(Z,S(Iv(In))),(In,S(Iv(Z)))]
    if haslog(In):
      tr=[tr[1],tr[0]]
    for Nm,Dn in tr:
      mk=len(LOG)
      lh=LH[0]
      say("0*inf form, write as")
      say(" ("+P(Nm)+")/("+P(Dn)+")")
      r=quot(Nm,Dn,dp)
      if r is not None and r[0] in ('v','i'):
        if not ov:
          return r
        c=S(prod(ov))
        if r[0]=='v':
          return ('v',S(Mu(c,r[1])))
        return ('i',r[1]*sv(('v',c)))
      del LOG[mk:]
      LH[0]=lh
    return None
  if zs:
    if nb:
      say("squeeze: bounded*0 -> 0")
    return ('v',ZERO)
  if ins:
    if nb:
      return ('d',)
    return ('i',sg)
  if nb:
    return ('b',)
  return ('v',S(prod(ov)))

def lsum(ts,dp):
  Rs=[]
  for t in ts:
    r=lim(t,dp+1)
    if r is None or r[0]=='d':
      return None
    Rs.append(r)
  ps=[r for r in Rs if r[0]=='i' and r[1]>0]
  ng=[r for r in Rs if r[0]=='i' and r[1]<0]
  if ps and ng:
    return infinf(ts,dp)
  if ps:
    return ('i',1)
  if ng:
    return ('i',-1)
  for r in Rs:
    if r[0]=='b':
      return ('b',)
  return ('v',S(('+',[r[1] for r in Rs])))

def infinf(ts,dp):
  a,s,inf=AP
  e=S(('+',ts))
  if inf:
    L=lead(e)
    if L:
      say("inf-inf, dominant term: "+P(mono(L)))
      if L[1]>0:
        return ('i',sv(('v',L[0])))
      if L[1]==0:
        return ('v',L[0])
      return ('v',ZERO)
  ns=[]
  ds=[]
  hd=False
  for t in ts:
    n,d,h=nd(t)
    ns.append(n)
    ds.append(d)
    hd=hd or h
  if hd:
    nu=[]
    for i in range(len(ts)):
      nu.append(Mu(ns[i],prod(ds[:i]+ds[i+1:])))
    Nm=tidy(S(('+',nu)))
    Dn=tidy(S(prod(ds)))
    say("inf-inf, combine fractions:")
    say(" ("+P(Nm)+")/("+P(Dn)+")")
    return quot(Nm,Dn,dp)
  for i in range(len(ts)):
    if hashalf(ts[i]):
      A_=ts[i]
      B=S(('+',ts[:i]+ts[i+1:]))
      Nm=tidy(expand(S(Sb(Pw(A_,N(2)),Pw(B,N(2))))))
      Dn=S(Sb(A_,B))
      say("inf-inf, multiply by conjugate:")
      say(" ("+P(Nm)+")/("+P(Dn)+")")
      return quot(Nm,Dn,dp)
  A_=ts[0]
  B=S(('+',ts[1:]))
  say("inf-inf, factor out "+P(A_))
  r=lim(S(Mu(B,Iv(A_))),dp+1)
  if r is not None and r[0]=='i':
    rb=lim(B,dp+1)
    if rb is not None and rb[0]=='i':
      return rb
    return None
  if r is None or r[0]!='v':
    return None
  c=1+ev(r[1])
  if abs(c)<1e-12:
    return None
  ra=lim(A_,dp+1)
  if ra is None or ra[0]!='i':
    return None
  if c>0:
    return ('i',ra[1])
  return ('i',-ra[1])

def expres(L):
  if L is None:
    return None
  if L[0]=='v':
    return ('v',S(F('exp',L[1])))
  if L[0]=='i':
    if L[1]>0:
      return ('i',1)
    return ('v',ZERO)
  if L[0]=='b':
    return ('b',)
  return None

def lpow(e,dp):
  b,k=e[1],e[2]
  if free(k):
    Lb=lim(b,dp+1)
    if Lb is None:
      return None
    kv=ev(k)
    if Lb[0]=='v':
      if isz(Lb):
        if kv>0:
          return ('v',ZERO)
        return ('i',sgn(e))
      try:
        v=S(Pw(Lb[1],k))
        ev(v)
        return ('v',v)
      except Exception:
        return None
    if Lb[0]=='i':
      if kv>0:
        return ('i',sgn(e))
      return ('v',ZERO)
    if Lb[0]=='b' and kv>0:
      return ('b',)
    return None
  ex=S(Mu(k,F('ln',b)))
  if not free(b):
    say("power form: write as e^("+P(ex)+")")
  return expres(lim(ex,dp+1))

def lfun(e,dp):
  n,g=e[1],e[2]
  L=lim(g,dp+1)
  if L is None:
    return None
  if L[0]=='v':
    if n=='ln' and isz(L):
      return ('i',-1)
    try:
      v=S(Sf(n,L[1]))
      w=ev(v)
      if abs(w)>1e12:
        return ('i',sgn(e))
      return ('v',v)
    except Exception:
      if n=='tan':
        return ('i',sgn(e))
      return None
  if L[0]=='i':
    p=L[1]>0
    if n=='exp':
      if p:
        return ('i',1)
      return ('v',ZERO)
    if n=='ln' and p:
      return ('i',1)
    if n=='atan':
      return ('v',S(Mu(N(R(L[1],2)),PI)))
    if n=='sin' or n=='cos':
      return ('b',)
    if n=='tan':
      return ('d',)
    return None
  if L[0]=='b' and n in ('sin','cos','exp','atan'):
    return ('b',)
  return None

def nev(H,Fv):
  Pv=Fv[:]
  n=len(H)
  for m in range(1,n):
    for i in range(n-m):
      Pv[i]=(-H[i+m]*Pv[i]+H[i]*Pv[i+1])/(H[i]-H[i+m])
  return Pv[0]

def fat(e,h):
  a,s,inf=AP
  if inf:
    return ev(e,1/h)
  return ev(e,ev(a)+s*h)

def numlim(e):
  a,s,inf=AP
  dv=[]
  ovf=False
  bad=0
  for k in range(2,17 if inf else 9):
    h=10.0**(-k/2) if inf else 10.0**(-k)
    try:
      dv.append(fat(e,h))
    except OverflowError:
      ovf=True
      break
    except Exception:
      bad+=1
  if ovf and dv:
    return ('i',1 if dv[-1]>0 else -1)
  if len(dv)<3:
    if bad:
      return ('u',)
    return ('d',)
  if abs(dv[-1])>1e7 and abs(dv[-1])>10*abs(dv[-3]):
    return ('i',1 if dv[-1]>0 else -1)
  H=[]
  Fv=[]
  for k in range(3,12):
    h=0.1/2**k
    try:
      Fv.append(fat(e,h))
      H.append(h)
    except Exception:
      break
  if len(H)>=7:
    v1=nev(H[-7:-1],Fv[-7:-1])
    v2=nev(H[-6:],Fv[-6:])
    if abs(v1-v2)<1e-7*(1+abs(v2)):
      return ('v',v2,True)
  d=[dv[i+1]-dv[i] for i in range(len(dv)-1)]
  same=True
  for x in d:
    if x*d[-1]<=0:
      same=False
  if same and abs(d[-1])>0.3*abs(d[-2]) and abs(dv[-1])>abs(dv[0]):
    return ('i',1 if d[-1]>0 else -1)
  if abs(dv[-1]-dv[-2])<1e-4*(1+abs(dv[-1])):
    return ('v',dv[-1],False)
  return ('d',)

def side(e,a,s,inf):
  AP[0]=a
  AP[1]=s
  AP[2]=inf
  del LOG[:]
  LH[0]=0
  try:
    e2=deabs(e)
    r=lim(e2)
  except MemoryError:
    r=None
    LOG.append("(out of memory)")
  except Exception:
    r=None
  lg=LOG[:]
  try:
    nm=numlim(e)
  except Exception:
    nm=('d',)
  if nm[0]=='u':
    return ('u',),lg
  if r is not None and r[0]=='v' and nm[0]=='v':
    if abs(ev(r[1])-nm[1])>(1e-6 if nm[2] else 1e-2)*(1+abs(nm[1])):
      r=None
  elif r is not None and r[0]=='v' and nm[0]=='i':
    r=None
  elif r is not None and r[0]=='i' and nm[0]=='i' and r[1]!=nm[1]:
    r=None
  if r is None:
    lg.append("(no exact method, numeric)")
    if nm[0]=='v':
      return ('n',nm[1]),lg
    if nm[0]=='i':
      return nm,lg
    return ('d',),lg
  return r,lg

def txt(r):
  t=r[0]
  if t=='v':
    return P(lncomb(r[1]))
  if t=='i':
    if r[1]>0:
      return "+inf"
    return "-inf"
  if t=='n':
    return "~"+str(r[1])
  if t=='b':
    return "DNE (oscillates)"
  if t=='u':
    return "undefined"
  return "DNE"

def val(r):
  if r[0]=='v':
    return ev(r[1])
  if r[0]=='n':
    return r[1]
  return None

def same(r1,r2):
  if r1[0]=='i' and r2[0]=='i':
    return r1[1]==r2[1]
  v1=val(r1)
  v2=val(r2)
  if v1 is None or v2 is None:
    return False
  return abs(v1-v2)<1e-7*(1+abs(v1))

def show(r):
  print("= "+txt(r))
  if r[0]=='v':
    print("~ "+str(ev(r[1])))
  if r[0]=='n':
    g=ident(r[1])
    if g:
      print("maybe = "+g)

def main():
  print("LIMIT finder")
  print("point: 2, pi/2, inf, -inf")
  print("add + or - for one side")
  while True:
    s=input("f(x)=")
    if s=="":
      break
    f=readf(s)
    if f is None:
      continue
    if hassym(f):
      print("Use x only here")
      print("(2 variables: SADDLE)")
      continue
    p=input("x -> ").replace(" ","").lower()
    if p=="":
      print("Enter a point: 2, pi/2, inf")
      continue
    try:
      sd=0
      if len(p)>1 and p[-1] in '+-':
        sd=1 if p[-1]=='+' else -1
        p=p[:-1]
      inf=False
      a=ZERO
      if p in ('inf','oo','+inf','+oo'):
        inf=True
        sd=1
      elif p in ('-inf','-oo'):
        inf=True
        sd=1
        f=S(subst(f,Ng(X)))
        print("x->-inf: put x=-t, t->inf")
        print("(steps shown with x as t)")
      else:
        a=S(parse(p))
        ev(a)
    except MemoryError:
      print("Out of memory")
      continue
    except Exception:
      print("Can't read the point")
      continue
    if gc:
      gc.collect()
    if sd!=0:
      r,lg=side(f,a,sd,inf)
      for l in lg:
        print(l)
      show(r)
      continue
    r1,l1=side(f,a,1,False)
    r2,l2=side(f,a,-1,False)
    if r1[0]=='u' and r2[0]=='u':
      print("f undefined near this point")
      continue
    if r1[0]=='u' or r2[0]=='u':
      if r1[0]=='u':
        r,lg,w="right","left",l2
        r1=r2
      else:
        r,lg,w="left","right",l1
      print("f undefined on "+r+" side,")
      print("using "+lg+"-hand limit")
      for l in w:
        print(l)
      show(r1)
      continue
    if same(r1,r2) or (r1[0] in 'bd' and r2[0] in 'bd'):
      for l in l1:
        print(l)
      if l2!=l1 and l2:
        print("(left side:)")
        for l in l2:
          print(l)
      show(r1)
      continue
    if l1:
      print("right side:")
      for l in l1:
        print(l)
    if l2:
      print("left side:")
      for l in l2:
        print(l)
    print("x->a+: "+txt(r1))
    print("x->a-: "+txt(r2))
    print("= DNE (one-sided limits differ)")

main()
