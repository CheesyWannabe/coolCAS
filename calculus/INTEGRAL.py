from INTR import *
if gc:
  gc.collect()

def single():
  s=input("f(x)=")
  f=readf(s)
  if f is None:
    return
  try:
    Fx=I(f)
    if Fx is not None and not chk(f,Fx):
      Fx=None
  except MemoryError:
    print("Out of memory")
    Fx=None
  except Exception:
    Fx=None
  if hassym(f):
    print("Use x only here")
    print("(modes 2/3 for y,z)")
    return
  if Fx is not None:
    print("F(x)="+P(lncomb(Fx))+"+C")
  else:
    print("No elementary F(x) found")
  sa=input("lower (blank=skip)=")
  if sa=="":
    return
  sb=input("upper=")
  try:
    a=S(parse(sa))
    b=S(parse(sb))
    ex,nq,st=definite(f,Fx,a,b)
  except Exception:
    print("Bad limits")
    return
  if ex is not None:
    print("= "+P(lncomb(ex)))
    print("~ "+str(ev(ex)))
  elif nq is None or st==0:
    print("Diverges or undefined")
    print("on this interval")
  else:
    if Fx is not None:
      print("(F jumps on this range,")
      print(" using numeric value)")
    print("~ "+str(nq))
    g=None
    if st==2:
      g=ident(nq)
    if g:
      print("maybe = "+g)

def multi(n):
  df="dydx"
  if n==3:
    df="dzdydx"
  s=input("f=")
  o=input("order (blank="+df+")=")
  if o=="":
    o=df
  od=[c for c in o.lower() if c in 'xyz']
  if len(od)!=n or len(set(od))!=n:
    print("Order must use "+str(n)+" of x,y,z")
    return
  f=readf(s)
  if f is None:
    return
  try:
    lm=[]
    for v in od:
      lo=S(parse(input(v+" from=")))
      hi=S(parse(input(v+" to=")))
      lm.append((lo,hi))
  except Exception:
    print("Can't read a limit")
    return
  for i in range(n):
    ok=od[i+1:]
    for v in vs(lm[i][0],[])+vs(lm[i][1],[]):
      if v in 'xyz' and v not in ok:
        print("Limits of "+od[i]+" can't use "+v)
        return
  Gp,ip=chainp(f,od,lm)
  G=None
  if ip==n and numc(Gp):
    G=Gp
  if G is None and n==2:
    sr=swapreg(od,lm)
    if sr and same(f,od,lm,sr[0],sr[1]):
      od2,lm2=sr
      print("Reversed order d"+od2[0]+"d"+od2[1]+":")
      for i in range(2):
        print(" "+od2[i]+": "+P(lm2[i][0])+" to "+P(lm2[i][1]))
      G=chain(f,od2,lm2)
  if G is not None:
    print("= "+P(lncomb(G)))
    print("~ "+str(ev(G)))
    return
  print("No exact form found")
  print("numeric, please wait...")
  val=numfall(f,od,lm,Gp,ip)
  if val is None:
    print("Can't evaluate")
    return
  print("~ "+str(val))
  g=ident(val)
  if g:
    print("maybe = "+g)

def main():
  print("INTEGRAL mini-CAS")
  print("1 single 2 double 3 triple")
  print("use x y z ^ sqrt pi e ln")
  while True:
    m=input("mode (blank=quit)=")
    if m=="":
      break
    if m=="1":
      single()
    elif m=="2" or m=="3":
      multi(int(m))

main()
