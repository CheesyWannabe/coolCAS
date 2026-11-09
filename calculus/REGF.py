from CORE import *
if gc:
  gc.collect()

def run():
  print("known: int_A^B f(t)dt = K")
  A_=S(parse(input("A=")))
  B_=S(parse(input("B=")))
  K=S(parse(input("K=")))
  out=ONE
  for v in ('x','y'):
    print(v+"-factor f(a*"+v+"+b), blank=none")
    s=input(" a=")
    b=ZERO
    if s!="":
      b=S(parse(input(" b=")))
    lo=S(parse(input(" "+v+" from=")))
    hi=S(parse(input(" "+v+" to=")))
    if s=="":
      w=S(Sb(hi,lo))
      print(" int d"+v+" = "+P(w))
      out=S(Mu(out,w))
      continue
    a=S(parse(s))
    t1=S(A(Mu(a,lo),b))
    t2=S(A(Mu(a,hi),b))
    print(" t=a"+v+"+b: t from "+P(t1)+" to "+P(t2))
    if abs(ev(t1)-ev(A_))<1e-12 and abs(ev(t2)-ev(B_))<1e-12:
      w=S(Mu(K,Iv(a)))
    elif abs(ev(t1)-ev(B_))<1e-12 and abs(ev(t2)-ev(A_))<1e-12:
      w=S(Ng(Mu(K,Iv(a))))
    else:
      print(" need int of f over that range")
      return
    print(" = (1/"+P(a)+")*int f dt = "+P(w))
    out=S(Mu(out,w))
  print("Fubini: product")
  print("= "+P(out))
