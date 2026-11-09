from COCV import *

def integ():
  f=readf(input("f(x,y,z)="))
  if f is None:
    return
  g0=conv(f)
  print("f = "+mp(g0))
  J=jac()
  print("dV = "+mp(S(J))+" "+("dr dt" if SYS[0]==1 else "dr dt dz" if SYS[0]==2 else "dr dp dt"))
  g=tidy(simp(S(Mu(g0,J))))
  print("integrand: "+mp(g))
  dfo=["drdt","dzdrdt","drdpdt"][SYS[0]-1]
  o=input("order (blank="+dfo+")=").replace(' ','')
  if o=="":
    o=dfo
  aliases()
  try:
    od=[]
    for c in o:
      if c in ALIAS:
        od.append(ALIAS[c])
      elif c=='z' and SYS[0]==2:
        od.append('z')
    n=2 if SYS[0]==1 else 3
    if len(od)!=n or len(set(od))!=n:
      print("bad order")
      return
    lm=[]
    for v in od:
      lo=S(parse(input(mv(v)+" from=")))
      hi=S(parse(input(mv(v)+" to=")))
      lm.append((lo,hi))
  finally:
    ALIAS.clear()
  ws=" ".join(["int["+mv(od[i])+"="+mp(lm[i][0])+".."+mp(lm[i][1])+"]" for i in range(n-1,-1,-1)])
  print(ws)
  print(" ("+mp(g)+") "+" ".join(["d"+mv(v) for v in od]))
  report(g,od,lm)

def run(s):
  setsys(s)
  try:
    integ()
  finally:
    MAP.clear()
    ALIAS.clear()
