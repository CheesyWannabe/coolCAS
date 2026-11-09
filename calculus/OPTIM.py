from NLIB import *
if gc:
  gc.collect()

def lam(fx,fy,gx,gy,p):
  x,y=p[0],p[1]
  for a,b in ((fx,gx),(fy,gy)):
    try:
      if abs(at(b,x,y))>1e-9:
        return val(S(Mu(a,Iv(b))),p)
    except Exception:
      pass
  return None,0.0

def run(m,f,g,c,R):
  fx=tidy(D(f))
  fy=tidy(dy(f))
  gx=tidy(D(g))
  gy=tidy(dy(g))
  print("Lagrange system (L=lambda):")
  print(" "+P(fx)+" = L*("+P(gx)+")")
  print(" "+P(fy)+" = L*("+P(gy)+")")
  print(" "+P(g)+" = "+P(c))
  F=(tidy(S(Sb(Mu(fx,gy),Mu(fy,gx)))),S(Sb(g,c)))
  print("eliminate L: fx*gy-fy*gx=0")
  print("searching...")
  cand=[]
  for p in solve2(F,R):
    Lw,Lv=lam(fx,fy,gx,gy,p)
    w,v=val(f,p)
    print(pstr(p)+" L="+vstr(Lw,Lv))
    print("  f="+vstr(w,v))
    cand.append((v,w,p))
  if m=='2':
    cv=ev(c)
    print("interior critical points:")
    k=0
    for p in solve2((fx,fy),R):
      try:
        if not at(g,p[0],p[1])<cv-1e-9:
          continue
      except Exception:
        continue
      w,v=val(f,p)
      print(pstr(p)+": f="+vstr(w,v))
      cand.append((v,w,p))
      k+=1
    if k==0:
      print(" none inside")
  if not cand:
    print("No candidates found")
    return
  if m=='1':
    print("(curve unbounded? then these")
    print(" are only candidates)")
  for nm,big in (("MAX",True),("MIN",False)):
    b=cand[0][0]
    for q in cand:
      if (big and q[0]>b) or (not big and q[0]<b):
        b=q[0]
    tie=[q for q in cand if abs(q[0]-b)<1e-9*(1+abs(b))]
    w=tie[0][1]
    for q in tie:
      if q[1] is not None:
        w=q[1]
    print(nm+" f = "+vstr(w,b))
    print("  at "+", ".join([pstr(q[2]) for q in tie]))

def main():
  print("OPTIM: Lagrange multipliers")
  print("1 max/min on curve g=c")
  print("2 on closed region g<=c")
  while True:
    m=input("mode (blank=quit)=")
    if m=="":
      break
    if m not in ('1','2'):
      continue
    f=readf(input("f(x,y)="))
    if f is None:
      continue
    g=readf(input("g(x,y)="))
    if g is None:
      continue
    c=readf(input("c="))
    if c is None:
      continue
    r=input("search |x|,|y|<= (blank=5)=")
    try:
      R=5.0
      if r!="":
        R=abs(ev(S(parse(r))))
    except Exception:
      print("Bad range")
      continue
    if gc:
      gc.collect()
    try:
      run(m,f,g,c,R)
    except MemoryError:
      print("Out of memory")
    except Exception as ex:
      print("Error: "+repr(ex))

main()
