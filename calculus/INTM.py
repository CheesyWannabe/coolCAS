from INTK import *
if gc:
  gc.collect()

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

def chainp(f,od,lm):
  G=f
  for i in range(len(od)):
    r=step(G,od[i],lm[i][0],lm[i][1])
    if r is None:
      return G,i
    G=r
    print(" "+mv(od[i])+"="+mp(lm[i][0])+".."+mp(lm[i][1])+": "+mp(lncomb(G)))
  return G,len(od)

def chain(f,od,lm):
  G,i=chainp(f,od,lm)
  if i==len(od) and numc(G):
    return G
  return None

def numfall(f,od,lm,G,i):
  if i>0:
    print("numeric on remaining "+str(len(od)-i)+"D")
    return numint(G,od[i:],lm[i:],min(5,7-len(od)+i))
  return numint(f,od,lm,6-len(od))

def numint(f,od,lm,m):
  sv={}
  for k in KV:
    sv[k]=KV[k]
  try:
    val=nint(f,od,lm,len(od)-1,nodes(m))
  except Exception:
    val=None
  for k in sv:
    KV[k]=sv[k]
  return val


def report(f,od,lm):
  G,i=chainp(f,od,lm)
  if i==len(od) and numc(G):
    print("= "+P(lncomb(G)))
    print("~ "+str(ev(G)))
    return G,ev(G)
  print("No exact form found")
  print("please wait...")
  v=numfall(f,od,lm,G,i)
  if v is None:
    print("Can't evaluate")
    return None,None
  print("~ "+str(v))
  g=ident(v)
  if g:
    print("maybe = "+g)
  return None,v
