from math import *
ln=log

def gcd(a,b):
  a,b=abs(a),abs(b)
  while b:
    a,b=b,a%b
  return a

def rat(x,t):
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
      return (s*h1,k1)
    f=y-a
    if f<=0:
      return None
    y=1/f
  return None

def term(c,s):
  if s=="":
    return str(c)
  if c==1:
    return s
  if c==-1:
    return "-"+s
  return str(c)+s

def show(p,q,s):
  r=term(p,s)
  if q!=1:
    r=r+"/"+str(q)
  return r

def sqf(n):
  k=1
  i=2
  while i*i<=n:
    while n%(i*i)==0:
      n=n//(i*i)
      k=k*i
    i+=1
  return k,n

def rad(x,t,u):
  r=rat(x*x,2*abs(x)*t+t*t)
  if not r:
    return None
  p,q=r
  if p<=0 or p*q>10**7:
    return None
  k,m=sqf(p*q)
  if m==1:
    return None
  g=gcd(k,q)
  k,q=k//g,q//g
  if x<0:
    k=-k
  return show(k,q,u+"sqrt("+str(m)+")")

SQ=[2,3,5,6,7,10,11,13,14,15,17,19,21,22,23,26,29,30]
RT=[sqrt(n) for n in SQ]

def mix(x,t):
  for d in range(1,D+1):
    for b in range(1,D+1):
      for j in range(len(SQ)):
        for sg in (1,-1):
          r=x*d-sg*b*RT[j]
          a=int(floor(r+0.5))
          if a!=0 and abs(r-a)<=d*t:
            g=gcd(gcd(a,b),d)
            a,bb,dd=a//g,b//g,d//g
            o="+"
            if sg<0:
              o="-"
            s=str(a)+o+term(bb,"sqrt("+str(SQ[j])+")")
            if dd>1:
              s="("+s+")/"+str(dd)
            return s
  return None

def solve(x,t):
  if abs(x)<=t:
    return "0"
  r=rat(x,t)
  if r:
    return show(r[0],r[1],"")
  r=rat(x/pi,t/pi)
  if r:
    return show(r[0],r[1],"pi")
  r=rad(x,t,"")
  if r:
    return r
  r=mix(x,t)
  if r:
    return r
  r=rad(x/pi,t/pi,"pi*")
  if r:
    return r
  r=mix(x*x,2*abs(x)*t)
  if r:
    if x<0:
      return "-sqrt("+r+")"
    return "sqrt("+r+")"
  r=rat(x/(pi*pi),t/9.8)
  if r:
    return show(r[0],r[1],"pi^2")
  r=rat(x/e,t/e)
  if r:
    return show(r[0],r[1],"e")
  if abs(x)<30:
    y=exp(x)
    r=rat(y,y*t)
    if r and r!=(1,1):
      return "ln("+show(r[0],r[1],"")+")"
  return None

ep=1.0
while 1.0+ep/2>1.0:
  ep=ep/2

def tol(s,x):
  m=max(1,abs(x))
  f=max(1e-12,100*ep)*m
  try:
    float(s)
  except:
    return f
  s=s.lower()
  x10=0
  if "e" in s:
    s,k=s.split("e")
    x10=int(k)
  if "." not in s:
    return f
  d=len(s.split(".")[1])-x10
  if d<4:
    return f
  return max(10**(-d),f)

print("EXACT VALUE FINDER")
print("Type a decimal or expr")
print("eg 0.8660254038")
print("or sin(pi/3)")
print("Blank line = quit")
while True:
  s=input("x=")
  if s=="":
    break
  s=s.replace("^","**")
  try:
    x=float(eval(s))
  except:
    print("Can't read that")
    continue
  t=tol(s,x)
  Q=int(min(1000,0.1/sqrt(t/max(1,abs(x)))))
  D=max(2,min(12,Q//10))
  r=solve(x,t)
  if r:
    print("= "+r)
  else:
    print("No simple form found")
