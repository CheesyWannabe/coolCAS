try:
  import gc
except ImportError:
  gc=None
import sys

def unload(me):
  for k in list(sys.modules.keys()):
    if k!=me and (k[:3] in ('INT','CHG','COC','COI','COS') or k in ('CORE','NLIB')):
      del sys.modules[k]
  if gc:
    gc.collect()

def main():
  print("COORDS: polar/cyl/spherical")
  print("1 polar  2 cylindrical")
  print("3 spherical")
  print("names: r (rho), t=theta, p=phi")
  while True:
    s=input("system (blank=quit)=")
    if s=="":
      break
    if s not in ('1','2','3'):
      continue
    m=input("1 integral  2 convert eq =")
    mod=None
    try:
      mod=__import__('COSF' if m=='2' else 'COIN')
      mod.run(int(s))
    except MemoryError:
      print("Out of memory")
    except Exception as ex:
      print("Error: "+repr(ex))
    mod=None
    unload('COORDS')

main()
