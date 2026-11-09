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
  print("CHGVAR: change of variables")
  print("u=u(x,y), v=v(x,y)")
  while True:
    s=input("f(x,y) (blank=quit)=")
    if s=="":
      break
    mod=None
    try:
      mod=__import__('CHGM')
      mod.main(s)
    except MemoryError:
      print("Out of memory")
    except Exception as ex:
      print("Error: "+repr(ex))
    mod=None
    unload('CHGVAR')

main()
