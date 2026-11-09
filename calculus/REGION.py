try:
  import gc
except ImportError:
  gc=None
import sys

def unload(me):
  for k in list(sys.modules.keys()):
    if k!=me and (k[:3] in ('INT','REG') or k in ('CORE','NLIB')):
      del sys.modules[k]
  if gc:
    gc.collect()

def main():
  print("REGION: double integrals")
  print("1 polygon (vertices)")
  print("2 polygon from lines")
  print("3 between two curves")
  print("4 f(ax+b)f(cy+d), known int")
  while True:
    m=input("mode (blank=quit)=")
    if m=="":
      break
    if m not in ('1','2','3','4'):
      continue
    mod=None
    try:
      if m=='4':
        mod=__import__('REGF')
        mod.run()
      elif m=='3':
        mod=__import__('REGB')
        mod.run()
      else:
        mod=__import__('REGP')
        mod.run(m)
    except MemoryError:
      st=''
      try:
        st=mod.STG[0]
      except Exception:
        pass
      print("Out of memory"+(" while "+st if st else ""))
    except Exception as ex:
      print("Error: "+repr(ex))
    mod=None
    unload('REGION')

main()
