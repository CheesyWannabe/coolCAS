import gc
import sys
def fr():
  gc.collect()
  try:
    return str(gc.mem_free())
  except Exception:
    return "?"
print("start free: "+fr())
for m in ("INTA","INTB","INTC","INTD","INTE","INTF","INTG","INTH","INTI","INTJ","INTK","INTM"):
  try:
    __import__(m)
    print(m+" ok, free "+fr())
  except Exception as e:
    print(m+" FAILED: "+repr(e))
    break
for m in ("INTR","REGP","REGB","CHGM","CHGF","COIN","COSF"):
  try:
    __import__(m)
    print("+"+m+" ok, free "+fr())
  except Exception as e:
    print("+"+m+" FAILED: "+repr(e))
  for k in ("INTR","REGP","REGB","CHGM","CHGF","COCV","COIN","COSF"):
    if k in sys.modules:
      del sys.modules[k]
