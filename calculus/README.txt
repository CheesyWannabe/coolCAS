TI-84 Evo Math Suite  (version 1.4.1, Nov 2026)
Python programs for the TI-84 Evo. Exact answers with working shown.

PROGRAMS (run these)
  EXACT     Decimal -> exact value (fractions, sqrt, pi, e, ln)
  INTEGRAL  Integrals: 1 single, 2 double, 3 triple. Shows each
              step; reverses order (dydx <-> dxdy) when stuck
  LIMIT     Limits with steps: L'Hopital, factor & cancel,
              conjugate, e^(g ln f), squeeze, one-sided, DNE
  SADDLE    Critical points of f(x,y): saddle / max / min
  OPTIM     Lagrange multipliers: prints the system, finds all
              candidates; max/min on a curve g=c or on a closed
              region g<=c (interior + boundary)
  REGION    Double integrals over: polygon (vertices), polygon
              from lines, region between two curves (finds the
              intersections), and f(ax+b)f(cy+d) with a known
              integral of f (Fubini + substitution)
  CHGVAR    Change of variables u(x,y), v(x,y): Jacobian, inverse,
              g(u,v), the new iterated integral, and its value
  COORDS    Polar / cylindrical / spherical: converts f and dV,
              sets up and evaluates the integral; also converts
              surfaces to bounds (z=2 -> r=2/cos(p))
  VECTOR    Gradient & geometry: line meets plane, tangent
              plane parallel to a plane, directional derivative,
              rate along a path r(t), chain rule z(u,v,w) (also
              with only partials known), tangent plane at a point
  POTENT    Find f(x,y) from fx and fy (+ f(x0,y0)=c)
  MEMTEST   Diagnostic: prints free memory while loading libs

LIBRARY FILES (must be on the calculator, do not run)
  INTA INTB INTC        core: numbers, simplifier, derivatives
  INTD .. INTK, INTM    integration engine
  INTR                  order reversal (used by INTEGRAL)
  CORE, NLIB            parser + equation solver (lighter)
  REGP REGB REGF        REGION's modes (loaded on demand)
  CHGM CHGF             CHGVAR's parts (loaded on demand)
  COCV COIN COSF        COORDS's parts (loaded on demand)

WHICH FILES EACH PROGRAM NEEDS
  EXACT              EXACT only
  INTEGRAL           INTA..INTK INTM INTR
  REGION             INTA..INTK INTM CORE REGP REGB REGF
  CHGVAR             INTA..INTK INTM CHGM CHGF
  COORDS             INTA..INTK INTM COCV COIN COSF
  POTENT             INTA..INTK INTM
  LIMIT SADDLE       INTA INTB INTC CORE
  OPTIM              INTA INTB INTC CORE NLIB
  VECTOR             INTA INTB INTC CORE
  Easiest: send all 34 .py files.

MEMORY
  REGION, CHGVAR and COORDS load only the mode you pick and
  free it afterwards. If you ever see "Out of memory", run
  MEMTEST and look at the free memory numbers it prints.

CLEAR BRACKETS
  Fractions are printed as (top)/(bottom) when either part
  has several factors, eg (cos(p)*sin(p)*r^3)/((r^2+2)^2).
  To turn this off, edit INTA: change CLEAR=[1] to CLEAR=[0].

TYPING TIPS
  ^ or ** for powers, implicit multiply ok (2x, xy)
  ln(x) needs brackets; log(x) also means ln
  sqrt(), cbrt(), sin cos tan, asin acos atan, abs(), e, pi
  LIMIT points: 2, pi/2, inf, -inf, add + or - for one side (0+)
  COORDS names: r (also rho), t = theta, p = phi
  CHGVAR limits use u and v

LIMITS OF THE TOOL
  Not a full CAS. If no exact form is found it says so and
  gives a numeric answer (with a "maybe =" guess).
  SADDLE/OPTIM search a box |x|,|y| <= 5 by default.

CHANGES
  1.4.1 clear brackets around fractions
  1.4  added VECTOR; OPTIM shows exact max/min values even
       when the point is messy; parser accepts u v w
  1.3.1 REGION intersection search uses far less memory
  1.3  REGION/CHGVAR/COORDS split into on-demand parts to
       fit in memory; order reversal moved to INTR
  1.2  added OPTIM, REGION, CHGVAR, COORDS; INTM and NLIB
       libraries; better numeric fallback when only the
       outer integral fails
  1.1  added POTENT
