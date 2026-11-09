# coolCAS

**A mini computer algebra system for the TI-84 Evo, written in Python.**

The TI-84 Evo has no CAS: integrals, limits and solvers give you decimals. coolCAS adds exact answers, with the working shown, in about 120 KB of Python that runs in the Evo's built-in Python app.

```
f(x)=1/(x^3+3)
F(x)=3^(1/3)*(ln|x+3^(1/3)|-ln(x^2-3^(1/3)*x+3^(2/3))/2
     +sqrt(3)*atan(3^(1/6)*(2x-3^(1/3))/3))/9+C
lower=1  upper=2
= ...   ~ 0.1617446618
```

---

## Programs

| Program | What it does |
|---|---|
| `EXACT` | Turns a decimal into an exact value: `0.8660254038` → `sqrt(3)/2` |
| `INTEGRAL` | Single, double and triple integrals, showing every antiderivative. Reverses the order of integration (dy dx ↔ dx dy) when it gets stuck |
| `LIMIT` | Limits with steps: L'Hôpital, factor & cancel, conjugates, `e^(g ln f)` for 1^∞ / 0⁰ / ∞⁰, the squeeze theorem, one-sided limits, DNE |
| `SADDLE` | Critical points of f(x,y), classified as saddle / max / min, including degenerate cases where D = 0 |
| `OPTIM` | Lagrange multipliers: prints the system, finds every candidate, gives max/min on a curve g = c or a closed region g ≤ c |
| `REGION` | Double integrals over polygons (from vertices or from lines), regions between two curves (finds the intersections), and separable `f(ax+b)f(cy+d)` with a known ∫f |
| `CHGVAR` | Change of variables u(x,y), v(x,y): Jacobian, inverse map, g(u,v), the new iterated integral and its value |
| `COORDS` | Polar, cylindrical and spherical integrals. Also converts surfaces to bounds (`z=2` → `ρ = 2/cos φ`) |
| `VECTOR` | Line ∩ plane, tangent plane parallel to a plane, directional derivative, rate along a path r(t), chain rule (even when only the partials are known), tangent plane at a point |
| `POTENT` | Finds f(x,y) from f_x and f_y (checks they're consistent first), with an optional condition f(x₀,y₀) = c |
| `MEMTEST` | Diagnostic: prints free memory as each library loads |

## Examples

**Reversing the order of integration**
```
f=1/sqrt(y^4+1)    y: x^(1/3) to 1    x: 0 to 1
dy: no elementary F
Reversed order dxdy:
 x: 0 to y^3
 y: 0 to 1
...
= sqrt(2)/2-1/2
```

**L'Hôpital, twice**
```
f(x)=((x-2)(ln(x)-ln(2)))/(x^2-5x+6+2ln(x)-ln(4))   x -> 2
0/0 -> L'Hopital: ...
0/0 -> L'Hopital: (2/x^2+1/x)/(2-2/x^2)
= 2/3
```

**Change of variables with no clean inverse**
```
f=(x+2y)(x^2+y^2)   u=x^2+y^2   v=y-2x
J=d(u,v)/d(x,y)=2x+4y
f/|J| = x^2/2+y^2/2
g(u,v) = u/2
= 175/4
```

**Surface to spherical bounds**
```
equation: z=sqrt(x^2+y^2)   (r cancels) p = pi/4
equation: z=2               r = 2/cos(p)
```

## Installation

1. Download the latest release and unzip it.
2. Send **all `.py` files** to the calculator with TI Connect CE.
3. On the calculator, open the **Python** app → **FILES** → pick a program → run.

Only the programs in the table above are meant to be run. The other files are libraries that the programs load themselves.

<details>
<summary>Which files each program needs</summary>

| Program | Files |
|---|---|
| EXACT | `EXACT` |
| INTEGRAL | `INTA`–`INTK`, `INTM`, `INTR` |
| REGION | `INTA`–`INTK`, `INTM`, `CORE`, `REGP`, `REGB`, `REGF` |
| CHGVAR | `INTA`–`INTK`, `INTM`, `CHGM`, `CHGF` |
| COORDS | `INTA`–`INTK`, `INTM`, `COCV`, `COIN`, `COSF` |
| POTENT | `INTA`–`INTK`, `INTM` |
| LIMIT, SADDLE, VECTOR | `INTA`, `INTB`, `INTC`, `CORE` |
| OPTIM | `INTA`, `INTB`, `INTC`, `CORE`, `NLIB` |

</details>

## Typing tips

- `^` or `**` for powers. Implicit multiplication works: `2x`, `3xy`, `2sin(x)`
- `ln(x)` needs brackets. `log(x)` also means ln
- Functions: `sqrt cbrt sin cos tan asin acos atan sec csc cot abs exp`, constants `e` and `pi`
- Fractional powers need brackets: `x^(1/3)`, not `x^1/3`
- LIMIT points: `2`, `pi/2`, `inf`, `-inf`. Add `+` or `-` for one side: `0+`
- COORDS variables: `r` (r or ρ), `t` (θ), `p` (φ). Translate limits before typing: y → `r*sin(t)`
- CHGVAR limits are typed in `u` and `v`. VECTOR's chain rule uses `u`, `v`, `w`

## How it works

coolCAS is a small symbolic engine built from scratch for MicroPython:

- **Exact arithmetic** with rationals, plus a simplifier that handles surds, π, e and logs
- **Symbolic differentiation**
- **Integration** using a table, linear and u-substitution, integration by parts (including tabular), trig powers and products, Weierstrass t = tan(x/2), root and exponential substitutions, and **full partial fractions**: repeated factors, irreducible quadratics, x³ + c, x⁴ + c
- **Every result is checked.** Antiderivatives are verified by differentiating them back, and definite values are compared with a tanh-sinh numerical integral. If an exact method fails, coolCAS says so and falls back to a numeric answer with a "maybe =" guess. It never presents a guess as exact.

### Memory

The Evo's Python heap is about **124 KB**, and the integration engine alone is about 47 KB of code. To fit:

- Programs with several modes (`REGION`, `CHGVAR`, `COORDS`) are small menus that **load each mode on demand and unload it afterwards**
- Lighter programs (`LIMIT`, `SADDLE`, `OPTIM`, `VECTOR`) use `CORE`, a slim parser and simplifier, instead of the full engine
- File names stay within the calculator's 8-character limit

If you see `Out of memory`, run `MEMTEST` and include its output in your bug report.

## Limitations

coolCAS is a *mini* CAS:

- Some elementary integrals aren't covered, for example x·√(quadratic). Those fall back to a numeric answer.
- Non-elementary integrals such as e^(−x²) and sin(x)/x are given numerically.
- `SADDLE` and `OPTIM` search a box (|x|, |y| ≤ 5 by default, adjustable). Points outside the box won't be found.
- `REGION` polygons must be convex.
- Results have been tested on desktop Python against a reference CAS. On-calculator testing is ongoing.

## Version history

| Version | Changes |
|---|---|
| [1.4.1](https://github.com/CheesyWannabe/coolCAS/releases/tag/v1.4.1) | Fixed [#1](https://github.com/CheesyWannabe/coolCAS/issues/1): results now show clear brackets around fractions |
| 1.4 | Added `VECTOR`. `OPTIM` shows exact max/min values even when coordinates are messy |
| 1.3.1 | `REGION`'s intersection search uses far less memory |
| 1.3 | `REGION` / `CHGVAR` / `COORDS` split into on-demand parts to fit in memory |
| 1.2 | Added `OPTIM`, `REGION`, `CHGVAR`, `COORDS` |
| 1.1 | Added `POTENT` |
| 1.0 | `EXACT`, `INTEGRAL` (single / double / triple), `LIMIT`, `SADDLE` |

## Contributing

Bug reports are welcome. Please include the exact input, what you expected, and a photo of the screen. For crashes or memory errors, also include the output of `MEMTEST`.
