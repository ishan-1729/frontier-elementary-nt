"""06_units.py: the units of R, by descent on T2 = sigma1^2 + 2Q (hypothesis_H.md, Section 4).

Notation of 04_ring.py, 1.7: l1 is the only real root of f; sigma1 and Q are multiplicative, sigma1(u) Q(u) = N(u),
Q >= 0, and T2 = sigma1^2 + 2Q is a positive definite quadratic form. For a nonzero rational triple y put
x = sigma1(y)^2 > 0 and q = Q(y) > 0, so x q^2 = N(y)^2 and T2(y) = x + 2q. All arithmetic is exact; real numbers
enter only through a rational interval that contains l1.

Checks:
  4.1 a rational interval of width < 2^-200 containing l1;
  4.2 e := sigma1(U) lies in (0, 1); E := Q(U) = 1/e > 1 (sigma1(U) Q(U) = N(U) = 1); sigma1(U^-1) = E, Q(U^-1) = e;
  4.3 descent lemma. With rho_lo = 2/(E(E+1)) and rho_hi = 2E^2/(E+1):
         T2(y U) - T2(y) = (E - 1)(2q - x(E+1)/E^2),    T2(y U^-1) - T2(y) = (E - 1)(x(E+1) - 2q/E),
      so T2(y U) < T2(y) iff x/q > rho_hi, and T2(y U^-1) < T2(y) iff x/q < rho_lo. The rational bound E < E_UP and
      the rounded constants RH_UP >= rho_hi, RL_DN <= rho_lo of the configuration are the ones hypothesis_H.md prints.
      If neither neighbour is smaller and y is a unit (x q^2 = 1, so x/q = q^-3), then q lies in [A, B] (rational
      bounds below) and T2(y) = h(q) := q^-2 + 2q <= max(h(A), h(B)) =: C, by two endpoint identities for h. Scaled
      form: for |N(y)| = n, T2(y) <= C n^(2/3) under the same condition;
  4.4 the units of small height: every triple with T2 <= C is listed (an exact enumeration of a rational lower form
      of T2, as in 05_cover.py 2.5.3), and the only units on the list are 1 and -1;
  4.5 conclusion: the units are exactly +-U^k (k unique); every unit is sigma U^l v^2 with sigma = +-1, l in {0,1},
      v a unit; the four classes +-1, +-U are distinct modulo squares of units.
Shared with the copy of this script in the other card that proves (H); only the configuration block differs.
"""
import time
from fractions import Fraction as Fr

import sympy as sp

import hring
from hring import mul, neg, norm, power, ONE, F_COEFFS, BASIS_POW, f
from hlib import check, bisect_root, real_forms, iv_eval, Iv, ldl3, coset_enumerate

# ---- configuration (the only card-specific part) ------------------------------------------------------------------
ROOT_BRACKET = (-2, -1)
UNIT = hring.EPS
UNIT_INV = hring.EPS_INV
UNIT_NAME = 'EPS'
E_UP = Fr(8137981, 10**6)            # E < E_UP (hypothesis_H.md 4.2)
RH_UP = Fr(1449483, 10**5)           # rho_hi < RH_UP (hypothesis_H.md 4.3)
RL_DN = Fr(2689446, 10**8)           # rho_lo > RL_DN (hypothesis_H.md 4.3)
C_TEXT = Fr(676515, 10**5)           # C <= C_TEXT, the rounded C of hypothesis_H.md (4.3)
# --------------------------------------------------------------------------------------------------------------------

T0 = time.time()
print("4.1 root interval")
lo, hi = bisect_root(f, ROOT_BRACKET[0], ROOT_BRACKET[1], 201)
check(f"l1 in [lo, hi] with f(lo) <= 0 <= f(hi) (exact) and width < 2^-200; l1 = {float(lo):.15f} (the only real"
      " root, 04_ring.py 1.7)", f(lo)*f(hi) <= 0 and hi - lo < Fr(1, 2**200))
LI = Iv(lo, hi)
L = sp.Symbol('L')
sigma1, Q, s_, p_ = real_forms(F_COEFFS, BASIS_POW, L)
fL = sp.Poly(f(L), L)


def red(expr):
    return sp.rem(sp.Poly(sp.expand(expr), L), fL).as_expr()


def sig_iv(u):
    return iv_eval(sigma1(u), L, LI)


def Q_iv(u):
    return iv_eval(red(Q(u)), L, LI)


def T2_iv(u):
    return iv_eval(red(sigma1(u)**2 + 2*Q(u)), L, LI)


print(f"4.2 the real values of {UNIT_NAME} = {UNIT}")
e_iv = sig_iv(UNIT)
E_iv = Q_iv(UNIT)
check(f"N({UNIT_NAME}) = 1 and {UNIT_NAME} * {UNIT_NAME}^-1 = 1 with {UNIT_NAME}^-1 = {UNIT_INV}",
      norm(UNIT) == 1 and mul(UNIT, UNIT_INV) == ONE)
check(f"e := sigma1({UNIT_NAME}) lies in {e_iv}, so 0 < e < 1", 0 < e_iv.lo and e_iv.hi < 1)
check(f"E := Q({UNIT_NAME}) lies in {E_iv}; sigma1 Q = N = 1 (identity modulo f(L)) and the product interval contains 1,"
      " so E = 1/e > 1", red(sigma1(UNIT)*Q(UNIT) - 1) == 0 and (e_iv*E_iv).lo <= 1 <= (e_iv*E_iv).hi and E_iv.lo > 1)
check(f"sigma1({UNIT_NAME}^-1) = E and Q({UNIT_NAME}^-1) = e (identities modulo f(L): sigma1(U) sigma1(U^-1) = 1 and"
      " Q(U) Q(U^-1) = 1)", red(sigma1(UNIT)*sigma1(UNIT_INV) - 1) == 0 and red(Q(UNIT)*Q(UNIT_INV) - 1) == 0)
E = E_iv

print("4.3 descent lemma")
x, q, Es, m, w = sp.symbols('x q E m w', positive=True)
T2_y = x + 2*q
T2_yU = x/Es**2 + 2*q*Es
T2_yUi = x*Es**2 + 2*q/Es
check("T2(y U) - T2(y) = (E - 1)(2q - x(E+1)/E^2) and T2(y U^-1) - T2(y) = (E - 1)(x(E+1) - 2q/E)",
      sp.simplify(T2_yU - T2_y - (Es - 1)*(2*q - x*(Es + 1)/Es**2)) == 0
      and sp.simplify(T2_yUi - T2_y - (Es - 1)*(x*(Es + 1) - 2*q/Es)) == 0)
rho_hi_s, rho_lo_s = 2*Es**2/(Es + 1), 2/(Es*(Es + 1))
check("as E > 1, x > 0, q > 0: T2(y U) < T2(y) iff x/q > rho_hi = 2E^2/(E+1), and T2(y U^-1) < T2(y) iff"
      " x/q < rho_lo = 2/(E(E+1))",
      sp.simplify(2*q - x*(Es + 1)/Es**2 - q*(Es + 1)/Es**2*(rho_hi_s - x/q)) == 0
      and sp.simplify(x*(Es + 1) - 2*q/Es - q*(Es + 1)*(x/q - rho_lo_s)) == 0)
check(f"rho_hi increases and rho_lo decreases with E: d/dE 2E^2/(E+1) = 2E(E+2)/(E+1)^2 > 0 and"
      " d/dE 2/(E(E+1)) = -2(2E+1)/(E^2(E+1)^2) < 0",
      sp.simplify(sp.diff(rho_hi_s, Es) - 2*Es*(Es + 2)/(Es + 1)**2) == 0
      and sp.simplify(sp.diff(rho_lo_s, Es) + 2*(2*Es + 1)/(Es**2*(Es + 1)**2)) == 0)
RH = 2*E_UP**2/(E_UP + 1)
RL = 2/(E_UP*(E_UP + 1))
check(f"E < E_UP = {E_UP} = {float(E_UP)}; so rho_hi < 2 E_UP^2/(E_UP + 1) = {float(RH):.9f} <= RH_UP = {float(RH_UP)}"
      f" and rho_lo > 2/(E_UP(E_UP + 1)) = {float(RL):.9e} >= RL_DN = {float(RL_DN):.7e}",
      E.hi < E_UP and RH <= RH_UP and RL >= RL_DN)
RH, RL = RH_UP, RL_DN


def cbrt_floor(v, scale):
    k = int(float(v)**(1/3)*scale) + 2
    while Fr(k, scale)**3 > v:
        k -= 1
    return Fr(k, scale)


def cbrt_ceil(v, scale):
    k = int(float(v)**(1/3)*scale) - 2
    while Fr(k, scale)**3 < v:
        k += 1
    return Fr(k, scale)


A = cbrt_floor(1/RH, 10**6)
B = cbrt_ceil(1/RL, 10**6)
check(f"A = {A} with A^3 RH_UP = {float(A**3*RH):.9f} <= 1, so q^3 >= 1/rho_hi > 1/RH_UP >= A^3 gives q >= A",
      A**3*RH <= 1 and A > 0)
check(f"B = {B} with B^3 RL_DN = {float(B**3*RL):.9f} >= 1, so q^3 <= 1/rho_lo < 1/RL_DN <= B^3 gives q <= B",
      B**3*RL >= 1 and A < B)


def hfun(zv):
    return 1/zv**2 + 2*zv


hs = 1/q**2 + 2*q
Asym, Bsym = sp.symbols('A B', positive=True)
check("for a unit (x q^2 = 1): x/q = q^-3 and T2 = h(q) = q^-2 + 2q",
      sp.simplify((x/q).subs(x, 1/q**2) - q**-3) == 0 and sp.simplify(T2_y.subs(x, 1/q**2) - hs) == 0)
check("h(q) - h(A) = (q - A)(2 - (q + A)/(q^2 A^2)) and h(q) - h(B) = (q - B)(2 - (q + B)/(q^2 B^2)); for"
      " 0 < A <= B, (q + A)/A^2 = q/A^2 + 1/A >= q/B^2 + 1/B = (q + B)/B^2",
      sp.simplify(hs - hfun(Asym) - (q - Asym)*(2 - (q + Asym)/(q**2*Asym**2))) == 0
      and sp.simplify(hs - hfun(Bsym) - (q - Bsym)*(2 - (q + Bsym)/(q**2*Bsym**2))) == 0)
print("    So for A <= q <= B: if 2 - (q + A)/(q^2 A^2) <= 0 then h(q) <= h(A); otherwise 2 - (q + B)/(q^2 B^2) > 0 too,")
print("    and h(q) <= h(B). Hence h(q) <= max(h(A), h(B)) on [A, B].")
C = max(hfun(A), hfun(B))
check(f"h(A) = {float(hfun(A)):.6f}..., h(B) = {float(hfun(B)):.6f}...; C := max(h(A), h(B)) <= C_TEXT = {float(C_TEXT)}",
      C <= C_TEXT)
check("scaled form: with n = m^3 and q = m^2 w, x = n^2/q^2 gives x/q = w^-3 and T2 = m^2 h(w); so if neither"
      " neighbour has smaller T2, T2(y) <= C |N(y)|^(2/3)",
      sp.simplify((m**6/(m**2*w)**2)/(m**2*w) - w**-3) == 0
      and sp.simplify(m**6/(m**2*w)**2 + 2*m**2*w - m**2*(1/w**2 + 2*w)) == 0)
print(f"    Lemma (units): if T2(u {UNIT_NAME}) >= T2(u) and T2(u {UNIT_NAME}^-1) >= T2(u), then T2(u) <= C <= {float(C_TEXT)}.")

print("4.4 the units of small height")
a, b, c = sp.symbols('a b c')
T2s = sp.expand(sigma1((a, b, c))**2 + 2*Q((a, b, c)))
syms = (a, b, c)
gram = [[red(sp.Rational(1, 2)*sp.diff(T2s, syms[k], syms[j])) for j in range(3)] for k in range(3)]
G_iv = [[iv_eval(gram[k][j], L, LI) for j in range(3)] for k in range(3)]
RHO = Fr(1, 2**30)
G_mid = [[Fr(round((G_iv[k][j].lo + G_iv[k][j].hi)/2*2**30), 2**30) for j in range(3)] for k in range(3)]
check("each entry of the Gram matrix of T2 lies within rho = 2^-30 of the rational matrix G'",
      all(G_mid[k][j] - RHO <= G_iv[k][j].lo and G_iv[k][j].hi <= G_mid[k][j] + RHO for k in range(3) for j in range(3)))
G_low = [[G_mid[k][j] - (3*RHO if k == j else 0) for j in range(3)] for k in range(3)]
dd = ldl3(G_low)
check(f"Glow := G' - 3 rho I has positive pivots {[f'{float(t):.6f}' for t in dd[:3]]}, and Glow(y) <= T2(y)",
      min(dd[:3]) > 0)
lst = coset_enumerate(G_low, (Fr(0), Fr(0), Fr(0)), C)
units_small = sorted(tuple(int(t) for t in u) for u in lst if norm(u) in (1, -1))
check(f"the {len(lst)} integer triples with Glow <= C contain every triple with T2 <= C; the units among them are exactly"
      f" (1,0,0) and (-1,0,0): {units_small}", units_small == [(-1, 0, 0), (1, 0, 0)])
if len(lst) <= 10:
    print("    the triples:", ", ".join(str(tuple(int(t) for t in u)) for u in sorted(lst)))

print("4.5 conclusion")
print(f"    Every unit u is +-{UNIT_NAME}^k: the set {{u {UNIT_NAME}^k : T2(u {UNIT_NAME}^k) <= T2(u)}} is finite (T2 is")
print("    positive definite) and nonempty; its element v of least T2 has T2(v U) >= T2(v) and T2(v U^-1) >= T2(v), so")
print("    T2(v) <= C (4.3) and v = +-1 (4.4).")
check(f"0 < sigma1({UNIT_NAME}) < 1, so sigma1({UNIT_NAME}^k) = sigma1({UNIT_NAME})^k = +-1 only for k = 0: k is unique",
      0 < e_iv.lo and e_iv.hi < 1)
check("signs: sigma1(-1) = -1 < 0 < sigma1(U) and sigma1(v^2) = sigma1(v)^2 > 0, so -1 and -U are not squares of units",
      e_iv.lo > 0)
sqrtE_hi = Fr(int(float(E.hi)**0.5*10**6) + 1, 10**6)
check(f"size: a unit v with v^2 = {UNIT_NAME} would have sigma1(v)^2 = e <= {float(e_iv.hi):.6f} and"
      f" Q(v) = sqrt(E) <= {sqrtE_hi}, so T2(v) <= {float(e_iv.hi + 2*sqrtE_hi):.6f} <= C; then v = +-1 by 4.4, but"
      f" (+-1)^2 = 1 != {UNIT_NAME}", sqrtE_hi**2 >= E.hi and e_iv.hi + 2*sqrtE_hi <= C and power(ONE, 2) != UNIT)
print("    So the classes +-1, +-U are distinct modulo squares of units, and every unit is sigma U^l v^2 with")
print("    sigma = +-1 (the sign of sigma1(u)), l in {0, 1} and v a unit.")
ok = True
for kexp in range(-3, 4):
    for sgn in (1, -1):
        u = power(UNIT, kexp) if sgn == 1 else neg(power(UNIT, kexp))
        steps = 0
        while True:
            t0, tp, tm = T2_iv(u), T2_iv(mul(u, UNIT)), T2_iv(mul(u, UNIT_INV))
            if tp.hi < t0.lo:
                u = mul(u, UNIT)
            elif tm.hi < t0.lo:
                u = mul(u, UNIT_INV)
            else:
                break
            steps += 1
        ok = ok and u == (sgn, 0, 0) and steps == abs(kexp)
check("illustration (not used): the greedy descent takes +-U^k to +-1 in |k| steps, |k| <= 3", ok)
print(f"    T2({UNIT_NAME}) in {T2_iv(UNIT)}, T2({UNIT_NAME}^-1) in {T2_iv(UNIT_INV)}")
print(f"all checks of 06_units.py passed ({time.time() - T0:.1f} s)")
