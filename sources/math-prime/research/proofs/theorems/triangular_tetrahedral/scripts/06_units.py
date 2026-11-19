"""06_units.py: the units of R = Z[l], by descent on T2 = sigma1^2 + 2Q (hypothesis_H.md, Section 4).

Notation of 04_ring.py, 1.6: l1 is the only real root of f(t) = t^3 - 4t + 6; sigma1(u) = a + b l1 + c l1^2;
Q(u) = sigma1(adj(u)) >= 0; sigma1 and Q are multiplicative and sigma1(u) Q(u) = N(u); T2 = sigma1^2 + 2Q is a
positive definite quadratic form. For a nonzero rational triple y put x = sigma1(y)^2 > 0 and q = Q(y) > 0, so
x q^2 = N(y)^2 and T2(y) = x + 2q. All arithmetic is exact (fractions.Fraction, sympy); real numbers enter only
through the rational interval [lo, hi] that contains l1.

Checks:
  4.1 a rational interval [lo, hi] of width < 2^-200 with f(lo) < 0 < f(hi), inside (-13/5, -5/2);
  4.2 sigma1(EPS) in (1/606, 1/605); E := Q(EPS) = sigma1(adj EPS) = sigma1(EPS_INV) in (605.08, 605.09);
      sigma1(EPS) Q(EPS) = N(EPS) = 1 (identity mod f, and the product of the two intervals contains 1);
      so sigma1(EPS) = 1/E, sigma1(EPS^-1) = E, Q(EPS^-1) = 1/E;
  4.3 descent lemma. With rho_lo = 2/(E(E+1)) and rho_hi = 2E^2/(E+1) (so rho_hi/rho_lo = E^3):
         T2(y EPS)    - T2(y) = (E - 1)(2q - x(E+1)/E^2),
         T2(y EPS^-1) - T2(y) = (E - 1)(x(E+1) - 2q/E),
      so T2(y EPS) < T2(y) iff x/q > rho_hi, and T2(y EPS^-1) < T2(y) iff x/q < rho_lo.
      Rational bounds: rho_hi < RH = 1209 and rho_lo > RL = 5453/10^9 (monotonicity in E, E < 60509/100).
      If neither neighbour is smaller and y is a unit (x q^2 = 1, so x/q = q^-3), then q^3 lies in
      [1/rho_hi, 1/rho_lo], hence q = Q(y) lies in [A, B] with A = 469/5000, B = 2841/50 (A^3 RH <= 1 <= B^3 RL),
      and T2(y) = h(q) := q^-2 + 2q <= max(h(A), h(B)) = h(A) < 113.85 < C = 114
      (h(q) <= max(h(A), h(B)) on [A, B]: identities below). Scaled form: for any nonzero rational y with
      |N(y)| = n, w = Q(y)/n^(2/3) plays the role of q and T2(y) = n^(2/3) h(w), so in that case
      T2(y) <= h(A) n^(2/3) < 114 n^(2/3).
  4.4 T2 = 3(a + 8c/3)^2 + d2 (b + k c)^2 + d3 c^2 (d2 = 3l1^2 - 8, k = 4l1/d2, d3 = 716/(3 d2); 04: 1.6), so
      T2 <= 114 forces |c| <= 2, |b| <= 5, |a| <= 11; among the 23*11*5 = 1265 triples of this box, T2 <= 114 is
      decided exactly by intervals (187 triples, 0 included); the units among them are exactly (1,0,0) and (-1,0,0);
  4.5 conclusion: |sigma1(EPS)| is in (0,1), so EPS^k = +-1 only for k = 0; every unit is +-EPS^k with k unique
      (least-T2 element of the orbit {u EPS^k}); every unit is sigma EPS^l v^2 (sigma = +-1, l in {0,1},
      v a unit); the classes +-1, +-EPS are distinct modulo squares of units (sign of sigma1, and a "size"
      argument: a square root of EPS would be a unit with T2 < 114 other than +-1). Illustration (not used):
      the greedy descent takes +-EPS^k to +-1 in |k| steps for |k| <= 4.
"""
from fractions import Fraction as Fr
import itertools
import sympy as sp
from cubic import mul, neg, norm, adj, power, check, ONE, EPS, EPS_INV

C = 114


def f(x):
    return x**3 - 4*x + 6


class Iv:
    """Closed interval with rational endpoints; exact interval arithmetic (each result encloses the true value)."""

    def __init__(self, lo_, hi_=None):
        self.lo = Fr(lo_)
        self.hi = Fr(lo_ if hi_ is None else hi_)
        assert self.lo <= self.hi

    @staticmethod
    def of(x):
        return x if isinstance(x, Iv) else Iv(x)

    def __add__(self, o):
        o = Iv.of(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + (-Iv.of(o))

    def __rsub__(self, o):
        return Iv.of(o) - self

    def __mul__(self, o):
        o = Iv.of(o)
        ps = [self.lo*o.lo, self.lo*o.hi, self.hi*o.lo, self.hi*o.hi]
        return Iv(min(ps), max(ps))

    __rmul__ = __mul__

    def sq(self):
        if self.lo >= 0:
            return Iv(self.lo**2, self.hi**2)
        if self.hi <= 0:
            return Iv(self.hi**2, self.lo**2)
        return Iv(0, max(self.lo**2, self.hi**2))

    def inv(self):
        assert self.lo > 0 or self.hi < 0
        return Iv(1/self.hi, 1/self.lo)

    def __truediv__(self, o):
        return self*Iv.of(o).inv()

    def __rtruediv__(self, o):
        return Iv.of(o)*self.inv()

    def contains(self, x):
        return self.lo <= x <= self.hi

    def __str__(self):
        return f"[{float(self.lo):.12g}, {float(self.hi):.12g}]"


print("4.1 root interval")
check("f(-13/5) = -147/125 < 0 < 3/8 = f(-5/2)", f(Fr(-13, 5)) == Fr(-147, 125) and f(Fr(-5, 2)) == Fr(3, 8))
lo, hi = Fr(-13, 5), Fr(-5, 2)
while hi - lo >= Fr(1, 2**200):
    mid = (lo + hi)/2
    assert f(mid) != 0
    if f(mid) < 0:
        lo = mid
    else:
        hi = mid
check(f"l1 in [lo, hi], f(lo) < 0 < f(hi), width < 2^-200, l1 = {float(lo):.15f} (the only real root: 04: 1.6)",
      f(lo) < 0 < f(hi) and hi - lo < Fr(1, 2**200))
LI = Iv(lo, hi)


def sigma1_iv(u):
    return u[0] + u[1]*LI + u[2]*(LI*LI)


def Q_iv(u):
    return sigma1_iv(adj(u))


def T2_iv(u):
    return sigma1_iv(u).sq() + 2*Q_iv(u)


print("4.2 the real values of EPS")
L = sp.symbols('L')
FL = L**3 - 4*L + 6


def red(expr):
    return sp.expand(sp.rem(sp.expand(expr), FL, L))


def sigma1_sym(u):
    return u[0] + u[1]*L + u[2]*L**2


sE = sigma1_iv(EPS)
E = Q_iv(EPS)
check("sigma1(EPS) = l1^2 - 5 l1 - 19 lies in (1/606, 1/605), in particular 0 < sigma1(EPS) < 1",
      Fr(1, 606) < sE.lo and sE.hi < Fr(1, 605))
check("adj(EPS) = EPS_INV, so E := Q(EPS) = sigma1(EPS_INV) = 40 l1^2 - 101 l1 + 95", adj(EPS) == EPS_INV)
print(f"    sigma1(EPS) in {sE}, E = Q(EPS) in {E}")
check("E in (605.08, 605.09)", Fr(60508, 100) < E.lo and E.hi < Fr(60509, 100))
check("sigma1(EPS) * Q(EPS) = N(EPS) = 1: exact identity mod f(L), and the product interval contains 1",
      red(sigma1_sym(EPS)*sigma1_sym(EPS_INV) - 1) == 0 and (sE*E).contains(1))
check("1/sigma1(EPS) lies in the interval of E as well (consistency)", (1/sE).lo <= E.hi and E.lo <= (1/sE).hi)
check("sigma1(EPS^-1) = sigma1(EPS_INV) = E and Q(EPS^-1) = sigma1(adj EPS_INV) = sigma1(EPS) = 1/E",
      adj(EPS_INV) == EPS and (sigma1_iv(EPS_INV).lo, sigma1_iv(EPS_INV).hi) == (E.lo, E.hi)
      and (Q_iv(EPS_INV).lo, Q_iv(EPS_INV).hi) == (sE.lo, sE.hi))

print("4.3 descent lemma")
x, q, Es, n, m, w, A_, B_ = sp.symbols('x q E n m w A B', positive=True)
T2_y = x + 2*q
T2_yE = x/Es**2 + 2*q*Es          # sigma1(y EPS)^2 + 2 Q(y EPS) with sigma1(EPS) = 1/E, Q(EPS) = E
T2_yEi = x*Es**2 + 2*q/Es         # sigma1(y EPS^-1)^2 + 2 Q(y EPS^-1) with sigma1(EPS^-1) = E, Q(EPS^-1) = 1/E
check("T2(y EPS) - T2(y) = (E - 1)(2q - x(E+1)/E^2)", sp.simplify(T2_yE - T2_y - (Es - 1)*(2*q - x*(Es + 1)/Es**2)) == 0)
check("T2(y EPS^-1) - T2(y) = (E - 1)(x(E+1) - 2q/E)", sp.simplify(T2_yEi - T2_y - (Es - 1)*(x*(Es + 1) - 2*q/Es)) == 0)
rho_hi_s, rho_lo_s = 2*Es**2/(Es + 1), 2/(Es*(Es + 1))
check("as E > 1, x > 0, q > 0: T2(y EPS) < T2(y) iff x/q > rho_hi = 2E^2/(E+1), T2(y EPS^-1) < T2(y) iff x/q < rho_lo = 2/(E(E+1))",
      sp.simplify(2*q - x*(Es + 1)/Es**2 - q*(Es + 1)/Es**2*(rho_hi_s - x/q)) == 0
      and sp.simplify(x*(Es + 1) - 2*q/Es - q*(Es + 1)*(x/q - rho_lo_s)) == 0 and E.lo > 1)
check("rho_hi / rho_lo = E^3", sp.simplify(rho_hi_s/rho_lo_s - Es**3) == 0)
rho_hi = 2*E*E/(E + 1)
rho_lo = 2/(E*(E + 1))
print(f"    rho_lo in {rho_lo}, rho_hi in {rho_hi}")
Eh = Fr(60509, 100)
RH, RL = Fr(1209), Fr(5453, 10**9)
check("monotone in E > 0: 2Eh^2/(Eh+1) - 2E^2/(E+1) = 2(Eh-E)(Eh E + Eh + E)/((Eh+1)(E+1)) and "
      "2/(E(E+1)) - 2/(Eh(Eh+1)) = 2(Eh-E)(Eh+E+1)/(E Eh (E+1)(Eh+1))",
      sp.simplify(2*B_**2/(B_ + 1) - 2*Es**2/(Es + 1) - 2*(B_ - Es)*(B_*Es + B_ + Es)/((B_ + 1)*(Es + 1))) == 0
      and sp.simplify(2/(Es*(Es + 1)) - 2/(B_*(B_ + 1)) - 2*(B_ - Es)*(B_ + Es + 1)/(Es*B_*(Es + 1)*(B_ + 1))) == 0)
check(f"so with E < Eh = 60509/100: rho_hi < 2Eh^2/(Eh+1) = {float(2*Eh**2/(Eh + 1)):.6f} < RH = 1209 and "
      f"rho_lo > 2/(Eh(Eh+1)) = {float(2/(Eh*(Eh + 1))):.6e} > RL = 5453/10^9",
      2*Eh**2/(Eh + 1) < RH and 2/(Eh*(Eh + 1)) > RL and rho_hi.hi < RH and rho_lo.lo > RL)
A, B = Fr(469, 5000), Fr(2841, 50)
check(f"A = 469/5000 = 0.0938: A^3 * 1209 = {float(A**3*RH):.6f} <= 1, so q^3 >= 1/rho_hi > 1/1209 >= A^3 gives q > A",
      A**3*RH <= 1)
check(f"B = 2841/50 = 56.82: B^3 * 5453/10^9 = {float(B**3*RL):.6f} >= 1, so q^3 <= 1/rho_lo < 10^9/5453 <= B^3 gives q < B",
      B**3*RL >= 1)


def hfun(z):
    return 1/z**2 + 2*z


hA, hB = hfun(A), hfun(B)
print(f"    h(A) = {hA} = {float(hA):.6f},  h(B) = {hB} = {float(hB):.6f}")
check("h(A) < 113.85 < 114 and h(B) < 113.65 < 114", hA < Fr(11385, 100) and hB < Fr(11365, 100) and max(hA, hB) < C)
hs = 1/q**2 + 2*q
check("for a unit (x q^2 = 1): x/q = q^-3 and T2 = x + 2q = h(q) = q^-2 + 2q",
      sp.simplify((x/q).subs(x, 1/q**2) - q**-3) == 0 and sp.simplify(T2_y.subs(x, 1/q**2) - hs) == 0)
check("endpoint lemma, identities: h(q) - h(A) = (q - A)(2 - (q + A)/(q^2 A^2)), h(q) - h(B) = (q - B)(2 - (q + B)/(q^2 B^2)),"
      " and (q + A)/A^2 - (q + B)/B^2 = q(1/A^2 - 1/B^2) + (1/A - 1/B) >= 0 for 0 < A <= B",
      sp.simplify(hs - (1/A_**2 + 2*A_) - (q - A_)*(2 - (q + A_)/(q**2*A_**2))) == 0
      and sp.simplify(hs - (1/B_**2 + 2*B_) - (q - B_)*(2 - (q + B_)/(q**2*B_**2))) == 0
      and sp.simplify((q + A_)/A_**2 - (q + B_)/B_**2 - (q*(1/A_**2 - 1/B_**2) + (1/A_ - 1/B_))) == 0)
print("    (so for A <= q <= B: if 2 - (q+A)/(q^2A^2) <= 0 then h(q) <= h(A); otherwise 2 - (q+B)/(q^2B^2) > 0 and")
print("    h(q) <= h(B). Equivalently: h is convex, h''(q) = 6/q^4 > 0.)")
check("h''(q) = 6/q^4", sp.simplify(sp.diff(hs, q, 2) - 6/q**4) == 0)
print("    Lemma (units): if T2(u EPS) >= T2(u) and T2(u EPS^-1) >= T2(u), then A < Q(u) < B and T2(u) <= h(A) < 114.")
print("    So T2(u) > 114 (indeed T2(u) >= 113.85) implies T2(u EPS) < T2(u) or T2(u EPS^-1) < T2(u).")
check("scaled form: with n = m^3 (m = n^(1/3) > 0) and q = m^2 w: x = n^2/q^2 gives x/q = w^-3 and T2 = m^2 h(w)",
      sp.simplify((m**6/(m**2*w)**2)/(m**2*w) - w**-3) == 0
      and sp.simplify(m**6/(m**2*w)**2 + 2*m**2*w - m**2*(1/w**2 + 2*w)) == 0)
print("    Lemma (scaled): for a nonzero rational triple y with |N(y)| = n, if neither T2(y EPS) nor T2(y EPS^-1) is")
print("    smaller than T2(y), then w = Q(y)/n^(2/3) lies in (A, B) and T2(y) = n^(2/3) h(w) <= h(A) n^(2/3) < 114 n^(2/3).")
# the exact threshold, for information: C* = h(q*) with q*^3 = 1/rho_hi = (E+1)/(2E^2), C*^3 = 4(E^2+E+1)^3/(E^2 (E+1)^2)
Rs = (Es + 1)/(2*Es**2)
Cstar3 = 4*(E*E + E + 1)*(E*E + E + 1)*(E*E + E + 1)/(E*E*(E + 1)*(E + 1))
check("remark (not used): the exact threshold is C* = h(q*) with q*^3 = 1/rho_hi = (E+1)/(2E^2); h(q)^3 = (1 + 2q^3)^3/q^6,"
      " so C*^3 = 4(E^2+E+1)^3/(E^2(E+1)^2); q^2 (h(Eq) - h(q)) = (E-1)(2q^3 - (E+1)/E^2), so also C* = h(E q*)"
      " (the two neighbours tie there); and 113.62 < C* < 113.63",
      sp.simplify(hs**3 - (1 + 2*q**3)**3/q**6) == 0
      and sp.simplify((1 + 2*Rs)**3/Rs**2 - 4*(Es**2 + Es + 1)**3/(Es**2*(Es + 1)**2)) == 0
      and sp.simplify(q**2*(hs.subs(q, Es*q) - hs) - (Es - 1)*(2*q**3 - (Es + 1)/Es**2)) == 0
      and Fr(11362, 100)**3 < Cstar3.lo and Cstar3.hi < Fr(11363, 100)**3)

print("4.4 the units of small height")
d2 = 3*LI*LI - 8
kk = 4*LI/d2
d3 = Fr(716, 3)/d2
print(f"    d2 = 3l1^2 - 8 in {d2}, k = 4l1/d2 in {kk}, d3 = 716/(3 d2) in {d3}")
a, b, c = sp.symbols('a b c')
T2L = 3*a**2 + 16*a*c + (3*L**2 - 8)*b**2 + 8*L*b*c + (32 - 18*L - 4*L**2)*c**2
check("mod f(L): (3L^2-8) T2 = 3(3L^2-8)(a + 8c/3)^2 + ((3L^2-8) b + 4L c)^2 + (716/3) c^2 (04: 1.6), "
      "and T2 = sigma1^2 + 2Q",
      red((3*L**2 - 8)*T2L - (3*(3*L**2 - 8)*(a + sp.Rational(8, 3)*c)**2 + ((3*L**2 - 8)*b + 4*L*c)**2
                               + sp.Rational(716, 3)*c**2)) == 0
      and red(sigma1_sym((a, b, c))**2 + 2*sigma1_sym(adj((a, b, c))) - T2L) == 0)
check("d2 > 0 and d3 > 0, so each of the three terms is >= 0", d2.lo > 0 and d3.lo > 0)
kabs = max(abs(kk.lo), abs(kk.hi))
check(f"T2 <= 114 forces |c| <= 2: 9 d3 >= {float(9*d3.lo):.4f} > 114", 9*d3.lo > C)
check(f"then |b| <= 5: for |b| >= 6, |c| <= 2: d2 (b + kc)^2 >= d2 (6 - 2|k|)^2 >= {float(d2.lo*(6 - 2*kabs)**2):.4f} > 114",
      6 - 2*kabs > 0 and d2.lo*(6 - 2*kabs)**2 > C)
check(f"then |a| <= 11: for |a| >= 12, |c| <= 2: 3(a + 8c/3)^2 >= 3(12 - 16/3)^2 = 400/3 > 114",
      3*(12 - Fr(16, 3))**2 == Fr(400, 3) and Fr(400, 3) > C)
print("    (in general T2 <= M gives c^2 <= M/d3, (b + kc)^2 <= M/d2, (a + 8c/3)^2 <= M/3: finitely many triples)")
box = list(itertools.product(range(-11, 12), range(-5, 6), range(-2, 3)))
below, undecided = [], []
for u in box:
    tv = T2_iv(u)
    if tv.hi <= C:
        below.append(u)
    elif tv.lo <= C:
        undecided.append(u)
check(f"box of {len(box)} triples; T2 <= 114 decided for every one (undecided: {len(undecided)}); "
      f"{len(below)} triples have T2 <= 114", len(box) == 1265 and not undecided)
units_small = [u for u in below if norm(u) in (1, -1)]
units_box = [u for u in box if norm(u) in (1, -1)]
check(f"the units with T2 <= 114 are exactly (1,0,0) and (-1,0,0) (T2 = 3); units anywhere in the box: {units_box}",
      sorted(units_small) == [(-1, 0, 0), (1, 0, 0)])
second = sorted(set(round(float(T2_iv(u).lo), 6) for u in below))[:6]
print(f"    the six smallest values of T2 on the box (all triples, not only units): {second}")

print("4.5 conclusion")
check("0 < sigma1(EPS) < 1, so sigma1(EPS^k) = sigma1(EPS)^k is < 1 for k > 0 and > 1 for k < 0, while"
      " |sigma1(+-1)| = 1: EPS^k = +-1 only for k = 0", 0 < sE.lo and sE.hi < 1)
print("    Every unit u is +-EPS^k: the set {u EPS^k : T2(u EPS^k) <= T2(u)} is finite (4.4) and nonempty; its element")
print("    v = u EPS^j of least T2 has T2(v EPS) >= T2(v) and T2(v EPS^-1) >= T2(v), so T2(v) < 114 (4.3) and v = +-1 (4.4).")
print("    k is unique since EPS^k = +-1 only for k = 0. Writing k = 2k1 + l: u = sigma EPS^l (EPS^k1)^2, l in {0,1}.")
check("signs: sigma1(-1) = -1 < 0 < sigma1(EPS), and sigma1(v^2) = sigma1(v)^2 > 0, so -1 and -EPS are not squares"
      " and sigma = +-1 is the sign of sigma1(u)", sE.lo > 0)
check("size: a unit v with v^2 = EPS would have Q(v)^2 = E and sigma1(v)^2 = 1/E, so T2(v) = 1/E + 2 sqrt(E) < 114"
      " (1/E < 1 and 4E < 113^2); then v = +-1 by 4.4, but (+-1)^2 = 1 != EPS: EPS is not a square",
      1/E.lo < 1 and 4*E.hi < 113**2 and power((1, 0, 0), 2) != EPS)
print("    So the four classes +-1, +-EPS are distinct modulo squares of units (pairwise quotients are -1, +-EPS^-1, -EPS,")
print("    none a square), and every unit is sigma EPS^l v^2 with sigma = +-1, l in {0,1}, v a unit.")
# illustration (not used by the proof): the greedy descent of 4.3 on +-EPS^k, |k| <= 4
ok = True
for kexp in range(-4, 5):
    for sgn in (1, -1):
        u = power(EPS, kexp) if sgn == 1 else neg(power(EPS, kexp))
        steps = 0
        while True:
            t0, tp, tm = T2_iv(u), T2_iv(mul(u, EPS)), T2_iv(mul(u, EPS_INV))
            if tp.hi < t0.lo:
                u = mul(u, EPS)
            elif tm.hi < t0.lo:
                u = mul(u, EPS_INV)
            else:
                break
            steps += 1
        ok = ok and u == (sgn, 0, 0) and steps == abs(kexp) and T2_iv(u).hi <= C
check("illustration: the greedy descent (multiply by EPS or EPS^-1 while T2 decreases) takes +-EPS^k to +-1 in |k| steps,"
      " |k| <= 4", ok)
print(f"    T2(EPS) in {T2_iv(EPS)}, T2(EPS^-1) in {T2_iv(EPS_INV)}  (both > 114)")
print("all checks of 06_units.py passed")
