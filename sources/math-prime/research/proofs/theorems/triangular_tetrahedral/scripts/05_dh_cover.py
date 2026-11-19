"""05_dh_cover.py: Dedekind-Hasse division in R = Z[l], l^3 = 4l - 6 (hypothesis_H.md, Section 2).

Lemma 2.1.
  (a) For every real triple xi in D = [-1/2,1/2] x [-1/2,1/2] x [0,1/2] (coordinates in the basis 1, l, l^2) there
      are s in {1,2} and an integer triple kappa with |N(s*xi - kappa)| < 1.
  (b) For each of the 7 classes xi in (1/2)R with xi not in R there is an integer triple kappa with
      0 < |N(xi - kappa)| < 1.
Every rational triple is brought into D by subtracting nearest integers and, when its third coordinate is negative,
by negating it (N(-x) = -N(x)). That reduction and Consequence 2.2 are written in hypothesis_H.md; 2.7 below runs them.

This script checks, in exact rational arithmetic (Fractions; sympy for the polynomial identities):
  2.1 the Taylor identity N(v+w) = N(v) + sum_i g_i(v) w_i + sum_{i<=j} q_ij(v) w_i w_j + N(w);
  2.2 the table dh_boxes.txt: each box lies in D, has dyadic ends, s in {1,2} and an integer kappa, and, with
      v = s*centre - kappa and half-widths h,
          B = |N(v)| + sum_i |g_i(v)| s h_i + sum_{i<=j} |q_ij(v)| s^2 h_i h_j + sum_m |n_m| (s h)^m  <  1
      (n_m: the coefficients of N), so |N(s*xi - kappa)| <= B < 1 on the whole box;
  2.3 the boxes cover D, by two arguments: they are exactly the leaves of a binary split tree of D (each internal
      node halved at the midpoint of one axis), replayed from the paths in the file; and each cell of the grid of
      step 1/64 that tiles D lies in exactly one box;
  2.4 part (b): the 7 classes of (1/2)R other than R, each with kappa and N(xi - kappa);
  2.5 remark, not used by (H): R is not norm-Euclidean. xi* = (79/109, 247/327, 142/327) has EPS*xi* - xi* in R and
      N(xi* - 2) = -121/109, and |N(y)| >= 121/109 for every y in xi* + R. The proof is a scaled descent on
      T2 = sigma1^2 + 2Q under EPS^(+-1), followed by an exact enumeration. Also: of the 654 classes of fixed points
      of EPS mod R, exactly 2 (those of xi* and -xi*) contain no element with |N| < 1;
  2.6 remark: Bernstein enclosures of N(s*xi - kappa) on each box. They are a second and sharper certificate, and
      they give the largest value of |N(s*xi - kappa)| on D with the table's choices;
  2.7 confirmation, not used by the proof: random Dedekind-Hasse divisions 0 < |N(s*alpha - kappa*beta)| < |N(beta)|,
      with s and kappa found by the reduction, the table and the list of 2.4.
The table was found by cover_search_K.py; the search is not part of the proof.
"""
import itertools
import math
import random
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

from cubic import mul, norm, adj, sub, check, EPS, EPS_INV, ONE

T_START = time.time()
HALF = Fr(1, 2)
D_LO = (-HALF, -HALF, Fr(0))
D_HI = (HALF, HALF, HALF)
AXES = 'xyz'


def grad(v):
    """The gradient (g_1, g_2, g_3) of N at v = (x, y, z)."""
    x, y, z = v
    return (3*x**2 + 16*x*z - 4*y**2 + 18*y*z + 16*z**2,
            -8*x*y + 18*x*z - 18*y**2 + 24*z**2,
            8*x**2 + 18*x*y + 32*x*z + 48*y*z + 108*z**2)


def quad(v):
    """The second-order Taylor coefficients q_ij (i <= j) of N at v = (x, y, z), checked in 2.1."""
    x, y, z = v
    return {(0, 0): 3*x + 8*z, (0, 1): -8*y + 18*z, (0, 2): 16*x + 18*y + 32*z,
            (1, 1): -4*x - 18*y, (1, 2): 18*x + 48*z, (2, 2): 16*x + 24*y + 108*z}


CUBIC = {(3, 0, 0): 1, (2, 0, 1): 8, (1, 2, 0): -4, (1, 1, 1): 18, (1, 0, 2): 16, (0, 3, 0): -6,
         (0, 1, 2): 24, (0, 0, 3): 36}

print("2.1 Taylor identity")
x, y, z, w1, w2, w3 = sp.symbols('x y z w1 w2 w3')
V, W = (x, y, z), (w1, w2, w3)
taylor = (norm(V) + sum(grad(V)[i]*W[i] for i in range(3))
          + sum(co*W[i]*W[j] for (i, j), co in quad(V).items()) + norm(W))
check("N(v+w) = N(v) + sum g_i(v) w_i + sum_(i<=j) q_ij(v) w_i w_j + N(w), with g and q as printed below",
      sp.expand(norm(tuple(V[i] + W[i] for i in range(3))) - taylor) == 0)
cubic_poly = sp.Poly(norm(W), w1, w2, w3)
check("the coefficients n_m of N(a,b,c) = a^3 + 8a^2c - 4ab^2 + 18abc + 16ac^2 - 6b^3 + 24bc^2 + 36c^3",
      dict(zip(cubic_poly.monoms(), cubic_poly.coeffs())) == CUBIC)
print("    g_1 =", grad(V)[0])
print("    g_2 =", grad(V)[1])
print("    g_3 =", grad(V)[2])
print("    q_ij:", ", ".join(f"q{i + 1}{j + 1} = {co}" for (i, j), co in quad(V).items()))
print("    sum of |n_m| =", sum(abs(c) for c in CUBIC.values()))


def taylor_bound(lo, hi, s, kap):
    """B for the box [lo, hi] with s and kappa (exact)."""
    c = [(lo[i] + hi[i])/2 for i in range(3)]
    H = [s*(hi[i] - lo[i])/2 for i in range(3)]
    v = tuple(s*c[i] - kap[i] for i in range(3))
    g = grad(v)
    return (abs(norm(v)) + sum(abs(g[i])*H[i] for i in range(3))
            + sum(abs(co)*H[i]*H[j] for (i, j), co in quad(v).items())
            + sum(abs(co)*H[0]**m[0]*H[1]**m[1]*H[2]**m[2] for m, co in CUBIC.items()))


def box_from_path(path):
    lo, hi = list(D_LO), list(D_HI)
    for ch in path:
        ax = AXES.index(ch.lower())
        mid = (lo[ax] + hi[ax])/2
        if ch.islower():
            hi[ax] = mid
        else:
            lo[ax] = mid
    return tuple(lo), tuple(hi)


print("2.2 the table dh_boxes.txt")
UNIT_EXP, ROWS = None, []
with open('dh_boxes.txt', encoding='ascii') as fh:
    for line in fh:
        t = line.split()
        if not t or t[0].startswith('#'):
            continue
        if t[0] == 'unit':
            UNIT_EXP = int(t[1].split('^-')[1])
            continue
        if len(t) != 11:
            check(f"line '{line.strip()}' has 11 fields", False)
        ROWS.append((t[0], tuple(int(u) for u in t[1:7]), int(t[7]), tuple(int(u) for u in t[8:11])))
check(f"the table has a unit line (2^-{UNIT_EXP}) and {len(ROWS)} boxes", UNIT_EXP is not None and len(ROWS) > 0)
UNIT = Fr(1, 2**UNIT_EXP)
BOXES = []
for n, (path, ends, s, kap) in enumerate(ROWS, 1):
    lo = tuple(e*UNIT for e in ends[:3])
    hi = tuple(e*UNIT for e in ends[3:])
    tag = f"box {n} ({path})"
    if not (path and all(ch in 'xyzXYZ' for ch in path)):
        check(f"{tag}: the path is a word in x, y, z, X, Y, Z", False)
    if box_from_path(path) != (lo, hi):
        check(f"{tag}: the listed corners are those of the path", False)
    if not all(D_LO[i] <= lo[i] < hi[i] <= D_HI[i] for i in range(3)):
        check(f"{tag}: lies in D", False)
    if not all((e.denominator & (e.denominator - 1)) == 0 for e in lo + hi):
        check(f"{tag}: dyadic ends", False)
    if s not in (1, 2):
        check(f"{tag}: s in {{1,2}}", False)
    B = taylor_bound(lo, hi, s, kap)
    if not B < 1:
        check(f"{tag}: B = {B} < 1", False)
    BOXES.append((path, lo, hi, s, kap, B))
check(f"all {len(BOXES)} boxes: the path is a word in x,y,z,X,Y,Z and gives the listed corners, the box lies in D,"
      f" its ends are dyadic (multiples of 2^-{UNIT_EXP}), s in {{1,2}}, kappa is an integer triple, and B < 1", True)
n1 = sum(1 for b in BOXES if b[3] == 1)
print(f"    boxes: {len(BOXES)}; with s = 1: {n1}; with s = 2: {len(BOXES) - n1}")
worst = max(BOXES, key=lambda b: b[5])
print(f"    largest B: {worst[5]} = {float(worst[5]):.9f}, on box {worst[0]} (s = {worst[3]}, kappa = {worst[4]})")
big = max(max(abs(t) for t in b[4]) for b in BOXES)
print(f"    largest |kappa_i|: {big}; boxes with all |kappa_i| <= 4: {sum(1 for b in BOXES if max(map(abs, b[4])) <= 4)}"
      f", with some |kappa_i| > 4: {sum(1 for b in BOXES if max(map(abs, b[4])) > 4)}")
print(f"    longest path: {max(len(b[0]) for b in BOXES)} splits; most splits along one axis (x, y, z):",
      [max(b[0].lower().count(a) for b in BOXES) for a in AXES])
freq = {}
for b in BOXES:
    freq[(b[3], b[4])] = freq.get((b[3], b[4]), 0) + 1
top = sorted(freq.items(), key=lambda kv: (-kv[1], kv[0]))
print(f"    {len(freq)} different pairs (s, kappa); the most frequent:",
      "; ".join(f"s={s} kappa={k}: {n}" for (s, k), n in top[:8]))
edges = [Fr(1, 2), Fr(9, 10), Fr(99, 100), Fr(999, 1000), Fr(1)]
counts, prev = [], Fr(0)
for e in edges:
    counts.append(f"[{prev}, {e}): {sum(1 for b in BOXES if prev <= b[5] < e)}")
    prev = e
print("    number of boxes by B:", "; ".join(counts))

print("2.3 the boxes cover D")
PATHS = [b[0] for b in BOXES]
check("no path occurs twice", len(set(PATHS)) == len(PATHS))


def replay(prefix, paths):
    """The leaves below the node `prefix`, in depth-first order, after checking that they form a split tree."""
    if paths == [prefix]:
        return [prefix]
    if not paths:
        check(f"node '{prefix}' carries at least one box", False)
    if prefix in paths:
        check(f"the box '{prefix}' has no other box inside it", False)
    axes = {p[len(prefix)].lower() for p in paths}
    if len(axes) != 1:
        check(f"node '{prefix}' is split along one axis only (found {sorted(axes)})", False)
    ax = axes.pop()
    lower = [p for p in paths if p[len(prefix)] == ax]
    upper = [p for p in paths if p[len(prefix)] == ax.upper()]
    return replay(prefix + ax, lower) + replay(prefix + ax.upper(), upper)


order = replay('', PATHS)
check(f"the {len(PATHS)} paths are exactly the leaves of a binary split tree of D: every node is either a box of the"
      " table or is halved at the midpoint of one axis into two nodes, and the file lists the leaves in"
      " depth-first order", order == PATHS)
print("    By induction from the leaves, the boxes below a node cover that node's box (a box that is halved is the")
print("    union of its two halves). So the boxes cover D, and two of them meet at most in boundary points.")
nx, nz = 2**UNIT_EXP, 2**(UNIT_EXP - 1)
cells = np.zeros((nx, nx, nz), dtype=np.int64)
off = (nx//2, nx//2, 0)
for path, ends, s, kap in ROWS:
    cells[ends[0] + off[0]:ends[3] + off[0], ends[1] + off[1]:ends[4] + off[1], ends[2]:ends[5]] += 1
check(f"second argument: each of the {nx}*{nx}*{nz} = {cells.size} cells of side 1/{nx} that tile D lies in exactly"
      " one box", bool((cells == 1).all()))
vol = sum((b[2][0] - b[1][0])*(b[2][1] - b[1][1])*(b[2][2] - b[1][2]) for b in BOXES)
check(f"the volumes of the boxes add up to vol(D) = 1/2 (sum = {vol})", vol == HALF)

print("2.4 part (b): the 7 classes xi in (1/2)R with xi not in R")
HALF_LIST = {}
for e in itertools.product((0, 1), repeat=3):
    if e == (0, 0, 0):
        continue
    xi = tuple(Fr(t, 2) for t in e)
    found = None
    for kap in itertools.product(range(-2, 3), repeat=3):
        n = norm(sub(xi, kap))
        if n != 0 and abs(n) < 1:
            key = (max(abs(t) for t in kap), sum(abs(t) for t in kap), abs(n))
            if found is None or key < found[0]:
                found = (key, kap, n)
    check(f"xi = {tuple(str(t) for t in xi)}: kappa = {found[1]}, N(xi - kappa) = N{tuple(int(2*t) for t in sub(xi, found[1]))}/8"
          f" = {found[2]}, and 0 < |N| < 1", found is not None and 0 < abs(found[2]) < 1)
    HALF_LIST[e] = found[1]

print("2.5 remark (not used by (H)): R is not norm-Euclidean")
XS = (Fr(79, 109), Fr(247, 327), Fr(142, 327))
d_plus = sub(mul(EPS, XS), XS)
d_minus = sub(mul(EPS_INV, XS), XS)
check(f"EPS*xi* - xi* = {tuple(str(t) for t in d_plus)} and EPS^-1*xi* - xi* = {tuple(str(t) for t in d_minus)} are"
      " integer triples, so EPS and EPS^-1 map the coset C = xi* + R into itself",
      all(t.denominator == 1 for t in d_plus + d_minus) and mul(EPS, EPS_INV) == ONE)
check("N(xi* - 2) = N(-139/109, 247/327, 142/327) = -121/109", norm(sub(XS, (2, 0, 0))) == Fr(-121, 109))
print("  2.5.1 the real root l1 of f(t) = t^3 - 4t + 6, and the two real forms sigma1 and Q")


def f(t):
    return t**3 - 4*t + 6


def dec(q, digits=24):
    """A rational as a decimal string with the given number of significant digits (for printing only)."""
    return str(sp.Float(sp.Rational(q.numerator, q.denominator), digits))


lo1, hi1 = Fr(-3), Fr(-2)
check("f(-3) = -9 < 0 < 6 = f(-2)", f(lo1) == -9 and f(hi1) == 6)
for _ in range(64):
    mid = (lo1 + hi1)/2
    if f(mid) <= 0:
        lo1 = mid
    else:
        hi1 = mid
check(f"l1 in [lo, hi] = [{dec(lo1)}, {dec(hi1)}]: f(lo) <= 0 < f(hi), width 2^-64",
      f(lo1) <= 0 < f(hi1) and hi1 - lo1 == Fr(1, 2**64))
L, t_ = sp.symbols('L t')
fL = sp.Poly(L**3 - 4*L + 6, L)
s0, p0 = -L, L**2 - 4


def zero_mod_f(expr):
    return sp.rem(sp.Poly(sp.expand(expr), L), fL).is_zero


def sigma1(u):
    return u[0] + u[1]*L + u[2]*L**2


def Qform(u):
    A, B = u[0] - u[2]*p0, u[1] + u[2]*s0
    return A**2 + s0*A*B + p0*B**2


a, b, c, d, e_, g = sp.symbols('a b c d e g')
U, U2 = (a, b, c), (d, e_, g)
check("f(t) = (t - L)(t^2 - s t + p) + f(L) with s = -L, p = L^2 - 4; so at L = l1 the other two roots are those"
      " of t^2 - s t + p, whose discriminant s^2 - 4p = 16 - 3 l1^2 is < 0 (l1^2 > hi^2 > 16/3)",
      sp.expand((t_ - L)*(t_**2 - s0*t_ + p0) + (L**3 - 4*L + 6) - (t_**3 - 4*t_ + 6)) == 0 and hi1**2 > Fr(16, 3))
print("    sigma1(a,b,c) = a + b l1 + c l1^2; Q(a,b,c) = A^2 + s A B + p B^2 with A = a - c p, B = b + c s")
print("    (Q is |value at a complex root|^2: a + b t + c t^2 = A + B t modulo t^2 - s t + p)")
check("sigma1(uv) = sigma1(u) sigma1(v) (identity modulo f(L))",
      zero_mod_f(sigma1(mul(U, U2)) - sigma1(U)*sigma1(U2)))
check("Q(uv) = Q(u) Q(v) (identity modulo f(L))", zero_mod_f(Qform(mul(U, U2)) - Qform(U)*Qform(U2)))
check("N(u) = sigma1(u) Q(u) (identity modulo f(L))", zero_mod_f(sigma1(U)*Qform(U) - norm(U)))
A_, B_ = sp.symbols('A B')
check("Q = (A + sB/2)^2 + (p - s^2/4) B^2 with p - s^2/4 = 3 l1^2/4 - 4 > 0, so Q >= 0, and Q(y) > 0 when N(y) != 0",
      sp.expand(A_**2 + s0*A_*B_ + p0*B_**2 - (A_ + s0*B_/2)**2 - (sp.Rational(3, 4)*L**2 - 4)*B_**2) == 0
      and Fr(3, 4)*hi1**2 - 4 > 0)
T2_gram = [[3, 0, 8], [0, 3*L**2 - 8, 4*L], [8, 4*L, -4*L**2 - 18*L + 32]]
check("T2 := sigma1^2 + 2Q = 3a^2 + 16ac + (3 l1^2 - 8) b^2 + 8 l1 bc + (-4 l1^2 - 18 l1 + 32) c^2 (modulo f(L))",
      zero_mod_f(sigma1(U)**2 + 2*Qform(U)
                 - sum(T2_gram[i][j]*U[i]*U[j] for i in range(3) for j in range(3))))
print("  2.5.2 the scaled descent lemma under EPS = (-19,-5,1) = l^2 - 5l - 19")
E_lo, E_hi = hi1**2 - 5*hi1 - 19, lo1**2 - 5*lo1 - 19
check(f"E := sigma1(EPS) = l1^2 - 5 l1 - 19 lies in [E-, E+] = [{dec(E_lo)}, {dec(E_hi)}] (t^2 - 5t - 19"
      " decreases for t < 5/2), and 0 < E < 1", 0 < E_lo <= E_hi < 1)
print("    N(EPS) = 1 = E Q(EPS) gives Q(EPS) = 1/E; EPS*EPS^-1 = 1 gives sigma1(EPS^-1) = 1/E, Q(EPS^-1) = E.")
print("    For y with r = sigma1(y)^2 and q = Q(y): T2(EPS y) = E^2 r + 2q/E and T2(EPS^-1 y) = r/E^2 + 2Eq.")
Es, r_, q_, tt = sp.symbols('E r q t', positive=True)
check("T2(EPS y) - T2(y) = (1 - E)(2q/E - (1 + E) r) and T2(EPS^-1 y) - T2(y) = (1 - E)((1 + E) r/E^2 - 2q)",
      sp.simplify(Es**2*r_ + 2*q_/Es - r_ - 2*q_ - (1 - Es)*(2*q_/Es - (1 + Es)*r_)) == 0
      and sp.simplify(r_/Es**2 + 2*Es*q_ - r_ - 2*q_ - (1 - Es)*((1 + Es)*r_/Es**2 - 2*q_)) == 0)
lam_s, mu_s = 2/(Es*(1 + Es)), 2*Es**2/(1 + Es)
check("with lambda = 2/(E(1+E)) and mu = 2E^2/(1+E): (1+E) lambda = 2/E and (1+E) mu/E^2 = 2; so r > lambda q gives"
      " T2(EPS y) < T2(y), and r < mu q gives T2(EPS^-1 y) < T2(y)",
      sp.simplify((1 + Es)*lam_s - 2/Es) == 0 and sp.simplify((1 + Es)*mu_s/Es**2 - 2) == 0)
psi = (tt + 2)**3/tt
check("psi(t) = (t+2)^3/t has psi'(t) = 2(t+2)^2(t-1)/t^2 (decreasing on (0,1], increasing on [1,oo)), and"
      " T2^3/N^2 = (r + 2q)^3/(r q^2) = psi(r/q)",
      sp.simplify(sp.diff(psi, tt) - 2*(tt + 2)**2*(tt - 1)/tt**2) == 0
      and sp.simplify((r_ + 2*q_)**3/(r_*q_**2) - psi.subs(tt, r_/q_)) == 0)
lam_hi = 2/(E_lo*(1 + E_lo))
mu_lo = 2*E_lo**2/(1 + E_lo)


def psi_q(t):
    return (t + 2)**3/t


K3 = max(psi_q(mu_lo), psi_q(lam_hi))
check(f"lambda <= lambda+ = 2/(E-(1+E-)) = {float(lam_hi):.6f} and mu >= mu- = 2E-^2/(1+E-) = {float(mu_lo):.6e}"
      " (lambda decreases and mu increases with E), and mu- < 1 < lambda+", mu_lo < 1 < lam_hi)


def cube_root_up(X, scale):
    """The least k/scale with (k/scale)^3 >= X (X > 0 rational)."""
    k = int(round(float(X)**(1/3)*scale)) - 2
    while Fr(k, scale)**3 < X:
        k += 1
    return Fr(k, scale)


C0_hi = cube_root_up(K3, 10**4)
print(f"    K3 := max(psi(mu-), psi(lambda+)) = {float(K3):.6f} (psi(mu-) = {float(psi_q(mu_lo)):.6f},"
      f" psi(lambda+) = {float(psi_q(lam_hi)):.6f}); C0 := K3^(1/3) <= {C0_hi} = {float(C0_hi)}")
print("    Lemma. For every triple y with N(y) != 0 (every rational y != 0, as N vanishes only at 0):")
print("    T2(EPS y) < T2(y), or T2(EPS^-1 y) < T2(y), or T2(y)^3 <= K3 N(y)^2.")
print("    Proof. N = sigma1 Q != 0 gives r > 0 and q > 0. r > lambda q and r < mu q are the first two cases.")
print("    Otherwise mu- <= mu <= r/q <= lambda <= lambda+, and T2(y)^3 = N(y)^2 psi(r/q) <= K3 N(y)^2, because psi")
print("    decreases and then increases, so it is largest at an end of [mu-, lambda+].")
print("    Orbit minimum. T2(EPS^k y) = E^(2k) r + 2 E^(-k) q tends to infinity as k -> +-infinity, so the orbit")
print("    {EPS^k y} has an element y0 of least T2; for it the lemma gives T2(y0)^3 <= K3 N(y0)^2 = K3 N(y)^2,")
print("    that is T2(y0) <= C0 |N(y)|^(2/3).")
print("  2.5.3 a rational lower bound for T2 and an exact enumeration of C")


class Iv:
    """A closed interval with rational ends; +, - and * enclose the exact results."""

    def __init__(self, lo, hi=None):
        self.lo, self.hi = Fr(lo), Fr(lo if hi is None else hi)

    @staticmethod
    def of(o):
        return o if isinstance(o, Iv) else Iv(o)

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
        p = (self.lo*o.lo, self.lo*o.hi, self.hi*o.lo, self.hi*o.hi)
        return Iv(min(p), max(p))

    __rmul__ = __mul__


l1 = Iv(lo1, hi1)
G_iv = [[Iv(3), Iv(0), Iv(8)], [Iv(0), 3*l1*l1 - 8, 4*l1], [Iv(8), 4*l1, -4*l1*l1 - 18*l1 + 32]]
RHO = Fr(1, 2**20)
G_mid = [[Fr(round((G_iv[i][j].lo + G_iv[i][j].hi)/2*2**20), 2**20) for j in range(3)] for i in range(3)]
check("each entry of the Gram matrix G of T2 lies within rho = 2^-20 of the rational matrix G' (entries rounded to"
      " multiples of 2^-20)",
      all(G_mid[i][j] - RHO <= G_iv[i][j].lo and G_iv[i][j].hi <= G_mid[i][j] + RHO for i in range(3) for j in range(3)))
G_low = [[G_mid[i][j] - (3*RHO if i == j else 0) for j in range(3)] for i in range(3)]
print("    so T2(y) - y^T G' y >= -rho (sum |y_i|)^2 >= -3 rho sum y_i^2, and T2(y) >= Glow(y) := y^T (G' - 3 rho I) y.")
print("    G' - 3 rho I =", [[str(t) for t in row] for row in G_low])
dd1 = G_low[0][0]
m12, m13 = G_low[0][1]/dd1, G_low[0][2]/dd1
dd2 = G_low[1][1] - dd1*m12**2
m23 = (G_low[1][2] - dd1*m12*m13)/dd2
dd3 = G_low[2][2] - dd1*m13**2 - dd2*m23**2
Y1, Y2, Y3 = sp.symbols('y1 y2 y3')
YY = (Y1, Y2, Y3)
check(f"Glow(y) = d1 (y1 + m12 y2 + m13 y3)^2 + d2 (y2 + m23 y3)^2 + d3 y3^2 with d1 = {float(dd1):.6f},"
      f" d2 = {float(dd2):.6f}, d3 = {float(dd3):.6f} all > 0 (exact), so Glow is positive definite",
      sp.expand(sum(G_low[i][j]*YY[i]*YY[j] for i in range(3) for j in range(3))
                - dd1*(Y1 + m12*Y2 + m13*Y3)**2 - dd2*(Y2 + m23*Y3)**2 - dd3*Y3**2) == 0
      and min(dd1, dd2, dd3) > 0)


def coset_range(x0, alpha, beta, gamma):
    """All x in x0 + Z with alpha*(x + beta)^2 <= gamma (alpha > 0). The float estimate of the range is certified:
    the scan stops at points that fail the test, one on each side of the centre -beta."""
    if gamma < 0:
        return []
    centre, rad = -beta, math.sqrt(float(gamma/alpha))
    n_lo = math.floor(float(centre - x0) - rad) - 1
    n_hi = math.ceil(float(centre - x0) + rad) + 1
    while not (x0 + n_lo < centre and alpha*(x0 + n_lo + beta)**2 > gamma):
        n_lo -= 1
    while not (x0 + n_hi > centre and alpha*(x0 + n_hi + beta)**2 > gamma):
        n_hi += 1
    return [x0 + n for n in range(n_lo + 1, n_hi) if alpha*(x0 + n + beta)**2 <= gamma]


def enumerate_coset(z0, T):
    """All y in z0 + Z^3 with Glow(y) <= T (exact)."""
    out = []
    for y3 in coset_range(z0[2], dd3, 0, T):
        r3 = T - dd3*y3**2
        for y2 in coset_range(z0[1], dd2, m23*y3, r3):
            r2 = r3 - dd2*(y2 + m23*y3)**2
            for yy1 in coset_range(z0[0], dd1, m12*y2 + m13*y3, r2):
                out.append((yy1, y2, y3))
    return out


T_e = cube_root_up(K3*Fr(121, 109)**2, 1000)
check(f"T_e = {T_e} = {float(T_e)} satisfies T_e^3 >= K3 (121/109)^2", T_e**3 >= K3*Fr(121, 109)**2)
listC = enumerate_coset(XS, T_e)
check("every listed y is in C and has Glow(y) <= T_e",
      all(all((yv - xv).denominator == 1 for yv, xv in zip(yy, XS)) for yy in listC)
      and all(sum(G_low[i][j]*yy[i]*yy[j] for i in range(3) for j in range(3)) <= T_e for yy in listC))
normsC = sorted((abs(norm(yy)), yy) for yy in listC)
minC = normsC[0][0]
check(f"the {len(listC)} elements y of C with Glow(y) <= T_e have min |N(y)| = {minC}, attained only at"
      f" {[tuple(str(t) for t in yy) for n, yy in normsC if n == minC]} = xi* - 2",
      minC == Fr(121, 109) and [yy for n, yy in normsC if n == minC] == [sub(XS, (2, 0, 0))])
print("    the next values of |N| on the list:", ", ".join(str(n) for n, _ in normsC[1:6]))
print("    Conclusion. Let y in C with |N(y)| < 121/109, and y0 the element of least T2 in its EPS-orbit (inside C).")
print("    Then Glow(y0)^3 <= T2(y0)^3 <= K3 N(y)^2 < K3 (121/109)^2 <= T_e^3, so y0 is on the list, where")
print("    |N| >= 121/109 > |N(y)| = |N(y0)|: impossible. So min over C of |N| = 121/109, attained at xi* - 2:")
print("    no kappa in R has |N(xi* - kappa)| < 1, and R is not norm-Euclidean.")
print("  2.5.4 the 654 classes of fixed points of EPS mod R (exploration)")
EM1 = sub(EPS, ONE)
NEM1 = norm(EM1)
check(f"EPS - 1 = {EM1} has norm {NEM1}; the fixed points of EPS mod R are the xi with (EPS - 1) xi in R, that is"
      f" xi in adj(EPS - 1) R / {NEM1}, adj(EPS - 1) = {adj(EM1)}", NEM1 == -654)


def reduce1(v):
    return tuple(t - math.floor(t) for t in v)


gens = [reduce1(tuple(Fr(t, NEM1) for t in mul(adj(EM1), ei))) for ei in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
classes = {(Fr(0), Fr(0), Fr(0))}
frontier = list(classes)
while frontier:
    new = []
    for v in frontier:
        for gv in gens:
            w = reduce1(tuple(p + q for p, q in zip(v, gv)))
            if w not in classes:
                classes.add(w)
                new.append(w)
    frontier = new
check(f"the group they generate modulo Z^3 has {len(classes)} = |N(EPS - 1)| classes, each fixed by EPS, and the class"
      " of xi* is one of them",
      len(classes) == 654 and all(all(t.denominator == 1 for t in mul(EM1, v)) for v in classes)
      and reduce1(XS) in classes)
T_1 = cube_root_up(K3, 1000)
check(f"T_1 = {T_1} satisfies T_1^3 >= K3", T_1**3 >= K3)
mins = []
for v in sorted(classes):
    if v == (0, 0, 0):
        continue
    lst = enumerate_coset(v, T_1)
    mins.append((min(abs(norm(yy)) for yy in lst), v, len(lst)))
mins.sort(reverse=True)
good = [m for m in mins if m[0] >= 1]
check(f"for each of the 653 classes zeta != 0, list y in zeta + R with Glow(y) <= T_1 ({min(m[2] for m in mins)} to"
      f" {max(m[2] for m in mins)} elements each): exactly {len(good)} classes have no listed y with |N(y)| < 1, the"
      " classes of xi* and -xi*",
      len(good) == 2 and {m[1] for m in good} == {reduce1(XS), reduce1(tuple(-t for t in XS))})
print("    By the orbit-minimum argument with 1 in place of 121/109, these 2 classes have no element with |N| < 1;")
print("    the other 651 contain one (listed). Largest values of min |N| over the listed elements:")
for m in mins[:6]:
    print(f"        {str(m[0]):>8s} = {float(m[0]):.6f}   class of {tuple(str(t) for t in m[1])}")

print("2.6 remark: Bernstein enclosures (a second, sharper certificate; not needed for the proof)")
uu = sp.symbols('u')


def bern(i, var):
    return sp.binomial(3, i)*var**i*(1 - var)**(3 - i)


TB = [[Fr(math.comb(I, i), math.comb(3, i)) if i <= I else Fr(0) for i in range(4)] for I in range(4)]
check("u^i = sum_(I >= i) C(I,i)/C(3,i) B_I(u) for i = 0..3, where B_I(u) = C(3,I) u^I (1-u)^(3-I) >= 0 on [0,1]"
      " and sum_I B_I(u) = 1",
      all(sp.expand(sum(sp.Rational(TB[I][i].numerator, TB[I][i].denominator)*bern(I, uu) for I in range(4))
                    - uu**i) == 0 for i in range(4))
      and sp.expand(sum(bern(I, uu) for I in range(4)) - 1) == 0)
print("    So on each box, with xi = lo + (hi - lo) u, u in [0,1]^3, the cubic P(u) = N(s*xi - kappa) is")
print("    sum b_IJK B_I(u1) B_J(u2) B_K(u3) with weights >= 0 adding up to 1: |N(s*xi - kappa)| <= max |b_IJK|,")
print("    and the corner coefficients (I, J, K in {0,3}) are the values at the corners. The monomial coefficients")
print("    of P come from the Taylor identity 2.1 at the corner v = s*lo - kappa, with w_i = s (hi_i - lo_i) u_i.")


def bernstein_coeffs(lo, hi, s, kap):
    A = tuple(s*lo[i] - kap[i] for i in range(3))
    Dl = tuple(s*(hi[i] - lo[i]) for i in range(3))
    mono = {(0, 0, 0): norm(A)}
    gA = grad(A)
    for i in range(3):
        ei = [0, 0, 0]
        ei[i] = 1
        mono[tuple(ei)] = gA[i]*Dl[i]
    for (i, j), co in quad(A).items():
        eij = [0, 0, 0]
        eij[i] += 1
        eij[j] += 1
        mono[tuple(eij)] = mono.get(tuple(eij), 0) + co*Dl[i]*Dl[j]
    for m, co in CUBIC.items():
        mono[m] = mono.get(m, 0) + co*Dl[0]**m[0]*Dl[1]**m[1]*Dl[2]**m[2]
    arr = [[[mono.get((i, j, k), 0) for k in range(4)] for j in range(4)] for i in range(4)]
    arr = [[[sum(TB[I][i]*arr[i][j][k] for i in range(I + 1)) for k in range(4)] for j in range(4)] for I in range(4)]
    arr = [[[sum(TB[J][j]*arr[I][j][k] for j in range(J + 1)) for k in range(4)] for J in range(4)] for I in range(4)]
    arr = [[[sum(TB[K][k]*arr[I][J][k] for k in range(K + 1)) for K in range(4)] for J in range(4)] for I in range(4)]
    return arr


overall = None
corner_ok = True
bern_max = []
for path, lo, hi, s, kap, B in BOXES:
    bc = bernstein_coeffs(lo, hi, s, kap)
    for e1, e2, e3 in itertools.product((0, 1), repeat=3):
        corner = (lo[0] if e1 == 0 else hi[0], lo[1] if e2 == 0 else hi[1], lo[2] if e3 == 0 else hi[2])
        if bc[3*e1][3*e2][3*e3] != norm(tuple(s*corner[i] - kap[i] for i in range(3))):
            corner_ok = False
    mx, where = max((abs(bc[I][J][K]), (I, J, K)) for I, J, K in itertools.product(range(4), repeat=3))
    bern_max.append(mx)
    if overall is None or mx > overall[0]:
        overall = (mx, where, path, lo, hi, s, kap)
check("on every box the corner coefficients equal the values of N(s*xi - kappa) at the corners", corner_ok)
check(f"on every box max |b_IJK| < 1; the largest over all boxes is {overall[0]} = {float(overall[0]):.6f}",
      all(m < 1 for m in bern_max))
mx, (I, J, K), path, lo, hi, s, kap = overall
if all(t in (0, 3) for t in (I, J, K)):
    corner = tuple(lo[i] if e == 0 else hi[i] for i, e in enumerate((I//3, J//3, K//3)))
    print(f"    it is the corner coefficient at xi = {tuple(str(t) for t in corner)} of box {path} (s = {s},"
          f" kappa = {kap}), so it is attained:")
    print(f"    with the table's choices, max over D of |N(s*xi - kappa)| = {mx} = {float(mx):.6f} exactly.")
else:
    print(f"    it is the coefficient {(I, J, K)} of box {path} (not a corner): with the table's choices,"
          f" |N(s*xi - kappa)| <= {mx} on D.")
ratio = max(b[5]/m for b, m in zip(BOXES, bern_max))
print(f"    Taylor bound / Bernstein bound on a box: at most {float(ratio):.2f};"
      f" largest Bernstein bound among the boxes with s = 1: {float(max(m for b, m in zip(BOXES, bern_max) if b[3] == 1)):.6f},"
      f" with s = 2: {float(max(m for b, m in zip(BOXES, bern_max) if b[3] == 2)):.6f}")

print("2.7 confirmation (not used by the proof): random Dedekind-Hasse divisions")
LEAF = {b[0]: (b[3], b[4]) for b in BOXES}
AXIS_OF = {}
for p in PATHS:
    for i in range(len(p)):
        AXIS_OF[p[:i]] = p[i].lower()


def locate(xi):
    node, lo, hi = '', list(D_LO), list(D_HI)
    while node not in LEAF:
        ch = AXIS_OF[node]
        ax = AXES.index(ch)
        mid = (lo[ax] + hi[ax])/2
        if xi[ax] <= mid:
            hi[ax], node = mid, node + ch
        else:
            lo[ax], node = mid, node + ch.upper()
    return LEAF[node]


def dh_step(alpha, beta):
    """s and kappa with 0 < |N(s alpha - kappa beta)| < |N(beta)|, as in Consequence 2.2; None if beta | alpha."""
    nb = norm(beta)
    xi = tuple(Fr(t, nb) for t in mul(alpha, adj(beta)))
    if all(t.denominator == 1 for t in xi):
        return None
    if all((2*t).denominator == 1 for t in xi):
        e = tuple(int(2*t) % 2 for t in xi)
        kb = HALF_LIST[e]
        return 1, tuple(int(xi[i] - Fr(e[i], 2)) + kb[i] for i in range(3)), 'b'
    r = tuple(math.floor(t + HALF) for t in xi)
    dlt = tuple(xi[i] - r[i] for i in range(3))
    sg = -1 if dlt[2] < 0 else 1
    xd = tuple(sg*t for t in dlt)
    s, kap = locate(xd)
    return s, tuple(s*r[i] + sg*kap[i] for i in range(3)), 'a'


rng = random.Random(20261015)
tests = []
for _ in range(3000):
    beta = tuple(rng.randint(-20, 20) for _ in range(3))
    alpha = tuple(rng.randint(-60, 60) for _ in range(3))
    if beta != (0, 0, 0):
        tests.append((alpha, beta))
for _ in range(300):
    bp = tuple(rng.randint(-10, 10) for _ in range(3))
    gam = tuple(rng.randint(-10, 10) for _ in range(3))
    e = rng.choice([v for v in itertools.product((0, 1), repeat=3) if v != (0, 0, 0)])
    if bp != (0, 0, 0):
        tests.append((mul(bp, tuple(2*gam[i] + e[i] for i in range(3))), tuple(2*t for t in bp)))
tests.append(((237, 247, 142), (327, 0, 0)))            # alpha/beta = xi*
tests.append(((-237, -247, -142), (327, 0, 0)))         # alpha/beta = -xi*
stats = {'a1': 0, 'a2': 0, 'b': 0, 'divides': 0}
ok = True
for alpha, beta in tests:
    res = dh_step(alpha, beta)
    if res is None:
        stats['divides'] += 1
        continue
    s, kappa, how = res
    rem_n = abs(norm(sub(tuple(s*t for t in alpha), mul(kappa, beta))))
    if not 0 < rem_n < abs(norm(beta)):
        ok = False
    stats['b' if how == 'b' else f'a{s}'] += 1
xs_s = dh_step((237, 247, 142), (327, 0, 0))[0]
check(f"{len(tests)} pairs: 0 < |N(s alpha - kappa beta)| < |N(beta)| whenever beta does not divide alpha"
      f" (s = 1 from the table: {stats['a1']}, s = 2 from the table: {stats['a2']}, from the list 2.4: {stats['b']},"
      f" beta | alpha: {stats['divides']}); alpha/beta = xi* uses s = {xs_s}", ok and xs_s == 2)
print(f"all checks of 05_dh_cover.py passed ({time.time() - T_START:.1f} s)")
