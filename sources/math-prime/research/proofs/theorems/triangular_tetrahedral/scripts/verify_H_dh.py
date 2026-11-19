"""verify_H_dh.py: an independent check of Section 2 of hypothesis_H.md (Lemma 2.1, Corollary 2.2, Division 2.3,
Remarks 2.4 and 2.5) for R = Z[l], l^3 = 4l - 6, where the triple (a, b, c) stands for a + b*l + c*l^2.

Written 2026-10-15 by Claude (Opus 5.5) as an independent verifier.

Independence. This script imports none of cubic.py, core.py, lifting2.py or the 0x_*.py scripts and copies no code
from them; the only project file it reads is the table dh_boxes.txt. Every decision is made in exact arithmetic
(Python int and Fraction). Floats appear only in printed approximations and in one cross-check labelled as a guide.

  0   own arithmetic: the product is polynomial multiplication followed by long division by f(t) = t^3 - 4t + 6;
      the norm is the determinant of the multiplication matrix. A polynomial identity is proved by evaluating it on
      a grid with more points per variable than its degree in that variable.
  1   the real root l1 of f: Newton's iteration in rationals, certified by an exact sign change; f has no other
      real root and no rational root.
  (a) the table: format, every box inside D = [-1/2,1/2] x [-1/2,1/2] x [0,1/2], s in {1,2}, kappa integral;
      the counts of Remark 2.4.
  (b) Lemma 2.1(a) box by box, by N = sigma1 * Q. sigma1 is linear and Q is a positive semidefinite quadratic form,
      so |sigma1| and Q are convex and are largest at a corner of a box: sup |N| <= max|sigma1| * max Q, evaluated
      with a rational enclosure of l1, and boxes are split dyadically where that product is >= 1 (not 05's Taylor
      bound). Also: a bracket for the largest value of |N(s xi - kappa)| over the table, and the text's own Taylor
      data recomputed from the formulas printed in Section 2.1.
  (c) the cover of D: exact volume sum and pairwise disjoint interiors; a sweep; random rational points of D.
  (d) Lemma 2.1(b), and Corollary 2.2 / Division 2.3 on 2000 pairs (alpha, beta), some with alpha/beta in xi* + R.
  (e) Remark 2.5: EPS xi* - xi* in R, N(xi* - 2) = -121/109, and no y in C = xi* + R has |N(y)| < 121/109, by a
      normalisation with powers of EPS (1 <= |sigma1(y)| < E) and an exact enumeration (not 05's T2 descent); also
      the 654 classes of fixed points of EPS.
  (f) negative controls: planted errors that the checks above must catch.
"""
import itertools
import math
import os
import random
import re
import sys
import time
from fractions import Fraction as Fr

import numpy as np

T_START = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
TABLE_PATH = os.path.join(HERE, 'dh_boxes.txt')
PASSED = [0]


def check(name, cond):
    """Print PASS or FAIL; stop at the first FAIL."""
    print(("PASS  " if cond else "FAIL  ") + name, flush=True)
    if not cond:
        print(f"STOPPED at the first FAIL ({PASSED[0]} checks had passed)", flush=True)
        sys.exit(1)
    PASSED[0] += 1


def note(*lines):
    for line in lines:
        print("    " + line)


def fl(x, digits=9):
    return f"{float(x):.{digits}g}"


def tstr(v):
    return "(" + ", ".join(str(t) for t in v) + ")"


# =================================================================================================================
print("0. own arithmetic in R_Q = Q[t]/(f), f(t) = t^3 - 4t + 6; a triple (a, b, c) is a + b t + c t^2")
FPOLY = (6, -4, 0, 1)                     # f = 6 - 4t + 0t^2 + t^3, monic
BASIS = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
ONE = (1, 0, 0)
ZERO = (0, 0, 0)
GRID4 = list(itertools.product(range(4), repeat=3))
GRID3 = list(itertools.product(range(3), repeat=3))


def pmul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def prem(p):
    """Remainder of the polynomial p (coefficients, lowest degree first) on long division by the monic f."""
    r = list(p) + [0] * max(0, 3 - len(p))
    for k in range(len(r) - 1, 2, -1):
        c = r[k]
        if c:
            for i in range(4):
                r[k - 3 + i] -= c * FPOLY[i]
    return tuple(r[:3])


def rmul(u, v):
    return prem(pmul(u, v))


def radd(u, v):
    return tuple(x + y for x, y in zip(u, v))


def rsub(u, v):
    return tuple(x - y for x, y in zip(u, v))


def rsc(r, u):
    return tuple(r * x for x in u)


def rpow(u, k):
    out = ONE
    for _ in range(k):
        out = rmul(out, u)
    return out


def mat(u):
    """Matrix of multiplication by u: column j is u * t^j (rows are coordinates)."""
    cols = [rmul(u, e) for e in BASIS]
    return [[cols[j][i] for j in range(3)] for i in range(3)]


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


# The formulas as printed in hypothesis_H.md (1.1), (1.2), (1.3), typed from the text:
def text_product(u, v):
    a, b, c = u
    d, e, f = v
    return (a * d - 6 * b * f - 6 * c * e, a * e + b * d + 4 * b * f + 4 * c * e - 6 * c * f,
            a * f + b * e + c * d + 4 * c * f)


def text_matrix(u):
    a, b, c = u
    return [[a, -6 * c, -6 * b], [b, a + 4 * c, 4 * b - 6 * c], [c, b, a + 4 * c]]


def ncub(v):
    """N(a,b,c) = a^3 + 8a^2c - 4ab^2 + 18abc + 16ac^2 - 6b^3 + 24bc^2 + 36c^3 (1.2)."""
    a, b, c = v
    return (a ** 3 + 8 * a * a * c - 4 * a * b * b + 18 * a * b * c + 16 * a * c * c - 6 * b ** 3
            + 24 * b * c * c + 36 * c ** 3)


def text_adj(v):
    a, b, c = v
    return (a * a + 8 * a * c - 4 * b * b + 6 * b * c + 16 * c * c, -a * b - 6 * c * c, b * b - a * c - 4 * c * c)


check("t*t^2 = t^3 = -6 + 4t and t^2*t^2 = t^4 = -6t + 4t^2 by long division by f",
      rmul(BASIS[1], BASIS[2]) == (-6, 4, 0) and rmul(BASIS[2], BASIS[2]) == (0, -6, 4))
check("the product printed in (1.1) equals polynomial product mod f on the 9 pairs of basis triples; both are"
      " bilinear, so they agree on all triples",
      all(text_product(e1, e2) == rmul(e1, e2) for e1 in BASIS for e2 in BASIS))
check("the product is commutative and associative on basis triples (9 pairs, 27 triples), hence everywhere"
      " (bi-/trilinearity), and (1,0,0) is its unit",
      all(rmul(e1, e2) == rmul(e2, e1) for e1 in BASIS for e2 in BASIS)
      and all(rmul(rmul(e1, e2), e3) == rmul(e1, rmul(e2, e3)) for e1 in BASIS for e2 in BASIS for e3 in BASIS)
      and all(rmul(ONE, e) == e for e in BASIS))
check("the matrix printed in (1.2) is the matrix of multiplication (columns u*1, u*l, u*l^2): equal on the basis,"
      " both linear in u", all(mat(e) == text_matrix(e) for e in BASIS))
check("det M(a,b,c) equals the printed cubic form N(a,b,c) on {0,1,2,3}^3; both are cubic forms (degree <= 3 in"
      " each variable), so they are the same polynomial", all(det3(mat(v)) == ncub(v) for v in GRID4))
check("u * adj(u) = (N(u), 0, 0) with adj as printed in (1.3), on {0,1,2,3}^3 (degree <= 3 in each variable), so it"
      " is an identity", all(rmul(v, text_adj(v)) == (ncub(v), 0, 0) for v in GRID4))
note("Multiplicativity: (uv)x = u(vx) for all x (associative, commutative), so M(uv) = M(u)M(v) and",
     "N(uv) = det M(uv) = N(u)N(v) for all rational triples; N(-u) = -N(u) since N is a cubic form.")
_rng0 = random.Random(1)
_rt = [tuple(Fr(_rng0.randint(-40, 40), _rng0.randint(1, 9)) for _ in range(3)) for _ in range(60)]
check("sanity: N(uv) = N(u)N(v) on 30 random pairs of rational triples",
      all(ncub(rmul(_rt[2 * i], _rt[2 * i + 1])) == ncub(_rt[2 * i]) * ncub(_rt[2 * i + 1]) for i in range(30)))

# =================================================================================================================
print("1. the real root l1 of f(t) = t^3 - 4t + 6")
P = 100
TWO_P = 1 << P


def f_at(t):
    return t ** 3 - 4 * t + 6


def f_scaled(L):
    """2^(3P) f(L / 2^P), an integer."""
    return L ** 3 - 4 * L * (1 << (2 * P)) + 6 * (1 << (3 * P))


_t = Fr(-5, 2)
for _ in range(10):
    _t = _t - f_at(_t) / (3 * _t * _t - 4)
    _t = Fr(round(_t * (1 << (2 * P))), 1 << (2 * P))
_L = math.floor(_t * TWO_P)
while f_scaled(_L) >= 0:
    _L -= 1
while f_scaled(_L + 1) < 0:
    _L += 1
LLO, LHI = _L, _L + 1
LAM_LO, LAM_HI = Fr(LLO, TWO_P), Fr(LHI, TWO_P)
check(f"Newton's iteration from -5/2 (rounded to 2^-{2 * P}), then an exact sign change: f(lo) < 0 < f(hi) for"
      f" [lo, hi] = [L, L+1]/2^{P}, so f has a root l1 in it (width 2^-{P}); l1 = {float(LAM_LO):.16f}",
      f_scaled(LLO) < 0 < f_scaled(LHI) and LHI < 0)
check("l1 is the only real root: f(t) = (t+2) t (t-2) + 6 (a cubic identity, checked at t = 0,1,2,3). On [-2,0] and"
      " [2,oo) the product is >= 0; on [0,2] it is >= -4 (t(t-2) = (t-1)^2 - 1 >= -1, 0 < t+2 <= 4); so f >= 2 on"
      " [-2,oo). On (-oo,-2], f' = 3t^2 - 4 >= 8 > 0, and f(-3) = -9 < 0 < 6 = f(-2); our interval lies in (-3,-2)",
      all((t + 2) * t * (t - 2) + 6 == f_at(t) for t in range(4)) and f_at(-3) == -9 and f_at(-2) == 6
      and -3 * TWO_P < LLO and LHI < -2 * TWO_P)
check("f has no rational root (a rational root of a monic integer polynomial is an integer dividing 6:"
      " f(r) != 0 for r = +-1, +-2, +-3, +-6), so the cubic f is irreducible over Q and 1, l1, l1^2 are linearly"
      " independent over Q: sigma1(y) = y0 + y1 l1 + y2 l1^2 != 0 for every rational triple y != 0",
      all(f_at(r) != 0 for r in (1, -1, 2, -2, 3, -3, 6, -6)))


def lam_iv(c0, c1, c2):
    """[lo, hi] containing c0 + c1 l1 + c2 l1^2 (rational c); l1 in [LAM_LO, LAM_HI] < 0, l1^2 in [LAM_HI^2, LAM_LO^2]."""
    t1 = sorted((c1 * LAM_LO, c1 * LAM_HI))
    t2 = sorted((c2 * LAM_HI ** 2, c2 * LAM_LO ** 2))
    return c0 + t1[0] + t2[0], c0 + t1[1] + t2[1]


def sig_iv_scaled(c0, c1, c2):
    """For integers c: integers lo <= 2^(2P) (c0 + c1 l1 + c2 l1^2) <= hi."""
    p1, p2 = (c1 * LLO) << P, (c1 * LHI) << P
    q1, q2 = c2 * LHI * LHI, c2 * LLO * LLO
    base = c0 << (2 * P)
    return base + min(p1, p2) + min(q1, q2), base + max(p1, p2) + max(q1, q2)


# =================================================================================================================
print("(a) the table dh_boxes.txt")
D_LO, D_HI = (-32, -32, 0), (32, 32, 32)       # D in units of 2^-6
VOL_D = (D_HI[0] - D_LO[0]) * (D_HI[1] - D_LO[1]) * (D_HI[2] - D_LO[2])


def parse_table(path):
    unit, rows, bad = None, [], []
    with open(path, encoding='ascii') as fh:
        for ln, line in enumerate(fh, 1):
            tok = line.split()
            if not tok or tok[0].startswith('#'):
                continue
            if tok[0] == 'unit':
                unit = ' '.join(tok[1:])
                continue
            if (len(tok) != 11 or not re.fullmatch(r'[xyzXYZ]+', tok[0])
                    or not all(re.fullmatch(r'[+-]?[0-9]+', x) for x in tok[1:])):
                bad.append(ln)
                continue
            v = [int(x) for x in tok[1:]]
            rows.append({'line': ln, 'path': tok[0], 'lo': tuple(v[0:3]), 'hi': tuple(v[3:6]), 's': v[6],
                         'kap': tuple(v[7:10])})
    return unit, rows, bad


def row_problem(r):
    if not all(D_LO[i] <= r['lo'][i] < r['hi'][i] <= D_HI[i] for i in range(3)):
        return "not a nondegenerate box inside D"
    if r['s'] not in (1, 2):
        return "s not in {1, 2}"
    if not all(isinstance(k, int) for k in r['kap']):
        return "kappa not integral"
    return None


def path_box(path):
    lo, hi = [Fr(x) for x in D_LO], [Fr(x) for x in D_HI]
    for ch in path:
        i = 'xyz'.index(ch.lower())
        mid = (lo[i] + hi[i]) / 2
        if ch.islower():
            hi[i] = mid
        else:
            lo[i] = mid
    return tuple(lo), tuple(hi)


UNIT, ROWS, BAD = parse_table(TABLE_PATH)
check(f"the header says 'unit 2^-6', and every data line has the form: path (letters x y z X Y Z), six integer"
      f" corner coordinates, s, three integers kappa ({len(ROWS)} lines, {len(BAD)} malformed)",
      UNIT == '2^-6' and not BAD and len(ROWS) > 0)
check(f"{len(ROWS)} boxes, as stated in Section 2.1 (398)", len(ROWS) == 398)
_probs = [(r['path'], row_problem(r)) for r in ROWS if row_problem(r)]
check("every box lies in D: -32 <= lo_a < hi_a <= 32, -32 <= lo_b < hi_b <= 32, 0 <= lo_c < hi_c <= 32 (units of"
      " 1/64); s in {1, 2}; kappa is an integer triple (parsed as integer literals)", not _probs)
check("file consistency (not used below): each path, read as midpoint splits of D, gives the listed corners;"
      " no path occurs twice",
      all(path_box(r['path']) == (r['lo'], r['hi']) for r in ROWS) and len({r['path'] for r in ROWS}) == len(ROWS))
_n1 = sum(1 for r in ROWS if r['s'] == 1)
_kmax = max(max(abs(k) for k in r['kap']) for r in ROWS)
_k4 = sum(1 for r in ROWS if max(abs(k) for k in r['kap']) <= 4)
_dmax = max(len(r['path']) for r in ROWS)
check(f"Remark 2.4 counts: s = 1 on {_n1} boxes and s = 2 on {len(ROWS) - _n1} (text: 263, 135); all |kappa_i| <= 4"
      f" on {_k4} boxes (text: 375); largest |kappa_i| = {_kmax} (text: 22); longest path {_dmax} splits (text: 15)",
      (_n1, len(ROWS) - _n1, _k4, _kmax, _dmax) == (263, 135, 375, 22, 15))

# =================================================================================================================
print("(b) Lemma 2.1(a) on every box, by N = sigma1 * Q (a different bound from the Taylor bound of 05)")
# Q = a^2 - L ab + (8 - L^2) ac + (L^2 - 4) b^2 + 6 bc + (16 - 6L - 4L^2) c^2: each coefficient as an element
# c0 + c1 L + c2 L^2 of Q[L]/(f).
QCOEF = {(2, 0, 0): (1, 0, 0), (1, 1, 0): (0, -1, 0), (1, 0, 1): (8, 0, -1),
         (0, 2, 0): (-4, 0, 1), (0, 1, 1): (6, 0, 0), (0, 0, 2): (16, -6, -4)}


def q_triple(v, coef=QCOEF):
    """Q(v) for a numeric triple v, as (A, B, C) = A + B L + C L^2; its real value is A + B l1 + C l1^2."""
    A = B = C = 0
    for (i, j, k), (c0, c1, c2) in coef.items():
        mono = v[0] ** i * v[1] ** j * v[2] ** k
        A += c0 * mono
        B += c1 * mono
        C += c2 * mono
    return (A, B, C)


def q_ring(v, coef=QCOEF):
    """Q(v) computed in Q[L]/(f) when the coordinates of v are themselves elements of Q[L]/(f)."""
    total = ZERO
    for (i, j, k), co in coef.items():
        term = co
        for x, e in zip(v, (i, j, k)):
            for _ in range(e):
                term = rmul(term, x)
        total = radd(total, term)
    return total


def nq_identity(coef):
    return all(rmul(v, q_triple(v, coef)) == (ncub(v), 0, 0) for v in GRID4)


check("N = sigma1 * Q in Q[L]/(f)[a,b,c], with sigma1 = a + bL + cL^2 and Q = a^2 - L ab + (8 - L^2) ac + (L^2 - 4)"
      " b^2 + 6 bc + (16 - 6L - 4L^2) c^2: (a,b,c) * Q(a,b,c) = (N(a,b,c), 0, 0) on {0,1,2,3}^3, degree <= 3 in each"
      " variable, so it is an identity; evaluating L at the root l1 gives N(v) = sigma1(v) Q(v) for all real v",
      nq_identity(QCOEF))


def sos_identity():
    for a, b, c in GRID3:
        w = (Fr(a + 4 * c), Fr(-b, 2), Fr(-c, 2))          # a - (L/2) b + ((8 - L^2)/2) c
        u = (b, -c, 0)                                     # b - L c
        k = (-4, 0, Fr(3, 4))                              # (3L^2 - 16)/4
        if radd(rmul(w, w), rmul(k, rmul(u, u))) != q_triple((a, b, c)):
            return False
    return True


check("Q = (a - (L/2) b + ((8 - L^2)/2) c)^2 + ((3L^2 - 16)/4) (b - L c)^2 in Q[L]/(f)[a,b,c]: equal on {0,1,2}^3,"
      " degree <= 2 in each variable", sos_identity())
K_LO = 3 * LAM_HI ** 2 - 16
check(f"3 l1^2 - 16 >= 3 hi^2 - 16 = {fl(K_LO)} > 0 (|l1| >= |hi|). So at L = l1, Q(x) >= 0 for every real x: Q is"
      f" positive semidefinite", K_LO > 0)
check("Q is only semidefinite, as (1.6) says: Q(L^2 - 4, L, 1) = 0 in Q[L]/(f) (and Q(1,0,0) = 1)",
      q_ring(((-4, 0, 1), (0, 1, 0), ONE)) == ZERO and q_triple(ONE) == ONE)
note("Convexity. For a positive semidefinite form, Q(tx + (1-t)y) = tQ(x) + (1-t)Q(y) - t(1-t)Q(x-y) <= tQ(x) +",
     "(1-t)Q(y): Q is convex. |sigma1| is the absolute value of a linear form: convex. xi -> s xi - kappa is affine,",
     "so on a box both xi -> |sigma1(s xi - kappa)| and xi -> Q(s xi - kappa) are convex; every point of a box is a",
     "convex combination of its 8 corners, so each is at most its largest corner value. With N = sigma1 Q, Q >= 0:",
     "  sup over the box of |N(s xi - kappa)| <= max_corners |sigma1(s xi - kappa)| * max_corners Q(s xi - kappa).",
     f"The corner values are enclosed with l1 in [L, L+1]/2^{P}, in exact integer arithmetic.")


def corners(lo, hi):
    return [(hi[0] if e0 else lo[0], hi[1] if e1 else lo[1], hi[2] if e2 else lo[2])
            for e0 in (0, 1) for e1 in (0, 1) for e2 in (0, 1)]


def box_data(lo, hi, m, s, kap):
    """For the box [lo, hi]/2^m and (s, kappa): (S, Qm, cmax, cX) with |sigma1(s xi - kappa)| <= S / 2^(m+2P) and
    Q(s xi - kappa) <= Qm / 2^(2m+2P) at every corner xi, and cmax / 2^(3m) = the largest exact |N| at a corner cX."""
    S = Qm = cmax = 0
    cX = None
    for X in corners(lo, hi):
        V = (s * X[0] - (kap[0] << m), s * X[1] - (kap[1] << m), s * X[2] - (kap[2] << m))
        lo1, hi1 = sig_iv_scaled(*V)
        S = max(S, -lo1, hi1)
        Qm = max(Qm, sig_iv_scaled(*q_triple(V))[1])
        n = abs(ncub(V))
        if cX is None or n > cmax:
            cmax, cX = n, X
    return S, Qm, cmax, cX


def certify_box(r, tau=Fr(1), max_depth=40):
    """Prove sup over box r of |N(s xi - kappa)| < tau by the corner product bound, splitting into 8 octants when
    the bound is >= tau. Fails at a corner with exact |N| >= tau, or at the depth limit."""
    s, kap = r['s'], r['kap']
    stack = [(r['lo'], r['hi'], 6, 0)]
    info = {'leaves': 0, 'depth': 0, 'worst': Fr(0), 'corner': (Fr(-1), None), 'vol': Fr(0)}
    while stack:
        lo, hi, m, d = stack.pop()
        S, Qm, cmax, cX = box_data(lo, hi, m, s, kap)
        cval = Fr(cmax, 1 << (3 * m))
        if cval > info['corner'][0]:
            info['corner'] = (cval, tuple(Fr(x, 1 << m) for x in cX))
        bound = Fr(S * Qm, 1 << (3 * m + 4 * P))
        if bound < tau:
            info['leaves'] += 1
            info['depth'] = max(info['depth'], d)
            info['worst'] = max(info['worst'], bound)
            info['vol'] += Fr((hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]), 1 << (3 * m))
            continue
        if cval >= tau:
            info['fail'] = (f"box {r['path']} (s = {s}, kappa = {kap}): at the corner xi = "
                            f"{tstr(Fr(x, 1 << m) for x in cX)}, |N(s xi - kappa)| = {cval} = {fl(cval)} >= {tau}")
            return False, info
        if d >= max_depth:
            info['fail'] = f"box {r['path']}: depth limit {max_depth} reached"
            return False, info
        lo2, hi2 = tuple(2 * x for x in lo), tuple(2 * x for x in hi)
        mid = tuple(lo[i] + hi[i] for i in range(3))
        for e in itertools.product((0, 1), repeat=3):
            stack.append((tuple(mid[i] if e[i] else lo2[i] for i in range(3)),
                          tuple(hi2[i] if e[i] else mid[i] for i in range(3)), m + 1, d + 1))
    box_vol = Fr((r['hi'][0] - r['lo'][0]) * (r['hi'][1] - r['lo'][1]) * (r['hi'][2] - r['lo'][2]), 1 << 18)
    info['vol_ok'] = info['vol'] == box_vol
    return True, info


def bound_all(rows, tau=Fr(1)):
    agg = {'direct': 0, 'split': 0, 'leaves': 0, 'depth': 0, 'worst': Fr(0), 'corner': (Fr(-1), None, None),
           'vol_ok': True}
    for r in rows:
        ok, info = certify_box(r, tau)
        if not ok:
            return False, info['fail']
        agg['direct' if info['leaves'] == 1 else 'split'] += 1
        agg['leaves'] += info['leaves']
        agg['depth'] = max(agg['depth'], info['depth'])
        agg['worst'] = max(agg['worst'], info['worst'])
        agg['vol_ok'] = agg['vol_ok'] and info['vol_ok']
        if info['corner'][0] > agg['corner'][0]:
            agg['corner'] = (info['corner'][0], info['corner'][1], r)
    return True, agg


_t0 = time.time()
OK_B, AGG = bound_all(ROWS)
check(f"Lemma 2.1(a): on every one of the {len(ROWS)} boxes, sup |N(s xi - kappa)| < 1 by the corner product bound."
      + (f" {AGG['direct']} boxes directly; {AGG['split']} boxes needed dyadic subdivision (octants); {AGG['leaves']}"
         f" sub-boxes in all, at most {AGG['depth']} splits deep; largest bound on a sub-box {fl(AGG['worst'])};"
         f" the sub-boxes of each box fill it (volumes add up) [{time.time() - _t0:.1f} s]"
         if OK_B else f" FAILED: {AGG}"), OK_B and AGG['vol_ok'])
_cv, _cx, _cr = AGG['corner']
check(f"the largest exact corner value met is |N(s xi - kappa)| = {_cv} = {fl(_cv)} at xi = {tstr(_cx)}, box"
      f" {_cr['path']} (s = {_cr['s']}, kappa = {_cr['kap']}): Remark 2.4 names 4093/4096 at (-5/16, -5/16, 1/4),"
      f" s = 1, kappa = (-1, 0, 0)",
      _cv == Fr(4093, 4096) and _cx == (Fr(-5, 16), Fr(-5, 16), Fr(1, 4)) and _cr['s'] == 1
      and _cr['kap'] == (-1, 0, 0))
TAU_SHARP = Fr(4093, 4096) + Fr(1, 2 ** 20)
_t0 = time.time()
OK_S, AGG_S = bound_all(ROWS, TAU_SHARP)
check(f"bracket for Remark 2.4: refining until every sub-box bound is < 4093/4096 + 2^-20 succeeds"
      + (f" ({AGG_S['leaves']} sub-boxes, at most {AGG_S['depth']} splits deep) [{time.time() - _t0:.1f} s]. So the"
         f" largest value of |N(s xi - kappa)| over the boxes, each with its own (s, kappa), lies in [4093/4096,"
         f" 4093/4096 + 2^-20): the stated sharp value 4093/4096 is confirmed to within 2^-20"
         if OK_S else f" FAILED: {AGG_S}"), OK_S)

print("  (b2) the text's own Taylor data, typed from Section 2.1 and recomputed (a check of the text, not of the lemma)")


def g_text(v):
    x, y, z = v
    return (3 * x * x + 16 * x * z - 4 * y * y + 18 * y * z + 16 * z * z,
            -8 * x * y + 18 * x * z - 18 * y * y + 24 * z * z,
            8 * x * x + 18 * x * y + 32 * x * z + 48 * y * z + 108 * z * z)


def qq_text(v):
    x, y, z = v
    return {(0, 0): 3 * x + 8 * z, (0, 1): -8 * y + 18 * z, (0, 2): 16 * x + 18 * y + 32 * z,
            (1, 1): -4 * x - 18 * y, (1, 2): 18 * x + 48 * z, (2, 2): 16 * x + 24 * y + 108 * z}


NCOEF = {(3, 0, 0): 1, (2, 0, 1): 8, (1, 2, 0): -4, (1, 1, 1): 18, (1, 0, 2): 16, (0, 3, 0): -6, (0, 1, 2): 24,
         (0, 0, 3): 36}


def taylor_ok():
    for v in GRID4:
        g, qq, nv = g_text(v), qq_text(v), ncub(v)
        for w in GRID4:
            rhs = (nv + sum(g[i] * w[i] for i in range(3)) + sum(c * w[i] * w[j] for (i, j), c in qq.items())
                   + ncub(w))
            if ncub(radd(v, w)) != rhs:
                return False
    return True


check("N(v + w) = N(v) + sum g_i(v) w_i + sum_{i<=j} q_ij(v) w_i w_j + N(w) with g_i, q_ij exactly as printed in"
      " Section 2.1: equal on {0,1,2,3}^6 (4096 points, degree <= 3 in each of the 6 variables), so an identity",
      taylor_ok())
check("the ten coefficients n_m of N are those of (1.2) (sum of |n_m| = 113)",
      all(sum(c * v[0] ** m[0] * v[1] ** m[1] * v[2] ** m[2] for m, c in NCOEF.items()) == ncub(v) for v in GRID4)
      and sum(abs(c) for c in NCOEF.values()) == 113)


def taylor_B(r):
    s, kap = r['s'], r['kap']
    c = tuple(Fr(r['lo'][i] + r['hi'][i], 128) for i in range(3))
    h = tuple(Fr(s * (r['hi'][i] - r['lo'][i]), 128) for i in range(3))
    v = tuple(s * c[i] - kap[i] for i in range(3))
    g, qq = g_text(v), qq_text(v)
    return (abs(ncub(v)) + sum(abs(g[i]) * h[i] for i in range(3))
            + sum(abs(co) * h[i] * h[j] for (i, j), co in qq.items())
            + sum(abs(co) * h[0] ** m[0] * h[1] ** m[1] * h[2] ** m[2] for m, co in NCOEF.items()))


_TB = [(taylor_B(r), r['path']) for r in ROWS]
_tbmax = max(_TB)
check(f"the text's bound B, recomputed from these formulas: B < 1 on all {len(ROWS)} boxes, largest B = {_tbmax[0]}"
      f" on box {_tbmax[1]} (text: 8191/8192)",
      all(b < 1 for b, _ in _TB) and _tbmax[0] == Fr(8191, 8192))

# =================================================================================================================
print("(c) the boxes cover D (not the split-tree replay nor the 1/64-cell count of 05)")


def volume_units(rows):
    return sum((r['hi'][0] - r['lo'][0]) * (r['hi'][1] - r['lo'][1]) * (r['hi'][2] - r['lo'][2]) for r in rows)


def first_overlap(rows):
    """The first pair of boxes whose interiors meet, or None (all pairs are tested)."""
    for i in range(len(rows)):
        a = rows[i]
        for j in range(i + 1, len(rows)):
            b = rows[j]
            if all(a['lo'][k] < b['hi'][k] and b['lo'][k] < a['hi'][k] for k in range(3)):
                return (a['path'], b['path'])
    return None


def sweep(rows):
    """Sweep in x, then in y: in every cell [xa,xb] x [ya,yb] between consecutive breakpoints, the z-intervals of the
    boxes that span the cell must chain from 0 to 1/2. Returns (first gap, first z-overlap, number of cells)."""
    xs = sorted({D_LO[0], D_HI[0]} | {r['lo'][0] for r in rows} | {r['hi'][0] for r in rows})
    xs = [x for x in xs if D_LO[0] <= x <= D_HI[0]]
    gap = overlap = None
    cells = 0
    for xa, xb in zip(xs, xs[1:]):
        act = [r for r in rows if r['lo'][0] <= xa and r['hi'][0] >= xb]
        ys = sorted({D_LO[1], D_HI[1]} | {r['lo'][1] for r in act} | {r['hi'][1] for r in act})
        ys = [y for y in ys if D_LO[1] <= y <= D_HI[1]]
        for ya, yb in zip(ys, ys[1:]):
            cells += 1
            ivs = sorted((r['lo'][2], r['hi'][2]) for r in act if r['lo'][1] <= ya and r['hi'][1] >= yb)
            reach = D_LO[2]
            for lz, hz in ivs:
                if lz > reach and gap is None:
                    gap = ((xa, xb), (ya, yb), (reach, lz))
                if lz < reach and overlap is None:
                    overlap = ((xa, xb), (ya, yb), (lz, reach))
                reach = max(reach, hz)
            if reach < D_HI[2] and gap is None:
                gap = ((xa, xb), (ya, yb), (reach, D_HI[2]))
    return gap, overlap, cells


def points_check(rows, pts, chunk=1000):
    """Every point (numerators n, denominators d) of D lies in some box, and |N(s xi - kappa)| < 1 exactly for every
    box containing it."""
    LO = np.array([r['lo'] for r in rows], dtype=np.int64)
    HI = np.array([r['hi'] for r in rows], dtype=np.int64)
    NUM = np.array([p[0] for p in pts], dtype=np.int64)
    DEN = np.array([p[1] for p in pts], dtype=np.int64)
    pairs = multi = 0
    worst = (Fr(0), None)
    for st in range(0, len(pts), chunk):
        x64 = 64 * NUM[st:st + chunk]
        den = DEN[st:st + chunk]
        inside = np.all((LO[None, :, :] * den[:, None, :] <= x64[:, None, :])
                        & (x64[:, None, :] <= HI[None, :, :] * den[:, None, :]), axis=2)
        for i in range(inside.shape[0]):
            idx = np.flatnonzero(inside[i])
            n, d = pts[st + i]
            if idx.size == 0:
                return False, f"the point {tstr(Fr(n[t], d[t]) for t in range(3))} lies in no box"
            multi += idx.size > 1
            Dc = math.lcm(*d)
            for j in idx:
                r = rows[int(j)]
                s, kap = r['s'], r['kap']
                V = tuple(s * n[t] * (Dc // d[t]) - kap[t] * Dc for t in range(3))
                nv = abs(ncub(V))
                pairs += 1
                if not nv < Dc ** 3:
                    return False, (f"at xi = {tstr(Fr(n[t], d[t]) for t in range(3))}, box {r['path']}:"
                                   f" |N(s xi - kappa)| = {fl(Fr(nv, Dc ** 3))} >= 1")
                if nv * worst[0].denominator > worst[0].numerator * Dc ** 3:
                    worst = (Fr(nv, Dc ** 3), r['path'])
    return True, {'pairs': pairs, 'multi': multi, 'worst': worst}


VOL = volume_units(ROWS)
check(f"the volumes of the boxes add up exactly to vol(D) = 1/2 (sum = {VOL}/2^18 = {Fr(VOL, 1 << 18)})", VOL == VOL_D)
_t0 = time.time()
_ov = first_overlap(ROWS)
check(f"the interiors of the boxes are pairwise disjoint: all {len(ROWS) * (len(ROWS) - 1) // 2} pairs tested"
      f" (some axis separates the two open boxes) [{time.time() - _t0:.1f} s]", _ov is None)
note("Hence the boxes cover D: their union U is closed and lies in D. If some p in D were not in U, a small closed",
     "cube K in D near p would miss U, and vol(D) >= vol(U) + vol(K) > vol(U) = sum of the volumes (interiors",
     "disjoint) = vol(D), a contradiction.")
_gap, _zov, _cells = sweep(ROWS)
check(f"second argument, a sweep: in each of the {_cells} cells [xa,xb] x [ya,yb] cut out by consecutive"
      f" x-breakpoints and then y-breakpoints, the z-intervals of the boxes spanning the cell, sorted, chain from 0 to"
      f" 1/2 with no gap (and no overlap: {_zov is None}). Every (x, y, z) in D lies in such a cell and z-interval,"
      f" so in a box", _gap is None and _zov is None)
_rng = random.Random(716_2026)
PTS = []
for _ in range(100000):
    _d = tuple(_rng.randint(1, 10 ** 6) for _ in range(3))
    PTS.append(((_rng.randint(-(_d[0] // 2), _d[0] // 2), _rng.randint(-(_d[1] // 2), _d[1] // 2),
                 _rng.randint(0, _d[2] // 2)), _d))
_t0 = time.time()
OK_P, INFO_P = points_check(ROWS, PTS)
check(f"100000 random rational points of D (each coordinate n/d with d uniform in [1, 10^6]): each lies in some"
      f" box, and |N(s xi - kappa)| < 1 exactly for every box containing it"
      + (f" ({INFO_P['pairs']} point-box pairs, {INFO_P['multi']} points in more than one box; largest value"
         f" {fl(INFO_P['worst'][0])}) [{time.time() - _t0:.1f} s]" if OK_P else f": {INFO_P}"), OK_P)
PTS_DY = []
for _ in range(20000):
    _d = tuple(1 << _rng.randint(1, 8) for _ in range(3))
    PTS_DY.append(((_rng.randint(-(_d[0] // 2), _d[0] // 2), _rng.randint(-(_d[1] // 2), _d[1] // 2),
                    _rng.randint(0, _d[2] // 2)), _d))
OK_PD, INFO_PD = points_check(ROWS, PTS_DY)
check(f"20000 random dyadic points of D (denominators 2^1..2^8, so many on faces, edges and corners of boxes): same"
      + (f" ({INFO_PD['pairs']} point-box pairs, {INFO_PD['multi']} points in more than one box; largest value"
         f" {fl(INFO_PD['worst'][0])})" if OK_PD else f": {INFO_PD}"), OK_PD)

# =================================================================================================================
print("(d) Lemma 2.1(b), Corollary 2.2 and Division 2.3")
HALF_TABLE = {   # the table in the proof of Lemma 2.1(b): x -> (kappa, N(x/2 - kappa))
    (0, 0, 1): ((-1, 1, 0), Fr(1, 2)), (0, 1, 0): ((0, 0, 0), Fr(-3, 4)), (0, 1, 1): ((2, 0, 0), Fr(-1, 4)),
    (1, 0, 0): ((0, 0, 0), Fr(1, 8)), (1, 0, 1): ((0, 1, 1), Fr(-7, 8)), (1, 1, 0): ((0, 1, 0), Fr(3, 8)),
    (1, 1, 1): ((-1, 2, 0), Fr(3, 8))}
check("the table of Lemma 2.1(b) lists the 7 vectors x in {0,1}^3 other than 0; every xi with 2xi in R and xi not in"
      " R is x/2 + r with r in R (take x = 2xi mod 2), and N(x/2 + r - (kappa + r)) = N(x/2 - kappa)",
      set(HALF_TABLE) == set(itertools.product((0, 1), repeat=3)) - {(0, 0, 0)})
for _x, (_k, _val) in HALF_TABLE.items():
    _xi = tuple(Fr(t, 2) for t in _x)
    _n = ncub(rsub(_xi, _k))
    _alt = sum(1 for kk in itertools.product(range(-3, 4), repeat=3) if 0 < abs(ncub(rsub(_xi, kk))) < 1)
    check(f"xi = {tstr(_xi)}, kappa = {_k}: N(xi - kappa) = {_n} (text: {_val}), 0 < |N| < 1"
          f" [{_alt} kappa in [-3,3]^3 work]", _n == _val and 0 < abs(_n) < 1)


def solve(alpha, beta):
    """The rational triple xi with xi * beta = alpha: Cramer's rule for M(beta) xi = alpha (own determinant)."""
    M = mat(beta)
    d = det3(M)
    out = []
    for j in range(3):
        Mj = [list(row) for row in M]
        for i in range(3):
            Mj[i][j] = alpha[i]
        out.append(Fr(det3(Mj), d))
    return tuple(out)


def integral(v):
    return all(Fr(t).denominator == 1 for t in v)


def reduce_D(xi):
    """r = nearest integers; xi' = xi - r, or r - xi when the third coordinate of xi - r is negative."""
    r = tuple(math.floor(t + Fr(1, 2)) for t in xi)
    d = rsub(xi, r)
    return (r, -1, rsc(-1, d)) if d[2] < 0 else (r, 1, d)


def boxes_at(rows, xi):
    return [r for r in rows if all(r['lo'][i] <= 64 * xi[i] <= r['hi'][i] for i in range(3))]


def divide(alpha, beta, rows, half):
    """All divisions (s, kappa, rho = s alpha - kappa beta, how) that Corollary 2.2 / Division 2.3 provide."""
    xi = solve(alpha, beta)
    if rmul(xi, beta) != tuple(Fr(t) for t in alpha):
        raise ValueError("Cramer's rule failed")
    if integral(xi):
        return 'divides', xi, []
    if integral(rsc(2, xi)):
        e = tuple(int(2 * t) % 2 for t in xi)
        r = tuple(int(t - Fr(e[i], 2)) for i, t in enumerate(xi))
        kap = radd(r, half[e][0])
        return 'half', xi, [(1, kap, rsub(alpha, rmul(kap, beta)))]
    r, sg, xid = reduce_D(xi)
    out = []
    for b in boxes_at(rows, xid):
        s = b['s']
        kap = radd(rsc(s, r), rsc(sg, b['kap']))
        out.append((s, kap, rsub(rsc(s, alpha), rmul(kap, beta)), b['path']))
    return 'table', xi, out


EPS, EPSI = (-19, -5, 1), (95, -101, 40)
XS = (Fr(79, 109), Fr(247, 327), Fr(142, 327))
XS_NUM = (237, 247, 142)                         # 327 xi*


def eps_pow(k):
    return rpow(EPS, k) if k >= 0 else rpow(EPSI, -k)


def gen_pairs(seed):
    rg = random.Random(seed)
    out = [(XS_NUM, (327, 0, 0), 'C'), (rsc(-1, XS_NUM), (327, 0, 0), 'C')]
    while len(out) < 300:                        # alpha/beta = +-(EPS^k xi* + r), in +-C
        k, sg = rg.randint(-2, 2), rg.choice((1, -1))
        r = tuple(rg.randint(-3, 3) for _ in range(3))
        gam = tuple(rg.randint(-5, 5) for _ in range(3))
        if gam != ZERO:
            num = radd(rmul(eps_pow(k), XS_NUM), rsc(327, r))
            out.append((rsc(sg, rmul(num, gam)), rmul((327, 0, 0), gam), 'C'))
    while len(out) < 500:                        # alpha/beta = x/2 + r
        bp = tuple(rg.randint(-12, 12) for _ in range(3))
        x = rg.choice(sorted(HALF_TABLE))
        r = tuple(rg.randint(-5, 5) for _ in range(3))
        if bp != ZERO:
            out.append((rmul(bp, radd(rsc(2, r), x)), rsc(2, bp), 'half'))
    while len(out) < 600:                        # beta | alpha
        b = tuple(rg.randint(-15, 15) for _ in range(3))
        if b != ZERO:
            out.append((rmul(b, tuple(rg.randint(-15, 15) for _ in range(3))), b, 'div'))
    while len(out) < 2000:                       # generic
        b = tuple(rg.randint(-25, 25) for _ in range(3))
        if b != ZERO:
            out.append((tuple(rg.randint(-400, 400) for _ in range(3)), b, 'gen'))
    return out


def division_test(pairs, rows, half):
    stats = {'divides': 0, 'half': 0, 's=1': 0, 's=2': 0}
    c_boxes, c_points = set(), set()
    for alpha, beta, kind in pairs:
        how, xi, divs = divide(alpha, beta, rows, half)
        if how == 'divides':
            stats['divides'] += 1
            continue
        if not divs:
            return False, f"alpha = {alpha}, beta = {beta}: the reduced point is in no box", stats, c_boxes
        nb = abs(ncub(beta))
        for d in divs:
            s, kap, rho = d[0], d[1], d[2]
            nr = ncub(rho)
            if not (0 < abs(nr) < nb and nr == ncub(rsub(rsc(s, xi), kap)) * ncub(beta)):
                return False, f"alpha = {alpha}, beta = {beta}, s = {s}, kappa = {kap}: N(rho) = {nr}", stats, c_boxes
            if kind == 'C':
                c_boxes.add((d[3], s))
        if kind == 'C':
            c_points.add(reduce_D(xi)[2])
        stats['half' if how == 'half' else f"s={divs[0][0]}"] += 1
    stats['C points'] = c_points
    return True, "", stats, c_boxes


PAIRS = gen_pairs(20261015)
_t0 = time.time()
OK_D, MSG_D, ST_D, CB_D = division_test(PAIRS, ROWS, HALF_TABLE)
check(f"{len(PAIRS)} pairs (alpha, beta), xi = alpha/beta by Cramer's rule: 0 < |N(s alpha - kappa beta)| < |N(beta)|"
      f" and N(s alpha - kappa beta) = N(s xi - kappa) N(beta), for every box containing the reduced point"
      f" (beta | alpha: {ST_D['divides']}; Lemma 2.1(b): {ST_D['half']}; table with s = 1: {ST_D['s=1']}, with"
      f" s = 2: {ST_D['s=2']}) [{time.time() - _t0:.1f} s]" + ("" if OK_D else f": {MSG_D}"), OK_D)
_xr = reduce_D(XS)[2]
check(f"the 300 pairs with alpha/beta in xi* + R or -xi* + R (EPS^k xi* + r, k = -2..2, and the pair 237 + 247l +"
      f" 142l^2 over 327 itself) all reduce to xi' = {tstr(_xr)}; the boxes containing it are {sorted(CB_D)}, all"
      f" with s = 2, as they must be by Remark 2.5", bool(CB_D) and all(s == 2 for _, s in CB_D)
      and ST_D['C points'] == {_xr} and reduce_D(rsc(-1, XS))[2] == _xr)

# =================================================================================================================
print("(e) Remark 2.5: R is not norm-Euclidean (a normalisation by EPS powers, not 05's T2 descent)")
check("EPS = (-19,-5,1) and EPS^-1 = (95,-101,40): EPS * EPS^-1 = 1, N(EPS) = N(EPS^-1) = 1",
      rmul(EPS, EPSI) == ONE and ncub(EPS) == 1 and ncub(EPSI) == 1)
_dp, _dm = rsub(rmul(EPS, XS), XS), rsub(rmul(EPSI, XS), XS)
check(f"EPS xi* - xi* = {tstr(_dp)} (text: (-6, -27, -10)) and EPS^-1 xi* - xi* = {tstr(_dm)} are in R; with EPS R ="
      f" R this gives EPS C = C for C = xi* + R", _dp == (-6, -27, -10) and integral(_dm))
check(f"N(xi* - 2) = {ncub(rsub(XS, (2, 0, 0)))} (text: -121/109), and 327 xi* = {tstr(rsc(327, XS))}",
      ncub(rsub(XS, (2, 0, 0))) == Fr(-121, 109) and rsc(327, XS) == XS_NUM)
E_LO, E_HI = lam_iv(*EPSI)
e_lo, e_hi = lam_iv(*EPS)
check(f"E := sigma1(EPS^-1) = 95 - 101 l1 + 40 l1^2 in [{fl(E_LO, 15)}, {fl(E_HI, 15)}], inside (605.08, 605.09);"
      f" sigma1(EPS) in (1/606, 1/605); sigma1(EPS) sigma1(EPS^-1) = sigma1(1) = 1, so sigma1(EPS) = 1/E",
      Fr(60508, 100) < E_LO <= E_HI < Fr(60509, 100) and Fr(1, 606) < e_lo <= e_hi < Fr(1, 605))


def sigma_identity():
    for v in BASIS:
        a, b, c = v
        w = (Fr(a + 4 * c), Fr(-b, 2), Fr(-c, 2))
        u = (b, -c, 0)
        rhs = radd(radd(w, rmul((0, Fr(3, 2), 0), u)), rmul((-4, 0, 3), (c, 0, 0)))
        if rhs != v:
            return False
    return True


check("sigma1 = w + (3L/2) u + (3L^2 - 4) c with w = a - (L/2) b + ((8 - L^2)/2) c and u = b - L c (linear in"
      " (a,b,c): checked on the basis); with the sum of squares of (b), Q = w^2 + ((3L^2 - 16)/4) u^2",
      sigma_identity())
note("Normalisation. Let y in C with |N(y)| < 121/109. y != 0 (xi* is not in R), so sigma1(y) != 0 (Section 1).",
     "E > 1, so some k in Z has E^k <= |sigma1(y)| < E^(k+1). y' = EPS^k y lies in C (EPS^+-1 C = C), N(y') = N(y),",
     "sigma1(y') = sigma1(EPS)^k sigma1(y) = E^-k sigma1(y), so 1 <= |sigma1(y')| < E; and Q(y') = |N(y')|/|sigma1(y')|",
     "<= |N(y)| < 121/109. So it suffices to list every y in C with |sigma1(y)| < E and Q(y) <= 121/109 and to show",
     "|N| >= 121/109 on that list.",
     "Enumeration. Q <= T gives |w| <= sqrt(T) and |u| <= sqrt(4T/(3 l1^2 - 16)); then",
     "|c| <= (E + |w| + (3|l1|/2)|u|)/(3 l1^2 - 4), b in l1 c + [-|u|, |u|], a in (l1/2) b - ((8 - l1^2)/2) c + [-|w|, |w|],",
     "all with rational outer bounds; candidates are kept exactly (a superset of the set to be listed).")
FP_LO = 3 * LAM_HI ** 2 - 4
ABS_L_HI = -LAM_LO


def sqrt_up(x, den=10 ** 6):
    k = math.isqrt(math.floor(x * den * den))
    while Fr(k, den) ** 2 < x:
        k += 1
    return Fr(k, den)


def coset_candidates(zeta, T, E_up):
    """A list containing every y in zeta + Z^3 with |sigma1(y)| <= E_up and Q(y) <= T."""
    rw, ru = sqrt_up(T), sqrt_up(4 * T / K_LO)
    Y2 = (E_up + rw + Fr(3, 2) * ABS_L_HI * ru) / FP_LO
    out = []
    for n2 in range(math.floor(-Y2 - zeta[2]), math.ceil(Y2 - zeta[2]) + 1):
        y2 = zeta[2] + n2
        if abs(y2) > Y2:
            continue
        lo1 = min(LAM_LO * y2, LAM_HI * y2) - ru
        hi1 = max(LAM_LO * y2, LAM_HI * y2) + ru
        for n1 in range(math.floor(lo1 - zeta[1]), math.ceil(hi1 - zeta[1]) + 1):
            y1 = zeta[1] + n1
            if not lo1 <= y1 <= hi1:
                continue
            c_lo, c_hi = lam_iv(-4 * y2, y1 / 2, y2 / 2)
            for n0 in range(math.floor(c_lo - rw - zeta[0]), math.ceil(c_hi + rw - zeta[0]) + 1):
                y0 = zeta[0] + n0
                if c_lo - rw <= y0 <= c_hi + rw:
                    out.append((y0, y1, y2))
    return out, (rw, ru, Y2)


T_XS = Fr(121, 109)
_t0 = time.time()
CANDS, (_rw, _ru, _Y2) = coset_candidates(XS, T_XS, E_HI)
_norms = sorted((abs(ncub(y)), y) for y in CANDS)
_mins = [y for n, y in _norms if n == _norms[0][0]]
check(f"with T = 121/109: |w| <= {_rw}, |u| <= {_ru}, |c| <= {fl(_Y2)}; {len(CANDS)} candidates, all in C; least |N|"
      f" among them = {_norms[0][0]} (no candidate below 121/109), attained only at {[tstr(y) for y in _mins]};"
      f" the next values {', '.join(str(n) for n, _ in _norms[1:6])} [{time.time() - _t0:.1f} s]",
      all(integral(rsub(y, XS)) for y in CANDS) and _norms[0][0] == T_XS and len(_mins) == 1)
_ks = [k for k in range(-4, 5) if rmul(eps_pow(k), rsub(XS, (2, 0, 0))) == _mins[0]]
_s2 = lam_iv(*rsub(XS, (2, 0, 0)))
check(f"that element is EPS^k (xi* - 2) with k = {_ks}; xi* - 2 itself has sigma1 in [{fl(_s2[0])}, {fl(_s2[1])}],"
      f" |sigma1| < 1, so it is not in the normalised range [1, E) and the normalised orbit point is this one",
      len(_ks) == 1 and max(abs(_s2[0]), abs(_s2[1])) < 1)
note("Conclusion: no y in C has |N(y)| < 121/109, and N(xi* - 2) = -121/109. So min over C of |N| = 121/109; no kappa",
     "in R has |N(xi* - kappa)| < 1, and 327 xi* = 237 + 247l + 142l^2 cannot be divided by 327 with a smaller",
     "remainder. As N(-y) = -N(y), the same holds on -C. Moreover, since the list has a single element of norm",
     "121/109 and every y in C with |N(y)| = 121/109 normalises into the list (Q(y') <= 121/109), the elements of C",
     "with |N| = 121/109 are exactly the EPS^k (xi* - 2), k in Z.")
# A cross-check with floats (a guide only): scan a large box of C and keep the points that floats put clearly inside.
_lam = float(LAM_LO)
_E = float(E_LO)
_cset = {tuple(int(y[i] - XS[i]) for i in range(3)) for y in CANDS}
_found, _edge = [], False
_N1, _N0 = np.meshgrid(np.arange(-130, 131), np.arange(-130, 131), indexing='ij')
for _n2 in range(-50, 51):
    _y2 = 142 / 327 + _n2
    _y1 = 247 / 327 + _N1
    _y0 = 79 / 109 + _N0
    _w = _y0 - _lam / 2 * _y1 + (8 - _lam ** 2) / 2 * _y2
    _u = _y1 - _lam * _y2
    _Qf = _w ** 2 + (3 * _lam ** 2 - 16) / 4 * _u ** 2
    _sf = _y0 + _lam * _y1 + _lam ** 2 * _y2
    _mask = (np.abs(_sf) < _E * (1 - 1e-7)) & (_Qf < float(T_XS) * (1 - 1e-7))
    for _i, _j in zip(*np.nonzero(_mask)):
        _found.append((int(_N0[_i, _j]), int(_N1[_i, _j]), _n2))
        _edge |= abs(int(_N0[_i, _j])) == 130 or abs(int(_N1[_i, _j])) == 130 or abs(_n2) == 50
check(f"cross-check with floats (a guide, not part of the proof): scanning |a|, |b| <= 130, |c| <= 50 in C finds"
      f" {len(_found)} points with |sigma1| < E and Q < 121/109 (none on the edge of the scan); all are candidates",
      not _edge and set(_found) <= _cset)

print("  (e2) the 654 classes of fixed points of EPS mod R (the last paragraph of Remark 2.5)")
EM1 = rsub(EPS, ONE)


def hnf_lower(M):
    """Unimodular integer column operations bring M to lower-triangular form with positive diagonal."""
    A = [list(row) for row in M]
    for i in range(3):
        for k in range(i + 1, 3):
            while A[i][k] != 0:
                qq = A[i][i] // A[i][k]
                for rr in range(3):
                    A[rr][i] -= qq * A[rr][k]
                for rr in range(3):
                    A[rr][i], A[rr][k] = A[rr][k], A[rr][i]
        if A[i][i] < 0:
            for rr in range(3):
                A[rr][i] = -A[rr][i]
    return A


def frac(v):
    return tuple(t - math.floor(t) for t in v)


_H = hnf_lower(mat(EM1))
_dg = (_H[0][0], _H[1][1], _H[2][2])
CLASSES = {frac(solve(r, EM1)) for r in itertools.product(*(range(d) for d in _dg))}
check(f"EPS - 1 = {EM1}, N(EPS - 1) = {ncub(EM1)}. A Hermite form of M(EPS - 1) (lower triangular, diagonal {_dg}) gives"
      f" representatives r of Z^3 / M(EPS - 1) Z^3; the fixed classes are (EPS - 1)^-1 r mod R: {len(CLASSES)} distinct"
      f" classes, among them those of xi*, -xi* and 2 xi*",
      _H[0][1] == _H[0][2] == _H[1][2] == 0 and _dg[0] * _dg[1] * _dg[2] == abs(ncub(EM1)) == 654
      and len(CLASSES) == 654 and {frac(XS), frac(rsc(-1, XS)), frac(rsc(2, XS))} <= CLASSES)
_t0 = time.time()
CLASS_MIN = []
for _z in sorted(CLASSES):
    if _z == ZERO:
        continue
    _cl, _ = coset_candidates(_z, Fr(1), E_HI)
    CLASS_MIN.append((min((abs(ncub(y)) for y in _cl), default=Fr(10 ** 9)), _z, len(_cl)))
CLASS_MIN.sort(reverse=True)
_big = [m for m in CLASS_MIN if m[0] >= 1]
_nxt = [m for m in CLASS_MIN if m[0] < 1]
check(f"for each of the 653 classes z != 0, the candidates with T = 1 ({min(m[2] for m in CLASS_MIN)} to"
      f" {max(m[2] for m in CLASS_MIN)} each): exactly {len(_big)} classes have no element with |N| < 1, those of xi*"
      f" and -xi*; for the other 651 the least |N| on the list is the exact minimum of |N| over the class (an element"
      f" with |N| = m < 1 normalises into the list), and the largest of these minima is {_nxt[0][0]}, for"
      f" {sum(1 for m in _nxt if m[0] == _nxt[0][0])} classes (text: 99/109) [{time.time() - _t0:.1f} s]",
      len(_big) == 2 and {m[1] for m in _big} == {frac(XS), frac(rsc(-1, XS))} and _nxt[0][0] == Fr(99, 109))
for _m in CLASS_MIN[:6]:
    note(f"min |N| = {str(_m[0]):>8s} = {fl(_m[0], 6)}   class of {tstr(_m[1])}")
_m2 = [m for m in CLASS_MIN if m[1] == frac(rsc(2, XS))][0]
check(f"the class of 2 xi* (used by the table near +-xi*) has minimum |N| = {_m2[0]} < 1", _m2[0] < 1)

# =================================================================================================================
print("(f) negative controls: planted errors that the checks must catch")


def run_all(rows, pts):
    res = {'table': (not [r for r in rows if row_problem(r)], None)}
    res['bound'] = bound_all(rows)
    res['volume'] = (volume_units(rows) == VOL_D, volume_units(rows))
    res['disjoint'] = (first_overlap(rows) is None, None)
    res['sweep'] = (sweep(rows)[0] is None, None)
    res['points'] = points_check(rows, pts)
    return res


def control(name, rows, must_fail):
    res = run_all(rows, PTS)
    failed = [k for k, v in res.items() if not v[0]]
    detail = f"; bound check: {res['bound'][1]}" if 'bound' in failed else ""
    check(f"control '{name}': the checker fails (failing: {', '.join(failed) or 'none'}; required: "
          f"{', '.join(must_fail)}){detail}", bool(failed) and set(must_fail) <= set(failed))


def modified(rows, i, **kw):
    out = [dict(r) for r in rows]
    out[i].update(kw)
    return out


_iw = [r['path'] for r in ROWS].index('xyxyYXzxy')
control(f"kappa of box xyxyYXzxy (the box of the maximum) changed from {ROWS[_iw]['kap']} to (0, 0, 0)",
        modified(ROWS, _iw, kap=(0, 0, 0)), ['bound'])
_kd = []
for _i in (_iw, 0, 199, len(ROWS) - 1):
    for _ax in range(3):
        for _sg in (1, -1):
            _k = list(ROWS[_i]['kap'])
            _k[_ax] += _sg
            _kd.append(not certify_box(dict(ROWS[_i], kap=tuple(_k)))[0])
check(f"control: each of the 6 unit changes of kappa on 4 boxes (xyxyYXzxy, the first, the 200th, the last) is"
      f" caught by the bound check: {sum(_kd)} of {len(_kd)}", all(_kd))
_r0 = ROWS[0]
control(f"a corner of box {_r0['path']} moved: hi_a {_r0['hi'][0]} -> {_r0['hi'][0] - 1} (leaves a gap)",
        modified(ROWS, 0, hi=(_r0['hi'][0] - 1,) + _r0['hi'][1:]), ['volume', 'sweep', 'points'])
control(f"a corner of box {_r0['path']} moved: hi_a {_r0['hi'][0]} -> {_r0['hi'][0] + 1} (overlaps a neighbour)",
        modified(ROWS, 0, hi=(_r0['hi'][0] + 1,) + _r0['hi'][1:]), ['volume', 'disjoint'])
control(f"box {_r0['path']} shifted by +1/64 along a (same volume)",
        modified(ROWS, 0, lo=(_r0['lo'][0] + 1,) + _r0['lo'][1:], hi=(_r0['hi'][0] + 1,) + _r0['hi'][1:]),
        ['disjoint', 'sweep'])
_ibig = max(range(len(ROWS)), key=lambda i: volume_units([ROWS[i]]))
_ismall = min(range(len(ROWS)), key=lambda i: volume_units([ROWS[i]]))
control(f"the largest box {ROWS[_ibig]['path']} dropped", ROWS[:_ibig] + ROWS[_ibig + 1:],
        ['volume', 'sweep', 'points'])
control(f"the smallest box {ROWS[_ismall]['path']} dropped (volume {Fr(volume_units([ROWS[_ismall]]), 1 << 18)})",
        ROWS[:_ismall] + ROWS[_ismall + 1:], ['volume', 'sweep'])
_qbad = dict(QCOEF)
_qbad[(0, 1, 1)] = (5, 0, 0)
check("control: the identity check N = sigma1 * Q fails when the coefficient 6 of bc in Q is changed to 5",
      not nq_identity(_qbad))
_hbad = dict(HALF_TABLE)
_hbad[(1, 0, 1)] = ((0, 0, 0), None)
_okb = division_test([p for p in PAIRS if p[2] == 'half'], ROWS, _hbad)[0]
check(f"control: with kappa for x = (1,0,1) in Lemma 2.1(b) changed to (0,0,0) (N(1/2, 0, 1/2) ="
      f" {ncub((Fr(1, 2), 0, Fr(1, 2)))}), the division test fails", not _okb)
_c2, _ = coset_candidates(frac(rsc(2, XS)), T_XS, E_HI)
_m2b = min(abs(ncub(y)) for y in _c2)
check(f"control: the same enumeration on the EPS-stable class of 2 xi* does find small norms (least |N| = {_m2b} < 1),"
      f" so the bound 121/109 on C is a property of C, not of the method", _m2b < 1)

print(f"SUMMARY: all {PASSED[0]} checks passed ({time.time() - T_START:.1f} s). Section 2 of hypothesis_H.md confirmed"
      f" independently: (a) the table is well formed; (b) Lemma 2.1(a) holds on every box by the factored bound"
      f" |N| <= max|sigma1| * max Q ({AGG['split']} of {len(ROWS)} boxes split, {AGG['leaves']} sub-boxes), and the"
      f" largest value is 4093/4096 to within 2^-20; (c) the boxes cover D (volume + disjointness, and a sweep);"
      f" (d) Lemma 2.1(b) and 2000 divisions; (e) min over xi* + R of |N| = 121/109, so R is not norm-Euclidean,"
      f" and of the 654 fixed classes only +-xi* have no element of norm below 1 (next 99/109); (f) every planted"
      f" error was caught.")
