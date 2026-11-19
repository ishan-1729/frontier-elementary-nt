"""verify_coverage.py -- independent re-check of 01_coverage.py (every element of norm c lies in the unit orbit of a
listed seed), for the three quartic fields of 1611.a1 and c in {1, 4, 9, 36}.

Independence.  Nothing is imported from 01_coverage.py, and nothing from core.py except, at the very end, its data
tables CLASSES / UNITS / KNOWN, only to compare them with the copies retyped below.  Everything else is re-done here
by other methods:
  * field arithmetic: polynomials modulo f in Fractions, back-substitution in the (triangular) integral basis,
    structure constants re-derived and checked integral; quotients by an exact Gaussian solve;
  * the norm form N(x0,..,x3) = Res_t(f(t), x(t)) computed symbolically (sympy), checked to have integer
    coefficients and to satisfy N(n - m alpha) = F(m,n) identically;
  * real embeddings: sympy's exact real-root isolation refined to width < 2^-220, each interval re-checked by an
    exact sign change of f, and the number of real roots by sympy's Sturm count (must be 2);
  * the complex embedding: a Rouche certificate with exact rational bounds for the disk |z - z0| < 2^-300 around a
    200-digit approximation z0 with Im z0 > 2^-300 (exactly one root there), cross-checked against sympy's exact
    complex isolating rectangle;
  * own outward-rounded interval arithmetic (Fraction endpoints on a 2^-360 grid); s3 = Re^2 + Im^2;
  * the polygon P+ by exact half-plane clipping of the simplex, cross-checked against the basic feasible
    solutions of all pairs of constraints;
  * the WEAKER linear bound of the lemma: kappa_v = min over the vertices of r_v and
        E(w) <= T'  with T' the least integer such that T'^4 >= c^2 / (kappa_1 kappa_2 kappa_3^2);
    (b_min and the product bound T of 01 are recomputed only for comparison);
  * a rational lower form A with x^T A x <= E(x) for all real x, from a Gershgorin-type diagonal shift
    (not 4 rho I), and an exact, float-free Fincke-Pohst enumeration in integers of all x with x^T A x <= T'_max;
  * orbits: the norm-one elements found are identified exactly as +-eps1^a eps2^b (table lookup); by the lemma
    with c = 1 this proves that every unit of norm +1 of O lies in <-1, eps1, eps2>.  Then two elements of norm c
    lie in one orbit iff their quotient lies in O (exact solve), which gives an exact orbit decomposition.
A failed check prints CHECK FAILED and the script exits with status 1.
"""
import sys, time, math, random
from fractions import Fraction as Fr
from math import isqrt, gcd
import sympy
from sympy import Poly, symbols, Rational, resultant, expand
import mpmath

T0 = time.time()
FAILS = []


def check(cond, msg):
    if cond:
        print(f'   [ok] {msg}')
    else:
        print(f'   CHECK FAILED: {msg}')
        FAILS.append(msg)
    return cond


def lcm(a, b):
    return a // gcd(a, b) * b


# ======================================================================================================================
# Data, retyped (compared with core.py's tables at the end)
# ======================================================================================================================
DATA = {
    '10': dict(label='(1,0)', F=(36, 0, 0, -6, 1),
               BASIS=[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, Fr(1, 6), 0], [0, 0, 0, Fr(1, 6)]],
               UNITS=((0, 1, -1, 0), (-86, -41, -117, 30)),
               KNOWN=[(0, 1), (1, 0), (1, 2), (1, 6)]),
    '01': dict(label='(0,1)', F=(-5, 14, 24, 2, 1),
               BASIS=[[1, 0, 0, 0], [0, 1, 0, 0], [Fr(1, 6), Fr(2, 3), Fr(1, 6), 0], [Fr(1, 3), Fr(1, 2), 0, Fr(1, 6)]],
               UNITS=((0, -1, 1, 0), (-21, 36, 16, 14)),
               KNOWN=[(-1, 1), (0, 1), (1, 1), (4, 1)]),
    '11': dict(label='(1,1)', F=(-8, 4, 24, 4, 1),
               BASIS=[[1, 0, 0, 0], [0, 1, 0, 0], [Fr(2, 3), Fr(1, 3), Fr(1, 6), 0], [Fr(1, 3), 0, 0, Fr(1, 12)]],
               UNITS=((-1, 0, 1, 0), (1, -3, -5, -2)),
               KNOWN=[(-1, 1), (0, 1), (2, 1)]),
}
RHS = (1, 4, 9, 36)
EXPONENTS = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)]

# What 01_coverage.out.txt printed (retyped from that file, for comparison only).
OUT01 = {
    '10': dict(T={1: 34, 4: 68, 9: 102, 36: 203}, pts={1: 74, 4: 278, 9: 638, 36: 2536},
               count={1: 8, 4: 4, 9: 2, 36: 6}, norbits={1: 1, 4: 1, 9: 1, 36: 1},
               seeds={1: [(-1, 0, 0, 0)], 4: [(2, -1, 0, 0)], 9: [(-3, 0, 6, -1)], 36: [(0, 0, 6, -1)]},
               verts=[(0.000041, 0.024017, 0.975942), (0.977009, 0.014872, 0.008120),
                      (0.000022, 0.763916, 0.236062), (0.525471, 0.472772, 0.001758)],
               bmin=7.676147e-07,
               f={(-1, -1): (1.19979e-05, 1688.14, 7.02656), (-1, 0): (0.537363, 31.8074, 0.241881),
                  (-1, 1): (24067.4, 0.599304, 0.00832649), (0, -1): (2.23274e-05, 53.0739, 29.0496),
                  (0, 1): (44788, 0.0188417, 0.0344239), (1, -1): (4.155e-05, 1.6686, 120.099),
                  (1, 0): (1.86094, 0.0314392, 4.13426), (1, 1): (83347.8, 0.000592368, 0.142317)}),
    '01': dict(T={1: 27, 4: 54, 9: 80, 36: 160}, pts={1: 44, 4: 168, 9: 388, 36: 1556},
               count={1: 4, 4: 2, 9: 4, 36: 4}, norbits={1: 1, 4: 1, 9: 1, 36: 1},
               seeds={1: [(-1, 0, 0, 0)], 4: [(-1, -1, 0, 0)], 9: [(-3, -2, -3, -1)], 36: [(1, -1, 0, 0)]},
               verts=[(0.316982, 0.000095, 0.682923), (0.861497, 0.000122, 0.138381),
                      (0.000244, 0.008369, 0.991387), (0.952237, 0.010693, 0.037071),
                      (0.000270, 0.734150, 0.265579), (0.000734, 0.945451, 0.053815)],
               bmin=2.010792e-06,
               f={(-1, -1): (0.000852461, 7756.77, 0.388886), (-1, 0): (3.00407, 112.973, 0.0542821),
                  (-1, 1): (10586.3, 1.6454, 0.0075769), (0, -1): (0.000283769, 68.6601, 7.16416),
                  (0, 1): (3524, 0.0145645, 0.139584), (1, -1): (9.44614e-05, 0.607755, 131.98),
                  (1, 0): (0.332882, 0.00885164, 18.4223), (1, 1): (1173.07, 0.00012892, 2.57145)}),
    '11': dict(T={1: 13, 4: 26, 9: 39, 36: 77}, pts={1: 18, 4: 86, 9: 180, 36: 730},
               count={1: 2, 4: 12, 9: 4, 36: 8}, norbits={1: 1, 4: 3, 9: 1, 36: 3},
               seeds={1: [(-1, 0, 0, 0)], 4: [(3, -1, -2, -1), (2, -1, -1, -1), (0, 0, -1, 0)],
                      9: [(-1, -1, 0, 0)], 36: [(6, -2, -3, -2), (-1, -2, -3, -1), (0, -1, -3, -1)]},
               verts=[(0.217232, 0.000505, 0.782262), (0.002379, 0.016638, 0.980984),
                      (0.683424, 0.000558, 0.316018), (0.003205, 0.838715, 0.158080),
                      (0.920947, 0.028128, 0.050924), (0.010084, 0.926055, 0.063861)],
               bmin=3.808329e-05,
               f={(-1, -1): (1218.1, 1.86667, 0.0209712), (-1, 0): (4.23946, 55.6599, 0.0650989),
                  (-1, 1): (0.0147549, 1659.65, 0.20208), (0, -1): (287.325, 0.0335371, 0.322144),
                  (0, 1): (0.00348038, 29.8177, 3.1042), (1, -1): (67.774, 0.000602537, 4.94853),
                  (1, 0): (0.235879, 0.0179663, 15.3612), (1, 1): (0.000820948, 0.535712, 47.6844)}),
}


def Fform(key, m, n):
    a0, a1, a2, a3, a4 = DATA[key]['F']
    return a0 * m ** 4 + a1 * m ** 3 * n + a2 * m ** 2 * n ** 2 + a3 * m * n ** 3 + a4 * n ** 4


# ======================================================================================================================
# Field arithmetic (own)
# ======================================================================================================================
class Order:
    """Q(alpha), f(alpha) = 0, f(t) = F(1,t) monic; O = Z-span of the rows of BASIS (power coordinates)."""

    def __init__(self, key):
        d = DATA[key]
        a0, a1, a2, a3, a4 = d['F']
        assert a4 == 1
        self.key, self.label = key, d['label']
        self.c = [Fr(a0), Fr(a1), Fr(a2), Fr(a3)]         # f = t^4 + c3 t^3 + c2 t^2 + c1 t + c0
        self.B = [[Fr(v) for v in row] for row in d['BASIS']]
        for i in range(4):
            assert self.B[i][i] != 0 and all(self.B[i][j] == 0 for j in range(i + 1, 4)), 'basis not triangular'
        self.st = {}
        for i in range(4):
            for j in range(4):
                co = self.coords(self.pmul(self.B[i], self.B[j]))
                if any(v.denominator != 1 for v in co):
                    raise ValueError('BASIS is not closed under multiplication')
                self.st[(i, j)] = tuple(int(v) for v in co)

    def pmul(self, p, q):
        """Product of power-coordinate vectors modulo f."""
        r = [Fr(0)] * 7
        for i in range(4):
            if p[i]:
                for j in range(4):
                    if q[j]:
                        r[i + j] += p[i] * q[j]
        c0, c1, c2, c3 = self.c
        for d in (6, 5, 4):                                   # t^d = -t^(d-4) (c3 t^3 + c2 t^2 + c1 t + c0)
            v = r[d]
            if v:
                r[d] = Fr(0)
                r[d - 1] -= v * c3
                r[d - 2] -= v * c2
                r[d - 3] -= v * c1
                r[d - 4] -= v * c0
        return r[:4]

    def power(self, x):
        return [sum(Fr(x[i]) * self.B[i][j] for i in range(4)) for j in range(4)]

    def coords(self, p):
        """Integral-basis coordinates of a power-coordinate vector (back-substitution, B lower triangular)."""
        x = [Fr(0)] * 4
        for j in (3, 2, 1, 0):
            x[j] = (Fr(p[j]) - sum(x[i] * self.B[i][j] for i in range(j + 1, 4))) / self.B[j][j]
        return x

    def mul(self, x, y):
        z = [0, 0, 0, 0]
        for i in range(4):
            if x[i]:
                for j in range(4):
                    if y[j]:
                        c = x[i] * y[j]
                        s = self.st[(i, j)]
                        z[0] += c * s[0]; z[1] += c * s[1]; z[2] += c * s[2]; z[3] += c * s[3]
        return tuple(z)

    def mulmatrix(self, x):
        """Columns = coordinates of x * omega_j."""
        cols = [self.mul(x, tuple(int(i == j) for i in range(4))) for j in range(4)]
        return [[cols[j][i] for j in range(4)] for i in range(4)]

    def quotient(self, x, y):
        """y / x as integral-basis coordinates (Fractions) by an exact Gaussian solve of M(x) q = y."""
        M = [[Fr(v) for v in row] + [Fr(y[i])] for i, row in enumerate(self.mulmatrix(x))]
        n = 4
        for col in range(n):
            piv = next(r for r in range(col, n) if M[r][col] != 0)
            M[col], M[piv] = M[piv], M[col]
            pv = M[col][col]
            M[col] = [v / pv for v in M[col]]
            for r in range(n):
                if r != col and M[r][col] != 0:
                    fct = M[r][col]
                    M[r] = [a - fct * b for a, b in zip(M[r], M[col])]
        return [M[i][n] for i in range(n)]

    def quotient_in_O(self, x, y):
        q = self.quotient(x, y)
        if all(v.denominator == 1 for v in q):
            return tuple(int(v) for v in q)
        return None


# ======================================================================================================================
# Interval arithmetic (own), Fraction endpoints rounded outward to a 2^-360 grid
# ======================================================================================================================
GRID = 2 ** 360


def _dn(q):
    return Fr(math.floor(q * GRID), GRID)


def _up(q):
    return Fr(math.ceil(q * GRID), GRID)


class I:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        lo = Fr(lo)
        hi = lo if hi is None else Fr(hi)
        if lo > hi:
            raise ValueError('empty interval')
        self.lo, self.hi = _dn(lo), _up(hi)

    def __add__(self, o):
        if not isinstance(o, I):
            o = Fr(o)
            return I(self.lo + o, self.hi + o)
        return I(self.lo + o.lo, self.hi + o.hi)
    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + (-o if isinstance(o, I) else -Fr(o))

    def __mul__(self, o):
        if not isinstance(o, I):
            o = Fr(o)
            a, b = self.lo * o, self.hi * o
            return I(min(a, b), max(a, b))
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return I(min(p), max(p))
    __rmul__ = __mul__

    def sq(self):
        if self.lo >= 0:
            return I(self.lo * self.lo, self.hi * self.hi)
        if self.hi <= 0:
            return I(self.hi * self.hi, self.lo * self.lo)
        return I(0, max(self.lo * self.lo, self.hi * self.hi))

    def contains(self, q):
        return self.lo <= q <= self.hi

    def mid(self):
        return (self.lo + self.hi) / 2

    def __float__(self):
        return float(self.mid())


def isum(terms):
    acc = I(0)
    for t in terms:
        acc = acc + t
    return acc


# ======================================================================================================================
# Rigorous embeddings (own)
# ======================================================================================================================
def sqrt_up(q, bits=480):
    """A rational >= sqrt(q) (q >= 0 rational), within 2^-bits."""
    q = Fr(q)
    s = isqrt(-((-q.numerator * 4 ** bits) // q.denominator))     # isqrt(ceil(q 4^bits))
    return Fr(s + 1, 2 ** bits)


def sqrt_dn(q, bits=480):
    q = Fr(q)
    s = isqrt((q.numerator * 4 ** bits) // q.denominator)          # isqrt(floor(q 4^bits))
    return Fr(s, 2 ** bits)


class Emb:
    """sigma_1 < sigma_2 (real roots of f) and sigma_3 (the root with Im > 0) of the order K, rigorously."""

    def __init__(self, K):
        self.K = K
        c0, c1, c2, c3 = K.c
        t = symbols('t')
        fs = Poly(t ** 4 + Rational(c3) * t ** 3 + Rational(c2) * t ** 2 + Rational(c1) * t + Rational(c0), t)
        fval = lambda q: q ** 4 + c3 * q ** 3 + c2 * q ** 2 + c1 * q + c0
        self.nreal = fs.count_roots()                          # Sturm, exact
        self.disc = sympy.discriminant(fs)
        real, cplx = fs.intervals(all=True, eps=Rational(1, 2 ** 40))     # refined isolating boxes
        self.real = []
        for (lo, hi), mult in real:
            lo, hi = fs.refine_root(lo, hi, eps=Rational(1, 2 ** 220))
            lo, hi = Fr(int(lo.p), int(lo.q)), Fr(int(hi.p), int(hi.q))
            self.real.append((lo, hi, fval(lo) * fval(hi) < 0, mult))
        self.real.sort()
        # complex root: Rouche disk around a 200-digit approximation
        mpmath.mp.dps = 200
        rts = mpmath.polyroots([1, int(c3), int(c2), int(c1), int(c0)], maxsteps=1000, extraprec=1200)
        zc = [z for z in rts if mpmath.im(z) > mpmath.mpf(10) ** -50]
        assert len(zc) == 1
        x0, y0 = Fr(mpmath.nstr(mpmath.re(zc[0]), 200)), Fr(mpmath.nstr(mpmath.im(zc[0]), 200))
        self.z0, self.rho = (x0, y0), Fr(1, 2 ** 300)
        self.rouche = self._rouche(x0, y0, self.rho)
        # sympy's exact isolating rectangle (corners a = lower left, b = upper right) lying in the closed upper
        # half-plane: must contain z0
        up_rects = []
        for ((a, b), mult) in cplx:
            ar, ai = sympy.re(a), sympy.im(a)
            br, bi = sympy.re(b), sympy.im(b)
            if bi > 0 and ai >= 0:
                up_rects.append((Fr(int(sympy.Rational(ar).p), int(sympy.Rational(ar).q)),
                                 Fr(int(sympy.Rational(ai).p), int(sympy.Rational(ai).q)),
                                 Fr(int(sympy.Rational(br).p), int(sympy.Rational(br).q)),
                                 Fr(int(sympy.Rational(bi).p), int(sympy.Rational(bi).q))))
        self.sympy_rect_ok = (len(up_rects) == 1 and up_rects[0][0] <= x0 <= up_rects[0][2]
                              and up_rects[0][1] <= y0 <= up_rects[0][3])
        # the enclosures
        (l1, h1, _, _), (l2, h2, _, _) = self.real[0], self.real[1]
        self.r1, self.r2 = I(l1, h1), I(l2, h2)
        self.zr, self.zi = I(x0 - self.rho, x0 + self.rho), I(y0 - self.rho, y0 + self.rho)
        # values of the basis elements
        self.V1, self.V2, self.R3, self.J3 = [], [], [], []
        for i in range(4):
            row = K.B[i]
            self.V1.append(self._real_horner(row, self.r1))
            self.V2.append(self._real_horner(row, self.r2))
            re, im = self._cx_horner(row, self.zr, self.zi)
            self.R3.append(re)
            self.J3.append(im)

    @staticmethod
    def _real_horner(row, r):
        v = I(row[3])
        for j in (2, 1, 0):
            v = v * r + row[j]
        return v

    @staticmethod
    def _cx_horner(row, zr, zi):
        vr, vi = I(row[3]), I(0)
        for j in (2, 1, 0):
            vr, vi = vr * zr - vi * zi + row[j], vr * zi + vi * zr
        return vr, vi

    def _rouche(self, x0, y0, rho):
        """|f(z0)| + |f''(z0)/2| rho^2 + |f'''(z0)/6| rho^3 + rho^4 < |f'(z0)| rho  (exact rational bounds),
        and Im z0 > rho.  Then f has exactly one root in the open disk |z - z0| < rho (Rouche against
        f'(z0)(z - z0)), and it is non-real with Im > 0."""
        c0, c1, c2, c3 = self.K.c

        def cmul(a, b):
            return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

        def cadd(a, b):
            return (a[0] + b[0], a[1] + b[1])

        def cpoly(coeffs, z):            # coeffs high -> low
            v = (Fr(coeffs[0]), Fr(0))
            for c in coeffs[1:]:
                v = cadd(cmul(v, z), (Fr(c), Fr(0)))
            return v
        z = (x0, y0)
        f0 = cpoly([1, c3, c2, c1, c0], z)
        f1 = cpoly([4, 3 * c3, 2 * c2, c1], z)
        f2 = cpoly([6, 3 * c3, c2], z)                       # f''/2
        f3 = cpoly([4, c3], z)                               # f'''/6
        ab2 = lambda w: w[0] * w[0] + w[1] * w[1]
        lhs = sqrt_up(ab2(f0)) + sqrt_up(ab2(f2)) * rho ** 2 + sqrt_up(ab2(f3)) * rho ** 3 + rho ** 4
        rhs = sqrt_dn(ab2(f1)) * rho
        return (lhs < rhs) and (y0 > rho)

    # element evaluations
    def sig1(self, x):
        return isum(self.V1[i] * x[i] for i in range(4) if x[i])

    def sig2(self, x):
        return isum(self.V2[i] * x[i] for i in range(4) if x[i])

    def sig3(self, x):
        return (isum(self.R3[i] * x[i] for i in range(4) if x[i]), isum(self.J3[i] * x[i] for i in range(4) if x[i]))

    def svec(self, x):
        re, im = self.sig3(x)
        return (self.sig1(x).sq(), self.sig2(x).sq(), re.sq() + im.sq())

    def E(self, x):
        s = self.svec(x)
        return s[0] + s[1] + s[2]

    def gram(self):
        G = [[None] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(i, 4):
                g = self.V1[i] * self.V1[j] + self.V2[i] * self.V2[j] + self.R3[i] * self.R3[j] + self.J3[i] * self.J3[j]
                G[i][j] = G[j][i] = g
        return G


# ======================================================================================================================
# Polygon P+ = {r >= 0, r1 + r2 + r3 = 1, g.r >= 0 for the cone constraints g} (exact)
# ======================================================================================================================
def clip_polygon(gs):
    """Sutherland-Hodgman clipping of the simplex by the half-planes g.r >= 0 (exact rationals)."""
    poly = [(Fr(1), Fr(0), Fr(0)), (Fr(0), Fr(1), Fr(0)), (Fr(0), Fr(0), Fr(1))]
    for g in gs:
        if not poly:
            break
        h = [sum(g[k] * P[k] for k in range(3)) for P in poly]
        new = []
        m = len(poly)
        for k in range(m):
            P, Q, hP, hQ = poly[k], poly[(k + 1) % m], h[k], h[(k + 1) % m]
            if hP >= 0:
                new.append(P)
            if (hP > 0 and hQ < 0) or (hP < 0 and hQ > 0):
                t = hP / (hP - hQ)
                new.append(tuple(P[i] + t * (Q[i] - P[i]) for i in range(3)))
        out = []
        for P in new:                                  # drop exact duplicates
            if not out or P != out[-1]:
                out.append(P)
        if len(out) > 1 and out[0] == out[-1]:
            out.pop()
        poly = out
    # drop points lying strictly inside a segment of their two neighbours (exact collinearity)
    changed = True
    while changed and len(poly) > 2:
        changed = False
        m = len(poly)
        for k in range(m):
            A, Bp, C = poly[k - 1], poly[k], poly[(k + 1) % m]
            u = [Bp[i] - A[i] for i in range(3)]
            w = [C[i] - A[i] for i in range(3)]
            cross = (u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0])
            if cross == (0, 0, 0):
                poly.pop(k)
                changed = True
                break
    return poly


def bfs_vertices(cons):
    """All feasible solutions of the nonsingular systems [g_i; g_j; (1,1,1)] r = (0,0,1) (cross-check)."""
    out = set()
    for i in range(len(cons)):
        for j in range(i + 1, len(cons)):
            M = [list(cons[i]) + [Fr(0)], list(cons[j]) + [Fr(0)], [Fr(1)] * 3 + [Fr(1)]]
            ok = True
            for col in range(3):
                piv = next((r for r in range(col, 3) if M[r][col] != 0), None)
                if piv is None:
                    ok = False
                    break
                M[col], M[piv] = M[piv], M[col]
                pv = M[col][col]
                M[col] = [v / pv for v in M[col]]
                for r in range(3):
                    if r != col and M[r][col] != 0:
                        fct = M[r][col]
                        M[r] = [a - fct * b for a, b in zip(M[r], M[col])]
            if not ok:
                continue
            r = tuple(M[k][3] for k in range(3))
            if all(sum(g[k] * r[k] for k in range(3)) >= 0 for g in cons):
                out.add(r)
    return out


def least_T(target):
    """Least integer T >= 1 with T^4 >= target (target a positive Fraction), exactly."""
    target = Fr(target)
    T = max(1, isqrt(isqrt(target.numerator // target.denominator)))
    while Fr(T) ** 4 >= target and T > 1:
        T -= 1
    while Fr(T) ** 4 < target:
        T += 1
    assert Fr(T) ** 4 >= target and (T == 1 or Fr(T - 1) ** 4 < target)
    return T


# ======================================================================================================================
# Exact integer Fincke-Pohst
# ======================================================================================================================
class Enumerator:
    """All x in Z^4 with x^T A x <= Bnd for a positive definite symmetric integer matrix A, in exact integers.

    A = U^T D U with U unit upper triangular (computed in Fractions and checked by exact reconstruction), so
    x^T A x = sum_i d_i (x_i + sum_{j>i} mu_ij x_j)^2.  Row i is scaled to integers: N_i = den_i x_i + sum_j M_ij x_j,
    and S * x^T A x = sum_i W_i N_i^2 with integers W_i > 0.  At level i the admissible x_i form the integer interval
    { x : W_i (den_i x + K_i)^2 <= b } = [ceil((-K_i - r)/den_i), floor((r - K_i)/den_i)], r = isqrt(b // W_i)."""

    def __init__(self, A):
        n = 4
        Q = [[Fr(A[i][j]) for j in range(n)] for i in range(n)]
        d = [Fr(0)] * n
        mu = [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
        for i in range(n):
            d[i] = Q[i][i] - sum(d[k] * mu[k][i] ** 2 for k in range(i))
            if d[i] <= 0:
                raise ValueError('not positive definite')
            for j in range(i + 1, n):
                mu[i][j] = (Q[i][j] - sum(d[k] * mu[k][i] * mu[k][j] for k in range(i))) / d[i]
        for i in range(n):                         # exact reconstruction A = U^T D U
            for j in range(n):
                assert sum(d[k] * mu[k][i] * mu[k][j] for k in range(n)) == Q[i][j]
        self.den = []
        self.M = [[0] * n for _ in range(n)]
        for i in range(n):
            dn = 1
            for j in range(i + 1, n):
                dn = lcm(dn, mu[i][j].denominator)
            self.den.append(dn)
            for j in range(i + 1, n):
                self.M[i][j] = int(mu[i][j] * dn)
        w = [d[i] / (self.den[i] ** 2) for i in range(n)]
        S = 1
        for v in w:
            S = lcm(S, v.denominator)
        self.S = S
        self.W = [int(v * S) for v in w]
        assert all(Fr(self.W[i], S) == w[i] and self.W[i] > 0 for i in range(n))
        self.d = d

    @staticmethod
    def rng(K, b, W, den):
        if b < 0:
            return range(0)
        r = isqrt(b // W)
        return range(-((K + r) // den), (r - K) // den + 1)

    def run(self, Bnd, leaf):
        """Calls leaf(x1, x2, x3, x0range, K0, b0) per level-1 node; the caller handles level 0."""
        W, den, M = self.W, self.den, self.M
        Btot = self.S * Bnd
        for x3 in self.rng(0, Btot, W[3], den[3]):
            N3 = den[3] * x3
            b2 = Btot - W[3] * N3 * N3
            K2 = M[2][3] * x3
            for x2 in self.rng(K2, b2, W[2], den[2]):
                N2 = den[2] * x2 + K2
                b1 = b2 - W[2] * N2 * N2
                K1 = M[1][2] * x2 + M[1][3] * x3
                node = None
                for x1 in self.rng(K1, b1, W[1], den[1]):
                    N1 = den[1] * x1 + K1
                    b0 = b1 - W[1] * N1 * N1
                    K0 = M[0][1] * x1 + M[0][2] * x2 + M[0][3] * x3
                    r0 = self.rng(K0, b0, W[0], den[0])
                    if len(r0):
                        node = leaf(x1, x2, x3, r0, K0, b0, Btot, node)


# ======================================================================================================================
# Main
# ======================================================================================================================
print('verify_coverage.py -- independent re-check of 01_coverage.py (exact; own arithmetic, own intervals)')
print()

X = symbols('x0:4')
tsym = symbols('t')
SUMMARY = {}
TOTAL_ORBITS = 0

for key in DATA:
    tk = time.time()
    d = DATA[key]
    K = Order(key)
    print('=' * 118)
    print(f'CLASS {d["label"]}   F = {d["F"]}   f(t) = t^4 + ({K.c[3]}) t^3 + ({K.c[2]}) t^2 + ({K.c[1]}) t + ({K.c[0]})')
    check(True, 'integral basis is triangular and closed under multiplication (structure constants integral)')

    # ---------------------------------------------------------------- norm form
    p_expr = sum(X[i] * sum(Rational(K.B[i][j].numerator, K.B[i][j].denominator) * tsym ** j for j in range(4))
                 for i in range(4))
    f_expr = tsym ** 4 + sum(Rational(K.c[j]) * tsym ** j for j in range(4))
    NP = Poly(expand(resultant(f_expr, p_expr, tsym)), *X)
    nterms = NP.terms()
    check(all(c.is_integer for _, c in nterms), f'norm form Res_t(f, x(t)) has {len(nterms)} terms, all integer')
    mm, nn = symbols('m n')
    ident = expand(NP.as_expr().subs({X[0]: nn, X[1]: -mm, X[2]: 0, X[3]: 0})
                   - sum(c * mm ** (4 - k) * nn ** k for k, c in enumerate(d['F'])))
    check(ident == 0, 'N(n - m alpha) = F(m,n) as a polynomial identity in m, n')
    NT = [(int(c), e) for e, c in nterms]

    def norm(x, NT=NT):
        s = 0
        for c, (e0, e1, e2, e3) in NT:
            s += c * x[0] ** e0 * x[1] ** e1 * x[2] ** e2 * x[3] ** e3
        return s
    # coefficient structure for the enumeration: N = sum_k C_k x0^k, C_k = sum_l D_kl(x2,x3) x1^l
    Dkl = {}
    for c, (e0, e1, e2, e3) in NT:
        Dkl.setdefault((e0, e1), []).append((c, e2, e3))
    check(Dkl.get((4, 0)) == [(1, 0, 0)], 'the coefficient of x0^4 in the norm form is 1')

    # ---------------------------------------------------------------- units
    e1, e2 = d['UNITS']
    check(norm(e1) == 1 and norm(e2) == 1, f'N(eps1) = N(eps2) = 1 for eps1 = {e1}, eps2 = {e2} (norm form)')
    one = (1, 0, 0, 0)
    inv1, inv2 = K.quotient_in_O(e1, one), K.quotient_in_O(e2, one)
    check(inv1 is not None and inv2 is not None and K.mul(e1, inv1) == one and K.mul(e2, inv2) == one,
          'eps1^-1, eps2^-1 lie in O (exact solve)')

    def upow(a, b):
        x = one
        for _ in range(abs(a)):
            x = K.mul(x, e1 if a > 0 else inv1)
        for _ in range(abs(b)):
            x = K.mul(x, e2 if b > 0 else inv2)
        return x

    # ---------------------------------------------------------------- embeddings
    Em = Emb(K)
    print(f'   disc(f) = {Em.disc};  real roots (Sturm count) = {Em.nreal};  real isolating intervals (width < 2^-220):')
    for (lo, hi, sc, mult) in Em.real:
        print(f'      [{float(lo):.15f}, {float(hi):.15f}]  exact sign change of f: {sc}')
    check(Em.nreal == 2 and len(Em.real) == 2 and all(r[2] for r in Em.real) and Em.real[0][1] < Em.real[1][0],
          'exactly two real roots, in two disjoint intervals with exact sign changes')
    check(Em.disc < 0, 'disc(f) < 0, consistent with signature (2,1)')
    check(Em.rouche, f'Rouche certificate: exactly one root in |z - z0| < 2^-300, z0 = '
                     f'{float(Em.z0[0]):.12f} + {float(Em.z0[1]):.12f} i (Im z0 > 2^-300)')
    check(Em.sympy_rect_ok, "z0 lies in sympy's exact isolating rectangle of the root with Im > 0")

    # ---------------------------------------------------------------- f-vectors of U and the polygon
    print('   U = {eps1^a eps2^b : a,b in {-1,0,1}, (a,b) != (0,0)};  f_v(u) = s_v(u) (rigorous intervals):')
    cons_cone = []
    fmatch = True
    for (a, b) in EXPONENTS:
        u = upow(a, b)
        f = Em.svec(u)
        prod = f[0] * f[1] * f[2].sq()
        okp = prod.contains(1) and norm(u) == 1
        ref = OUT01[key]['f'][(a, b)]
        rel = max(abs(float(f[v]) - ref[v]) / ref[v] for v in range(3))
        fmatch &= rel < 2e-5
        width = max(float(f[v].hi - f[v].lo) for v in range(3))
        print(f'      (a,b)=({a:2d},{b:2d})  f = ({float(f[0]):.6g}, {float(f[1]):.6g}, {float(f[2]):.6g})'
              f'   widths < {width:.1e}   N(u) = 1 and 1 in f1 f2 f3^2: {okp}   max rel. diff. to 01: {rel:.1e}')
        if not okp:
            check(False, f'norm/product check for u = eps1^{a} eps2^{b}')
        cons_cone.append(tuple(f[v].hi - 1 for v in range(3)))      # relaxed constraint (f_v^+ - 1).r >= 0
    check(fmatch, "the f-vectors agree with 01's printed values (6 significant digits)")
    verts = clip_polygon(cons_cone)
    cons_all = [(Fr(1), Fr(0), Fr(0)), (Fr(0), Fr(1), Fr(0)), (Fr(0), Fr(0), Fr(1))] + cons_cone
    check(set(verts) == bfs_vertices(cons_all),
          f'clipping gives {len(verts)} vertices = the basic feasible solutions of all pairs of constraints')
    allpos = all(all(v > 0 for v in P) for P in verts)
    check(allpos, 'every vertex of P+ has all coordinates > 0 (hypothesis of the lemma)')
    for P in verts:
        print(f'      vertex ({float(P[0]):.6f}, {float(P[1]):.6f}, {float(P[2]):.6f})')
    ref = OUT01[key]['verts']
    match = len(ref) == len(verts) and all(any(all(abs(float(P[k]) - R[k]) < 1.01e-6 for k in range(3)) for R in ref)
                                           for P in verts)
    check(match, "the vertex set agrees with 01's printed vertices (6 decimals)")
    kappa = [min(P[v] for P in verts) for v in range(3)]
    bmin = min(P[0] * P[1] * P[2] ** 2 for P in verts)
    kprod = kappa[0] * kappa[1] * kappa[2] ** 2
    print(f'   kappa = ({float(kappa[0]):.6e}, {float(kappa[1]):.6e}, {float(kappa[2]):.6e});  kappa1 kappa2 kappa3^2 = '
          f'{float(kprod):.6e};  b_min = {float(bmin):.6e}')
    check(abs(float(bmin) / OUT01[key]['bmin'] - 1) < 1e-6, "b_min agrees with 01's printed value")
    check(bmin >= kprod, 'b_min >= kappa1 kappa2 kappa3^2 (the linear bound is the weaker one)')
    Tlin = {c: least_T(Fr(c * c) / kprod) for c in RHS}
    Tprod = {c: least_T(Fr(c * c) / bmin) for c in RHS}
    print(f"   linear bound T'(c) (least T' with T'^4 >= c^2/(kappa1 kappa2 kappa3^2)): {Tlin}")
    print(f'   product bound T(c) (least T with T^4 >= c^2/b_min), for comparison:      {Tprod}')
    check(Tprod == OUT01[key]['T'], "the product bounds T(c) equal 01's T values")
    check(all(Tlin[c] >= Tprod[c] for c in RHS), "T'(c) >= T(c)")

    # ---------------------------------------------------------------- rational lower form A <= E
    KB = 40
    G = Em.gram()
    A0 = [[None] * 4 for _ in range(4)]
    dl = [[None] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            mid = G[i][j].mid()
            A0[i][j] = Fr(round(mid * 2 ** KB), 2 ** KB)
            dev = max(G[i][j].hi - A0[i][j], A0[i][j] - G[i][j].lo)
            dl[i][j] = Fr(math.ceil(dev * 2 ** KB), 2 ** KB)
    Aint = [[int((A0[i][j] - (sum(dl[i]) if i == j else 0)) * 2 ** KB) for j in range(4)] for i in range(4)]
    check(all(Aint[i][j] == Aint[j][i] for i in range(4) for j in range(4)), 'lower form is symmetric')
    shift = max(sum(dl[i]) for i in range(4))
    print(f'   lower form A = A0 - diag(sum_j delta_ij), entries in 2^-{KB} Z; largest diagonal shift {float(shift):.2e}')
    # random consistency test: x^T A x <= E(x) (rigorous interval) for random integer vectors
    rnd = random.Random(1)
    okA = True
    for _ in range(200):
        x = tuple(rnd.randint(-50, 50) for _ in range(4))
        qa = Fr(sum(Aint[i][j] * x[i] * x[j] for i in range(4) for j in range(4)), 2 ** KB)
        okA &= qa <= Em.E(x).lo
    check(okA, 'x^T A x <= E(x) on 200 random integer vectors (rigorous E), as the shift argument proves')
    En = Enumerator(Aint)
    check(True, f'A is positive definite (LDL pivots > 0) and A = U^T D U exactly; scale S has {En.S.bit_length()} bits')

    # ---------------------------------------------------------------- 01's lattice-point counts with E-bound T(c)
    def count_leaf_factory(store):
        def leaf(x1, x2, x3, r0, K0, b0, Btot, node):
            for x0 in r0:
                N0 = En.den[0] * x0 + K0
                store.append((Btot - b0 + En.W[0] * N0 * N0, (x0, x1, x2, x3)))
            return None
        return leaf
    small = []
    Tsmall = Tprod[36]
    En.run(Tsmall * 2 ** KB, count_leaf_factory(small))
    cnts = {c: sum(1 for (u, x) in small if any(x) and u <= En.S * 2 ** KB * Tprod[c]) for c in RHS}
    print(f'   nonzero lattice points with x^T A x <= T(c): {cnts}   (01 printed {OUT01[key]["pts"]})')
    check(cnts == OUT01[key]['pts'], "same numbers of lattice points as 01 at its bounds T(c)")

    # ---------------------------------------------------------------- the main enumeration at T'_max
    Tmax = max(Tlin.values())
    TSET = set(RHS) | {-1}          # -1 only for the by-product on units of norm -1 (below)
    found = []
    nleaf = [0]
    sampled = [0, 0]                # [leaves sampled, mismatches] for the fast norm evaluation

    def leaf(x1, x2, x3, r0, K0, b0, Btot, node):
        if node is None or node[0] != (x2, x3):
            p2 = [1, x2, x2 * x2, x2 ** 3, x2 ** 4]
            p3 = [1, x3, x3 * x3, x3 ** 3, x3 ** 4]
            Dv = {}
            for (k, l), lst in Dkl.items():
                Dv[(k, l)] = sum(c * p2[a] * p3[b] for (c, a, b) in lst)
            node = ((x2, x3), Dv)
        Dv = node[1]
        C = []
        for k in range(5):
            v = 0
            for l in range(4 - k, -1, -1):
                v = v * x1 + Dv.get((k, l), 0)
            C.append(v)
        C0, C1, C2, C3 = C[0], C[1], C[2], C[3]
        nleaf[0] += len(r0)
        if len(r0) and (x1 + 3 * x2 + 7 * x3) % 23 == 0:           # deterministic sample of all leaves
            x0s = r0[len(r0) // 2]
            nv = (((x0s + C3) * x0s + C2) * x0s + C1) * x0s + C0
            sampled[0] += 1
            sampled[1] += nv != norm((x0s, x1, x2, x3))
        for x0 in r0:
            nv = (((x0 + C3) * x0 + C2) * x0 + C1) * x0 + C0
            if nv in TSET:
                N0 = En.den[0] * x0 + K0
                found.append((nv, Btot - b0 + En.W[0] * N0 * N0, (x0, x1, x2, x3)))
        return node
    te = time.time()
    En.run(Tmax * 2 ** KB, leaf)
    print(f"   exact enumeration of all x with x^T A x <= T'_max = {Tmax}: {nleaf[0]} lattice points "
          f"(incl. 0), {len(found)} with N(x) in {{-1,1,4,9,36}}  [{time.time() - te:.1f} s]")
    okn = all(norm(x) == nv for (nv, _, x) in found)
    check(okn, 'the fast norm evaluation agrees with the norm form on every kept point')
    check(sampled[0] > 100 and sampled[1] == 0,
          f'the fast norm evaluation agrees with the norm form on {sampled[0]} sampled leaves, kept or not '
          f'(mismatches: {sampled[1]})')

    # ---------------------------------------------------------------- units of norm +1
    BU = 30
    table = {}
    row_a = {0: one}
    xa = one
    for a in range(1, BU + 1):
        xa = K.mul(xa, e1)
        row_a[a] = xa
    xa = one
    for a in range(1, BU + 1):
        xa = K.mul(xa, inv1)
        row_a[-a] = xa
    for a in range(-BU, BU + 1):
        y = row_a[a]
        table[y] = (1, a, 0)
        table[tuple(-v for v in y)] = (-1, a, 0)
        yp = y
        for b in range(1, BU + 1):
            yp = K.mul(yp, e2)
            table[yp] = (1, a, b)
            table[tuple(-v for v in yp)] = (-1, a, b)
        yp = y
        for b in range(1, BU + 1):
            yp = K.mul(yp, inv2)
            table[yp] = (1, a, -b)
            table[tuple(-v for v in yp)] = (-1, a, -b)
    check(len(table) == 2 * (2 * BU + 1) ** 2,
          f'the {2 * (2 * BU + 1) ** 2} units +-eps1^a eps2^b, |a|,|b| <= {BU}, are pairwise distinct')
    L = {c: [x for (nv, u, x) in found if nv == c and u <= En.S * 2 ** KB * Tlin[c]] for c in RHS}
    units_ok = all(x in table for x in L[1])
    print(f"   c = 1: {len(L[1])} elements of norm 1 with x^T A x <= T'(1) = {Tlin[1]}; exponents (sign, a, b):")
    print('      ' + ', '.join(str(table.get(x, 'NOT IN TABLE')) for x in sorted(L[1], key=lambda y: table.get(y, (0, 99, 99))[1:])))
    check(units_ok, f'every one of them is +-eps1^a eps2^b with |a|,|b| <= {BU} (exact table lookup)')
    print('      => by the lemma (c = 1): every unit of norm +1 of O lies in <-1, eps1, eps2>, so two elements of')
    print('         equal norm c != 0 lie in one orbit of <-1, eps1, eps2> iff their quotient lies in O.')
    Lm1 = [x for (nv, u, x) in found if nv == -1 and u <= En.S * 2 ** KB * Tlin[1]]
    print(f"   by-product: elements of norm -1 with x^T A x <= T'(1) = {Tlin[1]}: {len(Lm1)}"
          + ('  => (lemma with c = -1, same bound since c^2 = 1) O has no unit of norm -1, so the unit group of O'
             ' is exactly <-1, eps1, eps2> = {+-1} x <eps1> x <eps2> (eps1, eps2 are fundamental units of O)'
             if not Lm1 else f'  e.g. {Lm1[:3]}'))

    # ---------------------------------------------------------------- orbits for each c
    SUMMARY[key] = {}
    for c in RHS:
        elems = sorted(L[c], key=lambda x: (Em.E(x).lo, x))
        reps, members = [], []
        for x in elems:
            for k, w in enumerate(reps):
                if K.quotient_in_O(w, x) is not None:
                    members[k].append(x)
                    break
            else:
                reps.append(x)
                members.append([x])
        TOTAL_ORBITS += len(reps)
        # E-certified count at 01's bound T(c)
        Ecert = [Em.E(x) for x in elems]
        below = sum(1 for e in Ecert if e.hi <= OUT01[key]['T'][c])
        undecided = sum(1 for e in Ecert if e.lo <= OUT01[key]['T'][c] < e.hi)
        print(f"   c = {c:2d}: T'(c) = {Tlin[c]};  elements of norm c with x^T A x <= T'(c): {len(elems)};  "
              f'orbits: {len(reps)}   (01: {OUT01[key]["norbits"][c]})')
        print(f'      of these, E(x) <= T(c) = {OUT01[key]["T"][c]} certified for {below} (undecided {undecided}); '
              f'01 listed {OUT01[key]["count"][c]}')
        check(len(reps) == OUT01[key]['norbits'][c], f'orbit count for c = {c} agrees with 01')
        check(undecided == 0 and below == OUT01[key]['count'][c], f"count of norm-{c} elements with E <= T(c) agrees with 01")
        # 01's seeds
        seed_orbit = []
        for w in OUT01[key]['seeds'][c]:
            inlist = w in elems
            k = next((k for k, r in enumerate(reps) if K.quotient_in_O(r, w) is not None), None)
            seed_orbit.append(k)
            check(norm(w) == c and inlist and k is not None, f'01 seed {w}: N = {c}, listed here, orbit #{k}')
        check(sorted(seed_orbit) == list(range(len(reps))), "01's seeds lie in distinct orbits and meet every orbit")
        for k, w in enumerate(reps):
            ek = Em.E(w)
            print(f'      orbit #{k}: least-E member {w}  E = {float(ek):.4f}  ({len(members[k])} members listed)')
        # known solutions
        for (m_, n_) in d['KNOWN']:
            if Fform(key, m_, n_) != c:
                continue
            z = (n_, -m_, 0, 0)
            check(norm(z) == c, f'known solution (m,n) = ({m_},{n_}): N(n - m alpha) = F(m,n) = {c}')
            k = next((k for k, r in enumerate(reps) if K.quotient_in_O(r, z) is not None), None)
            w01 = OUT01[key]['seeds'][c][seed_orbit.index(k)] if k is not None else None
            q = K.quotient_in_O(w01, z) if w01 is not None else None
            rep = table.get(q, 'unit outside the table') if q is not None else None
            print(f'      known (m,n) = ({m_},{n_}): z = {z} lies in orbit #{k};  z = sign * eps1^a eps2^b * w with '
                  f'w = 01 seed {w01}: (sign, a, b) = {rep}')
            check(k is not None, f'known solution ({m_},{n_}) lies in a listed orbit')
        SUMMARY[key][c] = [OUT01[key]['seeds'][c][seed_orbit.index(k)] for k in range(len(reps))]
    print(f'   [class done in {time.time() - tk:.1f} s]')

# ======================================================================================================================
# Appendix: exact random tests of the product inequality used in the lemma (sanity only; the proof is in the text)
# ======================================================================================================================
print()
print('=' * 118)
print('APPENDIX: exact random tests of prod_i l_i(sum_j t_j V_j) >= min_j prod_i l_i(V_j) for affine l_i > 0 on the V_j')
rnd = random.Random(7)
bad = 0
for trial in range(2000):
    p = rnd.randint(2, 5)                       # number of points
    m = rnd.randint(1, 6)                       # number of factors
    V = [[Fr(rnd.randint(1, 50), rnd.randint(1, 9)) for _ in range(m)] for _ in range(p)]   # l_i(V_j) > 0
    tw = [Fr(rnd.randint(0, 9)) for _ in range(p)]
    if sum(tw) == 0:
        continue
    tw = [t / sum(tw) for t in tw]
    vals = [sum(tw[j] * V[j][i] for j in range(p)) for i in range(m)]   # l_i is affine: values combine convexly
    lhs = 1
    for v in vals:
        lhs *= v
    rhs = min(math.prod(V[j]) for j in range(p))
    bad += lhs < rhs
check(bad == 0, f'product inequality held in all random exact trials (violations: {bad})')

# ======================================================================================================================
# Data comparison with core.py (data tables only)
# ======================================================================================================================
print()
print('=' * 118)
print("DATA: the retyped tables against core.py's CLASSES / UNITS / KNOWN")
sys.path.insert(0, '.')
from core import CLASSES as CORE_CLASSES, UNITS as CORE_UNITS, KNOWN as CORE_KNOWN   # data only
for key in DATA:
    same = (tuple(CORE_CLASSES[key]['F']) == DATA[key]['F']
            and [[Fr(v) for v in row] for row in CORE_CLASSES[key]['BASIS']] == [[Fr(v) for v in row] for row in DATA[key]['BASIS']]
            and tuple(map(tuple, CORE_UNITS[key])) == DATA[key]['UNITS']
            and list(CORE_KNOWN[key]) == DATA[key]['KNOWN'])
    check(same, f'class {DATA[key]["label"]}: F, BASIS, UNITS, KNOWN identical to core.py')

print()
print(f'TOTAL ORBITS over the 12 pairs (class, c): {TOTAL_ORBITS}  (01: 16)')
check(TOTAL_ORBITS == 16, 'total orbit count 16')
print('SEEDS (01 seeds, one per orbit found here):')
for key in SUMMARY:
    print(f"    '{key}': " + '{' + ', '.join(f'{c}: {SUMMARY[key][c]}' for c in RHS) + '},')
print()
if FAILS:
    print(f'CHECK FAILED: {len(FAILS)} check(s) failed:')
    for m_ in FAILS:
        print('   - ' + m_)
    print(f'[done in {time.time() - T0:.1f} s]')
    sys.exit(1)
print('ALL CHECKS PASSED: the independent enumeration with the weaker linear bound finds exactly the orbits of 01.')
print(f'[done in {time.time() - T0:.1f} s]')
