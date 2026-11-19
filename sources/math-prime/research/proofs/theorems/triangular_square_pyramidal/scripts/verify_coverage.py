"""verify_coverage.py -- independent verification of the COVERAGE step (01_coverage.py) for 8604.a1
(triangular = square-pyramidal): every z in O with N(z) = c, c in {1, 4, 9, 36}, is +-eps1^a eps2^b w for one of
the seeds w printed by 01_coverage.out.txt, and those seeds lie in pairwise distinct orbits.

Independence.  Nothing is imported from core.py, lifting2.py, 01_coverage.py or 02_exclusion.py.  The field data
(forms G, bases, units, right sides, known solutions) are retyped below as literals.  The three files core.py,
01_coverage.out.txt and 02_exclusion.py / .out.txt are only READ AS TEXT at the end (core.py through the ast parser,
nothing is executed) to compare the retyped data and the printed results with what is recomputed here.
Methods used here, all different from 01/core:
  * the order: multiplication matrices M_i = BT^-1 (sum_k B[i][k] C^k) BT from the companion matrix C of g (checked
    integral, commutative, associative, M_0 = I); norm = determinant of the multiplication matrix (Leibniz formula);
    quotients by Cramer's rule; the symbolic norm form det(sum X_i M_i) equals Res_t(g(t), x(t)) (sympy, exact) and
    N(n - m beta) = G(m,n) identically;
  * real roots: an own Sturm chain over Q (exactly two real roots, g squarefree), isolation by Sturm counts from
    the Cauchy bound, refinement by exact sign changes to width 2^-300;
  * the complex pair: a Rouche certificate around a rational z0 (exactly one root in |z - z0| < 2^-300, Im z0 > 2^-300),
    and |sigma_3(x)|^2 = Re^2 + Im^2 by complex interval Horner evaluation (core's S,P formula is only cross-checked);
  * own outward-rounded interval arithmetic on a 2^-320 grid; f_v^+ = interval upper endpoint rounded UP to 2^-96;
  * the polygon P+ by exact half-plane clipping of the simplex, cross-checked against all basic feasible pairs;
  * BOTH bounds of the template's Lemma 3.2: the product bound b_min (vertex minimum of r1 r2 r3^2) and the weaker
    LINEAR bound of Remark 3.2(ii), kappa_v = least v-coordinate over the vertices,
        E(w)^4 <= c^2 / (kappa1 kappa2 kappa3^2),
    which uses only linearity (no log-concavity).  The enumeration is run to the linear bound T_lin(c) (feasible
    for all three classes, see the estimates printed), so the cross-check does not depend on step 4 of the lemma;
  * a rational lower form A = Ghat - eps I with eps >= Frobenius norm of the Gram enclosure error (not 4 rho I),
    and a float-free Fincke-Pohst enumeration in exact integers (ranges by isqrt); the innermost coordinate x0 is
    scanned with the norm written as a monic quartic in x0, evaluated in int64 only where a rigorous magnitude bound
    shows that no overflow is possible (else in Python integers); every kept point is re-checked by the determinant;
  * orbits: two elements of norm c are merged only when their exact quotient (Cramer) lies in O AND is identified
    exactly as +-eps1^a eps2^b (exponents guessed in floating point, then verified by exact multiplication).
A failed check prints CHECK FAILED and the script exits with status 1.  Run as
    python verify_coverage.py > verify_coverage.out.txt
"""
import sys, os, re, ast, time, math, random, itertools
from fractions import Fraction as Fr
from math import isqrt, lcm
import numpy as np
import sympy
import mpmath

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
FAILS = []
NOK = [0]


def check(cond, msg):
    if cond:
        NOK[0] += 1
        print(f'   [ok] {msg}')
    else:
        FAILS.append(msg)
        print(f'   CHECK FAILED: {msg}')
    sys.stdout.flush()
    return cond


# ======================================================================================================================
# 0. Data (retyped literals)
# ======================================================================================================================
# G(m,n) = a0 m^4 + a1 m^3 n + a2 m^2 n^2 + a3 m n^3 + n^4 with F = (a0, a1, a2, a3, 1); g(t) = G(1,t), beta a root of g.
# BASIS rows = coordinates of omega_0..omega_3 in 1, beta, beta^2, beta^3.  UNITS in omega-coordinates.
DATA = {
    '10': dict(label='(1,0)', F=(4, -5, -3, -2, 1),
               BASIS=[[1, 0, 0, 0], [0, 1, 0, 0], [Fr(1, 3), Fr(2, 3), Fr(1, 3), 0], [Fr(1, 3), 0, 0, Fr(1, 3)]],
               UNITS=((0, 1, -1, 0), (0, -2, 1, 0)),
               KNOWN=[(0, 1), (1, -1), (1, 0), (2, 1)]),
    '01': dict(label='(0,1)', F=(9, -72, 18, 0, 1),
               BASIS=[[1, 0, 0, 0], [0, 1, 0, 0], [Fr(1, 2), 0, Fr(1, 6), 0], [Fr(3, 4), Fr(3, 4), Fr(1, 12), Fr(1, 12)]],
               UNITS=((1, -1, 1, 0), (1, 0, -2, 0)),
               KNOWN=[(0, 1), (1, 0), (1, 3)]),
    '11': dict(label='(1,1)', F=(36, -72, 0, 0, 1),
               BASIS=[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, Fr(1, 6), 0], [0, Fr(1, 2), 0, Fr(1, 12)]],
               UNITS=((0, 0, -1, 0), (8, -4, 3, 0)),
               KNOWN=[(0, 1), (1, 0), (1, 4), (2, 1)]),
}
RHS = (1, 4, 9, 36)
EXPS = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)]
# claimed solutions of G(m,n) = c, one per pair +-(m,n) (from the task statement; compared with KNOWN)
CLAIMED = {
    '10': {1: [(0, 1)], 4: [(1, 0)], 9: [(1, -1), (2, 1)], 36: []},
    '01': {1: [(0, 1)], 4: [], 9: [(1, 0)], 36: [(1, 3)]},
    '11': {1: [(0, 1), (2, 1)], 4: [(1, 4)], 9: [], 36: [(1, 0)]},
}


def Gform(key, m, n):
    a0, a1, a2, a3, a4 = DATA[key]['F']
    return a0 * m ** 4 + a1 * m ** 3 * n + a2 * m ** 2 * n ** 2 + a3 * m * n ** 3 + a4 * n ** 4


# ======================================================================================================================
# 1. Exact linear algebra and the order O
# ======================================================================================================================
PERMS = [(p, (-1) ** sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j]))
         for p in itertools.permutations(range(4))]


def det4(M):
    """Leibniz formula (24 terms), exact for integers or Fractions."""
    s = 0
    for p, sg in PERMS:
        s += sg * M[0][p[0]] * M[1][p[1]] * M[2][p[2]] * M[3][p[3]]
    return s


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def matinv(A):
    n = len(A)
    M = [[Fr(A[i][j]) for j in range(n)] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [row[n:] for row in M]


class Order:
    """O = Z omega_0 + ... + Z omega_3 inside Q[t]/(g).  Elements are integer coordinate 4-tuples (columns)."""

    def __init__(self, key):
        d = DATA[key]
        self.key, self.label = key, d['label']
        a0, a1, a2, a3, a4 = d['F']
        assert a4 == 1
        self.g = [Fr(a0), Fr(a1), Fr(a2), Fr(a3), Fr(1)]           # g(t) = t^4 + a3 t^3 + a2 t^2 + a1 t + a0
        C = [[Fr(0)] * 4 for _ in range(4)]                         # multiplication by t on power coordinates
        for k in range(3):
            C[k + 1][k] = Fr(1)
        for k in range(4):
            C[k][3] = -self.g[k]
        Cp = [[[Fr(int(i == j)) for j in range(4)] for i in range(4)]]
        for _ in range(3):
            Cp.append(matmul(Cp[-1], C))
        self.B = [[Fr(v) for v in row] for row in d['BASIS']]
        BT = [[self.B[i][k] for i in range(4)] for k in range(4)]    # power coords = BT * omega coords
        BTi = matinv(BT)
        self.integral = True
        self.Mi = []
        for i in range(4):
            P = [[sum(self.B[i][k] * Cp[k][r][s] for k in range(4)) for s in range(4)] for r in range(4)]
            Mw = matmul(matmul(BTi, P), BT)
            if any(v.denominator != 1 for row in Mw for v in row):
                self.integral = False
            self.Mi.append([[int(v) for v in row] for row in Mw])

    def mat(self, x):
        Mi = self.Mi
        return [[x[0] * Mi[0][r][s] + x[1] * Mi[1][r][s] + x[2] * Mi[2][r][s] + x[3] * Mi[3][r][s]
                 for s in range(4)] for r in range(4)]

    def mul(self, x, y):
        M = self.mat(x)
        return tuple(M[r][0] * y[0] + M[r][1] * y[1] + M[r][2] * y[2] + M[r][3] * y[3] for r in range(4))

    def norm(self, x):
        return det4(self.mat(x))

    def quot(self, x, y):
        """y / x in omega-coordinates (Fractions), by Cramer's rule on M(x) q = y."""
        M = self.mat(x)
        D = det4(M)
        out = []
        for k in range(4):
            Mk = [row[:] for row in M]
            for r in range(4):
                Mk[r][k] = y[r]
            out.append(Fr(det4(Mk), D))
        return tuple(out)

    def quot_in_O(self, x, y):
        q = self.quot(x, y)
        if all(v.denominator == 1 for v in q):
            return tuple(int(v) for v in q)
        return None


def charpoly_coeffs(M):
    """Faddeev-LeVerrier: det(t I - M) = t^4 + c1 t^3 + c2 t^2 + c3 t + c4 (Fractions)."""
    n = 4
    I = [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
    Mk = [[Fr(0)] * n for _ in range(n)]
    cs = [Fr(1)]
    for k in range(1, n + 1):
        Mk = [[sum(M[i][l] * Mk[l][j] for l in range(n)) + cs[-1] * I[i][j] for j in range(n)] for i in range(n)]
        AM = matmul(M, Mk)
        ck = -sum(AM[i][i] for i in range(n)) / k
        cs.append(ck)
    return cs


# ======================================================================================================================
# 2. Polynomials over Q, Sturm chain, real roots, Rouche certificate for the complex root
# ======================================================================================================================
def ptrim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def prem(a, b):
    a = [Fr(v) for v in ptrim(a)]
    b = [Fr(v) for v in ptrim(b)]
    while len(a) >= len(b) and not (len(a) == 1 and a[0] == 0):
        k = len(a) - len(b)
        c = a[-1] / b[-1]
        for i in range(len(b)):
            a[i + k] -= c * b[i]
        a.pop()
        if not a:
            a = [Fr(0)]
        a = ptrim(a)
    return a


def peval(p, x):
    v = Fr(0)
    for c in reversed(p):
        v = v * x + c
    return v


def sturm_chain(p):
    p = ptrim([Fr(v) for v in p])
    dp = ptrim([k * p[k] for k in range(1, len(p))])
    ch = [p, dp]
    while True:
        r = prem(ch[-2], ch[-1])
        if len(r) == 1 and r[0] == 0:
            break
        ch.append([-v for v in r])
    return ch


def sign_changes(vals):
    s = [v for v in vals if v != 0]
    return sum(1 for i in range(len(s) - 1) if (s[i] > 0) != (s[i + 1] > 0))


def V_at(ch, x):
    return sign_changes([peval(p, x) for p in ch])


def V_inf(ch, sgn):
    return sign_changes([p[-1] * (1 if (sgn > 0 or (len(p) - 1) % 2 == 0) else -1) for p in ch])


def real_roots(g, wbits=300):
    """Exactly isolated real roots of g (Sturm), refined by exact sign changes to width < 2^-wbits."""
    ch = sturm_chain(g)
    R = 1 + max(abs(c) for c in g[:-1])          # Cauchy bound (g monic)
    nreal = V_inf(ch, -1) - V_inf(ch, +1)
    work, iso = [(Fr(-R), Fr(R))], []
    while work:
        a, b = work.pop()
        k = V_at(ch, a) - V_at(ch, b)
        if k == 0:
            continue
        if k == 1:
            iso.append((a, b))
            continue
        m = (a + b) / 2
        work += [(a, m), (m, b)]
    iso.sort()
    out = []
    for a, b in iso:
        ga = peval(g, a)
        assert ga * peval(g, b) < 0
        while b - a >= Fr(1, 1 << wbits):
            m = (a + b) / 2
            gm = peval(g, m)
            assert gm != 0                        # g has no rational root (irreducible, checked)
            if (gm < 0) == (ga < 0):
                a, ga = m, gm
            else:
                b = m
        out.append((a, b, V_at(ch, a) - V_at(ch, b), peval(g, a) * peval(g, b) < 0))
    return ch, nreal, out


def sqrt_hi(q, bits=400):
    q = Fr(q)
    n = -((-(q.numerator << (2 * bits))) // q.denominator)
    s = isqrt(n)
    if s * s < n:
        s += 1
    return Fr(s, 1 << bits)


def sqrt_lo(q, bits=400):
    q = Fr(q)
    return Fr(isqrt((q.numerator << (2 * bits)) // q.denominator), 1 << bits)


def ceval(coeffs_high_first, z):
    vr, vi = Fr(coeffs_high_first[0]), Fr(0)
    for c in coeffs_high_first[1:]:
        vr, vi = vr * z[0] - vi * z[1] + c, vr * z[1] + vi * z[0]
    return vr, vi


def complex_root(g, rho_bits=300):
    """Rational z0 (from a floating approximation) and a Rouche certificate: on |z - z0| = rho,
    |g(z) - g'(z0)(z - z0)| <= |g(z0)| + |g''(z0)/2| rho^2 + |g'''(z0)/6| rho^3 + rho^4 < |g'(z0)| rho,
    so g has exactly one root in |z - z0| < rho; with Im z0 > rho that root is non-real with Im > 0."""
    a0, a1, a2, a3 = [int(v) for v in g[:4]]
    mpmath.mp.dps = 160
    rts = mpmath.polyroots([1, a3, a2, a1, a0], maxsteps=500, extraprec=800)
    up = [z for z in rts if mpmath.im(z) > mpmath.mpf(10) ** -30]
    assert len(up) == 1
    S = 1 << 440
    x0 = Fr(int(mpmath.nint(mpmath.re(up[0]) * S)), S)
    y0 = Fr(int(mpmath.nint(mpmath.im(up[0]) * S)), S)
    rho = Fr(1, 1 << rho_bits)
    f0 = ceval([1, a3, a2, a1, a0], (x0, y0))
    f1 = ceval([4, 3 * a3, 2 * a2, a1], (x0, y0))
    f2 = ceval([6, 3 * a3, a2], (x0, y0))
    f3 = ceval([4, a3], (x0, y0))
    ab2 = lambda w: w[0] * w[0] + w[1] * w[1]
    lhs = sqrt_hi(ab2(f0)) + sqrt_hi(ab2(f2)) * rho ** 2 + sqrt_hi(ab2(f3)) * rho ** 3 + rho ** 4
    rhs = sqrt_lo(ab2(f1)) * rho
    return (x0, y0), rho, (lhs < rhs and y0 > rho), lhs / rhs


# ======================================================================================================================
# 3. Interval arithmetic (own), embeddings
# ======================================================================================================================
G_BITS = 320


def rdn(q):
    return Fr((q.numerator << G_BITS) // q.denominator, 1 << G_BITS)


def rup(q):
    return Fr(-((-(q.numerator << G_BITS)) // q.denominator), 1 << G_BITS)


class IV:
    """[lo, hi] with Fraction endpoints, every result rounded outward to the 2^-320 grid."""
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        lo = Fr(lo)
        hi = lo if hi is None else Fr(hi)
        if lo > hi:
            raise ValueError('empty interval')
        self.lo, self.hi = rdn(lo), rup(hi)

    @staticmethod
    def c(o):
        return o if isinstance(o, IV) else IV(o)

    def __add__(self, o):
        o = IV.c(o)
        return IV(self.lo + o.lo, self.hi + o.hi)
    __radd__ = __add__

    def __neg__(self):
        return IV(-self.hi, -self.lo)

    def __sub__(self, o):
        o = IV.c(o)
        return IV(self.lo - o.hi, self.hi - o.lo)

    def __rsub__(self, o):
        return IV.c(o) - self

    def __mul__(self, o):
        o = IV.c(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return IV(min(p), max(p))
    __rmul__ = __mul__

    def sq(self):
        if self.lo >= 0:
            return IV(self.lo * self.lo, self.hi * self.hi)
        if self.hi <= 0:
            return IV(self.hi * self.hi, self.lo * self.lo)
        return IV(0, max(self.lo * self.lo, self.hi * self.hi))

    def recip(self):
        if not (self.lo > 0 or self.hi < 0):
            raise ZeroDivisionError
        return IV(1 / self.hi, 1 / self.lo)

    def meets(self, o):
        o = IV.c(o)
        return self.lo <= o.hi and o.lo <= self.hi

    def contains(self, q):
        return self.lo <= q <= self.hi

    def mid(self):
        return (self.lo + self.hi) / 2

    def width(self):
        return self.hi - self.lo

    def __float__(self):
        return float(self.mid())


def isum(it):
    acc = IV(0)
    for t in it:
        acc = acc + t
    return acc


class Emb:
    """sigma_1 < sigma_2 (the real roots) and sigma_3 (the root with Im > 0), all rigorous."""

    def __init__(self, K):
        self.K = K
        g = K.g
        self.chain, self.nreal, rr = real_roots(g)
        self.rr = rr
        self.r1, self.r2 = IV(rr[0][0], rr[0][1]), IV(rr[1][0], rr[1][1])
        self.z0, self.rho, self.rouche, self.rouche_ratio = complex_root(g)
        x0, y0 = self.z0
        self.zr, self.zi = IV(x0 - self.rho, x0 + self.rho), IV(y0 - self.rho, y0 + self.rho)
        p1, p2, pc = [IV(1)], [IV(1)], [(IV(1), IV(0))]
        for _ in range(3):
            p1.append(p1[-1] * self.r1)
            p2.append(p2[-1] * self.r2)
            a, b = pc[-1]
            pc.append((a * self.zr - b * self.zi, a * self.zi + b * self.zr))
        self.V1, self.V2, self.VR, self.VI = [], [], [], []
        for i in range(4):
            row = K.B[i]
            self.V1.append(isum(p1[k] * row[k] for k in range(4) if row[k]))
            self.V2.append(isum(p2[k] * row[k] for k in range(4) if row[k]))
            self.VR.append(isum(pc[k][0] * row[k] for k in range(4) if row[k]))
            self.VI.append(isum(pc[k][1] * row[k] for k in range(4) if row[k]))
        self.fV1 = [float(v) for v in self.V1]
        self.fV2 = [float(v) for v in self.V2]

    def comps(self, x):
        s1 = isum(self.V1[i] * x[i] for i in range(4) if x[i])
        s2 = isum(self.V2[i] * x[i] for i in range(4) if x[i])
        re = isum(self.VR[i] * x[i] for i in range(4) if x[i])
        im = isum(self.VI[i] * x[i] for i in range(4) if x[i])
        return s1, s2, re, im

    def svec(self, x):
        s1, s2, re, im = self.comps(x)
        return (s1.sq(), s2.sq(), re.sq() + im.sq())

    def E(self, x):
        s = self.svec(x)
        return s[0] + s[1] + s[2]

    def gram(self):
        Gm = [[None] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(i, 4):
                Gm[i][j] = Gm[j][i] = (self.V1[i] * self.V1[j] + self.V2[i] * self.V2[j]
                                       + self.VR[i] * self.VR[j] + self.VI[i] * self.VI[j])
        return Gm

    def flog(self, x):
        """Floating-point (log|sigma_1 x|, log|sigma_2 x|), only to GUESS unit exponents."""
        a = sum(self.fV1[i] * x[i] for i in range(4))
        b = sum(self.fV2[i] * x[i] for i in range(4))
        return math.log(abs(a)), math.log(abs(b))


def s3_by_SP(Em, K, x):
    """core.py's formula: |sigma_3(x)|^2 = A^2 + S A B + P B^2, x = A + B t mod t^2 - S t + P,
    S = -a3 - sigma_1 - sigma_2, P = a0 / (sigma_1 sigma_2).  Used here only as a cross-check."""
    S = IV(-K.g[3]) - Em.r1 - Em.r2
    P = IV(K.g[0]) * (Em.r1 * Em.r2).recip()
    pc = [sum(Fr(x[i]) * K.B[i][k] for i in range(4)) for k in range(4)]
    A = IV(pc[0]) - P * pc[2] - S * P * pc[3]
    Bv = IV(pc[1]) + S * pc[2] + (S.sq() - P) * pc[3]
    return A.sq() + S * A * Bv + P * Bv.sq(), S, P


# ======================================================================================================================
# 4. The polygon P+ and the bounds
# ======================================================================================================================
FBITS = 96


def round_up(q, bits=FBITS):
    return Fr(-((-(q.numerator << bits)) // q.denominator), 1 << bits)


def dot3(g, r):
    return g[0] * r[0] + g[1] * r[1] + g[2] * r[2]


def clip_simplex(cons):
    """Vertices (cyclic order) of {r >= 0, r1+r2+r3 = 1, g.r >= 0 for g in cons}, by exact half-plane clipping."""
    poly = [(Fr(1), Fr(0), Fr(0)), (Fr(0), Fr(1), Fr(0)), (Fr(0), Fr(0), Fr(1))]
    for g in cons:
        h = [dot3(g, P) for P in poly]
        new = []
        m = len(poly)
        for k in range(m):
            P, Q, hP, hQ = poly[k], poly[(k + 1) % m], h[k], h[(k + 1) % m]
            if hP >= 0:
                new.append(P)
            if (hP > 0 > hQ) or (hP < 0 < hQ):
                t = hP / (hP - hQ)
                new.append(tuple(P[i] + t * (Q[i] - P[i]) for i in range(3)))
        out = []
        for P in new:
            if P not in out:
                out.append(P)
        poly = out
        if not poly:
            return []
    changed = True                      # drop points interior to an edge (exact collinearity)
    while changed and len(poly) > 2:
        changed = False
        for k in range(len(poly)):
            A, Bp, C = poly[k - 1], poly[k], poly[(k + 1) % len(poly)]
            u = [Bp[i] - A[i] for i in range(3)]
            w = [C[i] - A[i] for i in range(3)]
            if (u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]) == (0, 0, 0):
                poly.pop(k)
                changed = True
                break
    return poly


def basic_feasible(cons):
    """All feasible solutions of the nonsingular systems g_i.r = g_j.r = 0, r1+r2+r3 = 1 (incl. r_v >= 0)."""
    allc = [(Fr(1), Fr(0), Fr(0)), (Fr(0), Fr(1), Fr(0)), (Fr(0), Fr(0), Fr(1))] + list(cons)
    out = set()
    for i in range(len(allc)):
        for j in range(i + 1, len(allc)):
            p, q = allc[i], allc[j]
            # r is proportional to the cross product p x q (the line g_i.r = g_j.r = 0), normalized to sum 1
            cr = (p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0])
            s = cr[0] + cr[1] + cr[2]
            if s == 0:
                continue
            r = (cr[0] / s, cr[1] / s, cr[2] / s)
            if all(dot3(g, r) >= 0 for g in allc):
                out.add(r)
    return out


def least_T(target):
    """Least integer T >= 1 with T^4 >= target (exact)."""
    target = Fr(target)
    T = max(1, isqrt(isqrt(target.numerator // target.denominator)))
    while T > 1 and Fr(T - 1) ** 4 >= target:
        T -= 1
    while Fr(T) ** 4 < target:
        T += 1
    assert Fr(T) ** 4 >= target and (T == 1 or Fr(T - 1) ** 4 < target)
    return T


# ======================================================================================================================
# 5. Rational lower form and exact integer Fincke-Pohst
# ======================================================================================================================
KB = 40                     # Ghat on the 2^-40 grid


def lower_form(Gm):
    """Ghat = midpoints rounded to 2^-KB (symmetric), delta_ij >= |G*_ij - Ghat_ij|, eps >= sqrt(sum delta_ij^2)
    (a bound for the spectral norm of G* - Ghat), A = Ghat - eps I.  Then x^T A x <= x^T G* x = E(x) for all real x,
    and E(x) <= x^T A x + 2 eps |x|^2.  Returns 2^KB A as an integer matrix and eps."""
    Gh = [[None] * 4 for _ in range(4)]
    dl = [[None] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(i, 4):
            gh = Fr(round(Gm[i][j].mid() * (1 << KB)), 1 << KB)
            Gh[i][j] = Gh[j][i] = gh
            dl[i][j] = dl[j][i] = max(Gm[i][j].hi - gh, gh - Gm[i][j].lo)
    eps = round_up(sqrt_hi(sum(dl[i][j] ** 2 for i in range(4) for j in range(4))), KB)
    A = [[Gh[i][j] - (eps if i == j else 0) for j in range(4)] for i in range(4)]
    Aint = [[int(A[i][j] * (1 << KB)) for j in range(4)] for i in range(4)]
    assert all(Fr(Aint[i][j], 1 << KB) == A[i][j] for i in range(4) for j in range(4))
    return Aint, eps, Gh


def qform(Aint, x):
    return sum(Aint[i][j] * x[i] * x[j] for i in range(4) for j in range(4))


class Enumerator:
    """All x in Z^4 with x^T A x <= Bd (A symmetric positive definite, integer), exactly, without floats.
    x^T A x = sum_i d_i (x_i + sum_{j>i} mu_ij x_j)^2 (exact LDL^T over Q, verified by reconstruction).  With the
    integer forms L_i = den_i x_i + sum_{j>i} M_ij x_j and integers W_i = S d_i / den_i^2 > 0:
    S x^T A x = sum_i W_i L_i^2.  Given the outer coordinates and the remaining budget R_i, the admissible x_i are
    exactly those with |L_i| <= isqrt(R_i // W_i)."""

    def __init__(self, A):
        n = 4
        Q = [[Fr(v) for v in row] for row in A]
        R = [row[:] for row in Q]
        d = [None] * n
        mu = [[Fr(0)] * n for _ in range(n)]
        for i in range(n):
            d[i] = R[i][i]
            if d[i] <= 0:
                raise ValueError('lower form not positive definite')
            for j in range(i + 1, n):
                mu[i][j] = R[i][j] / d[i]
            for k in range(i + 1, n):
                for l in range(i + 1, n):
                    R[k][l] = R[k][l] - R[k][i] * R[i][l] / d[i]
        for a in range(n):                           # exact reconstruction A = U^T D U
            for b in range(n):
                s = Fr(0)
                for i in range(min(a, b) + 1):
                    ua = Fr(1) if a == i else mu[i][a]
                    ub = Fr(1) if b == i else mu[i][b]
                    s += d[i] * ua * ub
                assert s == Q[a][b]
        self.d = d
        self.den = []
        self.M = [[0] * n for _ in range(n)]
        for i in range(n):
            dn = 1
            for j in range(i + 1, n):
                dn = lcm(dn, mu[i][j].denominator)
            self.den.append(dn)
            for j in range(i + 1, n):
                self.M[i][j] = int(mu[i][j] * dn)
                assert Fr(self.M[i][j], dn) == mu[i][j]
        w = [d[i] / self.den[i] ** 2 for i in range(n)]
        S = 1
        for v in w:
            S = lcm(S, v.denominator)
        self.S = S
        self.W = [int(v * S) for v in w]
        assert all(Fr(self.W[i], S) == w[i] and self.W[i] > 0 for i in range(n))

    @staticmethod
    def rng(K, R, W, den):
        if R < 0:
            return 0, -1
        r = isqrt(R // W)
        return -((r + K) // den), (r - K) // den

    def points(self, Bd):
        """Every x with x^T A x <= Bd (pure Python)."""
        W, den, M, rng = self.W, self.den, self.M, self.rng
        Rt = self.S * Bd
        lo3, hi3 = rng(0, Rt, W[3], den[3])
        for x3 in range(lo3, hi3 + 1):
            L3 = den[3] * x3
            R2 = Rt - W[3] * L3 * L3
            K2 = M[2][3] * x3
            lo2, hi2 = rng(K2, R2, W[2], den[2])
            for x2 in range(lo2, hi2 + 1):
                L2 = den[2] * x2 + K2
                R1 = R2 - W[2] * L2 * L2
                K1 = M[1][2] * x2 + M[1][3] * x3
                lo1, hi1 = rng(K1, R1, W[1], den[1])
                for x1 in range(lo1, hi1 + 1):
                    L1 = den[1] * x1 + K1
                    R0 = R1 - W[1] * L1 * L1
                    K0 = M[0][1] * x1 + M[0][2] * x2 + M[0][3] * x3
                    lo0, hi0 = rng(K0, R0, W[0], den[0])
                    for x0 in range(lo0, hi0 + 1):
                        yield (x0, x1, x2, x3)

    def scan(self, Bd, Dkl, targets):
        """Every x with x^T A x <= Bd and N(x) in targets.  N(x) = x0^4 + C3 x0^3 + C2 x0^2 + C1 x0 + C0 with
        C_k = sum_l D_kl(x2,x3) x1^l (Dkl from the symbolic norm form).  The x0-range is evaluated in numpy int64
        only when |x0|^4 + |C3||x0|^3 + |C2||x0|^2 + |C1||x0| + |C0| < 2^62 on the whole range, which bounds every
        intermediate value of the Horner scheme; otherwise in Python integers."""
        W, den, M, rng = self.W, self.den, self.M, self.rng
        W0, W1, W2, W3 = W
        d0, d1, d2, d3 = den
        tset, tmin, tmax = set(targets), min(targets), max(targets)
        LIM = 1 << 62
        hits = []
        nnodes = npts = nslow = 0
        Rt = self.S * Bd
        lo3, hi3 = rng(0, Rt, W3, d3)
        for x3 in range(lo3, hi3 + 1):
            L3 = d3 * x3
            R2 = Rt - W3 * L3 * L3
            K2 = M[2][3] * x3
            lo2, hi2 = rng(K2, R2, W2, d2)
            p3 = (1, x3, x3 ** 2, x3 ** 3, x3 ** 4)
            for x2 in range(lo2, hi2 + 1):
                L2 = d2 * x2 + K2
                R1 = R2 - W2 * L2 * L2
                K1 = M[1][2] * x2 + M[1][3] * x3
                lo1, hi1 = rng(K1, R1, W1, d1)
                if lo1 > hi1:
                    continue
                p2 = (1, x2, x2 ** 2, x2 ** 3, x2 ** 4)
                Dv = {kl: sum(c * p2[a] * p3[b] for (c, a, b) in lst) for kl, lst in Dkl.items()}
                c3h = [Dv.get((3, l), 0) for l in (1, 0)]
                c2h = [Dv.get((2, l), 0) for l in (2, 1, 0)]
                c1h = [Dv.get((1, l), 0) for l in (3, 2, 1, 0)]
                c0h = [Dv.get((0, l), 0) for l in (4, 3, 2, 1, 0)]
                K0b = M[0][2] * x2 + M[0][3] * x3
                for x1 in range(lo1, hi1 + 1):
                    L1 = d1 * x1 + K1
                    R0 = R1 - W1 * L1 * L1
                    if R0 < 0:
                        continue
                    K0 = M[0][1] * x1 + K0b
                    r = isqrt(R0 // W0)
                    lo0, hi0 = -((r + K0) // d0), (r - K0) // d0
                    if lo0 > hi0:
                        continue
                    nnodes += 1
                    npts += hi0 - lo0 + 1
                    C3 = c3h[0] * x1 + c3h[1]
                    C2 = (c2h[0] * x1 + c2h[1]) * x1 + c2h[2]
                    C1 = ((c1h[0] * x1 + c1h[1]) * x1 + c1h[2]) * x1 + c1h[3]
                    C0 = (((c0h[0] * x1 + c0h[1]) * x1 + c0h[2]) * x1 + c0h[3]) * x1 + c0h[4]
                    if hi0 - lo0 >= 16:
                        X = max(-lo0, hi0, 1)
                        if X ** 4 + abs(C3) * X ** 3 + abs(C2) * X * X + abs(C1) * X + abs(C0) < LIM:
                            x0 = np.arange(lo0, hi0 + 1, dtype=np.int64)
                            v = (((x0 + C3) * x0 + C2) * x0 + C1) * x0 + C0
                            for i in np.flatnonzero((v >= tmin) & (v <= tmax)):
                                vi = int(v[i])
                                if vi in tset:
                                    hits.append((vi, (lo0 + int(i), x1, x2, x3)))
                            continue
                        nslow += 1
                    for x0 in range(lo0, hi0 + 1):
                        v = (((x0 + C3) * x0 + C2) * x0 + C1) * x0 + C0
                        if v in tset:
                            hits.append((v, (x0, x1, x2, x3)))
        return hits, nnodes, npts, nslow


# ======================================================================================================================
# 6. Reading the outputs of 01 and 02 and the data of core.py (as text only)
# ======================================================================================================================
def parse_01(path):
    txt = open(path, encoding='utf-8').read()
    res = {}
    for blk in re.split(r'\n=+\n', txt):
        m = re.match(r'CLASS \((\d),(\d)\)', blk.strip())
        if not m:
            continue
        key = m.group(1) + m.group(2)
        d = dict(f={}, verts=[], T={}, pts={}, normc={}, norbits={}, seeds={}, unitline=None)
        for mm in re.finditer(r'\(a,b\)=\(\s*(-?\d+),\s*(-?\d+)\)\s+f = \(([^)]*)\)', blk):
            d['f'][(int(mm.group(1)), int(mm.group(2)))] = tuple(float(v) for v in mm.group(3).split(','))
        for mm in re.finditer(r'vertex \(([^)]*)\)\s+r1 r2 r3\^2 = (\S+)', blk):
            d['verts'].append(tuple(float(v) for v in mm.group(1).split(',')))
        d['bmin'] = float(re.search(r'b_min = min r1 r2 r3\^2 over vertices = (\S+);', blk).group(1))
        d['nverts'] = int(re.search(r'relaxed polygon P\+ has (\d+) vertices', blk).group(1))
        for part in re.split(r'\n   c = ', blk)[1:]:
            mc = re.match(r'\s*(\d+): T = (\d+) \(T\^4 >= c\^2/b_min\);\s+lattice points with Gt <= T: (\d+);\s+'
                          r'of norm c: (\d+);\s+orbits of <-1,eps1,eps2>: (\d+)', part)
            c = int(mc.group(1))
            d['T'][c], d['pts'][c], d['normc'][c], d['norbits'][c] = (int(mc.group(k)) for k in range(2, 6))
            d['seeds'][c] = [(tuple(int(v) for v in s.group(1).split(',')), int(s.group(2)))
                             for s in re.finditer(r'seed w = \(([^)]*)\)\s+E\(w\) ~ \S+\s+\((\d+) listed members', part)]
            mu = re.search(r'containing \+-1: (True|False); elements of norm -1 with E <= T: (\d+)', part)
            if mu:
                d['unitline'] = (mu.group(1) == 'True', int(mu.group(2)))
            d.setdefault('notes', 0)
            d['notes'] += len(re.findall(r'NOTE:', part))
        res[key] = d
    ms = re.search(r'\nSEEDS = (\{.*?\n\})', txt, re.S)
    return res, ast.literal_eval(ms.group(1))


def parse_02_seeds(path):
    src = open(path, encoding='utf-8').read()
    ms = re.search(r'\nSEEDS = \{[^\n]*\n(.*?)\n\}', src, re.S)
    body = '\n'.join(re.sub(r'#.*$', '', line) for line in ms.group(1).split('\n'))
    return ast.literal_eval('{' + body + '}')


def parse_02_exponents(path):
    out, key = {}, None
    for line in open(path, encoding='utf-8'):
        m = re.match(r'CLASS \((\d),(\d)\)', line)
        if m:
            key = m.group(1) + m.group(2)
        m = re.search(r'--- orbit c = (\d+), seed w = \(([^)]*)\); known solutions in it: (\[.*?\]) \(exponents (\[.*?\])\)',
                      line)
        if m:
            out[(key, int(m.group(1)))] = (tuple(int(v) for v in m.group(2).split(',')),
                                           ast.literal_eval(m.group(3)), ast.literal_eval(m.group(4)))
    return out


def safe_eval(node):
    """Evaluate a literal-like expression of core.py (constants, tuples, lists, dicts, dict(...), Fr(...))."""
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Tuple):
        return tuple(safe_eval(e) for e in node.elts)
    if isinstance(node, ast.List):
        return [safe_eval(e) for e in node.elts]
    if isinstance(node, ast.Dict):
        return {safe_eval(k): safe_eval(v) for k, v in zip(node.keys, node.values)}
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -safe_eval(node.operand)
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
        if node.func.id == 'dict' and not node.args:
            return {kw.arg: safe_eval(kw.value) for kw in node.keywords}
        if node.func.id in ('Fr', 'Fraction') and not node.keywords:
            return Fr(*[safe_eval(a) for a in node.args])
    raise ValueError('unsupported expression in core.py: ' + ast.dump(node)[:80])


def parse_core_data(path):
    tree = ast.parse(open(path, encoding='utf-8').read())
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if name in ('CLASSES', 'UNITS', 'KNOWN', 'RHS'):
                out[name] = safe_eval(node.value)
    return out


# ======================================================================================================================
# 7. Main
# ======================================================================================================================
print('verify_coverage.py (8604.a1) -- independent verification of the coverage step (exact; own arithmetic)')
print('No import of core.py, lifting2.py, 01_coverage.py, 02_exclusion.py; their files are only read as text at the end.')
print()

OUT01, SEEDS01 = parse_01(os.path.join(HERE, '01_coverage.out.txt'))
SUMMARY = {}
TOTAL_ORBITS = 0
STATS = {}

for key in DATA:
    tk = time.time()
    d = DATA[key]
    ref = OUT01[key]
    print('=' * 118)
    a0, a1, a2, a3, _ = d['F']
    print(f'CLASS {d["label"]}   G = {d["F"]}   g(t) = G(1,t) = t^4 + ({a3}) t^3 + ({a2}) t^2 + ({a1}) t + ({a0})')
    print('-- 1. the order O')
    K = Order(key)
    check(K.integral, 'multiplication matrices BT^-1 (sum_k B[i][k] C^k) BT of omega_0..omega_3 are integral (O is closed '
                      'under multiplication)')
    I4 = [[int(i == j) for j in range(4)] for i in range(4)]
    check(K.Mi[0] == I4, 'omega_0 = 1 (M_0 = identity)')
    unit = [tuple(int(i == j) for i in range(4)) for j in range(4)]
    comm = all(K.mul(unit[i], unit[j]) == K.mul(unit[j], unit[i]) for i in range(4) for j in range(4))
    assoc = all(K.mat(K.mul(unit[i], unit[j])) == matmul(K.Mi[i], K.Mi[j]) for i in range(4) for j in range(4))
    check(comm and assoc, 'structure constants commutative and associative (M(omega_i omega_j) = M_i M_j)')
    B = K.B
    shape = (B[0] == [1, 0, 0, 0] and B[1] == [0, 1, 0, 0]
             and B[2][2] * B[3][3] - B[2][3] * B[3][2] != 0)
    check(shape, 'omega_1 = beta and the (beta^2, beta^3)-block of omega_2, omega_3 is invertible (marker shape)')
    # symbolic norm form
    X = sympy.symbols('X0:4')
    tt = sympy.Symbol('t')
    Msym = sympy.zeros(4, 4)
    for i in range(4):
        Msym += X[i] * sympy.Matrix(K.Mi[i])
    Ndet = sympy.Poly(sympy.expand(Msym.det(method='berkowitz')), *X)
    xt = sum(X[i] * sum(sympy.Rational(B[i][k].numerator, B[i][k].denominator) * tt ** k for k in range(4))
             for i in range(4))
    gt = sum(sympy.Integer(int(K.g[k])) * tt ** k for k in range(5))
    Nres = sympy.Poly(sympy.expand(sympy.resultant(gt, xt, tt)), *X)
    check(Ndet == Nres, f'norm form det(sum X_i M_i) = Res_t(g(t), x(t)) as polynomials ({len(Ndet.terms())} terms)')
    check(all(c.is_integer for _, c in Ndet.terms()), 'the norm form has integer coefficients')
    mm, nn = sympy.symbols('m n')
    idG = sympy.expand(Ndet.as_expr().subs({X[0]: nn, X[1]: -mm, X[2]: 0, X[3]: 0})
                       - sum(int(d['F'][k]) * mm ** (4 - k) * nn ** k for k in range(5)))
    check(idG == 0, 'N(n - m beta) = G(m,n) identically')
    Dkl = {}
    for (e0, e1, e2, e3), c in Ndet.terms():
        Dkl.setdefault((e0, e1), []).append((int(c), e2, e3))
    check(Dkl.get((4, 0)) == [(1, 0, 0)] and all(k[0] <= 4 for k in Dkl),
          'the norm form is a monic quartic in X0 (coefficient of X0^4 is 1)')
    disc = det4([[sum(matmul(K.Mi[i], K.Mi[j])[k][k] for k in range(4)) for j in range(4)] for i in range(4)])
    print(f'   disc(O) = det(Tr(omega_i omega_j)) = {disc} = {sympy.factorint(disc)}')
    gs = sympy.Poly(gt, tt)
    irr = len(sympy.factor_list(gs)[1]) == 1 and sympy.factor_list(gs)[1][0][1] == 1
    check(irr, 'g is irreducible over Q (sympy factor_list)')
    # p-maximality (information only; maximality is not used by the coverage step)
    pmax = {}
    for p in sorted(q for q, e in sympy.factorint(abs(disc)).items() if e >= 2):
        integral_found = []
        for a in itertools.product(range(p), repeat=4):
            if not any(a):
                continue
            cs = charpoly_coeffs([[Fr(v, p) for v in row] for row in K.mat(a)])
            if all(cv.denominator == 1 for cv in cs):
                integral_found.append(a)
        pmax[p] = integral_found
    print('   p-maximality at the primes p with p^2 | disc(O) (info; not used): '
          + '; '.join(f'p = {p}: ' + ('maximal' if not v else f'NOT maximal, e.g. ({v[0]})/p is integral')
                      for p, v in pmax.items()))
    # units
    e1, e2 = d['UNITS']
    one = (1, 0, 0, 0)
    check(K.norm(e1) == 1 and K.norm(e2) == 1, f'N(eps1) = N(eps2) = 1 (Leibniz determinant), eps1 = {e1}, eps2 = {e2}')
    i1, i2 = K.quot_in_O(e1, one), K.quot_in_O(e2, one)
    check(i1 is not None and i2 is not None and K.mul(e1, i1) == one and K.mul(e2, i2) == one,
          f'eps1^-1 = {i1}, eps2^-1 = {i2} lie in O (Cramer)')
    PW = {}

    def epow(a, b):
        if (a, b) not in PW:
            x = one
            for _ in range(abs(a)):
                x = K.mul(x, e1 if a > 0 else i1)
            for _ in range(abs(b)):
                x = K.mul(x, e2 if b > 0 else i2)
            PW[(a, b)] = x
        return PW[(a, b)]

    # ---------------------------------------------------------------------------------------------------- embeddings
    print('-- 2. embeddings (rigorous)')
    Em = Emb(K)
    cnt_sympy = sympy.Poly(gt, tt).count_roots()
    check(Em.nreal == 2 and len(Em.chain[-1]) == 1 and Em.chain[-1][0] != 0,
          f'own Sturm chain (degrees {[len(p) - 1 for p in Em.chain]}): exactly 2 real roots; last term a nonzero constant '
          f'(g squarefree);  sympy count_roots = {cnt_sympy}')
    check(cnt_sympy == 2, 'sympy count_roots agrees (2 real roots), so the other two roots are a complex-conjugate pair')
    mpmath.mp.dps = 60
    for v, (lo, hi, k, sc) in enumerate(Em.rr):
        print(f'      sigma_{v + 1} = {mpmath.nstr(mpmath.mpf((lo + hi).numerator) / (lo + hi).denominator / 2, 45)} '
              f'+- {float((hi - lo) / 2):.1e}  (rational endpoints, exact sign change of g)')
    check(len(Em.rr) == 2 and all(k == 1 and sc for (_, _, k, sc) in Em.rr) and Em.rr[0][1] < Em.rr[1][0],
          'two disjoint isolating intervals, each with Sturm count 1 and an exact sign change of g')
    x0, y0 = Em.z0
    mpmath.mp.dps = 60
    print(f'      sigma_3 = beta_3 with |beta_3 - z0| < 2^-300, z0 = {mpmath.nstr(mpmath.mpf(x0.numerator) / x0.denominator, 45)} '
          f'+ {mpmath.nstr(mpmath.mpf(y0.numerator) / y0.denominator, 45)} i')
    check(Em.rouche, f'Rouche certificate (exactly one root in the disk, Im z0 > radius; lhs/rhs = {float(Em.rouche_ratio):.1e})')
    # core's S, P route, as a cross-check of the formula
    S_iv = IV(-K.g[3]) - Em.r1 - Em.r2
    P_iv = IV(K.g[0]) * (Em.r1 * Em.r2).recip()
    check(S_iv.meets(2 * Em.zr) and P_iv.meets(Em.zr.sq() + Em.zi.sq()),
          "core's S = -a3 - sigma1 - sigma2 and P = a0/(sigma1 sigma2) agree with 2 Re beta_3 and |beta_3|^2 (Rouche disk)")
    Gm = Em.gram()
    detG = isum(sg * Gm[0][p[0]] * Gm[1][p[1]] * Gm[2][p[2]] * Gm[3][p[3]] for p, sg in PERMS)
    check(detG.contains(Fr(abs(disc), 4)),
          f'det of the interval Gram matrix of E contains |disc(O)|/4 = {Fr(abs(disc), 4)} (the four embeddings are right)')

    # ---------------------------------------------------------------------------------------------------- f-vectors
    print('-- 3. size factors f_v(u) = s_v(u), u in U = {eps1^a eps2^b : a,b in {-1,0,1}, (a,b) != (0,0)}')
    cons, fvals = [], {}
    fmatch, spmatch, normok = True, True, True
    for (a, b) in EXPS:
        u = epow(a, b)
        normok &= K.norm(u) == 1
        f = Em.svec(u)
        pr = f[0] * f[1] * f[2].sq()
        normok &= pr.contains(1)
        s3sp, _, _ = s3_by_SP(Em, K, u)
        spmatch &= s3sp.meets(f[2])
        fvals[(a, b)] = f
        fp = [round_up(f[v].hi) for v in range(3)]
        cons.append(tuple(fp[v] - 1 for v in range(3)))
        rf = ref['f'][(a, b)]
        rel = max(abs(float(f[v]) / rf[v] - 1) for v in range(3))
        fmatch &= rel < 5e-6
        print(f'      (a,b)=({a:2d},{b:2d})  u = {u}')
        print(f'          f = ({mpmath.nstr(mpmath.mpf(f[0].lo.numerator) / f[0].lo.denominator, 15)}, '
              f'{mpmath.nstr(mpmath.mpf(f[1].lo.numerator) / f[1].lo.denominator, 15)}, '
              f'{mpmath.nstr(mpmath.mpf(f[2].lo.numerator) / f[2].lo.denominator, 15)})  widths < '
              f'{float(max(x.width() for x in f)):.1e};  01 printed {rf}')
    check(normok, 'every u in U has N(u) = 1 (determinant) and 1 lies in the interval f1 f2 f3^2')
    check(spmatch, "core's S,P formula A^2 + S A B + P B^2 for |sigma_3(u)|^2 meets the direct Re^2 + Im^2 enclosure for all u in U")
    check(fmatch, "the f-vectors agree with 01's printed values (6 significant digits)")
    # independence of eps1, eps2 directly: sign patterns of f(eps) - 1
    sg = {}
    for ab in ((1, 0), (0, 1)):
        f = fvals[ab]
        sg[ab] = tuple(1 if f[v].lo > 1 else (-1 if f[v].hi < 1 else 0) for v in range(3))
    indep = (0 not in sg[(1, 0)] and 0 not in sg[(0, 1)] and sg[(1, 0)] != sg[(0, 1)]
             and sg[(1, 0)] != tuple(-v for v in sg[(0, 1)]))
    check(indep, f'sign patterns of log f(eps1) = {sg[(1, 0)]} and log f(eps2) = {sg[(0, 1)]} are neither equal nor '
                 'opposite, so l(eps1), l(eps2) are independent: eps1, eps2 multiplicatively independent, -1 not in <eps1,eps2>')

    # ---------------------------------------------------------------------------------------------------- polygon
    print('-- 4. the relaxed polygon P+ and the bounds (f_v^+ = upper endpoint rounded up to 2^-96)')
    verts = clip_simplex(cons)
    bfs = basic_feasible(cons)
    check(len(verts) > 0 and set(verts) == bfs,
          f'half-plane clipping gives {len(verts)} vertices = the basic feasible solutions of all pairs of the '
          f'{len(cons) + 3} constraints (01 printed {ref["nverts"]} vertices)')
    allpos = all(all(v > 0 for v in P) for P in verts)
    check(allpos, 'every vertex of P+ has all three coordinates > 0 (exact)')
    for P in verts:
        print(f'      vertex ({float(P[0]):.6e}, {float(P[1]):.6e}, {float(P[2]):.6e})   r1 r2 r3^2 = '
              f'{float(P[0] * P[1] * P[2] ** 2):.6e}')
    vm = len(ref['verts']) == len(verts) and all(
        any(all(abs(float(P[k]) - R[k]) <= 1.01e-6 for k in range(3)) for R in ref['verts']) for P in verts)
    check(vm, "the vertices agree with 01's printed vertices (6 decimals), as sets")
    bmin = min(P[0] * P[1] * P[2] ** 2 for P in verts)
    kap = [min(P[v] for P in verts) for v in range(3)]
    kprod = kap[0] * kap[1] * kap[2] ** 2
    print(f'   b_min = {float(bmin):.10e}   (01: {ref["bmin"]:.6e});   kappa = ({float(kap[0]):.6e}, {float(kap[1]):.6e}, '
          f'{float(kap[2]):.6e}),  kappa1 kappa2 kappa3^2 = {float(kprod):.6e}')
    check(abs(float(bmin) / ref['bmin'] - 1) < 1e-6, "b_min agrees with 01's printed value")
    check(kprod <= bmin, 'kappa1 kappa2 kappa3^2 <= b_min (the linear bound is the weaker one)')
    Tprod = {c: least_T(Fr(c * c) / bmin) for c in RHS}
    Tlin = {c: least_T(Fr(c * c) / kprod) for c in RHS}
    print(f'   product bound T(c)  = least T with T^4 >= c^2/b_min:                   {Tprod}')
    print(f'   linear bound T_lin(c) = least T with T^4 >= c^2/(kappa1 kappa2 kappa3^2): {Tlin}')
    check(Tprod == ref['T'], f"the product bounds T(c) equal 01's T values {ref['T']}")
    # negative control: a dependent set U' = {eps1^a : a = +-1, +-2} must violate the positivity hypothesis
    cons_dep = []
    for a in (-2, -1, 1, 2):
        f = Em.svec(epow(a, 0))
        cons_dep.append(tuple(round_up(f[v].hi) - 1 for v in range(3)))
    vdep = clip_simplex(cons_dep)
    check(len(vdep) > 0 and not all(all(v > 0 for v in P) for P in vdep),
          'negative control: for the dependent set U\' = {eps1^a : a = +-1, +-2} the polygon has a vertex with a zero '
          'coordinate (the positivity test is not vacuous)')

    # ---------------------------------------------------------------------------------------------------- lower form
    print('-- 5. rational lower form and exact enumeration')
    Aint, eps, Gh = lower_form(Gm)
    En = Enumerator(Aint)
    check(True, f'A = Ghat - eps I, eps = {float(eps):.3e} >= Frobenius norm of the Gram enclosure error; A positive '
                f'definite (exact LDL^T pivots > 0, A = U^T D U reconstructed exactly)')
    rnd = random.Random(20261015)
    okA = True
    for _ in range(100):
        x = tuple(rnd.randint(-60, 60) for _ in range(4))
        okA &= Fr(qform(Aint, x), 1 << KB) <= Em.E(x).lo
    check(okA, 'x^T A x <= E(x) (rigorous interval) on 100 random integer vectors (sanity check of the shift argument)')
    # validation of the enumerator against a brute-force box for a small bound
    Tb = 120
    pts_enum = set(En.points(Tb << KB))
    Ai = matinv([[Fr(v, 1 << KB) for v in row] for row in Aint])
    box = [isqrt(int(Tb * Ai[i][i]) + 1) + 1 for i in range(4)]       # |x_i| <= sqrt(T (A^-1)_ii) on the ellipsoid
    pts_box = set()
    for x in itertools.product(*[range(-b, b + 1) for b in box]):
        if qform(Aint, x) <= Tb << KB:
            pts_box.add(x)
    check(pts_enum == pts_box, f'enumerator = brute force over the box |x_i| <= {box} for x^T A x <= {Tb} '
                               f'({len(pts_enum)} points)')

    def E_le(x, T):
        """Certified decision of E(x) <= T: True / False, or None if undecided."""
        q = Fr(qform(Aint, x), 1 << KB)
        if q > T:
            return False
        if q + 2 * eps * sum(v * v for v in x) <= T:
            return True
        e = Em.E(x)
        return True if e.hi <= T else (False if e.lo > T else None)

    # lattice-point counts at 01's bounds
    Tsmall = max(Tprod.values())
    small = [x for x in En.points(Tsmall << KB) if any(x)]
    cnt_ok = True
    undec_total = 0
    for c in RHS:
        dec = [E_le(x, Tprod[c]) for x in small]
        undec = sum(1 for v in dec if v is None)
        undec_total += undec
        n_le = sum(1 for v in dec if v)
        n_nc = sum(1 for x, v in zip(small, dec) if v and K.norm(x) == c)
        cnt_ok &= (n_le == ref['pts'][c] and n_nc == ref['normc'][c] and undec == 0)
        print(f'      T(c) = {Tprod[c]:5d}:  nonzero x with E(x) <= T certified: {n_le} (undecided {undec});  of norm {c}: '
              f'{n_nc}   [01: {ref["pts"][c]} points, {ref["normc"][c]} of norm c]')
    check(cnt_ok, "same numbers of lattice points (E <= T(c)) and of norm-c points as 01 at each of its bounds")

    # the main enumeration at the linear bound
    vol = lambda T: math.pi ** 2 / 2 * T * T / math.sqrt(float(disc if disc > 0 else -disc) / 4)
    Tbig = max(Tlin.values())
    est = vol(Tbig)
    print(f'   feasibility: estimated number of lattice points with E <= T_lin(36) = {Tbig}: '
          f'(pi^2/2) T^2 / sqrt(det G) = {est:.3e}  -> enumerate to the LINEAR bound T_lin(c) for every c')
    te = time.time()
    hits, nnodes, npts, nslow = En.scan(Tbig << KB, Dkl, (-1, 1, 4, 9, 36))
    print(f'   enumerated every x with x^T A x <= {Tbig}: {npts} lattice points (incl. 0), {nnodes} (x1,x2,x3)-nodes '
          f'({nslow} long nodes in Python integers), {len(hits)} with N(x) in {{-1,1,4,9,36}}  [{time.time() - te:.1f} s]')
    detok = all(K.norm(x) == v for v, x in hits)
    check(detok, 'every kept point has that norm by the Leibniz determinant of its multiplication matrix')
    rs = random.Random(7)
    samp_bad, samp_n = 0, 0
    for _ in range(3000):                         # the fast norm evaluation on random lattice points, kept or not
        x = tuple(rs.randint(-400, 400) for _ in range(4))
        p2, p3 = (1, x[2], x[2] ** 2, x[2] ** 3, x[2] ** 4), (1, x[3], x[3] ** 2, x[3] ** 3, x[3] ** 4)
        val = sum(c * x[0] ** k * x[1] ** l * p2[a] * p3[b] for (k, l), lst in Dkl.items() for (c, a, b) in lst)
        samp_n += 1
        samp_bad += val != K.norm(x)
    check(samp_bad == 0, f'the coefficient form of the norm agrees with the determinant on {samp_n} random points')
    STATS[key] = (Tbig, npts, nnodes, round(time.time() - te, 1))
    Q = {x: qform(Aint, x) for _, x in hits}
    H = {c: sorted({x for v, x in hits if v == c}, key=lambda y: (Q[y], y)) for c in (-1, 1, 4, 9, 36)}

    # ---------------------------------------------------------------------------------------------------- units
    print('-- 6. the unit group (Remark 3.2(iii)): c = 1 and c = -1 with the linear bound T_lin(1)')
    l1, l2 = Em.flog(e1), Em.flog(e2)
    det_l = l1[0] * l2[1] - l1[1] * l2[0]

    def ident(q):
        """(sign, a, b) with q = sign * eps1^a eps2^b, verified exactly; None if no such representation is found."""
        lq = Em.flog(q)
        ga = (lq[0] * l2[1] - lq[1] * l2[0]) / det_l
        gb = (l1[0] * lq[1] - l1[1] * lq[0]) / det_l
        a0_, b0_ = int(round(ga)), int(round(gb))
        for da in (0, -1, 1):
            for db in (0, -1, 1):
                y = epow(a0_ + da, b0_ + db)
                if y == q:
                    return (1, a0_ + da, b0_ + db)
                if y == tuple(-v for v in q):
                    return (-1, a0_ + da, b0_ + db)
        return None

    U1 = [x for x in H[1] if Q[x] <= Tlin[1] << KB]
    Um1 = [x for x in H[-1] if Q[x] <= Tlin[1] << KB]
    ids = [ident(x) for x in U1]
    check(all(r is not None for r in ids),
          f'all {len(U1)} elements of norm 1 with x^T A x <= T_lin(1) = {Tlin[1]} are +-eps1^a eps2^b (exact)')
    exps = sorted(r for r in ids if r is not None)
    print(f'      their exponents (sign, a, b): {exps}')
    check(len(set(exps)) == len(exps), 'they are pairwise distinct as (sign, a, b)')
    check(len(Um1) == 0, f'no element of norm -1 with x^T A x <= T_lin(1) (found {len(Um1)})')
    print('      => by Lemma 3.2 with the linear bound (c = 1 and c = -1): every unit of O is +-eps1^a eps2^b, and with the')
    print('         independence above O^x = {+-1} x <eps1> x <eps2> (eps1, eps2 fundamental; not used by the proof)')
    check(ref['unitline'] == (True, 0), "01 printed the same unit-group by-product (one orbit containing +-1, no norm -1)")

    # ---------------------------------------------------------------------------------------------------- orbits
    print('-- 7. orbits of <-1, eps1, eps2> on the elements of norm c, with the linear bound')
    SUMMARY[key] = {}
    for c in RHS:
        Lc = [x for x in H[c] if Q[x] <= Tlin[c] << KB]
        reps, members, links = [], [], []
        merged_ok = True
        for x in Lc:
            for k, w in enumerate(reps):
                qt = K.quot_in_O(w, x)
                if qt is not None:
                    r = ident(qt)
                    if r is None:
                        merged_ok = False          # a unit quotient outside <-1,eps1,eps2>: keep separate
                        continue
                    members[k].append(x)
                    links.append((x, k, r))
                    break
            else:
                reps.append(x)
                members.append([x])
        TOTAL_ORBITS += len(reps)
        print(f'   c = {c:2d}: T_lin(c) = {Tlin[c]}; elements of norm {c} with x^T A x <= T_lin(c): {len(Lc)};  '
              f'orbits: {len(reps)}   (01: {ref["norbits"][c]} orbit(s) at T(c) = {Tprod[c]})')
        check(merged_ok, f'c = {c}: every merge is by an exact quotient identified as +-eps1^a eps2^b')
        check(len(reps) == ref['norbits'][c], f'c = {c}: orbit count {len(reps)} agrees with 01')
        # 01's seeds
        seed_orb = []
        for (w, nmem) in ref['seeds'][c]:
            k = None
            if w in Lc:
                k = next((kk for kk, mem in enumerate(members) if w in mem), None)
            seed_orb.append(k)
            check(K.norm(w) == c and k is not None, f'c = {c}: 01 seed w = {w} has N(w) = {c} and lies in orbit #{k} here')
        check(sorted(s for s in seed_orb if s is not None) == list(range(len(reps))),
              f"c = {c}: 01's seeds lie in pairwise distinct orbits and meet every orbit found here")
        # the members within 01's bound T(c)
        inT = [x for x in Lc if E_le(x, Tprod[c])]
        print(f'      members with E <= T(c) = {Tprod[c]} (certified): {len(inT)}; 01 listed '
              f'{sum(nm for _, nm in ref["seeds"][c])} members')
        check(len(inT) == sum(nm for _, nm in ref['seeds'][c]), f'c = {c}: same number of listed members as 01')
        # wider check: every element of norm c with x^T A x <= T_lin(36) lies in a seed orbit
        wide = [x for x in H[c] if x not in Lc]
        wide_ok = all(any(K.quot_in_O(s, x) is not None and ident(K.quot_in_O(s, x)) is not None
                          for (s, _) in ref['seeds'][c]) for x in wide)
        check(wide_ok, f'c = {c}: the {len(wide)} further elements of norm {c} with T_lin({c}) < x^T A x <= {Tbig} also lie '
                       f'in the seed orbits (exact quotients)')
        for k, w in enumerate(reps):
            print(f'      orbit #{k}: {len(members[k])} members; least x^T A x member {w}  E = {float(Em.E(w)):.4f}')
            sk = next((s for (s, _) in ref['seeds'][c] if s in members[k]), None)
            if sk is not None:
                rel = sorted(ident(K.quot_in_O(sk, x)) for x in members[k])
                print(f'         members = (sign, a, b) * eps1^a eps2^b * seed {sk}: {rel}')
        # known solutions: exponents relative to 01's seed, compared with 02's printed exponents
        SUMMARY[key][c] = [w for (w, _) in ref['seeds'][c]]
        sols = [(m_, n_) for (m_, n_) in d['KNOWN'] if Gform(key, m_, n_) == c]
        check(sorted(sols) == sorted(CLAIMED[key][c]), f'c = {c}: KNOWN solutions {sols} = the claimed list')
        for (m_, n_) in sols:
            z = (n_, -m_, 0, 0)
            w = ref['seeds'][c][0][0]
            qt = K.quot_in_O(w, z)
            r = ident(qt) if qt is not None else None
            print(f'      known (m,n) = ({m_},{n_}): z = {z} = (sign, a, b) = {r} times eps1^a eps2^b * seed {w}')
            check(K.norm(z) == c and r is not None, f'c = {c}: known solution ({m_},{n_}) lies in the seed orbit')
        # stress test of Lemma 3.2: a U-local minimum w (E(uw) >= E(w) for all u in U) has s(w)/E(w) in P+, hence
        # E(w) <= T(c).  Greedy descent by U from random elements eps1^a eps2^b * seed (moves decided in floating
        # point on log sizes, which has no cancellation); the end point is then rebuilt exactly and E <= T(c) certified.
        w = ref['seeds'][c][0][0]
        lsw = [math.log(float(s)) for s in Em.svec(w)]
        lf1 = [math.log(float(v)) for v in fvals[(1, 0)]]
        lf2 = [math.log(float(v)) for v in fvals[(0, 1)]]
        Ef = lambda a, b: sum(math.exp(lsw[v] + a * lf1[v] + b * lf2[v]) for v in range(3))
        rg = random.Random(1000 * int(key) + c)
        ends, desc_ok = set(), True
        for _ in range(40):
            a, b = rg.randint(-12, 12), rg.randint(-12, 12)
            while True:
                cur = Ef(a, b)
                nb = min(((Ef(a + da, b + db), da, db) for (da, db) in EXPS))
                if nb[0] < cur * (1 - 1e-12):
                    a, b = a + nb[1], b + nb[2]
                else:
                    break
            ends.add((a, b))
        for (a, b) in sorted(ends):
            x = K.mul(epow(a, b), w)
            desc_ok &= (E_le(x, Tprod[c]) is True) and (x in Lc or tuple(-v for v in x) in Lc)
        check(desc_ok, f'c = {c}: greedy descent by U from 40 random orbit elements ends at {sorted(ends)} '
                       f'(exponents relative to the seed); every end point has E <= T(c) = {Tprod[c]} (certified) and is listed')
    print(f'   [class done in {time.time() - tk:.1f} s]')

# ======================================================================================================================
# 8. Comparisons with the files of the folder (read as text)
# ======================================================================================================================
print()
print('=' * 118)
print('COMPARISON WITH core.py (data, via ast; nothing executed), 01_coverage.out.txt and 02_exclusion.py/.out.txt')
CORE = parse_core_data(os.path.join(HERE, 'core.py'))
for key in DATA:
    cc = CORE['CLASSES'][key]
    same = (tuple(cc['F']) == DATA[key]['F']
            and [[Fr(v) for v in row] for row in cc['BASIS']] == [[Fr(v) for v in row] for row in DATA[key]['BASIS']]
            and tuple(cc['RHS']) == RHS
            and tuple(map(tuple, CORE['UNITS'][key])) == DATA[key]['UNITS']
            and list(CORE['KNOWN'][key]) == DATA[key]['KNOWN'])
    check(same, f'class {DATA[key]["label"]}: F, BASIS, RHS, UNITS, KNOWN retyped here are identical to core.py')
check(tuple(CORE['RHS']) == RHS, 'RHS identical to core.py')
check(SEEDS01 == SUMMARY, "the SEEDS block printed by 01 equals its per-c seed lines, each placed above in its own orbit")
SEEDS02 = parse_02_seeds(os.path.join(HERE, '02_exclusion.py'))
check(SEEDS02 == SEEDS01, "02_exclusion.py's SEEDS table is the one printed by 01")
check(all(OUT01[k].get('notes', 0) == 0 for k in OUT01), '01 printed no NOTE about unit quotients outside <-1,eps1,eps2>')
EXP02 = parse_02_exponents(os.path.join(HERE, '02_exclusion.out.txt'))
ok02 = len(EXP02) == 12
for key in DATA:
    K = Order(key)
    e1, e2 = DATA[key]['UNITS']
    for c in RHS:
        w, sols, exps = EXP02[(key, c)]
        ok02 &= [w] == SEEDS01[key][c]
        ok02 &= sorted(sols) == sorted(CLAIMED[key][c])
        for (m_, n_), (a, b) in zip(sols, exps):
            z = (n_, -m_, 0, 0)
            # z = +- eps1^a eps2^b w, checked exactly
            x = w
            inv = lambda e: K.quot_in_O(e, (1, 0, 0, 0))
            for _ in range(abs(a)):
                x = K.mul(x, e1 if a > 0 else inv(e1))
            for _ in range(abs(b)):
                x = K.mul(x, e2 if b > 0 else inv(e2))
            ok02 &= x == z or x == tuple(-v for v in z)
check(ok02, "02's printed seeds, known solutions and exponents: every known z equals +-eps1^a eps2^b w exactly")

print()
print('RUN STATISTICS (class: T_lin(36), lattice points enumerated, nodes, seconds):')
for key, st in STATS.items():
    print(f'   {DATA[key]["label"]}: {st}')
print(f'TOTAL ORBITS over the 12 pairs (class, c): {TOTAL_ORBITS}  (01: {sum(sum(OUT01[k]["norbits"].values()) for k in OUT01)})')
check(TOTAL_ORBITS == 12 and TOTAL_ORBITS == sum(sum(OUT01[k]['norbits'].values()) for k in OUT01), 'total orbit count 12')
print('SEEDS (01 seeds, one per orbit found here):')
for key in SUMMARY:
    print(f"    '{key}': " + '{' + ', '.join(f'{c}: {SUMMARY[key][c]}' for c in RHS) + '},')
print()
print(f'{NOK[0]} checks passed, {len(FAILS)} failed.')
if FAILS:
    print('CHECK FAILED:')
    for m_ in FAILS:
        print('   - ' + m_)
    print(f'[done in {time.time() - T0:.1f} s]')
    sys.exit(1)
print('ALL CHECKS PASSED: with the linear bound of Remark 3.2(ii) (no log-concavity) and an independent exact')
print('enumeration, every element of norm c in {1,4,9,36} lies in the orbit of exactly one 01 seed, for all three classes.')
print(f'[done in {time.time() - T0:.1f} s]')
