"""Exact arithmetic for the three quartic Thue fields of 8604.a1 (triangular = square-pyramidal).

Shared by the scripts of this folder. Nothing here prints; the numbered scripts check and report.
The generic part (Field, exact linear algebra, intervals, Embeddings) is that of
../../triangular_tetrahedral/scripts/core.py; only the data (CLASSES, KNOWN, UNITS, descent_maps) is new.

Conventions
-----------
* The 2-descent of 8604.a1 (y^2 = x^3 - 9x + 81, already monic: X = x, Y = y) gives, for the three non-trivial
  classes, binary quartic CLASS FORMS F_cl(m,n), and every integral point of the class gives coprime (m,n) with
  F_cl(m,n) in {1,4,9,36} (research/proofs/targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt).
* For each class the Thue equation actually solved is that of a FIELD FORM
      G(m,n) = a0 m^4 + a1 m^3 n + a2 m^2 n^2 + a3 m n^3 + a4 n^4   (a4 = 1),
  with g(t) := G(1,t) = t^4 + a3 t^3 + a2 t^2 + a1 t + a0 monic and beta a root of g, chosen so that the maximal
  order O satisfies O cap (Q + Q beta) = Z + Z beta.  For classes (0,1) and (1,1), G = F_cl and beta = alpha.
  For class (1,0), (1 + alpha)/2 is integral, so beta = (1 + alpha)/2 and F_cl(m,n) = G(2m, n + m)
  (to_class(key, m, n) maps a solution of G to one of F_cl, or to None).
* G(m,n) = m^4 g(n/m) = N(n - m beta), so a solution of G(m,n) = c is an element z = n - m beta of norm c whose
  beta^2- and beta^3-coordinates vanish.
* Elements of the order O are integer 4-tuples x = (x0,x1,x2,x3), the coordinates in the integral basis
  omega_0..omega_3 listed in BASIS (rows = coefficients of 1, beta, beta^2, beta^3).  In all three bases
  omega_0 = 1 and omega_1 = beta, and omega_2, omega_3 are the only basis elements involving beta^2, beta^3.
  Hence z = n - m beta has coordinates (n, -m, 0, 0), and the two MARKERS are the coordinates x2, x3.
* Embeddings: sigma_1 < sigma_2 are the two real roots of g, beta_3 the complex root with Im > 0.
  For an element x, s_1 = sigma_1(x)^2, s_2 = sigma_2(x)^2, s_3 = |sigma_3(x)|^2, and N(x)^2 = s_1 s_2 s_3^2.
  E(x) = s_1 + s_2 + s_3 is a positive definite quadratic form on O (x) R.
* Rigour: real numbers are handled only through intervals with rational (dyadic) endpoints (class Iv), built from
  exact sign changes of g at rationals.  Every comparison used in a proof is between interval endpoints.
"""
from fractions import Fraction as Fr
from math import gcd
import itertools

CURVE = '8604.a1'
CURVE_EQ = 'y^2 = x^3 - 9x + 81'
GA, GB = -9, 81                  # monic model Y^2 = X^3 + GA X + GB (here equal to the LMFDB model)

# ----------------------------------------------------------------------------------------------------------
# The three classes
# ----------------------------------------------------------------------------------------------------------
# key: class (i mod 2, j mod 2) of the 2-descent; mu = X0 - theta of the representative point.
# F: the field form G (Thue equation solved); CLASS_F: the class form of the 2-descent; RHS: right sides needed.
CLASSES = {
    '10': dict(label='(1,0)', mu='3 - theta', F=(4, -5, -3, -2, 1), CLASS_F=(9, -72, -18, 0, 1),
               SUB='F_cl(m,n) = G(2m, n + m), beta = (1 + alpha)/2',
               BASIS=[[1, 0, 0, 0], [0, 1, 0, 0], [Fr(1, 3), Fr(2, 3), Fr(1, 3), 0], [Fr(1, 3), 0, 0, Fr(1, 3)]],
               RHS=(1, 4, 9, 36)),
    '01': dict(label='(0,1)', mu='-3 - theta', F=(9, -72, 18, 0, 1), CLASS_F=(9, -72, 18, 0, 1),
               SUB='G = F_cl, beta = alpha',
               BASIS=[[1, 0, 0, 0], [0, 1, 0, 0], [Fr(1, 2), 0, Fr(1, 6), 0], [Fr(3, 4), Fr(3, 4), Fr(1, 12), Fr(1, 12)]],
               RHS=(1, 4, 9, 36)),
    '11': dict(label='(1,1)', mu='-theta', F=(36, -72, 0, 0, 1), CLASS_F=(36, -72, 0, 0, 1),
               SUB='G = F_cl, beta = alpha',
               BASIS=[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, Fr(1, 6), 0], [0, Fr(1, 2), 0, Fr(1, 12)]],
               RHS=(1, 4, 9, 36)),
}
RHS = (1, 4, 9, 36)               # union of the right sides


def to_class(key, m, n):
    """A solution (m,n) of G(m,n) = c -> the solution (m', n') of F_cl(m',n') = c it gives, or None."""
    if key == '10':
        if m % 2:
            return None
        return (m // 2, n - m // 2)
    return (m, n)


def from_class(key, m, n):
    """A solution (m,n) of F_cl(m,n) = c -> the solution of G = c (inverse of to_class)."""
    if key == '10':
        return (2 * m, n + m)
    return (m, n)


# Known coprime solutions (m,n) of G(m,n) in RHS, one of each pair +-(m,n), from a box search |m|,|n| <= 2000
# (re-checked in 00_fields.py).
KNOWN = {
    '10': [(0, 1), (1, -1), (1, 0), (2, 1)],
    '01': [(0, 1), (1, 0), (1, 3)],
    '11': [(0, 1), (1, 0), (1, 4), (2, 1)],
}


def Fval(key, m, n):
    """The field form G of the class."""
    a0, a1, a2, a3, a4 = CLASSES[key]['F']
    return a0 * m**4 + a1 * m**3 * n + a2 * m**2 * n**2 + a3 * m * n**3 + a4 * n**4


def Fclass(key, m, n):
    """The class form F_cl of the 2-descent."""
    a0, a1, a2, a3, a4 = CLASSES[key]['CLASS_F']
    return a0 * m**4 + a1 * m**3 * n + a2 * m**2 * n**2 + a3 * m * n**3 + a4 * n**4


def fcoeffs(key):
    """Coefficients of the monic g(t) = G(1,t), low degree first: [c0, c1, c2, c3, 1]."""
    a0, a1, a2, a3, a4 = CLASSES[key]['F']
    assert a4 == 1
    return [a0, a1, a2, a3, 1]


# ----------------------------------------------------------------------------------------------------------
# Polynomial arithmetic in Q[t]/(g) on power coordinates (length-4 lists of Fractions)
# ----------------------------------------------------------------------------------------------------------
def pmulmod(u, v, fc):
    """Product of two power-coordinate vectors modulo the monic quartic with coefficients fc (low first)."""
    prod = [Fr(0)] * 7
    for i in range(4):
        if u[i] == 0:
            continue
        for j in range(4):
            prod[i + j] += u[i] * v[j]
    for d in range(6, 3, -1):            # t^d = t^(d-4) * t^4 = -t^(d-4) * sum_{k<4} c_k t^k
        c = prod[d]
        if c:
            prod[d] = Fr(0)
            for k in range(4):
                prod[d - 4 + k] -= c * fc[k]
    return prod[:4]


class Field:
    """The order O with integral basis BASIS in Q(beta), g(beta) = 0."""

    def __init__(self, key):
        self.key = key
        self.label = CLASSES[key]['label']
        self.fc = fcoeffs(key)
        self.B = [[Fr(c) for c in row] for row in CLASSES[key]['BASIS']]   # basis -> power coords
        self.Binv = mat_inv(self.B)                                          # power coords -> basis coords
        # multiplication table: mt[i][j] = coordinates of omega_i * omega_j (must be integers)
        self.mt = [[None] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                pc = pmulmod(self.B[i], self.B[j], self.fc)
                co = vec_mat(pc, self.Binv)
                assert all(c.denominator == 1 for c in co), 'BASIS does not span an order'
                self.mt[i][j] = [int(c) for c in co]

    # --- conversions -----------------------------------------------------------------------------------
    def to_power(self, x):
        return [sum(Fr(x[i]) * self.B[i][j] for i in range(4)) for j in range(4)]

    def from_power(self, pc):
        co = vec_mat([Fr(c) for c in pc], self.Binv)
        return co

    def from_mn(self, m, n):
        """Coordinates of n - m*beta."""
        return (n, -m, 0, 0)

    # --- ring operations -------------------------------------------------------------------------------
    def mul(self, x, y):
        z = [0, 0, 0, 0]
        mt = self.mt
        for i in range(4):
            xi = x[i]
            if not xi:
                continue
            for j in range(4):
                yj = y[j]
                if not yj:
                    continue
                c = xi * yj
                r = mt[i][j]
                z[0] += c * r[0]; z[1] += c * r[1]; z[2] += c * r[2]; z[3] += c * r[3]
        return tuple(z)

    def mulmod(self, x, y, q):
        return tuple(v % q for v in self.mul(x, y))

    def pow(self, x, e):
        """x^e for e >= 0, or e < 0 when x is a unit (uses the adjugate)."""
        if e < 0:
            n = self.norm(x)
            assert abs(n) == 1, 'negative power of a non-unit'
            x = tuple(n * c for c in self.adj(x))     # x^{-1} = N(x) * x^# when N(x) = +-1
            e = -e
        r = (1, 0, 0, 0)
        b = tuple(x)
        while e:
            if e & 1:
                r = self.mul(r, b)
            b = self.mul(b, b)
            e >>= 1
        return r

    def powmod(self, x, e, q):
        r = (1, 0, 0, 0)
        b = tuple(v % q for v in x)
        while e:
            if e & 1:
                r = self.mulmod(r, b, q)
            b = self.mulmod(b, b, q)
            e >>= 1
        return r

    def matrix(self, x):
        """Integer matrix of multiplication by x: column j = coordinates of x * omega_j."""
        M = [[0] * 4 for _ in range(4)]
        for j in range(4):
            col = [0, 0, 0, 0]
            for i in range(4):
                if x[i]:
                    for k in range(4):
                        col[k] += x[i] * self.mt[i][j][k]
            for k in range(4):
                M[k][j] = col[k]
        return M

    def norm(self, x):
        return det_int(self.matrix(x))

    def adj(self, x):
        """x^# with x * x^# = N(x): the first column of adj(M(x)) (an integer vector)."""
        A = adj_int(self.matrix(x))
        return tuple(A[k][0] for k in range(4))

    def divides_exactly(self, x, y):
        """Return y/x as integer coordinates if it lies in O, else None (exact)."""
        n = self.norm(x)
        num = self.mul(y, self.adj(x))
        if all(c % n == 0 for c in num):
            return tuple(c // n for c in num)
        return None


# ----------------------------------------------------------------------------------------------------------
# Exact linear algebra
# ----------------------------------------------------------------------------------------------------------
def mat_inv(A):
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


def vec_mat(v, A):
    """Row vector v times matrix A."""
    return [sum(v[i] * A[i][j] for i in range(len(v))) for j in range(len(A[0]))]


def det_int(M):
    """Exact integer determinant (Bareiss)."""
    n = len(M)
    A = [list(map(int, row)) for row in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((r for r in range(k + 1, n) if A[r][k] != 0), None)
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def adj_int(M):
    """Adjugate of an integer matrix (cofactor transpose), exact."""
    n = len(M)
    C = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            minor = [[M[r][c] for c in range(n) if c != j] for r in range(n) if r != i]
            C[i][j] = (-1) ** (i + j) * det_int(minor)
    return [[C[j][i] for j in range(n)] for i in range(n)]


# ----------------------------------------------------------------------------------------------------------
# Rigorous real arithmetic: intervals with dyadic endpoints
# ----------------------------------------------------------------------------------------------------------
PREC = 2 ** 240            # endpoints are rounded outward to multiples of 1/PREC


def _down(q):
    q = Fr(q)
    return Fr((q.numerator * PREC) // q.denominator, PREC)


def _up(q):
    q = Fr(q)
    return Fr(-((-q.numerator * PREC) // q.denominator), PREC)


class Iv:
    """Closed interval [lo, hi] with rational endpoints; operations are outward rounded, hence rigorous."""
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        if hi is None:
            hi = lo
        lo, hi = Fr(lo), Fr(hi)
        assert lo <= hi
        self.lo, self.hi = _down(lo), _up(hi)

    @staticmethod
    def _c(x):
        return x if isinstance(x, Iv) else Iv(x)

    def __add__(self, o):
        o = Iv._c(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)
    __radd__ = __add__

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + (-Iv._c(o))

    def __rsub__(self, o):
        return Iv._c(o) - self

    def __mul__(self, o):
        o = Iv._c(o)
        p = [self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi]
        return Iv(min(p), max(p))
    __rmul__ = __mul__

    def inv(self):
        assert self.lo > 0 or self.hi < 0, 'division by an interval containing 0'
        return Iv(1 / self.hi, 1 / self.lo)

    def __truediv__(self, o):
        return self * Iv._c(o).inv()

    def __rtruediv__(self, o):
        return Iv._c(o) * self.inv()

    def sq(self):
        if self.lo >= 0:
            return Iv(self.lo ** 2, self.hi ** 2)
        if self.hi <= 0:
            return Iv(self.hi ** 2, self.lo ** 2)
        return Iv(0, max(self.lo ** 2, self.hi ** 2))

    def width(self):
        return self.hi - self.lo

    def mid(self):
        return (self.lo + self.hi) / 2

    def pos(self):
        return self.lo > 0

    def neg(self):
        return self.hi < 0

    def __float__(self):
        return float(self.mid())

    def __repr__(self):
        return f'[{float(self.lo):.12g}, {float(self.hi):.12g}]'


def _feval(fc, x):
    return sum(Fr(c) * x ** k for k, c in enumerate(fc))


def real_root_intervals(fc, approx, width=Fr(1, 2 ** 200)):
    """Isolating intervals for the real roots of the monic quartic fc near the float approximations `approx`.

    For each approximation r, find rationals a < r < b with g(a) g(b) < 0 (exact), then bisect exactly.
    The number of real roots is certified separately (Sturm) in Embeddings and in 00_fields.py.
    """
    out = []
    for r in approx:
        d = Fr(1, 10 ** 6)
        a, b = Fr(r) - d, Fr(r) + d
        for _ in range(60):
            if _feval(fc, a) * _feval(fc, b) < 0:
                break
            d *= 2
            a, b = Fr(r) - d, Fr(r) + d
        fa = _feval(fc, a)
        assert fa * _feval(fc, b) < 0, 'no sign change'
        while b - a > width:
            c = (a + b) / 2
            fcv = _feval(fc, c)
            if fcv == 0:
                a = b = c
                break
            if (fcv < 0) == (fa < 0):
                a, fa = c, fcv
            else:
                b = c
        out.append(Iv(a, b))
    return out


class Embeddings:
    """Rigorous embeddings of a Field: sigma_1 < sigma_2 real, and |sigma_3|^2 through the real quadratic factor
    t^2 - S t + P of g carrying the complex pair (S = -c3 - sigma_1 - sigma_2, P = c0/(sigma_1 sigma_2))."""

    def __init__(self, K):
        import numpy as np
        import sympy
        self.K = K
        fc = K.fc
        # Sturm's theorem (sympy count_roots is exact): g has exactly two real roots, so the remaining pair is
        # complex conjugate and t^2 - S t + P below is the real quadratic factor carrying it.
        tt = sympy.Symbol('t')
        assert sympy.Poly(sum(c * tt ** k for k, c in enumerate(fc)), tt).count_roots() == 2, 'signature is not (2,1)'
        rts = np.roots([float(c) for c in reversed(fc)])
        real = sorted(r.real for r in rts if abs(r.imag) < 1e-9)
        assert len(real) == 2, 'expected signature (2,1)'
        self.r1, self.r2 = real_root_intervals(fc, real)
        assert self.r1.hi < self.r2.lo
        self.S = Iv(-fc[3]) - self.r1 - self.r2
        self.P = Iv(fc[0]) / (self.r1 * self.r2)
        # values of the power basis at the real roots
        self.pw1 = [Iv(1), self.r1, self.r1.sq(), self.r1.sq() * self.r1]
        self.pw2 = [Iv(1), self.r2, self.r2.sq(), self.r2.sq() * self.r2]
        # linear forms (on integral-basis coordinates) for sigma_1, sigma_2, and A, B with sigma_3 = A + B beta_3
        B = K.B
        self.L1 = [sum((self.pw1[j] * B[i][j] for j in range(4)), Iv(0)) for i in range(4)]
        self.L2 = [sum((self.pw2[j] * B[i][j] for j in range(4)), Iv(0)) for i in range(4)]
        S, P = self.S, self.P
        self.LA, self.LB = [], []
        for i in range(4):
            c0, c1, c2, c3 = B[i]
            self.LA.append(Iv(c0) - P * c2 - S * P * c3)
            self.LB.append(Iv(c1) + S * c2 + (S.sq() - P) * c3)

    def sig1(self, x):
        return sum((self.L1[i] * x[i] for i in range(4) if x[i]), Iv(0))

    def sig2(self, x):
        return sum((self.L2[i] * x[i] for i in range(4) if x[i]), Iv(0))

    def s3(self, x):
        A = sum((self.LA[i] * x[i] for i in range(4) if x[i]), Iv(0))
        Bv = sum((self.LB[i] * x[i] for i in range(4) if x[i]), Iv(0))
        return A.sq() + self.S * A * Bv + self.P * Bv.sq()

    def svec(self, x):
        """(s_1, s_2, s_3) = (sigma_1^2, sigma_2^2, |sigma_3|^2) as intervals."""
        return (self.sig1(x).sq(), self.sig2(x).sq(), self.s3(x))

    def E(self, x):
        s = self.svec(x)
        return s[0] + s[1] + s[2]

    def gram(self):
        """Interval Gram matrix G with E(x) = x^T G x."""
        G = [[None] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                g = self.L1[i] * self.L1[j] + self.L2[i] * self.L2[j]
                g = g + self.LA[i] * self.LA[j] + (self.S * (self.LA[i] * self.LB[j] + self.LB[i] * self.LA[j])) / 2 \
                    + self.P * self.LB[i] * self.LB[j]
                G[i][j] = g
        return G


# ----------------------------------------------------------------------------------------------------------
# Units used for the orbits (found by a lattice search, verified in 00_fields.py)
# ----------------------------------------------------------------------------------------------------------
# Two units of norm +1 per field (coordinates in the integral basis), found by a lattice search (LLL on
# weighted embeddings, 2401 weightings per field) and Lagrange-reduced in log space.  Every unit met by that
# search lies in <-1, eps1, eps2>; whether the pair is fundamental is NOT used anywhere.  Their norms are re-checked,
# and their independence follows from the cone certificate of 01_coverage.py.
UNITS = {
    '10': ((0, 1, -1, 0), (0, -2, 1, 0)),
    '01': ((1, -1, 1, 0), (1, 0, -2, 0)),
    '11': ((0, 0, -1, 0), (8, -4, 3, 0)),
}


def descent_maps(key):
    """The maps of the 2-descent (weierstrass_2descent.out.txt, section CURVE 8604.a1) from a solution (m,n) of the
    CLASS form F_cl to the point: X = X4(m,n)/g^2, Y = +- Y6(m,n)/g^3, with g = gcd of the conic coordinates
    (u,v,w), and the conic map (u,v,w)."""
    if key == '10':
        uvw = lambda m, n: (3*m*m + n*n, 6*m*m, 6*m*n)
        X4 = lambda m, n: 351*m**4 + 108*m**3*n + 18*m*m*n*n + 36*m*n**3 + 3*n**4
        Y6 = lambda m, n: (-6561*m**6 - 3240*m**5*n - 405*m**4*n*n + 1620*m**3*n**3 + 405*m*m*n**4 + 36*m*n**5
                           + 9*n**6)
    elif key == '01':
        uvw = lambda m, n: (-3*m*m + n*n, 6*m*m - 6*m*n, 6*m*n)
        X4 = lambda m, n: 297*m**4 - 108*m**3*n + 18*m*m*n*n + 36*m*n**3 - 3*n**4
        Y6 = lambda m, n: (-5103*m**6 + 2592*m**5*n - 405*m**4*n*n + 1620*m**3*n**3 - 405*m*m*n**4 + 36*m*n**5
                           + 9*n**6)
    elif key == '11':
        uvw = lambda m, n: (-6*m*n + n*n, 6*m*m, 6*m*n)
        X4 = lambda m, n: 324*m**4 - 36*m*m*n*n + 36*m*n**3
        Y6 = lambda m, n: (-5832*m**6 + 648*m**5*n - 1620*m**4*n*n + 1620*m**3*n**3 - 18*m*n**5 + 9*n**6)
    else:
        raise KeyError(key)
    return uvw, X4, Y6
