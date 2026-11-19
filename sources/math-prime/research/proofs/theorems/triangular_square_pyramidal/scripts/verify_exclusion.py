#!/usr/bin/env python3
"""verify_exclusion.py (8604.a1, triangular = square-pyramidal) -- independent re-implementation of the
exclusion step 02_exclusion.py.

Independent verifier (Claude subagent), 2026-10-15.  Nothing is imported from core.py, lifting2.py,
01_coverage.py or 02_exclusion.py.  The data are copied as literals: field forms, integral bases, units and
known solutions from core.py; the seeds from 01_coverage.out.txt; every count and every Proposition 17
instance that 02_exclusion.out.txt prints.  All arithmetic is written afresh, by routes that differ from
02_exclusion.py / lifting2.py:

  * The order O through its regular representation: the matrix of beta in the omega-basis is
    B^-T C_g B^T (C_g = companion matrix of g), M(omega_i) = omega_i(M(beta)), and these must be integral
    (the basis spans a ring).  Cross-checked against sympy products omega_i omega_j mod g.  Norms both as
    determinants and as resultants Res(g, x(t)).
  * G_q and Lambda_q.  h1 = ord(eps1), h2 = ord(eps2) mod qO; ALL relations eps1^a eps2^b == 1 with
    0 <= a < h1, 0 <= b < h2 are read off a table of the powers of eps1 against the powers of eps2^-1, and
    Lambda_q is the Hermite form of (h1,0), (0,h2) and these relations.  CERTIFIED: every Hermite basis vector
    maps to 1 mod qO (modular exponentiation; negative exponents through the exact inverse), and the index
    equals |G_q| counted independently as the closure (BFS) of {1} under eps1, eps2 mod qO.
  * Intersections of lattices: the integer kernel of [L ; -M] by unimodular row reduction, CERTIFIED by
    containment and [Z^2 : L cap M][Z^2 : L + M] = [Z^2 : L][Z^2 : M]; and separately as Lambda_lcm(q_i)
    (the kernel modulo the lcm of the moduli), recomputed from scratch.  Routines self-tested on random lattices.
  * Proposition 17 with exact integers (no truncation) in four bases of Lambda_q (Hermite; 02's basis
    (h,0),(-j,k) of proof.md 4.1; a unimodular change; a Lagrange-reduced basis).  Hypotheses checked as stated
    in lifting.binomial.multi/proof.md sections 1.1, 1.5 and 8, with Lambda(2) computed from its definition, and
    the strict inequalities e + delta < Lambda(2) and e + lambda < Lambda(2) of the proof.  The proof's
    mechanism (stable reductions of f1 and f2' = f2 - rho f1) is recomputed at each instance up to degree 4.
  * The sieve of proof.md 4.4 (certified rule; sibling rule only at the q where Proposition 17 holds) in 02's
    order, compared count by count; a direct check of every class modulo the final lattice; the sieve again in
    reversed order and with other modulus sets (greedy over odd moduli only, and over powers of 2 only).
  * Negative controls: planted false certificates; h1 Z x h2 Z in place of Lambda_q; -1 and phantom classes;
    basis dependence of nu(J) where nu(J) >= e.
  * A box scan |a|,|b| <= 300 modulo two primes, then exact (a sanity check only; not part of the proof).
Every discrepancy prints CHECK FAILED, and the exit status is 1 if there is any.
Run from this folder:   python verify_exclusion.py > verify_exclusion.out.txt
"""
import sys
import time
import random
from fractions import Fraction as Fr
from math import gcd, factorial

import numpy as np
import sympy as sp

T0 = time.time()
FAILS = []
NCHECK = [0]
INF = 10 ** 9                      # stands for nu(0) = +infinity


def check(cond, msg):
    NCHECK[0] += 1
    if not cond:
        FAILS.append(msg)
        print('   CHECK FAILED: ' + msg)
    return bool(cond)


def fnu(v):
    return 'inf' if v >= INF else str(v)


# ======================================================================================================
# DATA (copied as literals)
# ======================================================================================================
# field forms G(m,n) = a0 m^4 + a1 m^3 n + a2 m^2 n^2 + a3 m n^3 + a4 n^4 ; g(t) = G(1,t) ; core.py CLASSES[*]['F']
FORMS = {'10': (4, -5, -3, -2, 1), '01': (9, -72, 18, 0, 1), '11': (36, -72, 0, 0, 1)}
LABEL = {'10': '(1,0)', '01': '(0,1)', '11': '(1,1)'}
# rows = omega_0..omega_3 in the power basis 1, beta, beta^2, beta^3 (core.py BASIS)
BASES = {
    '10': [[1, 0, 0, 0], [0, 1, 0, 0], [Fr(1, 3), Fr(2, 3), Fr(1, 3), 0], [Fr(1, 3), 0, 0, Fr(1, 3)]],
    '01': [[1, 0, 0, 0], [0, 1, 0, 0], [Fr(1, 2), 0, Fr(1, 6), 0], [Fr(3, 4), Fr(3, 4), Fr(1, 12), Fr(1, 12)]],
    '11': [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, Fr(1, 6), 0], [0, Fr(1, 2), 0, Fr(1, 12)]],
}
UNITS = {'10': ((0, 1, -1, 0), (0, -2, 1, 0)), '01': ((1, -1, 1, 0), (1, 0, -2, 0)),
         '11': ((0, 0, -1, 0), (8, -4, 3, 0))}
KNOWN = {'10': [(0, 1), (1, -1), (1, 0), (2, 1)], '01': [(0, 1), (1, 0), (1, 3)],
         '11': [(0, 1), (1, 0), (1, 4), (2, 1)]}
# the claimed solution sets (all integer solutions of G = c up to sign), as stated for this R1 run
CLAIMED = {
    '10': {1: [(0, 1)], 4: [(1, 0)], 9: [(1, -1), (2, 1)], 36: []},
    '01': {1: [(0, 1)], 4: [], 9: [(1, 0)], 36: [(1, 3)]},
    '11': {1: [(0, 1), (2, 1)], 4: [(1, 4)], 9: [], 36: [(1, 0)]},
}
# seeds, as printed by 01_coverage.out.txt (SEEDS = ...)
SEEDS_01 = {
    '10': {1: [(-1, 0, 0, 0)], 4: [(1, 0, 2, -1)], 9: [(1, -2, 3, -1)], 36: [(-2, -2, 3, -1)]},
    '01': {1: [(-1, 0, 0, 0)], 4: [(1, 1, -2, 0)], 9: [(0, -1, 0, 0)], 36: [(3, -1, 0, 0)]},
    '11': {1: [(-1, 0, 0, 0)], 4: [(4, -1, 0, 0)], 9: [(0, 2, -3, 0)], 36: [(0, -1, 0, 0)]},
}
# Everything 02_exclusion.out.txt prints per orbit.
#   steps: (q, candidate classes mod the new lattice, certified, sibling-empty, left)
#   cert:  (printed class, solution (m,n), q, lambda, nu(J), e)        ["class ... contains the solution ..."]
#   inst:  (q, p, e, solution (m,n), lambda, nu(J), delta, bound)      ["Prop. 17 instance used: ..."]
#   G:     q -> (|G_q|, candidate classes mod Lambda_q)
ORBITS = [
    dict(key='10', c=1, w=(-1, 0, 0, 0), known=[(0, 1)], exps=[(0, 0)], used=[8, 16], index=224,
         steps=[(8, 2, 1, 0, 1), (16, 0, 0, 0, 0)], cert=[((0, 0), (0, 1), 8, 0, 1, 3)],
         inst=[(8, 2, 3, (0, 1), 0, 1, 1, 1)], G={8: (56, 2), 16: (224, 2)}),
    dict(key='10', c=4, w=(1, 0, 2, -1), known=[(1, 0)], exps=[(1, 0)], used=[8], index=56,
         steps=[(8, 1, 1, 0, 0)], cert=[((1, 0), (1, 0), 8, 0, 0, 3)],
         inst=[(8, 2, 3, (1, 0), 0, 0, 0, 1)], G={8: (56, 1)}),
    dict(key='10', c=9, w=(1, -2, 3, -1), known=[(1, -1), (2, 1)], exps=[(1, -1), (2, 0)], used=[8, 16], index=224,
         steps=[(8, 3, 2, 0, 1), (16, 0, 0, 0, 0)],
         cert=[((7, 1), (1, -1), 8, 0, 0, 3), ((2, 0), (2, 1), 8, 0, 1, 3)],
         inst=[(8, 2, 3, (1, -1), 0, 0, 0, 1), (8, 2, 3, (2, 1), 0, 1, 1, 1)], G={8: (56, 3), 16: (224, 3)}),
    dict(key='10', c=36, w=(-2, -2, 3, -1), known=[], exps=[], used=[4], index=28,
         steps=[(4, 0, 0, 0, 0)], cert=[], inst=[], G={4: (28, 0)}),
    dict(key='01', c=1, w=(-1, 0, 0, 0), known=[(0, 1)], exps=[(0, 0)], used=[8, 16], index=128,
         steps=[(8, 2, 1, 0, 1), (16, 0, 0, 0, 0)], cert=[((0, 0), (0, 1), 8, 0, 1, 3)],
         inst=[(8, 2, 3, (0, 1), 0, 1, 1, 1)], G={8: (32, 2), 16: (128, 2)}),
    dict(key='01', c=4, w=(1, 1, -2, 0), known=[], exps=[], used=[4], index=8,
         steps=[(4, 0, 0, 0, 0)], cert=[], inst=[], G={4: (8, 0)}),
    dict(key='01', c=9, w=(0, -1, 0, 0), known=[(1, 0)], exps=[(0, 0)], used=[8, 16], index=128,
         steps=[(8, 2, 1, 0, 1), (16, 0, 0, 0, 0)], cert=[((0, 0), (1, 0), 8, 0, 1, 3)],
         inst=[(8, 2, 3, (1, 0), 0, 1, 1, 1)], G={8: (32, 2), 16: (128, 2)}),
    dict(key='01', c=36, w=(3, -1, 0, 0), known=[(1, 3)], exps=[(0, 0)], used=[43, 13], index=77616,
         steps=[(43, 3, 1, 0, 2), (13, 0, 0, 0, 0)], cert=[((0, 0), (1, 3), 43, 0, 0, 1)],
         inst=[(43, 43, 1, (1, 3), 0, 0, 0, 0)], G={43: (3234, 3), 13: (336, 4)}),
    dict(key='11', c=1, w=(-1, 0, 0, 0), known=[(0, 1), (2, 1)], exps=[(0, 0), (2, 0)], used=[4, 7, 23],
         index=11616, steps=[(4, 2, 0, 0, 2), (7, 4, 0, 0, 4), (23, 2, 2, 0, 0)],
         cert=[((0, 0), (0, 1), 23, 0, 0, 1), ((2, 0), (2, 1), 23, 0, 0, 1)],
         inst=[(23, 23, 1, (0, 1), 0, 0, 0, 0), (23, 23, 1, (2, 1), 0, 0, 0, 0)],
         G={4: (4, 2), 7: (24, 2), 23: (2904, 7)}),
    dict(key='11', c=4, w=(4, -1, 0, 0), known=[(1, 4)], exps=[(0, 0)], used=[4, 8, 7, 23], index=23232,
         steps=[(4, 2, 0, 0, 2), (8, 4, 0, 0, 4), (7, 4, 0, 0, 4), (23, 2, 1, 1, 0)],
         cert=[((0, 0), (1, 4), 23, 0, 0, 1)], inst=[(23, 23, 1, (1, 4), 0, 0, 0, 0)],
         G={4: (4, 2), 8: (16, 4), 7: (24, 2), 23: (2904, 7)}),
    dict(key='11', c=9, w=(0, 2, -3, 0), known=[], exps=[], used=[8], index=16,
         steps=[(8, 0, 0, 0, 0)], cert=[], inst=[], G={8: (16, 0)}),
    dict(key='11', c=36, w=(0, -1, 0, 0), known=[(1, 0)], exps=[(0, 0)], used=[49], index=168,
         steps=[(49, 1, 1, 0, 0)], cert=[((0, 0), (1, 0), 49, 0, 0, 2)], inst=[(49, 7, 2, (1, 0), 0, 0, 0, 1)],
         G={49: (168, 1)}),
]
# "moduli available (|G_q|)" lines of 02_exclusion.out.txt (the moduli with |G_q| <= 5000)
AVAILABLE = {
    '10': {4: 28, 8: 56, 16: 224, 32: 896, 64: 3584, 3: 72, 9: 216, 27: 1944, 5: 124, 25: 3100, 7: 200, 11: 665,
           13: 1190, 17: 4912},
    '01': {4: 8, 8: 32, 16: 128, 32: 512, 64: 2048, 3: 36, 9: 108, 27: 972, 5: 124, 25: 3100, 7: 200, 11: 1330,
           13: 336, 17: 4912, 19: 1620, 43: 3234},
    '11': {4: 4, 8: 16, 16: 64, 32: 256, 64: 1024, 3: 12, 9: 108, 27: 972, 5: 124, 25: 3100, 7: 24, 49: 168,
           11: 1330, 13: 1190, 17: 4912, 19: 3620, 23: 2904, 31: 4965},
}
NEG02 = {'10': 8, '01': 4, '11': 4}        # "checked: -1 is not a residue of <eps1,eps2> mod q"
MODULI_02 = [4, 8, 16, 32, 64, 3, 9, 27, 81, 5, 25, 7, 49, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
MAXG = 5000
BOX_R = 300
BOX_PRIMES = (1000000007, 998244353)


# ======================================================================================================
# exact linear algebra over Q
# ======================================================================================================
def mmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def transp(A):
    return [list(r) for r in zip(*A)]


def q_inverse(A):
    """Inverse of a square rational matrix (Gauss-Jordan over Fractions)."""
    n = len(A)
    M = [[Fr(A[i][j]) for j in range(n)] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [row[n:] for row in M]


def det_mod(M, P):
    """Determinant of an integer matrix modulo the prime P."""
    n = len(M)
    A = [[v % P for v in row] for row in M]
    d = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c]), None)
        if piv is None:
            return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            d = -d
        d = d * A[c][c] % P
        iv = pow(A[c][c], -1, P)
        for r in range(c + 1, n):
            if A[r][c]:
                f = A[r][c] * iv % P
                A[r] = [(a - f * b) % P for a, b in zip(A[r], A[c])]
    return d % P


# ======================================================================================================
# The order O by its regular representation
# ======================================================================================================
class Order:
    def __init__(self, key):
        self.key = key
        a0, a1, a2, a3, a4 = FORMS[key]
        assert a4 == 1
        self.gc = [a0, a1, a2, a3, 1]                        # g(t) = sum gc[k] t^k, monic
        self.Bm = [[Fr(c) for c in row] for row in BASES[key]]
        # companion matrix: column j = power coordinates of beta * beta^j
        Cg = [[0] * 4 for _ in range(4)]
        for j in range(3):
            Cg[j + 1][j] = 1
        for k in range(4):
            Cg[k][3] = -self.gc[k]
        Bt = transp(self.Bm)                                  # power coords = Bt * omega coords
        self.Bti = q_inverse(Bt)
        Cp = [[[Fr(int(i == j)) for j in range(4)] for i in range(4)]]
        for _ in range(3):
            Cp.append(mmul(Cp[-1], Cg))
        self.Mom = []
        self.integral = True
        for i in range(4):
            Mp = [[sum(self.Bm[i][j] * Cp[j][r][c] for j in range(4)) for c in range(4)] for r in range(4)]
            Mo = mmul(mmul(self.Bti, Mp), Bt)
            if not all(x.denominator == 1 for row in Mo for x in row):
                self.integral = False
            self.Mom.append([[int(x) for x in row] for row in Mo])
        # structure constants: coordinates of omega_i omega_j = column j of M(omega_i)
        self.C = [[tuple(self.Mom[i][k][j] for k in range(4)) for j in range(4)] for i in range(4)]
        self.e1, self.e2 = UNITS[key]
        self._inv = {}
        self._pow = {}

    def mul(self, x, y):
        z0 = z1 = z2 = z3 = 0
        C = self.C
        for i in range(4):
            xi = x[i]
            if xi:
                Ci = C[i]
                for j in range(4):
                    yj = y[j]
                    if yj:
                        t = xi * yj
                        c = Ci[j]
                        z0 += t * c[0]
                        z1 += t * c[1]
                        z2 += t * c[2]
                        z3 += t * c[3]
        return (z0, z1, z2, z3)

    def mulq(self, x, y, q):
        z = self.mul(x, y)
        return (z[0] % q, z[1] % q, z[2] % q, z[3] % q)

    def mat(self, x):
        """M(x) with coords(x y) = M(x) coords(y)."""
        return [[sum(x[i] * self.Mom[i][r][c] for i in range(4)) for c in range(4)] for r in range(4)]

    def inverse(self, x):
        x = tuple(x)
        if x not in self._inv:
            Mi = q_inverse(self.mat(x))
            col = [Mi[r][0] for r in range(4)]
            assert all(c.denominator == 1 for c in col), 'not a unit of O'
            y = tuple(int(c) for c in col)
            assert self.mul(x, y) == (1, 0, 0, 0)
            self._inv[x] = y
        return self._inv[x]

    def pw(self, x, n):
        if n < 0:
            x, n = self.inverse(x), -n
        r, b = (1, 0, 0, 0), tuple(x)
        while n:
            if n & 1:
                r = self.mul(r, b)
            n >>= 1
            if n:
                b = self.mul(b, b)
        return r

    def pwq(self, x, n, q):
        if n < 0:
            x, n = self.inverse(x), -n
        r, b = (1 % q, 0, 0, 0), tuple(v % q for v in x)
        while n:
            if n & 1:
                r = self.mulq(r, b, q)
            n >>= 1
            if n:
                b = self.mulq(b, b, q)
        return r

    def unit(self, a, b):
        """eps1^a eps2^b, exact (cached)."""
        k = (a, b)
        if k not in self._pow:
            self._pow[k] = self.mul(self.pw(self.e1, a), self.pw(self.e2, b))
        return self._pow[k]

    def unit_mod(self, a, b, q):
        return self.mulq(self.pwq(self.e1, a, q), self.pwq(self.e2, b, q), q)

    def norm(self, x):
        return int(sp.Matrix(self.mat(x)).det())

    # --- sympy cross-checks --------------------------------------------------------------------------
    def sympy_poly(self, x):
        t = sp.Symbol('t')
        return sp.Poly(sum(sp.Rational(int(x[i])) * sum(sp.Rational(self.Bm[i][j].numerator, self.Bm[i][j].denominator)
                                                         * t ** j for j in range(4)) for i in range(4)), t, domain='QQ')

    def g_poly(self):
        t = sp.Symbol('t')
        return sp.Poly(sum(c * t ** k for k, c in enumerate(self.gc)), t, domain='QQ')

    def norm_res(self, x):
        return sp.resultant(self.g_poly(), self.sympy_poly(x))

    def sympy_structure_ok(self):
        t = sp.Symbol('t')
        g = self.g_poly()
        om = [sp.Poly(sum(sp.Rational(c.numerator, c.denominator) * t ** j for j, c in enumerate(row)), t, domain='QQ')
              for row in self.Bm]
        Bti = sp.Matrix([[sp.Rational(c.numerator, c.denominator) for c in row] for row in self.Bti])
        for i in range(4):
            for j in range(4):
                r = (om[i] * om[j]).rem(g)
                pc = sp.Matrix([r.coeff_monomial(t ** k) for k in range(4)])
                co = Bti * pc
                if not all(v.is_integer for v in co) or tuple(int(v) for v in co) != self.C[i][j]:
                    return False
        return True


def fval(key, m, n):
    a0, a1, a2, a3, a4 = FORMS[key]
    return a0 * m ** 4 + a1 * m ** 3 * n + a2 * m ** 2 * n ** 2 + a3 * m * n ** 3 + a4 * n ** 4


# ======================================================================================================
# Lattices of Z^2: lower Hermite form H = (A, B, C), rows (A, 0), (B, C), A, C > 0, 0 <= B < A
# ======================================================================================================
def xgcd(a, b):
    """(g, s, t) with s a + t b = g = gcd(a, b) >= 0."""
    s0, s1, t0, t1 = 1, 0, 0, 1
    while b:
        k = a // b
        a, b = b, a - k * b
        s0, s1 = s1, s0 - k * s1
        t0, t1 = t1, t0 - k * t1
    if a < 0:
        a, s0, t0 = -a, -s0, -t0
    return a, s0, t0


def hnf(vecs):
    """Hermite form of the lattice spanned by vecs (must have rank 2).  Adds the vectors one at a time:
    the pivot row (B,C) and a new (x,y) are replaced, by the unimodular matrix [[s,t],[y/g,-C/g]], with
    (sB+tx, g) and ((y/g)B - (C/g)x, 0); vectors with second coordinate 0 are absorbed into A by gcd."""
    A = 0
    piv = None
    for v in vecs:
        x, y = int(v[0]), int(v[1])
        if y == 0:
            A = gcd(A, x)
            continue
        if piv is None:
            piv = (x, y) if y > 0 else (-x, -y)
            continue
        B, C = piv
        g, s, t = xgcd(C, y)
        A = gcd(A, (y // g) * B - (C // g) * x)
        piv = (s * B + t * x, g)
    assert piv is not None and A > 0, 'lattice of rank < 2'
    return (A, piv[0] % A, piv[1])


def lred(H, v):
    """Canonical representative of v mod H: (a, b) with 0 <= a < A, 0 <= b < C."""
    A, B, C = H
    k = v[1] // C
    return ((v[0] - k * B) % A, v[1] - k * C)


def lhas(H, v):
    return lred(H, v) == (0, 0)


def lidx(H):
    return H[0] * H[2]


def lbasis(H):
    return [(H[0], 0), (H[1], H[2])]


def box_reps(H):
    return [(a, b) for b in range(H[2]) for a in range(H[0])]


def int_kernel(Z):
    """Rows w of a unimodular U with w Z = 0 spanning the integer left kernel of the integer matrix Z."""
    n, m = len(Z), len(Z[0])
    Z = [list(r) for r in Z]
    U = [[int(i == j) for j in range(n)] for i in range(n)]
    r = 0
    for c in range(m):
        while True:
            nz = [i for i in range(r, n) if Z[i][c] != 0]
            if not nz:
                break
            piv = min(nz, key=lambda i: abs(Z[i][c]))
            Z[r], Z[piv] = Z[piv], Z[r]
            U[r], U[piv] = U[piv], U[r]
            clean = True
            for i in range(r + 1, n):
                if Z[i][c]:
                    k = Z[i][c] // Z[r][c]
                    Z[i] = [a - k * b for a, b in zip(Z[i], Z[r])]
                    U[i] = [a - k * b for a, b in zip(U[i], U[r])]
                    if Z[i][c]:
                        clean = False
            if clean:
                break
        if any(Z[i][c] for i in range(r, n)):
            r += 1
    assert all(v == 0 for i in range(r, n) for v in Z[i])
    return U[r:]


def intersect(L, M, certify=True):
    l1, l2 = lbasis(L)
    m1, m2 = lbasis(M)
    K = int_kernel([[l1[0], l1[1]], [l2[0], l2[1]], [-m1[0], -m1[1]], [-m2[0], -m2[1]]])
    gens = [(u[0] * l1[0] + u[1] * l2[0], u[0] * l1[1] + u[1] * l2[1]) for u in K]
    I = hnf(gens)
    if certify:
        S = hnf(lbasis(L) + lbasis(M))
        ok = all(lhas(L, v) and lhas(M, v) for v in lbasis(I))
        ok = ok and lidx(I) * lidx(S) == lidx(L) * lidx(M)
        check(ok, f'intersection {L} cap {M} = {I} not certified')
    return I


def coset_reps(Lold, Lnew):
    """Vectors of Lold, one in each class of Lold / Lnew."""
    A, B, C = Lold
    S = []
    for v in lbasis(Lnew):
        t, r1 = divmod(v[1], C)
        s, r2 = divmod(v[0] - t * B, A)
        assert r1 == 0 and r2 == 0, 'Lnew is not inside Lold'
        S.append((s, t))
    HS = hnf(S)
    reps = [(i * A + j * B, j * C) for j in range(HS[2]) for i in range(HS[0])]
    assert len(reps) * lidx(Lold) == lidx(Lnew)
    return reps


def lagrange(v1, v2):
    """Lagrange-Gauss reduction of a basis of a lattice in Z^2 (unimodular change)."""
    n = lambda v: v[0] * v[0] + v[1] * v[1]
    while True:
        if n(v2) < n(v1):
            v1, v2 = v2, v1
        k = Fr(v1[0] * v2[0] + v1[1] * v2[1], n(v1))
        k = round(k)
        if k == 0:
            return v1, v2
        v2 = (v2[0] - k * v1[0], v2[1] - k * v1[1])


def selftest_lattices():
    rng = random.Random(8604)
    ok = True
    for _ in range(300):
        A, C = rng.randint(1, 12), rng.randint(1, 12)
        H0 = (A, rng.randrange(A), C)
        b1, b2 = lbasis(H0)
        gens = []
        for _ in range(rng.randint(2, 5)):
            i, j = rng.randint(-5, 5), rng.randint(-5, 5)
            gens.append((i * b1[0] + j * b2[0], i * b1[1] + j * b2[1]))
        gens += [b1, (b2[0] + 3 * b1[0], b2[1])]
        rng.shuffle(gens)
        ok &= hnf(gens) == H0
        A2, C2 = rng.randint(1, 12), rng.randint(1, 12)
        H1 = (A2, rng.randrange(A2), C2)
        I = intersect(H0, H1)
        for x in range(-30, 31):
            for y in range(-30, 31):
                ok &= lhas(I, (x, y)) == (lhas(H0, (x, y)) and lhas(H1, (x, y)))
        reps = coset_reps(H0, I)
        ok &= len({lred(I, r) for r in reps}) == len(reps) and all(lhas(H0, r) for r in reps)
    return ok


# ======================================================================================================
# G_q, Lambda_q, candidates
# ======================================================================================================
def prime_power(q):
    p = next(d for d in range(2, q + 1) if q % d == 0)
    e, t = 0, q
    while t % p == 0:
        t //= p
        e += 1
    return p, e, t == 1


def closure(O, q, cap=None):
    """G_q = <eps1, eps2> mod qO as a set, by breadth-first closure of {1}; None if larger than cap."""
    one = (1 % q, 0, 0, 0)
    u1 = tuple(v % q for v in O.e1)
    u2 = tuple(v % q for v in O.e2)
    seen = {one}
    frontier = [one]
    while frontier:
        nxt = []
        for g in frontier:
            for u in (u1, u2):
                h = O.mulq(g, u, q)
                if h not in seen:
                    seen.add(h)
                    nxt.append(h)
        if cap is not None and len(seen) > cap:
            return None
        frontier = nxt
    return seen


def powers_table(O, x, q):
    one = (1 % q, 0, 0, 0)
    xq = tuple(v % q for v in x)
    tab = {one: 0}
    y, k = xq, 1
    while y != one:
        tab[y] = k
        y = O.mulq(y, xq, q)
        k += 1
    return k, tab


class KernelLattice:
    """Lambda_q = ker(Z^2 -> (O/qO)^*, (a,b) -> eps1^a eps2^b), certified."""

    def __init__(self, O, q):
        self.q = q
        self.h1, tab1 = powers_table(O, O.e1, q)
        self.h2, _ = powers_table(O, O.e2, q)
        i2 = tuple(v % q for v in O.inverse(O.e2))
        rels = []
        y = (1 % q, 0, 0, 0)
        for b in range(self.h2):            # y = eps2^-b
            a = tab1.get(y)
            if a is not None:
                rels.append((a, b))
            y = O.mulq(y, i2, q)
        self.nrel = len(rels)
        self.H = hnf([(self.h1, 0), (0, self.h2)] + rels)
        G = closure(O, q)
        self.order = len(G)
        self.minus_one = (q - 1, 0, 0, 0) in G if q > 2 else True
        one = (1 % q, 0, 0, 0)
        ok = lidx(self.H) == self.order and self.h1 * self.h2 == self.nrel * self.order and self.H[0] == self.h1
        ok = ok and all(O.unit_mod(v[0], v[1], q) == one for v in lbasis(self.H))
        check(ok, f'Lambda_{q} not certified (index {lidx(self.H)}, |G_q| {self.order})')
        # -1 = eps^star mod q, if any (search over the box representatives)
        self.star = None
        if self.minus_one and q > 2:
            minus = (q - 1, 0, 0, 0)
            A, B, C = self.H
            y = (1, 0, 0, 0)
            u1 = tuple(v % q for v in O.e1)
            for b in range(C):
                x = y
                for a in range(A):
                    if x == minus:
                        self.star = (a, b)
                        break
                    x = O.mulq(x, u1, q)
                if self.star:
                    break
                y = O.mulq(y, tuple(v % q for v in O.e2), q)
            check(self.star is not None and O.unit_mod(self.star[0], self.star[1], q) == minus, f'-1 mod {q}')

    def candidates(self, O, w):
        """Box representatives (a,b) of Z^2 / Lambda_q with both markers of eps^(a,b) w == 0 mod q."""
        q = self.q
        A, B, C = self.H
        u1 = tuple(v % q for v in O.e1)
        u2 = tuple(v % q for v in O.e2)
        out = set()
        y = tuple(v % q for v in w)
        for b in range(C):
            x = y
            for a in range(A):
                if x[2] == 0 and x[3] == 0:
                    out.add((a, b))
                x = O.mulq(x, u1, q)
            y = O.mulq(y, u2, q)
        return out


# ======================================================================================================
# Proposition 17 (lifting.binomial.multi, section 8), exact
# ======================================================================================================
def nu(x, p):
    if x == 0:
        return INF
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def LAMBDA(n, p, e):
    """Lambda(n) = ceil((n theta + 1)/(p - 1)), theta = e(p - 1) - 1 (proof.md 1.5)."""
    th = e * (p - 1) - 1
    return -((-(n * th + 1)) // (p - 1))


def p17_data(O, q, z, basis):
    """x_k = l_k K_1 z, y_k = l_k K_2 z for T_j = M(eta_j), eta_j = eps^(basis_j), K_j = (T_j - I)/q (exact)."""
    p, e, _ = prime_power(q)
    if z[2] != 0 or z[3] != 0:
        return None                                   # the seed is not a common zero: no instance
    cols = []
    for vec in basis:
        eta = O.unit(vec[0], vec[1])
        if (eta[0] - 1) % q or any(c % q for c in eta[1:]):
            return None                               # eta is not == 1 mod qO
        ez = O.mul(eta, z)
        d = [ez[k] - z[k] for k in range(4)]
        cols.append((d[2] // q, d[3] // q))
    (x1, x2), (y1, y2) = cols
    lam = min(nu(x1, p), nu(x2, p), nu(y1, p), nu(y2, p))
    J = x1 * y2 - x2 * y1
    nJ = nu(J, p)
    return dict(p=p, e=e, x=(x1, x2), y=(y1, y2), lam=lam, nJ=nJ)


def p17_hypotheses(p, e, lam, nJ):
    """Proposition 17 exactly as stated: q = p^e (e >= 2 if p = 2), iota = [p = 2], lambda <= e-1-iota and
    delta = nu(J) - lambda <= e-1-iota.  Also the strict inequalities used in its proof."""
    iota = 1 if p == 2 else 0
    setting = bool(sp.isprime(p)) and e >= (2 if p == 2 else 1)
    delta = INF if (lam >= INF or nJ >= INF) else nJ - lam
    bound = e - 1 - iota
    holds = setting and lam <= bound and delta <= bound
    L2 = LAMBDA(2, p, e)
    strict = (e + delta < L2) and (e + lam < L2)
    return dict(iota=iota, setting=setting, delta=delta, bound=bound, holds=holds, L2=L2,
                L2ok=(L2 == 2 * e - iota), strict=strict)


def vfact(n, p):
    v, t = 0, p
    while t <= n:
        v += n // t
        t *= p
    return v


def p17_mechanism(O, q, z, basis, D=4):
    """Steps 1-3 of the proof of Proposition 17 at this instance: b^(k)_ij = l_k K1^i K2^j z (i+j <= D) mod p^M;
    after relabelling, f1 has b00 = 0 and minimizers in {(1,0),(0,1)} containing (1,0), with s* = e + lambda;
    f2' = f2 - rho f1 has E' = {(0,1)} and s*' = e + delta; both s* < Lambda(2) <= val(i,j) for i+j >= 2."""
    p, e, _ = prime_power(q)
    M = e * (D + 2) + 30
    PM = p ** M
    kap = []
    for vec in basis:
        eta = O.unit(vec[0], vec[1])
        kap.append(tuple(((eta[k] - (1 if k == 0 else 0)) // q) % PM for k in range(4)))
    b = {}
    row = tuple(c % PM for c in z)
    for i in range(D + 1):
        cur = row
        for j in range(D + 1 - i):
            b[(i, j)] = (cur[2], cur[3])
            cur = O.mulq(cur, kap[1], PM)
        row = O.mulq(row, kap[0], PM)
    nuM = lambda x: min(nu(x % PM, p), M)
    x, y = b[(1, 0)], b[(0, 1)]
    lam = min(nuM(v) for v in x + y)
    if lam >= M - 5:
        return False, 'lambda beyond the working precision'
    swap_i = min(nuM(x[0]), nuM(x[1])) > lam
    if swap_i:
        b = {(j, i): v for (i, j), v in b.items()}
        x, y = y, x
    swap_m = nuM(x[0]) > lam
    if swap_m:
        b = {k: (v[1], v[0]) for k, v in b.items()}
        x, y = (x[1], x[0]), (y[1], y[0])
    Pl, PM2 = p ** lam, p ** (M - lam)
    rho = ((x[1] // Pl) * pow((x[0] // Pl) % PM2, -1, PM2)) % PM2
    val1, val2 = {}, {}
    for (i, j), (u1, u2) in b.items():
        for val, bv, prec in ((val1, u1 % PM, M), (val2, (u2 - rho * u1) % PM2, M - lam)):
            vv = nu(bv, p)
            val[(i, j)] = INF if vv >= prec else (i + j) * e + vv - vfact(i, p) - vfact(j, p)
    L2 = LAMBDA(2, p, e)
    ok = val1[(0, 0)] == INF and val2[(0, 0)] == INF and val2[(1, 0)] == INF
    for val in (val1, val2):
        for (i, j), v in val.items():
            if i + j >= 1 and v < INF:
                ok = ok and v >= LAMBDA(i + j, p, e)
    m1 = min(val1.values())
    E1 = sorted(k for k, v in val1.items() if v == m1)
    m2 = min(val2.values())
    E2 = sorted(k for k, v in val2.items() if v == m2)
    nJ = nu((x[0] * y[1] - x[1] * y[0]) % PM2, p)
    ok = ok and m1 == e + lam and m1 < L2 and (1, 0) in E1 and set(E1) <= {(1, 0), (0, 1)}
    ok = ok and E2 == [(0, 1)] and m2 == e + (nJ - lam) and m2 < L2
    txt = (f'relabel(idx={int(swap_i)},mk={int(swap_m)}); f1: s*={m1}=e+lam, E={E1}; f2\': s*={m2}=e+delta, E={E2}; '
           f'Lambda(2)={L2}; val>=Lambda(i+j) for 1<=i+j<={D}')
    return ok, txt


def four_bases(H):
    A, B, C = H
    hb = [(A, 0), (B, C)]
    o2 = [(A, 0), ((B - A) if B else 0, C)]          # 02 / proof.md 4.1: (h, 0), (-j, k), 0 <= j < h
    uni = [(2 * A + 3 * B, 3 * C), (3 * A + 5 * B, 5 * C)]   # [[2,3],[3,5]] times the Hermite basis
    red = list(lagrange((A, 0), (B, C)))
    return {'hermite': hb, '02': o2, 'unimod': uni, 'reduced': red}


def basis_ok(H, bs):
    (a, b), (c, d) = bs
    return lhas(H, bs[0]) and lhas(H, bs[1]) and abs(a * d - b * c) == lidx(H)


# ======================================================================================================
# The sieve (proof.md 4.4)
# ======================================================================================================
def sieve(seq, KER, CAND, zeros, P17, sibling='same', record=None):
    """KER[q]: Hermite form of Lambda_q; CAND[q]: candidate box representatives mod Lambda_q; zeros: list of
    exponents; P17[(q, i)]: Proposition 17 holds at zero i for q.
    sibling = 'same' (Lemma 4.4(b): inside the Lambda_q-class of a zero at which Prop. 17 holds for that SAME q),
              'none', or 'any' (WRONG on purpose: Lambda_q-class of a zero for one q, Prop. 17 at any used q)."""
    L = (1, 0, 1)
    live = [(0, 0)]
    used, steps, sibs = [], [], []
    cert = {}
    for q in seq:
        Lq = KER[q]
        Ln = intersect(L, Lq)
        reps = coset_reps(L, Ln)
        newlive = []
        for r in live:
            for x in reps:
                v = (r[0] + x[0], r[1] + x[1])
                if lred(Lq, v) in CAND[q]:
                    newlive.append(lred(Ln, v))
        check(len(set(newlive)) == len(newlive), 'duplicate live classes')
        mods = used + [q]
        nc = ns = 0
        rest = []
        rec = dict(q=q, L=Ln, cand=list(newlive), cert=[], sib=[])
        for v in newlive:
            hits = [i for i, ex in enumerate(zeros) if lred(Ln, ex) == v]
            if hits:
                check(len(hits) == 1, f'two known zeros in one class {v}')
                i = hits[0]
                okq = [qq for qq in mods if P17.get((qq, i))]
                if okq:
                    check(i not in cert, 'a zero certified twice')
                    cert[i] = (okq[0], v, Ln)
                    rec['cert'].append(v)
                    nc += 1
                else:
                    rest.append(v)
                continue
            if sibling == 'same':
                sib = [(qq, i) for qq in mods for i, ex in enumerate(zeros)
                       if P17.get((qq, i)) and lred(KER[qq], v) == lred(KER[qq], ex)]
            elif sibling == 'any':
                sib = [(qq, i) for qq in mods for i, ex in enumerate(zeros)
                       if any(P17.get((q3, i)) for q3 in mods) and lred(KER[qq], v) == lred(KER[qq], ex)]
            else:
                sib = []
            if sib:
                ns += 1
                sibs.append((q, v, sib))
                rec['sib'].append(v)
                continue
            rest.append(v)
        steps.append((q, len(newlive), nc, ns, len(rest)))
        rec['left'] = list(rest)
        if record is not None:
            record.append(rec)
        used.append(q)
        L, live = Ln, rest
    return dict(steps=steps, live=live, cert=cert, L=L, sibs=sibs, used=used)


def direct_check(L, used, KER, CAND, zeros, P17):
    """Every class of Z^2 / L: excluded (non-candidate mod some Lambda_q) or inside the Lambda_q-class of a zero
    at which Prop. 17 holds for that q.  Returns (#excluded, #covered, list of unaccounted classes)."""
    zc = [(q, lred(KER[q], ex)) for q in used for i, ex in enumerate(zeros) if P17.get((q, i))]
    excl = cov = 0
    bad = []
    for v in box_reps(L):
        if any(lred(KER[q], v) not in CAND[q] for q in used):
            excl += 1
        elif any(lred(KER[q], v) == zr for q, zr in zc):
            cov += 1
        else:
            bad.append(v)
    return excl, cov, bad


# ======================================================================================================
# MAIN
# ======================================================================================================
print('verify_exclusion.py (8604.a1) -- independent re-implementation of the exclusion step (exact)')
print('nothing imported from core.py, lifting2.py, 01_coverage.py or 02_exclusion.py; data copied as literals;')
print('every number is compared with 02_exclusion.out.txt and 01_coverage.out.txt')
print()


def literals_vs_files():
    """The literals above against the CURRENT core.py, 01_coverage.out.txt and 02_exclusion.out.txt of this folder,
    read as text (never imported or executed)."""
    import ast
    import os
    import re
    here = os.path.dirname(os.path.abspath(__file__))
    out = {}
    try:
        core = open(os.path.join(here, 'core.py'), encoding='utf-8').read()
        o01 = open(os.path.join(here, '01_coverage.out.txt'), encoding='utf-8').read()
        o02 = open(os.path.join(here, '02_exclusion.out.txt'), encoding='utf-8').read()
    except OSError as exc:
        print(f'   (files not readable, comparison skipped: {exc})')
        return None
    fr = lambda c: f'Fr({c.numerator}, {c.denominator})' if Fr(c).denominator != 1 else str(int(c))
    ok_core = True
    for key in FORMS:
        bas = '[' + ', '.join('[' + ', '.join(fr(Fr(c)) for c in row) + ']' for row in BASES[key]) + ']'
        ok_core &= f'F={FORMS[key]}' in core and f'BASIS={bas}' in core
        ok_core &= f"'{key}': {UNITS[key]}" in core and f"'{key}': {KNOWN[key]}" in core
    out['core.py (forms, bases, units, known)'] = ok_core
    m = re.search(r'SEEDS = (\{.*?\n\})', o01, re.S)
    out['01_coverage.out.txt (SEEDS)'] = bool(m) and ast.literal_eval(m.group(1)) == SEEDS_01
    # 02_exclusion.out.txt
    parsed, avail, neg = [], {}, {}
    key = None
    cur = None
    for line in o02.splitlines():
        s = line.strip()
        mm = re.match(r'CLASS \((\d),(\d)\)', s)
        if mm:
            key = mm.group(1) + mm.group(2)
            continue
        mm = re.match(r'moduli available \(\|G_q\|\): (.*)$', s)
        if mm:
            avail[key] = {int(a): int(b) for a, b in re.findall(r'(\d+)\((\d+)\)', mm.group(1))}
            continue
        mm = re.match(r'checked: -1 is not a residue of <eps1,eps2> mod (\d+)', s)
        if mm:
            neg[key] = int(mm.group(1))
            continue
        mm = re.match(r'--- orbit c = (\d+), seed w = (\(.*?\)); known solutions in it: (\[.*?\]) \(exponents (\[.*?\])\)', s)
        if mm:
            cur = dict(key=key, c=int(mm.group(1)), w=ast.literal_eval(mm.group(2)), known=ast.literal_eval(mm.group(3)),
                       exps=ast.literal_eval(mm.group(4)), steps=[], cert=[], inst=[], G={})
            parsed.append(cur)
            continue
        mm = re.match(r'moduli used: (\[.*?\]);\s+final lattice index (\d+)', s)
        if mm:
            cur['used'], cur['index'] = ast.literal_eval(mm.group(1)), int(mm.group(2))
            continue
        mm = re.match(r'after q =\s*(\d+): (\d+) candidate .*?; (\d+) certified .*?, (\d+) empty .*?, (\d+) left', s)
        if mm:
            cur['steps'].append(tuple(int(g) for g in mm.groups()))
            continue
        mm = re.match(r'class (\(.*?\))(?: mod L = \(.*?\))? contains the solution (\(.*?\)): Proposition 17 at q = (\d+) \(lambda = (\d+), '
                      r'nu\(J\) = (\d+), e = (\d+)\)', s)
        if mm:
            cur['cert'].append((ast.literal_eval(mm.group(1)), ast.literal_eval(mm.group(2)))
                               + tuple(int(g) for g in mm.groups()[2:]))
            continue
        mm = re.match(r'Prop\. 17 instance used: q = (\d+) \(p = (\d+), e = (\d+)\), zero (\(.*?\)): lambda = (\d+), '
                      r'nu\(J\) = (\d+), delta = (\d+) <= e-1-iota = (-?\d+): True', s)
        if mm:
            g = mm.groups()
            cur['inst'].append((int(g[0]), int(g[1]), int(g[2]), ast.literal_eval(g[3]), int(g[4]), int(g[5]), int(g[6]),
                                int(g[7])))
            continue
        mm = re.match(r'q = (\d+): \|G_q\| = (\d+), (?:Lambda_q = .*?, )?candidate classes mod Lambda_q: (\d+)', s)
        if mm and cur is not None:
            cur['G'][int(mm.group(1))] = (int(mm.group(2)), int(mm.group(3)))
    fields = ('key', 'c', 'w', 'known', 'exps', 'used', 'index', 'steps', 'cert', 'inst', 'G')
    same = len(parsed) == len(ORBITS) and all(all(p[f] == o[f] for f in fields) for p, o in zip(parsed, ORBITS))
    out['02_exclusion.out.txt (12 orbits: seeds, zeros, moduli, steps, certified classes, Prop. 17 instances, |G_q|)'] = same
    out['02_exclusion.out.txt (moduli available)'] = avail == AVAILABLE
    out['02_exclusion.out.txt (-1 non-residue moduli)'] = neg == NEG02
    out['02_exclusion.out.txt ends with ALL ORBITS CERTIFIED'] = 'ALL ORBITS CERTIFIED' in o02
    return out


cmp_files = literals_vs_files()
if cmp_files is not None:
    for k, v in cmp_files.items():
        check(v, f'literals differ from {k}')
        print(f'   literals == current {k}: {v}')
    print()
st = selftest_lattices()
check(st, 'lattice self-test')
print(f'self-test of the lattice routines (Hermite form, intersection, coset representatives) on 300 random '
      f'lattices against brute-force membership in [-30,30]^2: {st}')

FIELDS, KERS = {}, {}


def ker(key, q):
    if (key, q) not in KERS:
        KERS[(key, q)] = KernelLattice(FIELDS[key], q)
    return KERS[(key, q)]


NEGTAB = {}
for key in ('10', '01', '11'):
    print('=' * 112)
    O = Order(key)
    FIELDS[key] = O
    g = O.g_poly()
    print(f'CLASS {LABEL[key]}  G = {FORMS[key]}  g(t) = {g.as_expr()}')
    check(O.integral, 'the basis does not span an order (non-integral structure constants)')
    sok = O.sympy_structure_ok()
    check(sok, 'structure constants disagree with sympy')
    check(O.Bm[0] == [1, 0, 0, 0] and O.Bm[1] == [0, 1, 0, 0], 'omega_0 = 1, omega_1 = beta')
    blk = O.Bm[2][2] * O.Bm[3][3] - O.Bm[2][3] * O.Bm[3][2]
    check(blk != 0, '(beta^2, beta^3)-block singular')
    irr = g.is_irreducible
    check(irr, 'g reducible')
    print(f'   regular representation integral: {O.integral}; equals sympy products omega_i omega_j mod g: {sok}; '
          f'omega_0 = 1, omega_1 = beta, det of the (beta^2,beta^3)-block = {blk}; g irreducible over Q: {irr}')
    print('   => an element of O lies in Z + Z beta iff its markers x2, x3 vanish, and n - m beta = (n, -m, 0, 0)')
    for nm, u in (('eps1', O.e1), ('eps2', O.e2)):
        nd, nr = O.norm(u), O.norm_res(u)
        ui = O.inverse(u)
        check(nd == 1 and nr == 1, f'{nm}: norm {nd} / {nr}')
        print(f'   {nm} = {u}: N = {nd} (det) = {nr} (resultant); inverse in O = {ui}')
    # sanity: G(m,n) = N(n - m beta) on a few pairs
    rr = random.Random(int(key, 2))
    okn = all(O.norm((n, -m, 0, 0)) == fval(key, m, n) for m, n in [(rr.randint(-9, 9), rr.randint(-9, 9))
                                                                      for _ in range(12)])
    check(okn, 'G(m,n) != N(n - m beta)')
    print(f'   G(m,n) = N(n - m beta) on 12 random pairs: {okn}')
    # |G_q| for every modulus of 02's list, and the skip rule |G_q| <= 5000
    got, neg = {}, {}
    for q in MODULI_02:
        p, e, pp = prime_power(q)
        assert pp
        S = closure(O, q, cap=MAXG)
        if S is None:
            continue
        got[q] = len(S)
        neg[q] = (q - 1, 0, 0, 0) in S
    agree = got == AVAILABLE[key]
    check(agree, f'moduli available: {got} vs printed {AVAILABLE[key]}')
    print(f'   |G_q| <= {MAXG} (BFS closure): ' + ', '.join(f'{q}({n})' for q, n in got.items()))
    print(f'   agrees with 02\'s "moduli available" line (values and skipped moduli): {agree}')
    NEGTAB[key] = neg
    first = next(q for q in got if not neg[q])
    check(first == NEG02[key], f'first modulus with -1 not a residue: {first} vs printed {NEG02[key]}')
    print('   -1 in G_q ?  ' + ', '.join(f'{q}:{"yes" if neg[q] else "NO"}' for q in got))
    print(f'   first modulus (02\'s order) with -1 not in G_q: {first} (02 prints {NEG02[key]}).  Hence -1 is not in '
          f'<eps1,eps2>: eps1^a eps2^b = -1 would give -1 in G_{first}.')
print()

# seeds of 02 == seeds of 01
s02 = {}
for orb in ORBITS:
    s02.setdefault(orb['key'], {}).setdefault(orb['c'], []).append(orb['w'])
check(s02 == SEEDS_01, 'seeds in 02_exclusion.out.txt differ from 01_coverage.out.txt')
print(f'seeds printed by 02_exclusion.out.txt == SEEDS of 01_coverage.out.txt: {s02 == SEEDS_01}')
print()

CTX = []
SUMMARY = []
INSTANCES = []
BASISSTAT = dict(pairs=0, agree=0, holds=0)
for orb in ORBITS:
    key, c, w = orb['key'], orb['c'], orb['w']
    O = FIELDS[key]
    print('=' * 112)
    print(f'CLASS {LABEL[key]}  c = {c}  seed w = {w}')
    nf0 = len(FAILS)
    # ---- (0) the seed and the known zeros ----------------------------------------------------------
    check(O.norm(w) == c and O.norm_res(w) == c, f'N(w) != {c}')
    known_c = sorted((m, n) for (m, n) in KNOWN[key] if fval(key, m, n) == c)
    check(known_c == sorted(CLAIMED[key][c]), f'KNOWN with G = c {known_c} vs claimed {CLAIMED[key][c]}')
    check(known_c == sorted(orb['known']), f'known {known_c} vs printed {orb["known"]}')
    zeros, zinfo = [], []
    R0 = 12
    for (m, n), ex in zip(orb['known'], orb['exps']):
        target = (n, -m, 0, 0)
        hits = []
        for a in range(-R0, R0 + 1):
            for b in range(-R0, R0 + 1):
                x = O.mul(O.unit(a, b), w)
                if x == target:
                    hits.append((a, b, '+'))
                elif x == tuple(-t for t in target):
                    hits.append((a, b, '-'))
        check(len(hits) == 1 and hits[0][:2] == ex, f'exponent of {(m, n)}: {hits} vs printed {ex}')
        zf = O.mul(O.unit(*ex), w)
        check(zf[2] == 0 and zf[3] == 0 and O.norm(target) == c, f'{(m, n)} is not a zero')
        zeros.append(ex)
        zinfo.append(dict(mn=(m, n), exp=ex, z=zf))
        print(f'   known solution (m,n) = {(m, n)}: eps1^{ex[0]} eps2^{ex[1]} w = {hits[0][2]}({n} - ({m}) beta) exactly; '
              f'the only exponent in |a|,|b| <= {R0}; N = {O.norm(target)}')
    if not zeros:
        print('   no known solution in this orbit')
    # ---- (1) Lambda_q, |G_q| and candidates for each modulus used --------------------------------------
    check(all(prime_power(q)[2] and q >= 3 for q in orb['used']), 'a modulus used is not a prime power >= 3')
    KER, CAND, PE = {}, {}, {}
    for q in orb['used']:
        kq = ker(key, q)
        KER[q] = kq.H
        CAND[q] = kq.candidates(O, w)
        PE[q] = prime_power(q)[:2]
        pg, pc = orb['G'][q]
        check(kq.order == pg and len(CAND[q]) == pc,
              f'q = {q}: |G_q| = {kq.order}, candidates = {len(CAND[q])} vs printed {pg}, {pc}')
        for ex in zeros:
            check(lred(kq.H, ex) in CAND[q], f'known zero {ex} is not a candidate mod {q}')
        print(f'   q = {q:3d} = {PE[q][0]}^{PE[q][1]}: |G_q| = {kq.order} (printed {pg}); Lambda_q = Hermite {kq.H} '
              f'[certified: basis -> 1 mod q, index = BFS |G_q|]; h1 = {kq.h1}, h2 = {kq.h2}; candidate classes '
              f'= {len(CAND[q])} (printed {pc}): {sorted(CAND[q])[:8]}{" ..." if len(CAND[q]) > 8 else ""}')
    # ---- (2) Proposition 17 at every (q, zero), four bases ---------------------------------------------
    P17 = {}
    for q in orb['used']:
        p, e = PE[q]
        bs = four_bases(KER[q])
        for nm, bb in bs.items():
            check(basis_ok(KER[q], bb), f'{nm} basis of Lambda_{q} is not a basis')
        for i, zz in enumerate(zinfo):
            res = {nm: p17_data(O, q, zz['z'], bb) for nm, bb in bs.items()}
            if not check(all(r is not None for r in res.values()), f'no coset system at q = {q}, zero {zz["mn"]}'):
                P17[(q, i)] = False
                continue
            lams = {r['lam'] for r in res.values()}
            nJs = {r['nJ'] for r in res.values()}
            hyp = {nm: p17_hypotheses(p, e, r['lam'], r['nJ']) for nm, r in res.items()}
            verdicts = {h['holds'] for h in hyp.values()}
            check(len(verdicts) == 1, f'Prop. 17 verdict depends on the basis at q = {q}, zero {zz["mn"]}')
            if min(lams) < e:
                check(len(lams) == 1, f'lambda < e depends on the basis at q = {q}, zero {zz["mn"]}')
            if min(nJs) < e:
                check(len(nJs) == 1, f'nu(J) < e depends on the basis at q = {q}, zero {zz["mn"]}')
            BASISSTAT['pairs'] += 1
            BASISSTAT['agree'] += int(len(lams) == 1 and len(nJs) == 1)
            BASISSTAT['holds'] += int(hyp['02']['holds'])
            h = hyp['02']
            P17[(q, i)] = h['holds']
            vals = '; '.join(f'{nm}: lam={fnu(r["lam"])} nuJ={fnu(r["nJ"])}' for nm, r in res.items())
            dep = '' if (len(lams) == 1 and len(nJs) == 1) else '  [basis-dependent values, all >= e: never certify]'
            print(f'   Prop. 17 at {zz["mn"]}, q = {q} (p = {p}, e = {e}, iota = {h["iota"]}, bound e-1-iota = {h["bound"]}): '
                  f'{vals} -> {"HOLDS" if h["holds"] else "fails"}{dep}')
            if h['holds']:
                for nm in ('02', 'unimod'):
                    okm, txt = p17_mechanism(O, q, zz['z'], bs[nm])
                    check(okm, f'mechanism of Prop. 17 at q = {q}, zero {zz["mn"]}, basis {nm}: {txt}')
                    if nm == '02':
                        print(f'      proof mechanism ({nm} basis): {txt}: {okm}')
                # T_j in SL_4(Z): det = N(eta_j) = 1, also checked mod a prime; T_1 T_2 = T_2 T_1 mod the prime
                Pp = 2 ** 61 - 1
                Ts = [[[v % Pp for v in row] for row in O.mat(O.unit(*vec))] for vec in bs['02']]
                dets = [det_mod(T, Pp) for T in Ts]
                comm = [[v % Pp for v in r] for r in mmul(Ts[0], Ts[1])] == \
                    [[v % Pp for v in r] for r in mmul(Ts[1], Ts[0])]
                check(dets == [1, 1] and comm, f'T_j not in SL_4 or not commuting at q = {q}')
    # compare with every printed Prop. 17 instance
    for (q, p_, e_, mn, plam, pnJ, pdel, pbound) in orb['inst']:
        i = next(k for k, zz in enumerate(zinfo) if zz['mn'] == mn)
        r = p17_data(O, q, zinfo[i]['z'], four_bases(KER[q])['02'])
        h = p17_hypotheses(PE[q][0], PE[q][1], r['lam'], r['nJ'])
        same = (PE[q] == (p_, e_) and (r['lam'], r['nJ'], h['delta'], h['bound']) == (plam, pnJ, pdel, pbound))
        check(same and h['holds'], f'printed instance {(q, mn)}: lam {plam} nuJ {pnJ} delta {pdel} bound {pbound}')
        check(h['strict'] and h['L2ok'], f'strict inequalities at {(q, mn)}')
        INSTANCES.append((key, c, q, p_, e_, mn, r['lam'], r['nJ'], h['delta'], h['bound'], h['L2'],
                          e_ + h['delta'], e_ + r['lam']))
        print(f'   printed instance q = {q} (p = {p_}, e = {e_}), zero {mn}: lambda = {plam}, nu(J) = {pnJ}, delta = '
              f'{pdel}, bound = {pbound}: recomputed equal: {same}; hypotheses hold: {h["holds"]}; strict: e+delta = '
              f'{e_ + h["delta"]} < Lambda(2) = {h["L2"]} (= 2e-iota: {h["L2ok"]}) and e+lambda = {e_ + r["lam"]} < '
              f'{h["L2"]}: {h["strict"]}')
    # ---- (3) the sieve in 02's order -----------------------------------------------------------------
    sv = sieve(orb['used'], KER, CAND, zeros, P17)
    check(sv['steps'] == orb['steps'], f'sieve counts {sv["steps"]} vs printed {orb["steps"]}')
    check(lidx(sv['L']) == orb['index'], f'final index {lidx(sv["L"])} vs printed {orb["index"]}')
    for (q, nl, nc, ns, nr) in sv['steps']:
        print(f'   sieve after q = {q:3d}: {nl} candidate classes mod the new lattice, {nc} certified, {ns} sibling-empty, '
              f'{nr} left')
    for (q, v, lst) in sv['sibs']:
        print(f'      sibling-empty class {v} at step q = {q}: inside the Lambda_qq-class of a zero with Prop. 17 at the '
              f'same qq, (qq, zero) in {[(qq, zinfo[i]["mn"]) for qq, i in lst]}')
    check(not sv['live'], 'live classes remain')
    check(sorted(sv['cert']) == list(range(len(zeros))), 'not every known zero certified exactly once')
    for (rep, mn, pq, plam, pnJ, pe) in orb['cert']:
        i = next(k for k, zz in enumerate(zinfo) if zz['mn'] == mn)
        q0, v, Lc = sv['cert'][i]
        okc = q0 == pq and lhas(Lc, (rep[0] - zeros[i][0], rep[1] - zeros[i][1])) and PE[pq][1] == pe
        check(okc, f'printed certified class {rep} of {mn} at q = {pq}')
        print(f'   printed "class {rep} contains {mn}: Prop. 17 at q = {pq}": rep == exponent {zeros[i]} mod the step '
              f'lattice {Lc}: {okc}')
    # final lattice two ways
    Lf = sv['L']
    lcm = 1
    for q in orb['used']:
        lcm = lcm * q // gcd(lcm, q)
    kl = ker(key, lcm)
    check(kl.H == Lf, f'final lattice {Lf} vs Lambda_lcm {kl.H}')
    print(f'   final lattice {Lf}, index {lidx(Lf)} (printed {orb["index"]}); recomputed from scratch as Lambda_{lcm} '
          f'(lcm of the moduli; BFS |G| = {kl.order}): {kl.H == Lf}; live classes left: {len(sv["live"])}')
    # ---- (3') direct check of all classes modulo the final lattice ------------------------------------
    excl, cov, bad = direct_check(Lf, orb['used'], KER, CAND, zeros, P17)
    check(not bad, f'direct check: {len(bad)} classes unaccounted')
    print(f'   direct check over all {lidx(Lf)} classes mod the final lattice: {excl} excluded by a candidate test, '
          f'{cov} inside the Lambda_q-class of a zero with Prop. 17 at that q, {len(bad)} unaccounted')
    # ---- (3'') reversed order -------------------------------------------------------------------------
    if len(orb['used']) > 1:
        sr = sieve(list(reversed(orb['used'])), KER, CAND, zeros, P17)
        okr = not sr['live'] and sorted(sr['cert']) == list(range(len(zeros)))
        check(okr, 'reversed order does not terminate')
        print(f'   reversed order {list(reversed(orb["used"]))}: steps {sr["steps"]}; terminates with no live class: {okr}; '
              f'certified at q = {[sr["cert"][i][0] for i in sorted(sr["cert"])]}; sibling-empty classes: '
              f'{[(q, v, [(qq, zinfo[i]["mn"]) for qq, i in lst]) for (q, v, lst) in sr["sibs"]]}')
    ok = len(FAILS) == nf0
    SUMMARY.append((key, c, w, len(zeros), ok))
    CTX.append(dict(orb=orb, KER=KER, CAND=CAND, zeros=zeros, zinfo=zinfo, P17=P17, sv=sv))
    print(f'   ORBIT {"CONFIRMED" if ok else "NOT CONFIRMED"}: the zeros of the markers are exactly the exponents of '
          f'the known solutions')
print()

# ======================================================================================================
# Different modulus sets (greedy), and the phantom classes
# ======================================================================================================
print('=' * 112)
print('OTHER MODULUS SETS: greedy sieve (fewest uncertified live classes, then smallest index), at most 10 moduli,')
print('over (i) the odd moduli of 02\'s list and (ii) the powers of 2 (4..64), each with |G_q| <= 5000; Prop. 17 in the')
print('Hermite basis; certified and same-q sibling rules')


def greedy(key, w, zinfo, pool, max_steps=8, max_index=10 ** 8):
    O = FIELDS[key]
    zeros = [zz['exp'] for zz in zinfo]
    KER, CAND, P17 = {}, {}, {}
    for q in pool:
        kq = ker(key, q)
        KER[q] = kq.H
        CAND[q] = kq.candidates(O, w)
        p, e, _ = prime_power(q)
        for i, zz in enumerate(zinfo):
            r = p17_data(O, q, zz['z'], four_bases(kq.H)['hermite'])
            P17[(q, i)] = p17_hypotheses(p, e, r['lam'], r['nJ'])['holds']
    L, live, used, cert = (1, 0, 1), [(0, 0)], [], {}
    while live and len(used) < max_steps:
        best = None
        for q in pool:
            if q in used:
                continue
            Ln = intersect(L, KER[q], certify=False)
            if lidx(Ln) > max_index:
                continue
            reps = coset_reps(L, Ln)
            mods = used + [q]
            rest, nc = [], {}
            for r in live:
                for x in reps:
                    v = (r[0] + x[0], r[1] + x[1])
                    if lred(KER[q], v) not in CAND[q]:
                        continue
                    v = lred(Ln, v)
                    hits = [i for i, ex in enumerate(zeros) if lred(Ln, ex) == v]
                    if hits:
                        okq = [qq for qq in mods if P17[(qq, hits[0])]]
                        if okq:
                            nc[hits[0]] = okq[0]
                        else:
                            rest.append(v)
                        continue
                    if any(P17[(qq, i)] and lred(KER[qq], v) == lred(KER[qq], ex)
                           for qq in mods for i, ex in enumerate(zeros)):
                        continue
                    rest.append(v)
            score = (len(rest), lidx(Ln))
            if best is None or score < best[0]:
                best = (score, q, Ln, rest, nc)
        if best is None:
            break
        _, q, L, live, nc = best
        used.append(q)
        cert.update(nc)
    done = not live and sorted(cert) == list(range(len(zeros)))
    return done, used, len(live), lidx(L), {zinfo[i]['mn']: qq for i, qq in cert.items()}, P17


EXTRA_ODD = [19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61]      # odd prime powers beyond 02's pool
ALT = []
FAILED_ALT = []
for cx in CTX:
    orb = cx['orb']
    key, c, w = orb['key'], orb['c'], orb['w']
    O = FIELDS[key]
    pools = {
        'odd': [q for q in AVAILABLE[key] if q % 2],
        'pow2': [q for q in AVAILABLE[key] if q % 2 == 0],
    }
    for nm, pool in pools.items():
        done, used, nlive, idx, certq, P17g = greedy(key, w, cx['zinfo'], pool, max_steps=10, max_index=10 ** 40)
        p17q = {zz['mn']: [q for q in pool if P17g[(q, i)]] for i, zz in enumerate(cx['zinfo'])}
        # re-run the chosen sequence through the main sieve: same verdict
        KERa = {q: ker(key, q).H for q in used}
        CANDa = {q: ker(key, q).candidates(O, w) for q in used}
        sva = sieve(used, KERa, CANDa, cx['zeros'], {(q, i): P17g[(q, i)] for q in used for i in range(len(cx['zeros']))})
        check((not sva['live']) == done, f'greedy verdict vs sieve for {key} c={c} {nm}')
        ALT.append((key, c, nm, done, used, nlive, idx))
        print(f'   class {LABEL[key]} c = {c:2d}, {nm:4s}: moduli {used}; final index {idx}; live left {nlive}; '
              f'TERMINATES: {done}; certified at {certq}; Prop. 17 holds (Hermite basis) at {p17q}')
        if not done:
            FAILED_ALT.append((cx, nm, used))
print()
print('   The non-terminating runs, examined:')
for cx, nm, used in FAILED_ALT:
    orb = cx['orb']
    key, c, w = orb['key'], orb['c'], orb['w']
    O = FIELDS[key]
    if nm == 'odd':
        # append odd prime powers outside 02's pool (|G_q| <= 30000), in order of |G_q|, until the sieve terminates
        ext = []
        for q in EXTRA_ODD:
            if q in AVAILABLE[key]:
                continue
            S = closure(O, q, cap=30000)
            if S is not None:
                ext.append((len(S), q))
        ext.sort()
        seq = list(used)
        verdict = None
        for g_, q in ext:
            seq.append(q)
            KERa = {qq: ker(key, qq).H for qq in seq}
            CANDa = {qq: ker(key, qq).candidates(O, w) for qq in seq}
            P = {}
            for qq in seq:
                pq, eq, _ = prime_power(qq)
                for i, zz in enumerate(cx['zinfo']):
                    r = p17_data(O, qq, zz['z'], four_bases(KERa[qq])['reduced'])
                    P[(qq, i)] = p17_hypotheses(pq, eq, r['lam'], r['nJ'])['holds']
            sva = sieve(seq, KERa, CANDa, cx['zeros'], P)
            if not sva['live']:
                verdict = (q, g_, sva['steps'][len(used):])
                break
        check(verdict is not None, f'odd-only sieve for {key} c={c} does not terminate even with extra moduli')
        if verdict:
            print(f'   class {LABEL[key]} c = {c}, odd: stalls inside 02\'s pool (|G_q| <= 5000); appending odd moduli outside the '
                  f'pool in order of |G_q| ({[q for _, q in ext]}), it terminates once q = {verdict[0]} (|G_q| = {verdict[1]}) '
                  f'is added; steps after the pool: {verdict[2]}.  A pool-size effect, not a p-adic obstruction.')
    else:
        # powers of 2 only: the 2-adic tower of candidate classes
        zi = cx['zinfo']
        z16 = [lred(ker(key, 16).H, zz['exp']) for zz in zi]
        cnt = []
        for k in range(2, 10):
            kq = ker(key, 2 ** k)
            cand = kq.candidates(O, w)
            outside = [v for v in cand if k < 4 or lred(ker(key, 16).H, v) not in z16]
            cnt.append((2 ** k, kq.H, len(cand), len(outside)))
        print(f'   class {LABEL[key]} c = {c}, pow2: candidate classes mod Lambda_(2^k) (Hermite form; all; outside the '
              f'Lambda_16-class of a known zero): {cnt}')
        # lift the candidates mod Lambda_64 along (a + 2^k i, b + 2^(k-1) j) to 2^K (shape of Lambda_(2^k) seen above)
        shape = all(ker(key, 2 ** k).H == (2 ** k, 0, 2 ** (k - 1)) for k in range(2, 10))
        if not shape:
            print('      (the lifting below assumes Lambda_(2^k) = 2^k Z x 2^(k-1) Z, which fails here: skipped)')
            continue
        K = 48
        cl = sorted(ker(key, 64).candidates(O, w))
        for k in range(6, K):
            q2 = 2 ** (k + 1)
            cl = [(a + (2 ** k) * i, b + (2 ** (k - 1)) * j) for (a, b) in cl for i in (0, 1) for j in (0, 1)
                  if O.mulq(O.unit_mod(a + (2 ** k) * i, b + (2 ** (k - 1)) * j, q2), w, q2)[2:] == (0, 0)]
        clusters = {}
        for (a, b) in cl:
            clusters.setdefault((a % 64, b % 32), []).append((a, b))
        print(f'      Lambda_(2^k) = 2^k Z x 2^(k-1) Z for k = 2..9: {shape}.  Lifting the 16 candidates mod Lambda_64 '
              f'level by level to 2^{K}: {len(cl)} candidates remain, in {len(clusters)} clusters '
              f'(exponents mod (64,32)): {[(k2, len(v)) for k2, v in clusters.items()]}')
        qK = 2 ** K
        for k2, v in clusters.items():
            a, b = v[0]
            x = O.mulq(O.unit_mod(a, b, qK), w, qK)
            n0, m0 = x[0], (-x[1]) % qK
            t = n0 * pow(m0, -1, qK) % qK if m0 % 2 else None
            sv1 = lagrange((1, t), (0, qK))[0] if t is not None else None
            known_here = [zz['mn'] for zz in zi if (zz['exp'][0] % 64, zz['exp'][1] % 32) == k2]
            print(f'      cluster {k2}: exponent == {(a % 2 ** 24, b % 2 ** 24)} mod 2^24; the element is n - m beta with '
                  f'n/m == {t % 2 ** 24 if t is not None else None} mod 2^24; shortest (m, n) with n == (n/m) m mod 2^{K}: '
                  f'{sv1}; known zero here: {known_here or "none"}')
        # (round 2) the lifting uses classes mod 2^k Z x 2^(k-1) Z; check that this lattice lies inside Lambda_(2^k)
        sub = all(O.unit_mod(2 ** k, 0, 2 ** k) == (1, 0, 0, 0) and O.unit_mod(0, 2 ** (k - 1), 2 ** k) == (1, 0, 0, 0)
                  for k in range(2, K + 1))
        check(sub, 'eps1^(2^k) == eps2^(2^(k-1)) == 1 mod 2^k fails for some k')
        z4 = [lred(ker(key, 4).H, zz['exp']) for zz in zi]
        bare = [v for k2, v in clusters.items() if not any((zz['exp'][0] % 64, zz['exp'][1] % 32) == k2 for zz in zi)]
        apart = all(lred(ker(key, 4).H, e) not in z4 for v in bare for e in v)
        check(bare and apart, 'the cluster without a known zero is not separated from the known zeros mod Lambda_4')
        print(f'      eps1^(2^k) == eps2^(2^(k-1)) == 1 mod 2^k for k = 2..{K}: {sub}, so each class mod 2^k Z x 2^(k-1) Z above')
        print(f'      lies in one class mod Lambda_(2^k).  The candidates without a known zero lie outside the Lambda_4-class of')
        print(f'      the known zero: {apart}.')
        print(f'      => no sieve with moduli 2^k, k <= {K}, finishes this orbit; odd moduli do (02 uses 43 and 13).  The')
        print('         candidates converge to a 2-adic zero of the markers that is not an integer zero (n/m has no small')
        print('         rational reconstruction): check_2adic_zero.py gives the Hensel certificate (round 2).')
print()

print('=' * 112)
print('NEGATIVE CONTROL (c): -1 and the phantom classes')
for key in ('10', '01', '11'):
    neg = NEGTAB[key]
    print(f'   class {LABEL[key]}: -1 is a residue of <eps1,eps2> exactly at the moduli '
          f'{[q for q in neg if neg[q]]}; not at {[q for q in neg if not neg[q]]}')
nophantom = []
for cx in CTX:
    orb, KER, CAND, zeros, zinfo, P17 = cx['orb'], cx['KER'], cx['CAND'], cx['zeros'], cx['zinfo'], cx['P17']
    key, c = orb['key'], orb['c']
    O = FIELDS[key]
    resid = [q for q in orb['used'] if ker(key, q).minus_one]
    if not resid:
        nophantom.append(f'{LABEL[key]} c={c} {orb["used"]}')
        continue
    print(f'   orbit class {LABEL[key]} c = {c}: moduli used at which -1 IS a residue: {resid}')
    for q in resid:
        s = ker(key, q).star
        if not zeros:
            print(f'      q = {q}: -1 == eps^{s} mod {q}, but the orbit has no known zero, so there is no phantom class '
                  f'(candidate classes mod Lambda_{q}: {len(CAND[q])})')
            continue
        for i, ex in enumerate(zeros):
            ph = lred(KER[q], (ex[0] + s[0], ex[1] + s[1]))
            isc = ph in CAND[q]
            own = lred(KER[q], ex)
            other = [zinfo[j]['mn'] for j, ex2 in enumerate(zeros) if lred(KER[q], ex2) == ph]
            check(isc and ph != own, f'phantom of {zinfo[i]["mn"]} at q = {q}')
            # fate of the phantom in the direct check: all classes mod the final lattice above it
            Lf = cx['sv']['L']
            fates = {}
            for v in box_reps(Lf):
                if lred(KER[q], v) != ph:
                    continue
                why = None
                for q2 in orb['used']:
                    if lred(KER[q2], v) not in CAND[q2]:
                        why = f'not a candidate mod {q2}'
                        break
                if why is None:
                    for q2 in orb['used']:
                        for j, ex2 in enumerate(zeros):
                            if P17[(q2, j)] and lred(KER[q2], v) == lred(KER[q2], ex2):
                                why = f'inside the Lambda_{q2}-class of {zinfo[j]["mn"]} (Prop. 17 at {q2})'
                if why is None:
                    why = 'UNACCOUNTED'
                fates[why] = fates.get(why, 0) + 1
            check('UNACCOUNTED' not in fates, 'phantom survives')
            print(f'      q = {q}: -1 == eps^{s} mod {q}; phantom of {zinfo[i]["mn"]} = class {ph} mod Lambda_{q} '
                  f'(the class of -z), a candidate: {isc}; it contains a known zero: {other or "no"}')
            print(f'         its {sum(fates.values())} subclasses mod the final lattice {Lf}: {fates}')
            # the phantom in the incremental sieve, in 02's order and reversed
            for nm, seq in (('02\'s order', orb['used']), ('reversed order', list(reversed(orb['used'])))):
                rec = []
                sieve(seq, KER, CAND, zeros, P17, record=rec)
                trace = []
                seen_q = False
                for r_ in rec:
                    seen_q = seen_q or r_['q'] == q
                    if not seen_q:
                        trace.append(f'after {r_["q"]}: (q = {q} not applied yet)')
                        continue
                    above = lambda lst: sum(1 for v in lst if lred(KER[q], v) == ph)
                    trace.append(f'after {r_["q"]}: {above(r_["cand"])} candidate class(es) above the phantom, '
                                 f'{above(r_["cert"])} certified, {above(r_["sib"])} sibling-empty, {above(r_["left"])} live')
                print(f'         {nm} {seq}: ' + '; '.join(trace))
            # is -1 a residue modulo the lcm of the moduli used?
            lcm = 1
            for qq in orb['used']:
                lcm = lcm * qq // gcd(lcm, qq)
            print(f'         -1 is a residue modulo lcm = {lcm} of the moduli used: {ker(key, lcm).minus_one}')
print(f'   the other {len(nophantom)} orbits use no modulus at which -1 is a residue, so no phantom class arises: '
      f'{"; ".join(nophantom)}')
print()

# ======================================================================================================
# NEGATIVE CONTROL (b): h1 Z x h2 Z in place of Lambda_q
# ======================================================================================================
print('=' * 112)
print('NEGATIVE CONTROL (b): h1 Z x h2 Z (h_i = order of eps_i mod q) in place of Lambda_q')
for cx in CTX:
    orb, zinfo = cx['orb'], cx['zinfo']
    key, c, w = orb['key'], orb['c'], orb['w']
    O = FIELDS[key]
    for q in orb['used']:
        kq = ker(key, q)
        h1, h2 = kq.h1, kq.h2
        ratio = h1 * h2 // kq.order
        A, B, C = kq.H
        # a vector of Lambda_q outside h1 Z x h2 Z
        wit = (B, C)
        inprod = (wit[0] % h1 == 0 and wit[1] % h2 == 0)
        okw = O.unit_mod(B, C, q) == (1, 0, 0, 0) and (not inprod) == (ratio > 1)
        check(okw and h1 * h2 == ratio * kq.order, f'h1 Z x h2 Z control at q = {q}')
        # candidates counted over the box [0,h1) x [0,h2): one per class of h1 Z x h2 Z
        ncand_box = 0
        u1 = tuple(v % q for v in O.e1)
        u2 = tuple(v % q for v in O.e2)
        y = tuple(v % q for v in w)
        zcount = [0] * len(zinfo)
        zcls = [lred(kq.H, zz['exp']) for zz in zinfo]       # class mod Lambda_q of each known zero
        for b in range(h2):
            x = y
            for a in range(h1):
                if x[2] == 0 and x[3] == 0:
                    ncand_box += 1
                cl = lred(kq.H, (a, b))
                for i, zc in enumerate(zcls):
                    if cl == zc:
                        zcount[i] += 1
                x = O.mulq(x, u1, q)
            y = O.mulq(y, u2, q)
        okc = ncand_box == ratio * len(cx['CAND'][q]) and all(n == ratio for n in zcount)
        check(okc, f'overcount at q = {q}')
        print(f'   class {LABEL[key]} c = {c:2d}, q = {q:3d}: h1 = {h1}, h2 = {h2}, h1*h2 = {h1 * h2} vs |G_q| = {kq.order} '
              f'(index [Lambda_q : h1 Z x h2 Z] = {ratio}); witness {wit} in Lambda_q (eps^{wit} == 1 mod {q}) '
              f'{"outside" if not inprod else "inside"} h1 Z x h2 Z; candidate pairs in [0,h1)x[0,h2): {ncand_box} = '
              f'{ratio} x {len(cx["CAND"][q])} classes mod Lambda_q; each known zero counted {zcount} times')
print('   With h1 Z x h2 Z, |G_q| would be overcounted by the index above, every candidate class (and every known')
print('   zero) would be counted index-many times, and a Prop. 17 coset system on eta_1 = eps1^h1, eta_2 = eps2^h2')
print('   would cover only one of the index-many subclasses of a Lambda_q-class, so its conclusion could not be used')
print('   for the whole Lambda_q-class (the certified and sibling rules of 4.4 would then be unsound).')
print()

# ======================================================================================================
# NEGATIVE CONTROL (a): planted false certificates
# ======================================================================================================
print('=' * 112)
print('NEGATIVE CONTROL (a): planted false certificates must be caught')
na_ok = True
# (a1) claim a non-zero is a zero: a candidate class without a known zero, at the first modulus used
for cx in CTX:
    orb, zinfo, zeros = cx['orb'], cx['zinfo'], cx['zeros']
    key, c, w = orb['key'], orb['c'], orb['w']
    O = FIELDS[key]
    q = orb['used'][0]
    fake = sorted(v for v in cx['CAND'][q] if all(lred(cx['KER'][q], ex) != v for ex in zeros))
    if not fake:
        continue
    v = fake[0]
    x = O.mul(O.unit(*v), w)
    caught = (x[2] % q == 0 and x[3] % q == 0) and (x[2], x[3]) != (0, 0)
    refused = p17_data(O, q, x, four_bases(cx['KER'][q])['02']) is None
    na_ok &= caught and refused
    show = lambda n: str(n) if abs(n) < 10 ** 15 else f'<{len(str(abs(n)))}-digit integer, == {n % q} mod {q}>'
    print(f'   (a1) class {LABEL[key]} c = {c}: plant exponent {v} (a candidate mod {q}: markers == 0 mod {q}) as a '
          f'"known zero": exact markers = ({show(x[2])}, {show(x[3])}) != 0 -> rejected: {caught}; Prop. 17 refuses this '
          f'seed: {refused}')
# (a2) drop a known zero: the sieve must fail
for cx in CTX:
    orb, zinfo, zeros, P17 = cx['orb'], cx['zinfo'], cx['zeros'], cx['P17']
    if not zeros:
        continue
    for i in range(len(zeros)):
        z2 = [ex for k, ex in enumerate(zeros) if k != i]
        P2 = {}
        for (q, k), v in P17.items():
            if k != i:
                P2[(q, k - (1 if k > i else 0))] = v
        sv = sieve(orb['used'], cx['KER'], cx['CAND'], z2, P2)
        _, _, bad = direct_check(sv['L'], orb['used'], cx['KER'], cx['CAND'], z2, P2)
        dead_zero_live = lred(sv['L'], zeros[i]) in sv['live']
        caught = len(sv['live']) >= 1 and len(bad) >= 1 and dead_zero_live
        na_ok &= caught
        print(f'   (a2) class {LABEL[orb["key"]]} c = {orb["c"]}: drop the known zero {zinfo[i]["mn"]}: sieve ends with '
              f'{len(sv["live"])} live class(es) (the dropped zero\'s class among them: {dead_zero_live}); direct check: '
              f'{len(bad)} unaccounted -> caught: {caught}')
# (a3) forged Proposition 17 instances: every (q, zero) at which the hypotheses fail, declared to hold
for cx in CTX:
    orb, zinfo, zeros, P17 = cx['orb'], cx['zinfo'], cx['zeros'], cx['P17']
    failing = [(q, i) for (q, i), v in P17.items() if not v]
    if not failing:
        continue
    O = FIELDS[orb['key']]
    for (q, i) in failing:
        p, e = prime_power(q)[:2]
        r = p17_data(O, q, zinfo[i]['z'], four_bases(cx['KER'][q])['02'])
        h = p17_hypotheses(p, e, r['lam'], r['nJ'])
        okm, _ = p17_mechanism(O, q, zinfo[i]['z'], four_bases(cx['KER'][q])['02'])
        na_ok &= (not h['holds']) and (not okm)
        print(f'   (a3) class {LABEL[orb["key"]]} c = {orb["c"]}, q = {q}, zero {zinfo[i]["mn"]}: lambda = {fnu(r["lam"])}, '
              f'nu(J) = {fnu(r["nJ"])}, delta = {fnu(h["delta"])}, bound {h["bound"]}: a forged "Prop. 17 holds" is '
              f'rejected by the hypothesis check: {not h["holds"]}, and by the mechanism check: {not okm}')
    # what a forged certificate would buy: the shortest prefix of 02's order that "terminates" with all forged
    Pf = {k: True for k in P17}
    for k in range(1, len(orb['used'])):
        pre = orb['used'][:k]
        sv = sieve(pre, cx['KER'], cx['CAND'], zeros, Pf)
        if not sv['live']:
            print(f'      with every instance forged, the prefix {pre} of 02\'s order would already "certify" the orbit '
                  f'(its conclusion happens to be true, the certificate is not): caught by the hypothesis check above')
            break
# (a4) the sibling rule at a q where Prop. 17 does not hold: class (1,1), c = 1, zero (2,1) dropped, q = 2 added
cx = next(cx for cx in CTX if cx['orb']['key'] == '11' and cx['orb']['c'] == 1)
O = FIELDS['11']
zinfo, zeros = cx['zinfo'], cx['zeros']
k2 = ker('11', 2)
d = (zeros[1][0] - zeros[0][0], zeros[1][1] - zeros[0][1])
u = O.unit(*d)
gd = 0
for t in ((u[0] - 1,) + u[1:]):
    gd = gcd(gd, t)
KER = dict(cx['KER'])
CAND = dict(cx['CAND'])
KER[2] = k2.H
CAND[2] = k2.candidates(O, cx['orb']['w'])
z1 = [zeros[0]]
P1 = {(q, 0): cx['P17'][(q, 0)] for q in cx['orb']['used']}
P1[(2, 0)] = False
seq = [2] + cx['orb']['used']
s_wrong = sieve(seq, KER, CAND, z1, P1, sibling='any')
s_right = sieve(seq, KER, CAND, z1, P1, sibling='same')
bite = (not s_wrong['live']) and len(s_right['live']) >= 1 and gd == 2
na_ok &= bite
print(f'   (a4) class (1,1) c = 1: the zeros {zinfo[0]["mn"]} (exponent {zeros[0]}) and {zinfo[1]["mn"]} (exponent {zeros[1]}) '
      f'differ by {d}, and eps^{d} - 1 = {tuple((u[0] - 1,) + u[1:])} has content {gd}: they share one Lambda_2-class '
      f'(Lambda_2 = {k2.H}).')
print(f'        Drop {zinfo[1]["mn"]}; moduli {seq}.  WRONG sibling rule (Lambda_q-class for one q, Prop. 17 at another): '
      f'live left {len(s_wrong["live"])} -> would falsely certify "zeros = {{{zinfo[0]["mn"]}}}"; correct same-q rule: live '
      f'left {len(s_right["live"])} -> caught.  The control bites: {bite}')
# (a5) the same with a forged Prop. 17 instance at q = 2 (p = 2, e = 1 violates e >= 2) and the CORRECT same-q rule
r2 = p17_data(O, 2, zinfo[0]['z'], four_bases(k2.H)['02'])
h2 = p17_hypotheses(2, 1, r2['lam'], r2['nJ'])
P5 = dict(P1)
P5[(2, 0)] = True                                    # forged
s_forged = sieve(seq, KER, CAND, z1, P5, sibling='same')
bite5 = (not s_forged['live']) and (not h2['holds']) and (not h2['setting'])
na_ok &= bite5
print(f'   (a5) same orbit, (2,1) dropped, same-q sibling rule, but Prop. 17 FORGED at q = 2 for {zinfo[0]["mn"]} (recomputed: '
      f'lambda = {fnu(r2["lam"])}, nu(J) = {fnu(r2["nJ"])}, bound e-1-iota = {h2["bound"]}, setting e >= 2 for p = 2: '
      f'{h2["setting"]}): live left {len(s_forged["live"])} -> a false certificate with a FALSE conclusion (the genuine '
      f'solution (2,1) is declared absent); the hypothesis check rejects the instance: {not h2["holds"]}.  Bites: {bite5}')
check(na_ok, 'negative controls (a)')
print(f'   all planted false certificates caught: {na_ok}')
print()

# ======================================================================================================
# NEGATIVE CONTROL (d): basis dependence where nu(J) >= e
# ======================================================================================================
print('=' * 112)
print('NEGATIVE CONTROL (d): lambda and nu(J) in the four bases at every (q, zero) with q among the moduli used')
dep_found = []
for cx in CTX:
    orb, zinfo = cx['orb'], cx['zinfo']
    O = FIELDS[orb['key']]
    for q in orb['used']:
        p, e = prime_power(q)[:2]
        for i, zz in enumerate(zinfo):
            res = {nm: p17_data(O, q, zz['z'], bb) for nm, bb in four_bases(cx['KER'][q]).items()}
            lams = sorted({r['lam'] for r in res.values()})
            nJs = sorted({r['nJ'] for r in res.values()})
            if len(lams) > 1 or len(nJs) > 1:
                dep_found.append((orb['key'], orb['c'], q, zz['mn'], lams, nJs, e))
for (key, c, q, mn, lams, nJs, e) in dep_found:
    print(f'   class {LABEL[key]} c = {c}, q = {q} (e = {e}), zero {mn}: lambda in {lams}, nu(J) in {nJs} -- all >= e, '
          f'and Prop. 17 fails there in every basis')
    check(min(lams) >= e or len(lams) == 1, 'basis dependence below e (lambda)')
    check(min(nJs) >= e, 'basis dependence below e (nu(J))')
if not dep_found:
    print('   no basis dependence at all')
print()

# ======================================================================================================
# Box scan (sanity only)
# ======================================================================================================
print('=' * 112)
print(f'BOX SCAN (sanity only, not part of the proof): |a|,|b| <= {BOX_R}, both markers == 0 mod {BOX_PRIMES[0]} and mod '
      f'{BOX_PRIMES[1]}, then exact')
for cx in CTX:
    orb = cx['orb']
    key, c, w = orb['key'], orb['c'], orb['w']
    O = FIELDS[key]
    R = BOX_R
    hits = None
    for P in BOX_PRIMES:
        E1 = np.zeros((2 * R + 1, 4), dtype=np.int64)
        cur = O.pwq(O.e1, -R, P)
        for i in range(2 * R + 1):
            E1[i] = cur
            cur = O.mulq(cur, O.e1, P)
        y = O.mulq(O.pwq(O.e2, -R, P), w, P)
        S = set()
        for j in range(2 * R + 1):
            My = np.array([[v % P for v in row] for row in O.mat(y)], dtype=np.int64)
            Z = (E1 @ My.T) % P
            for i in np.nonzero((Z[:, 2] == 0) & (Z[:, 3] == 0))[0]:
                S.add((int(i) - R, j - R))
            y = O.mulq(y, O.e2, P)
        hits = S if hits is None else hits & S
    exact = sorted(v for v in hits if O.mul(O.unit(*v), w)[2:] == (0, 0))
    check(exact == sorted(cx['zeros']), f'box scan {exact} vs {cx["zeros"]}')
    print(f'   class {LABEL[key]} c = {c:2d}: common zeros mod both primes in the box: {sorted(hits)}; exact: {exact}; '
          f'known exponents: {sorted(cx["zeros"])}')
print()

RC = 40
print('=' * 112)
print(f'CLASS SCAN (sanity only): for every Prop. 17 instance used by 02, the coset system T1^i T2^j z covering the '
      f'Lambda_q-class of the zero (02 basis), |i|,|j| <= {RC}; both markers == 0 mod {BOX_PRIMES[0]} and mod '
      f'{BOX_PRIMES[1]}')
for cx in CTX:
    orb = cx['orb']
    O = FIELDS[orb['key']]
    for (q, p_, e_, mn, *_r) in orb['inst']:
        i0 = next(k for k, zz in enumerate(cx['zinfo']) if zz['mn'] == mn)
        z = cx['zinfo'][i0]['z']
        (a1, b1), (a2, b2) = four_bases(cx['KER'][q])['02']
        hits = None
        for P in BOX_PRIMES:
            t1, t2 = O.unit_mod(a1, b1, P), O.unit_mod(a2, b2, P)
            x = O.mulq(O.unit_mod(-RC * a1, -RC * b1, P), z, P)
            t2m = O.unit_mod(-RC * a2, -RC * b2, P)
            S = set()
            for i in range(-RC, RC + 1):
                y = O.mulq(x, t2m, P)
                for j in range(-RC, RC + 1):
                    if y[2] == 0 and y[3] == 0:
                        S.add((i, j))
                    y = O.mulq(y, t2, P)
                x = O.mulq(x, t1, P)
            hits = S if hits is None else hits & S
        okcs = hits == {(0, 0)}
        check(okcs, f'class scan at q = {q}, zero {mn}: {sorted(hits)}')
        print(f'   class {LABEL[orb["key"]]} c = {orb["c"]:2d}, q = {q:2d}, zero {mn}: {(2 * RC + 1) ** 2} points of the class '
              f'(|a| up to {RC * (abs(a1) + abs(a2))}); common zeros mod both primes at (i,j) = {sorted(hits)} only: {okcs}')
print()

# ======================================================================================================
# Summary
# ======================================================================================================
print('=' * 112)
print('PROPOSITION 17 INSTANCES USED BY 02 (all recomputed exactly; values in the 02 basis, equal in all four bases)')
print('   class  c   q   p  e  zero      lambda nu(J) delta bound(e-1-iota)  Lambda(2)  e+delta  e+lambda  margin(delta)')
for (key, c, q, p, e, mn, lam, nJ, dl, bd, L2, ed, el) in INSTANCES:
    print(f'   {LABEL[key]} {c:2d} {q:3d} {p:3d} {e:2d}  {str(mn):9s} {lam:5d} {nJ:5d} {dl:5d} {bd:8d}        {L2:5d}     '
          f'{ed:5d}    {el:5d}     {bd - dl:5d}')
nb = sum(1 for t in INSTANCES if t[9] == t[8])
lam_vals = sorted({t[6] for t in INSTANCES})
strict_all = all(t[11] < t[10] and t[12] < t[10] for t in INSTANCES)
check(strict_all, 'strict inequalities')
print(f'   {len(INSTANCES)} instances; lambda takes the values {lam_vals} (an instance with lambda = 1 is used: '
      f'{1 in lam_vals}); delta = e-1-iota exactly in {nb} of them, slack in {len(INSTANCES) - nb}; strict inequalities '
      f'e+delta < Lambda(2) and e+lambda < Lambda(2) hold in all: {strict_all}')
print(f'   four bases (Hermite, 02, unimodular, reduced) at all {BASISSTAT["pairs"]} (q, zero) pairs of the moduli used: '
      f'lambda and nu(J) identical in {BASISSTAT["agree"]}; Prop. 17 holds at {BASISSTAT["holds"]} pairs (verdict '
      f'basis-independent at every pair); values differ only where they are >= e (control (d))')
print()
print('SUMMARY')
for key, c, w, nz, ok in SUMMARY:
    print(f'   class {LABEL[key]} c = {c:2d} seed {w}: {nz} known solution pair(s); independently confirmed: {ok}')
nt = [(LABEL[k], c, nm) for (k, c, nm, d, u, nl, ix) in ALT if not d]
print(f'   reversed order of 02\'s moduli: terminates for every orbit with more than one modulus')
print(f'   greedy over other sets (pool of 02, |G_q| <= 5000): {sum(1 for t in ALT if t[3])} of {len(ALT)} runs terminate; '
      f'not: {nt} (examined above)')
print(f'checks run: {NCHECK[0]}')
if FAILS:
    print(f'CHECK FAILED: {len(FAILS)} failed check(s)')
else:
    print('ALL 12 ORBITS INDEPENDENTLY CONFIRMED; no discrepancy with 02_exclusion.out.txt or 01_coverage.out.txt')
print(f'[done in {time.time() - T0:.1f} s]')
sys.exit(1 if FAILS else 0)
