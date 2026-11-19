#!/usr/bin/env python3
"""verify_exclusion.py -- independent re-implementation of the exclusion step (02_exclusion.py).

Independent verifier (Claude subagent), 2026-10-15.  Nothing is imported from core.py, lifting2.py or
02_exclusion.py.  The data are retyped: forms, integral bases and units from core.py; seeds, known exponents,
moduli sequences and every printed count from 02_exclusion.out.txt (seeds also from 01_coverage.out.txt).
All arithmetic is written afresh, by different routes from the audited scripts:
  * the order O: structure constants from sympy polynomial arithmetic in Q[t]/(f) and sympy's inverse of the
    basis matrix, checked integral; norms both as determinants (sympy) and as resultants Res(f, x(t));
  * G_q: breadth-first closure of {1} under the matrices M(eps1), M(eps2) acting on (O/qO) coordinates;
    Lambda_q: Hermite form of the relations met in that search, CERTIFIED by (i) each basis vector maps to 1
    mod qO (direct modular exponentiation, negative exponents through the exact inverse) and (ii) its index
    equals |G_q| (so it is the whole kernel of Z^2 -> G_q);
  * lattice intersections through dual lattices, L cap M = (L* + M*)*, CERTIFIED by containment and the index
    identity [Z^2 : L cap M] [Z^2 : L + M] = [Z^2 : L] [Z^2 : M] (L + M also by sympy's Hermite form);
  * Proposition 17 with exact integers (no truncation), in three bases of Lambda_q: own Hermite basis, the
    script's basis (h1,0),(-a1,b0), and a random unimodular change;
  * the sieve in the printed order of moduli (certified and sibling rules), compared count by count;
  * an independent direct check of every class of Z^2 modulo the final lattice (numpy, vectorized);
  * a brute-force scan of the box |a|,|b| <= R (two primes, then exact): a sanity check only, not a proof.
Every discrepancy prints CHECK FAILED; the exit status is 1 if there is any.
"""
import sys
import time
import random
from fractions import Fraction as Fr
from math import gcd

import numpy as np
from sympy import Matrix, Poly, Rational, symbols, QQ, isprime, real_roots
from sympy.matrices.normalforms import hermite_normal_form

T0 = time.time()
FAILS = []
INF = 10 ** 9          # stands for nu(0) = +infinity


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print('      CHECK FAILED: ' + msg)
    return bool(cond)


# ======================================================================================================
# DATA (retyped)
# ======================================================================================================
# F(m,n) = a0 m^4 + a1 m^3 n + a2 m^2 n^2 + a3 m n^3 + a4 n^4 ; f(t) = F(1,t)
FORMS = {'10': (36, 0, 0, -6, 1), '01': (-5, 14, 24, 2, 1), '11': (-8, 4, 24, 4, 1)}
LABEL = {'10': '(1,0)', '01': '(0,1)', '11': '(1,1)'}
# rows = omega_i in the power basis 1, alpha, alpha^2, alpha^3
BASES = {
    '10': [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, Fr(1, 6), 0], [0, 0, 0, Fr(1, 6)]],
    '01': [[1, 0, 0, 0], [0, 1, 0, 0], [Fr(1, 6), Fr(2, 3), Fr(1, 6), 0], [Fr(1, 3), Fr(1, 2), 0, Fr(1, 6)]],
    '11': [[1, 0, 0, 0], [0, 1, 0, 0], [Fr(2, 3), Fr(1, 3), Fr(1, 6), 0], [Fr(1, 3), 0, 0, Fr(1, 12)]],
}
UNITS = {
    '10': ((0, 1, -1, 0), (-86, -41, -117, 30)),
    '01': ((0, -1, 1, 0), (-21, 36, 16, 14)),
    '11': ((-1, 0, 1, 0), (1, -3, -5, -2)),
}
KNOWN = {'10': [(0, 1), (1, 0), (1, 2), (1, 6)], '01': [(-1, 1), (0, 1), (1, 1), (4, 1)], '11': [(-1, 1), (0, 1), (2, 1)]}
SEEDS_01 = {   # as printed by 01_coverage.out.txt
    '10': {1: [(-1, 0, 0, 0)], 4: [(2, -1, 0, 0)], 9: [(-3, 0, 6, -1)], 36: [(0, 0, 6, -1)]},
    '01': {1: [(-1, 0, 0, 0)], 4: [(-1, -1, 0, 0)], 9: [(-3, -2, -3, -1)], 36: [(1, -1, 0, 0)]},
    '11': {1: [(-1, 0, 0, 0)], 4: [(3, -1, -2, -1), (2, -1, -1, -1), (0, 0, -1, 0)], 9: [(-1, -1, 0, 0)],
           36: [(6, -2, -3, -2), (-1, -2, -3, -1), (0, -1, -3, -1)]},
}
# Everything 02_exclusion.out.txt prints per orbit.
#   steps: (q, candidate classes mod new lattice, certified, sibling-empty, left)
#   cert:  (printed class representative, solution (m,n), q, lambda, nu(J), e)
#   G:     q -> (|G_q|, candidate classes mod Lambda_q)
ORBITS = [
    dict(key='10', c=1, w=(-1, 0, 0, 0), known=[(0, 1)], exps=[(0, 0)], used=[4], index=8,
         steps=[(4, 1, 1, 0, 0)], cert=[((0, 0), (0, 1), 4, 0, 0, 2)], G={4: (8, 1)}),
    dict(key='10', c=4, w=(2, -1, 0, 0), known=[(1, 2)], exps=[(0, 0)], used=[4, 31, 7], index=36000,
         steps=[(4, 2, 0, 0, 2), (31, 4, 1, 1, 2), (7, 0, 0, 0, 0)], cert=[((0, 0), (1, 2), 31, 0, 0, 1)],
         G={4: (8, 2), 31: (450, 4), 7: (200, 8)}),
    dict(key='10', c=9, w=(-3, 0, 6, -1), known=[], exps=[], used=[4], index=8,
         steps=[(4, 0, 0, 0, 0)], cert=[], G={4: (8, 0)}),
    dict(key='10', c=36, w=(0, 0, 6, -1), known=[(1, 0), (1, 6)], exps=[(-1, 0), (1, 0)], used=[4, 31], index=1800,
         steps=[(4, 2, 0, 0, 2), (31, 4, 2, 2, 0)],
         cert=[((1, 0), (1, 6), 31, 0, 0, 1), ((119, 0), (1, 0), 31, 0, 0, 1)], G={4: (8, 2), 31: (450, 4)}),
    dict(key='01', c=1, w=(-1, 0, 0, 0), known=[(0, 1)], exps=[(0, 0)], used=[4, 7, 3], index=576,
         steps=[(4, 1, 0, 0, 1), (7, 2, 1, 0, 1), (3, 0, 0, 0, 0)], cert=[((0, 0), (0, 1), 7, 0, 0, 1)],
         G={4: (8, 1), 7: (144, 4), 3: (36, 6)}),
    dict(key='01', c=4, w=(-1, -1, 0, 0), known=[(-1, 1)], exps=[(0, 0)], used=[4, 31, 8], index=3840,
         steps=[(4, 2, 0, 0, 2), (31, 2, 1, 0, 1), (8, 0, 0, 0, 0)], cert=[((0, 0), (-1, 1), 31, 0, 0, 1)],
         G={4: (8, 2), 31: (3840, 4), 8: (32, 4)}),
    dict(key='01', c=9, w=(-3, -2, -3, -1), known=[(4, 1)], exps=[(2, 1)], used=[4, 8, 31, 3], index=11520,
         steps=[(4, 1, 0, 0, 1), (8, 2, 0, 0, 2), (31, 2, 1, 0, 1), (3, 0, 0, 0, 0)],
         cert=[((2, 1), (4, 1), 31, 0, 0, 1)], G={4: (8, 1), 8: (32, 2), 31: (3840, 6), 3: (36, 18)}),
    dict(key='01', c=36, w=(1, -1, 0, 0), known=[(1, 1)], exps=[(0, 0)], used=[4, 31, 8], index=3840,
         steps=[(4, 2, 0, 0, 2), (31, 2, 1, 0, 1), (8, 0, 0, 0, 0)], cert=[((0, 0), (1, 1), 31, 0, 0, 1)],
         G={4: (8, 2), 31: (3840, 6), 8: (32, 4)}),
    dict(key='11', c=1, w=(-1, 0, 0, 0), known=[(0, 1)], exps=[(0, 0)], used=[9, 3, 27, 13], index=118584,
         steps=[(9, 2, 1, 0, 1), (3, 1, 0, 0, 1), (27, 1, 0, 0, 1), (13, 1, 0, 1, 0)],
         cert=[((0, 0), (0, 1), 9, 0, 0, 2)], G={9: (216, 2), 3: (72, 6), 27: (1944, 2), 13: (2196, 20)}),
    dict(key='11', c=4, w=(3, -1, -2, -1), known=[], exps=[], used=[4], index=16,
         steps=[(4, 0, 0, 0, 0)], cert=[], G={4: (16, 0)}),
    dict(key='11', c=4, w=(2, -1, -1, -1), known=[], exps=[], used=[4], index=16,
         steps=[(4, 0, 0, 0, 0)], cert=[], G={4: (16, 0)}),
    dict(key='11', c=4, w=(0, 0, -1, 0), known=[], exps=[], used=[4], index=16,
         steps=[(4, 0, 0, 0, 0)], cert=[], G={4: (16, 0)}),
    dict(key='11', c=9, w=(-1, -1, 0, 0), known=[(-1, 1), (2, 1)], exps=[(0, 0), (1, -1)], used=[4, 8, 16, 41, 25],
         index=4166400,
         steps=[(4, 4, 0, 0, 4), (8, 4, 0, 0, 4), (16, 8, 1, 0, 7), (41, 10, 1, 1, 8), (25, 4, 0, 4, 0)],
         cert=[((0, 0), (-1, 1), 16, 0, 2, 4), ((2281, 39), (2, 1), 41, 0, 0, 1)],
         G={4: (16, 4), 8: (32, 4), 16: (128, 8), 41: (4200, 6), 25: (3100, 28)}),
    dict(key='11', c=36, w=(6, -2, -3, -2), known=[], exps=[], used=[4], index=16,
         steps=[(4, 0, 0, 0, 0)], cert=[], G={4: (16, 0)}),
    dict(key='11', c=36, w=(-1, -2, -3, -1), known=[], exps=[], used=[4], index=16,
         steps=[(4, 0, 0, 0, 0)], cert=[], G={4: (16, 0)}),
    dict(key='11', c=36, w=(0, -1, -3, -1), known=[], exps=[], used=[4], index=16,
         steps=[(4, 0, 0, 0, 0)], cert=[], G={4: (16, 0)}),
]
# "moduli available (|G_q|)" lines of 02_exclusion.out.txt (moduli with |G_q| <= 5000)
AVAILABLE = {
    '10': {4: 8, 8: 32, 16: 128, 32: 512, 64: 2048, 3: 36, 9: 108, 27: 972, 5: 124, 25: 3100, 7: 200, 49: 1400,
           11: 732, 13: 2196, 17: 4912, 19: 3429, 31: 450},
    '01': {4: 8, 8: 32, 16: 128, 32: 512, 64: 2048, 3: 36, 9: 108, 27: 972, 5: 124, 25: 3100, 7: 144, 11: 300,
           13: 2196, 17: 4912, 19: 3429, 31: 3840},
    '11': {4: 16, 8: 32, 16: 128, 32: 512, 64: 2048, 3: 72, 9: 216, 27: 1944, 5: 124, 25: 3100, 7: 200, 11: 1464,
           13: 2196, 17: 4912, 41: 4200},
}
MODULI_ALL = [4, 8, 16, 32, 64, 3, 9, 27, 81, 5, 25, 7, 49, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
MAXG = 5000
BOX_R = 400                         # brute-force box |a|,|b| <= BOX_R (sanity check only)
BOX_PRIMES = (1073741789, 1073741783)


# ======================================================================================================
# The order O (sympy for the field arithmetic; integers for everything after that)
# ======================================================================================================
TT = symbols('t')


class Order:
    def __init__(self, key):
        self.key = key
        a0, a1, a2, a3, a4 = FORMS[key]
        self.F = FORMS[key]
        self.f = Poly(a4 * TT ** 4 + a3 * TT ** 3 + a2 * TT ** 2 + a1 * TT + a0, TT, domain=QQ)
        self.Bm = Matrix([[Rational(Fr(c).numerator, Fr(c).denominator) for c in row] for row in BASES[key]])
        self.Bi = self.Bm.inv()
        self.om = [Poly(sum(self.Bm[i, j] * TT ** j for j in range(4)), TT, domain=QQ) for i in range(4)]
        self.C = [[None] * 4 for _ in range(4)]
        self.integral = True
        for i in range(4):
            for j in range(4):
                r = (self.om[i] * self.om[j]).rem(self.f)
                co = Matrix([self._pc(r)]) * self.Bi
                if not all(x.is_integer for x in co):
                    self.integral = False
                self.C[i][j] = tuple(int(x) for x in co)
        self.e1, self.e2 = UNITS[key]
        self._inv = {}
        self._unit = {}

    @staticmethod
    def _pc(r):
        cs = [Rational(c) for c in r.all_coeffs()[::-1]]
        return cs + [Rational(0)] * (4 - len(cs))

    def poly(self, x):
        return sum((self.om[i] * int(x[i]) for i in range(4)), Poly(0, TT, domain=QQ))

    def mul(self, x, y):
        z = [0, 0, 0, 0]
        C = self.C
        for i in range(4):
            xi = x[i]
            if xi:
                for j in range(4):
                    yj = y[j]
                    if yj:
                        p = xi * yj
                        c = C[i][j]
                        z[0] += p * c[0]
                        z[1] += p * c[1]
                        z[2] += p * c[2]
                        z[3] += p * c[3]
        return tuple(z)

    def mulmod(self, x, y, m):
        return tuple(v % m for v in self.mul(x, y))

    def mat(self, x):
        """M(x)[k][j] = k-th coordinate of x * omega_j, so coords(x*y) = M(x) coords(y)."""
        return [[sum(x[i] * self.C[i][j][k] for i in range(4)) for j in range(4)] for k in range(4)]

    def norm_det(self, x):
        return int(Matrix(self.mat(x)).det())

    def norm_res(self, x):
        return self.f.resultant(self.poly(x))

    def inv(self, x):
        x = tuple(x)
        if x not in self._inv:
            y = Matrix(self.mat(x)).LUsolve(Matrix([1, 0, 0, 0]))
            assert all(c.is_integer for c in y), 'inverse not in O'
            self._inv[x] = tuple(int(c) for c in y)
        return self._inv[x]

    def pw(self, x, n):
        if n < 0:
            x, n = self.inv(x), -n
        r, b = (1, 0, 0, 0), tuple(x)
        while n:
            if n & 1:
                r = self.mul(r, b)
            n >>= 1
            if n:
                b = self.mul(b, b)
        return r

    def pw_mod(self, x, n, m):
        if n < 0:
            x, n = self.inv(x), -n
        r, b = (1 % m, 0, 0, 0), tuple(v % m for v in x)
        while n:
            if n & 1:
                r = self.mulmod(r, b, m)
            n >>= 1
            if n:
                b = self.mulmod(b, b, m)
        return r

    def unit(self, a, b):
        """eps1^a eps2^b, exact."""
        if (a, b) not in self._unit:
            self._unit[(a, b)] = self.mul(self.pw(self.e1, a), self.pw(self.e2, b))
        return self._unit[(a, b)]

    def unit_mod(self, a, b, m):
        return self.mulmod(self.pw_mod(self.e1, a, m), self.pw_mod(self.e2, b, m), m)


def matmod(M, m):
    return [[v % m for v in row] for row in M]


# ======================================================================================================
# Lattices of Z^2: Hermite rows (A, B), (0, C) with A, C > 0 and 0 <= B < C
# ======================================================================================================
def xgcd(a, b):
    r0, r1, s0, s1, t0, t1 = a, b, 1, 0, 0, 1
    while r1:
        k = r0 // r1
        r0, r1 = r1, r0 - k * r1
        s0, s1 = s1, s0 - k * s1
        t0, t1 = t1, t0 - k * t1
    if r0 < 0:
        r0, s0, t0 = -r0, -s0, -t0
    return r0, s0, t0


def _hnf_add(H, v):
    A, B, C = H
    a, b = v
    g, s, t = xgcd(A, a)                    # s A + t a = g ;  [[s, t], [a/g, -A/g]] is unimodular
    d = (a // g) * B - (A // g) * b
    C2 = gcd(C, d)
    return (g, (s * B + t * b) % C2, C2)


def hnf(vectors):
    vs = [(int(v[0]), int(v[1])) for v in vectors if (v[0], v[1]) != (0, 0)]
    D = 0
    for i in range(len(vs)):
        for j in range(i + 1, len(vs)):
            D = abs(vs[i][0] * vs[j][1] - vs[i][1] * vs[j][0])
            if D:
                break
        if D:
            break
    assert D > 0, 'not of full rank'
    H = (D, 0, D)                           # D Z^2 lies in the lattice
    for v in vs:
        H = _hnf_add(H, v)
    return H


class Lat:
    def __init__(self, H):
        self.A, self.B, self.C = H
        assert self.A > 0 and self.C > 0 and 0 <= self.B < self.C

    def red(self, v):
        k = v[0] // self.A
        return (v[0] - k * self.A, (v[1] - k * self.B) % self.C)

    def has(self, v):
        return self.red(v) == (0, 0)

    def index(self):
        return self.A * self.C

    def basis(self):
        return [(self.A, self.B), (0, self.C)]

    def __repr__(self):
        return f'<({self.A},{self.B}),(0,{self.C})>'


def dual_rows(rows):
    (p, q), (r, s) = rows
    det = Fr(p) * s - Fr(q) * r
    return [(Fr(s) / det, Fr(-r) / det), (Fr(-q) / det, Fr(p) / det)]


def intersect(L, M):
    rows = dual_rows(L.basis()) + dual_rows(M.basis())
    den = 1
    for x, y in rows:
        den = den * x.denominator // gcd(den, x.denominator)
        den = den * y.denominator // gcd(den, y.denominator)
    S = Lat(hnf([(int(x * den), int(y * den)) for x, y in rows]))         # den (L* + M*)
    inter = dual_rows([(Fr(x, den), Fr(y, den)) for x, y in S.basis()])
    assert all(x.denominator == 1 and y.denominator == 1 for x, y in inter)
    I = Lat(hnf([(int(x), int(y)) for x, y in inter]))
    ok = all(L.has(v) and M.has(v) for v in I.basis())
    Lsum = Lat(hnf(L.basis() + M.basis()))
    ok = ok and I.index() * Lsum.index() == L.index() * M.index()
    Hs = hermite_normal_form(Matrix([[L.A, 0, M.A, 0], [L.B, L.C, M.B, M.C]]))     # columns = generators
    ok = ok and abs(Hs.det()) == Lsum.index()
    check(ok, f'lattice intersection {L} cap {M}')
    return I


def coset_reps(Lold, Lnew):
    """Vectors of Lold, one in each class of Lold / Lnew."""
    S = []
    for v in Lnew.basis():
        c1, r1 = divmod(v[0], Lold.A)
        c2, r2 = divmod(v[1] - c1 * Lold.B, Lold.C)
        assert r1 == 0 and r2 == 0, 'Lnew not inside Lold'
        S.append((c1, c2))
    HS = Lat(hnf(S))
    (u1a, u1b), (u2a, u2b) = Lold.basis()
    reps = [(i * u1a + j * u2a, i * u1b + j * u2b) for i in range(HS.A) for j in range(HS.C)]
    check(len(reps) * Lold.index() == Lnew.index() and len({Lnew.red(x) for x in reps}) == len(reps),
          'coset representatives')
    return reps


# ======================================================================================================
# G_q, Lambda_q, candidates
# ======================================================================================================
def group_bfs(O, q, cap=None):
    M1, M2 = matmod(O.mat(O.e1), q), matmod(O.mat(O.e2), q)
    one = (1, 0, 0, 0)
    ev = {one: (0, 0)}
    queue = [one]
    rels = []
    i = 0
    while i < len(queue):
        g = queue[i]
        i += 1
        a, b = ev[g]
        for M, (da, db) in ((M1, (1, 0)), (M2, (0, 1))):
            h = tuple((M[k][0] * g[0] + M[k][1] * g[1] + M[k][2] * g[2] + M[k][3] * g[3]) % q for k in range(4))
            nv = (a + da, b + db)
            o = ev.get(h)
            if o is None:
                ev[h] = nv
                queue.append(h)
                if cap is not None and len(ev) > cap:
                    return None, None
            elif nv != o:
                rels.append((nv[0] - o[0], nv[1] - o[1]))
    return ev, rels


class ModData:
    """Lambda_q (certified) and the group, for one field and one modulus q = p^e."""

    def __init__(self, O, q):
        self.q = q
        self.p = next(d for d in range(2, q + 1) if q % d == 0)
        e, t = 0, q
        while t % self.p == 0:
            t //= self.p
            e += 1
        assert t == 1, 'not a prime power'
        self.e = e
        self.ev, rels = group_bfs(O, q)
        self.order = len(self.ev)
        self.L = Lat(hnf(rels))
        one = (1, 0, 0, 0)
        ok = self.L.index() == self.order
        for v in self.L.basis():
            ok = ok and O.unit_mod(v[0], v[1], q) == one
        check(ok, f'Lambda_q for q = {q} (index {self.L.index()} vs |G_q| {self.order})')
        # bijection classes <-> group elements
        cls = {self.L.red(e_): g for g, e_ in self.ev.items()}
        check(len(cls) == self.order, f'classes vs group elements, q = {q}')
        self.cls = cls

    def candidates(self, O, w):
        Mw = matmod(O.mat(w), self.q)
        q = self.q
        out = set()
        for c, g in self.cls.items():
            r2 = (Mw[2][0] * g[0] + Mw[2][1] * g[1] + Mw[2][2] * g[2] + Mw[2][3] * g[3]) % q
            r3 = (Mw[3][0] * g[0] + Mw[3][1] * g[1] + Mw[3][2] * g[2] + Mw[3][3] * g[3]) % q
            if r2 == 0 and r3 == 0:
                out.add(c)
        return out


# ======================================================================================================
# Proposition 17, exact
# ======================================================================================================
def nu(x, p):
    if x == 0:
        return INF
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def prop17(O, p, e, z, basis):
    q = p ** e
    cols = []
    for lam in basis:
        eta = O.unit(*lam)
        check((eta[0] - 1) % q == 0 and all(c % q == 0 for c in eta[1:]), f'eta^{lam} not == 1 mod {q}')
        ez = O.mul(eta, z)
        d = [ez[k] - z[k] for k in range(4)]
        check(all(c % q == 0 for c in d), 'K z not integral')
        cols.append((d[2] // q, d[3] // q))          # (l_1 K_j z, l_2 K_j z), l = markers x2, x3
    (x1, x2), (y1, y2) = cols
    lam = min(nu(v, p) for v in (x1, x2, y1, y2))
    J = x1 * y2 - x2 * y1
    nJ = nu(J, p)
    iota = 1 if p == 2 else 0
    delta = INF if lam >= INF or nJ >= INF else nJ - lam
    ok = lam <= e - 1 - iota and delta <= e - 1 - iota and (p != 2 or e >= 2)
    return dict(lam=lam, nJ=nJ, delta=delta, ok=ok)


def Lam(n, p, e):
    """Lambda(n) = ceil((n theta + 1)/(p - 1)), theta = e(p - 1) - 1 (proof.md 1.5)."""
    th = e * (p - 1) - 1
    return -((-(n * th + 1)) // (p - 1))


def vfac(n, p):
    v, t = 0, p
    while t <= n:
        v += n // t
        t *= p
    return v


def p17_mechanism(O, p, e, z, basis, D=6):
    """Instance check of steps 1-3 of the proof of Proposition 17 at the seed z.

    K_j = (T_j - I)/q is multiplication by kappa_j = (eta_j - 1)/q in O, so b^(k)_ij = l_k K1^i K2^j z is the
    k-th marker of kappa1^i kappa2^j z.  Everything mod p^M; valuations below the precision are exact.
    Returns (ok, text): f1 has E within {(1,0),(0,1)} containing (1,0) after relabelling, f2' = f2 - rho f1 has
    E = {(0,1)} with s* = e + delta, and val >= Lambda(i+j) on every computed (i,j) with i+j >= 1."""
    q = p ** e
    M = e * (D + 2) + 40
    PM = p ** M
    kap = []
    for lam in basis:
        eta = O.unit(*lam)
        d = (eta[0] - 1, eta[1], eta[2], eta[3])
        kap.append(tuple((c // q) % PM for c in d))
    zz = tuple(c % PM for c in z)
    b = {}
    pw1 = zz
    for i in range(D + 1):
        cur = pw1
        for j in range(D + 1 - i):
            b[(i, j)] = (cur[2], cur[3])
            cur = O.mulmod(cur, kap[1], PM)
        pw1 = O.mulmod(pw1, kap[0], PM)
    nuM = lambda x: min(nu(x % PM, p), M)
    x, y = b[(1, 0)], b[(0, 1)]
    lam = min(nuM(v) for v in x + y)
    if lam >= M:
        return False, 'lambda beyond the precision'
    # relabel: indices first (so that the minimum is attained by an x), then markers (so that nu(x1) = lam)
    swap_idx = min(nuM(x[0]), nuM(x[1])) > lam
    if swap_idx:
        b = {(j, i): v for (i, j), v in b.items()}
        x, y = y, x
    swap_mk = nuM(x[0]) > lam
    if swap_mk:
        b = {k: (v[1], v[0]) for k, v in b.items()}
        x, y = (x[1], x[0]), (y[1], y[0])
    Pl = p ** lam
    PM2 = p ** (M - lam)
    rho = ((x[1] // Pl) * pow(x[0] // Pl, -1, PM2)) % PM2
    ok = True
    val1, val2 = {}, {}
    for (i, j), (u1, u2) in b.items():
        bb2 = (u2 - rho * u1) % PM2
        for val, bv, prec in ((val1, u1 % PM, M), (val2, bb2, M - lam)):
            vv = nu(bv, p)
            val[(i, j)] = INF if vv >= prec else (i + j) * e + vv - vfac(i, p) - vfac(j, p)
    for val in (val1, val2):
        for (i, j), v in val.items():
            if i + j >= 1 and v < INF:
                ok = ok and v >= Lam(i + j, p, e)
    ok = ok and val1[(0, 0)] == INF and val2[(0, 0)] == INF and val2[(1, 0)] == INF
    m1 = min(val1.values())
    E1 = sorted(k for k, v in val1.items() if v == m1)
    m2 = min(val2.values())
    E2 = sorted(k for k, v in val2.items() if v == m2)
    nJ = nu((x[0] * y[1] - x[1] * y[0]) % PM2, p)
    delta = nJ - lam
    ok = ok and m1 == e + lam and m1 < Lam(2, p, e) and (1, 0) in E1 and set(E1) <= {(1, 0), (0, 1)}
    ok = ok and E2 == [(0, 1)] and m2 == e + delta and m2 < Lam(2, p, e)
    txt = (f'relabel (indices {swap_idx}, markers {swap_mk}); f1: s* = {m1}, E = {E1}; f2\': s* = {m2} = e + delta, '
           f'E = {E2}; Lambda(2) = {Lam(2, p, e)}; val >= Lambda(i+j) on i+j <= {D}')
    return ok, txt


def three_bases(L, rng):
    A, B, C = L.A, L.B, L.C
    own = [(A, B), (0, C)]
    g, s, t = xgcd(B, C)
    h1, b0 = A * C // g, g
    a1 = (-(s * A)) % h1
    script = [(h1, 0), (-a1, b0)]
    check(L.has(script[0]) and L.has(script[1]) and h1 * b0 == L.index(), 'script-style basis')
    while True:
        U = [[rng.randint(-2, 2) for _ in range(2)] for _ in range(2)]
        if abs(U[0][0] * U[1][1] - U[0][1] * U[1][0]) == 1 and U not in ([[1, 0], [0, 1]], [[0, 1], [1, 0]]):
            break
    uni = [(U[i][0] * A, U[i][0] * B + U[i][1] * C) for i in range(2)]
    return {'own': own, 'script': script, 'unimod': uni}


def fmt_nu(v):
    return 'inf' if v >= INF else str(v)


# ======================================================================================================
# The sieve (printed order of moduli) and the direct check of all classes
# ======================================================================================================
def run_sieve(seq, data, zeros, p17, use_sib=True):
    L = Lat((1, 0, 1))
    live = [(0, 0)]
    used, steps, sib_uses = [], [], []
    certified = {}
    for q in seq:
        md, cands = data[q]
        Lnew = intersect(L, md.L)
        reps = coset_reps(L, Lnew)
        newlive = []
        for r in live:
            for x in reps:
                v = (r[0] + x[0], r[1] + x[1])
                if md.L.red(v) in cands:
                    newlive.append(Lnew.red(v))
        check(len(set(newlive)) == len(newlive), 'duplicate classes')
        mods = used + [q]
        ncert = nsib = 0
        rest = []
        for v in newlive:
            hits = [zi for zi, zz in enumerate(zeros) if Lnew.has((v[0] - zz['exp'][0], v[1] - zz['exp'][1]))]
            if hits:
                check(len(hits) == 1, f'two known zeros in one class {v}')
                zi = hits[0]
                okq = [qq for qq in mods if p17[(qq, zi)]]
                if okq:
                    ncert += 1
                    check(zi not in certified, 'zero certified twice')
                    certified[zi] = (okq[0], v, Lnew, okq)
                else:
                    rest.append(v)
                continue
            sib = [(qq, zi) for qq in mods for zi, zz in enumerate(zeros)
                   if p17[(qq, zi)] and data[qq][0].L.has((v[0] - zz['exp'][0], v[1] - zz['exp'][1]))]
            if sib and use_sib:
                nsib += 1
                sib_uses.append((q, v, [(qq, zeros[zi]['mn']) for qq, zi in sib]))
                continue
            rest.append(v)
        steps.append((q, len(newlive), ncert, nsib, len(rest)))
        used.append(q)
        L = Lnew
        live = rest
    return dict(steps=steps, live=live, certified=certified, L=L, sib=sib_uses, used=used)


def direct_check(L, used, data, zeros, p17):
    """Every class of Z^2 / L: excluded (non-candidate for some q) or inside the Lambda_q-class of a zero at which
    Proposition 17 holds for that q.  Returns (unaccounted, excluded, covered)."""
    A, Cc = L.A, L.C
    masks, zcls = {}, []
    for q in used:
        md, cands = data[q]
        mk = np.zeros(md.L.A * md.L.C, dtype=bool)
        for (x, y) in cands:
            mk[x * md.L.C + y] = True
        masks[q] = mk
        for zi, zz in enumerate(zeros):
            if p17[(q, zi)]:
                rr = md.L.red(zz['exp'])
                zcls.append((q, rr[0] * md.L.C + rr[1]))
    alive = killed_cand = covered = 0
    chunk = max(1, (1 << 20) // Cc)
    for a0 in range(0, A, chunk):
        a1 = min(A, a0 + chunk)
        aa = np.repeat(np.arange(a0, a1, dtype=np.int64), Cc)
        bb = np.tile(np.arange(Cc, dtype=np.int64), a1 - a0)
        dead = np.zeros(aa.shape, dtype=bool)
        cov = np.zeros(aa.shape, dtype=bool)
        for q in used:
            Lq = data[q][0].L
            k = aa // Lq.A
            idx = (aa - k * Lq.A) * Lq.C + np.mod(bb - k * Lq.B, Lq.C)
            dead |= ~masks[q][idx]
            for (qq, zix) in zcls:
                if qq == q:
                    cov |= (idx == zix)
        alive += int(np.count_nonzero(~(dead | cov)))
        killed_cand += int(np.count_nonzero(dead))
        covered += int(np.count_nonzero(cov & ~dead))
    return alive, killed_cand, covered


CTX = []


# ======================================================================================================
# MAIN
# ======================================================================================================
print('verify_exclusion.py -- independent re-implementation of the exclusion step (exact)')
print('(nothing imported from core.py, lifting2.py, 02_exclusion.py; data retyped; compared with 02_exclusion.out.txt)')
print()
rng = random.Random(20261015)
FIELDS = {}
MOD = {}
for key in ('10', '01', '11'):
    print('=' * 110)
    O = Order(key)
    FIELDS[key] = O
    print(f'CLASS {LABEL[key]}  F = {FORMS[key]}  f(t) = {O.f.as_expr()}')
    # ---- (0) data checks -------------------------------------------------------------------------
    check(O.integral, 'structure constants not integral: the basis does not span an order')
    check(list(O.Bm.row(0)) == [1, 0, 0, 0] and list(O.Bm.row(1)) == [0, 1, 0, 0], 'omega_0 = 1, omega_1 = alpha')
    blk = Matrix([[O.Bm[2, 2], O.Bm[2, 3]], [O.Bm[3, 2], O.Bm[3, 3]]])
    check(blk.det() != 0, 'alpha^2, alpha^3 block singular')
    check(O.f.is_irreducible, 'f reducible')
    nreal = len(real_roots(O.f))
    check(nreal == 2, f'signature: {nreal} real roots')
    disc = O.f.discriminant()
    print(f'   order: structure constants integral = {O.integral}; f irreducible = {O.f.is_irreducible}; '
          f'real roots = {nreal}; disc(f) = {disc}; [O : Z[alpha]] = {abs(1 / O.Bm.det())}; '
          f'disc(O) = {disc * O.Bm.det() ** 2}')
    # markers: x in Z + Z alpha  <=>  x2 = x3 = 0   (block of omega_2, omega_3 on alpha^2, alpha^3 is invertible)
    for nm, u in (('eps1', O.e1), ('eps2', O.e2)):
        nd, nr = O.norm_det(u), O.norm_res(u)
        ui = O.inv(u)
        check(nd == 1 and nr == 1 and O.mul(u, ui) == (1, 0, 0, 0), f'{nm} not a norm-1 unit of O')
        print(f'   {nm} = {u}: N = {nd} (det) = {nr} (resultant); inverse in O = {ui}')
    # multiplicative independence (floating point; information only, the coverage step certifies it)
    rts = sorted(r.real for r in np.roots([float(c) for c in O.f.all_coeffs()]) if abs(r.imag) < 1e-9)
    def logs(x):
        out = []
        for r in rts:
            v = sum(float(O.Bm[i, j]) * x[i] * r ** j for i in range(4) for j in range(4))
            out.append(np.log(abs(v)))
        return out
    l1, l2 = logs(O.e1), logs(O.e2)
    print(f'   log|sigma_1|, log|sigma_2| of eps1 = ({l1[0]:.6f}, {l1[1]:.6f}), eps2 = ({l2[0]:.6f}, {l2[1]:.6f}); '
          f'det = {l1[0] * l2[1] - l1[1] * l2[0]:.6f} (float, information only)')
    # ---- |G_q| for every printed available modulus, and the skip decisions --------------------------
    agree = True
    got = {}
    for q in MODULI_ALL:
        ev, _ = group_bfs(O, q, cap=MAXG)
        if ev is None:
            agree &= check(q not in AVAILABLE[key], f'q = {q}: |G_q| > {MAXG} but printed as available')
        else:
            got[q] = len(ev)
            agree &= check(AVAILABLE[key].get(q) == len(ev),
                           f'q = {q}: |G_q| = {len(ev)} vs printed {AVAILABLE[key].get(q)}')
    print(f'   |G_q| <= {MAXG} recomputed (BFS): ' + ', '.join(f'{q}({n})' for q, n in got.items()))
    print(f'   agrees with the printed "moduli available" line (values and skipped moduli): {agree}')
print()

# seeds: 02 vs 01
s02 = {}
for orb in ORBITS:
    s02.setdefault(orb['key'], {}).setdefault(orb['c'], []).append(orb['w'])
check(s02 == SEEDS_01, 'seeds of 02_exclusion.out.txt differ from 01_coverage.out.txt')
print(f'seeds of 02_exclusion.out.txt == seeds of 01_coverage.out.txt: {s02 == SEEDS_01}')
print()

SUMMARY = []
for orb in ORBITS:
    key, c, w = orb['key'], orb['c'], orb['w']
    O = FIELDS[key]
    print('=' * 110)
    print(f'CLASS {LABEL[key]}  c = {c}  seed w = {w}')
    nfail0 = len(FAILS)
    # ---- (0) seed and known zeros --------------------------------------------------------------------
    check(O.norm_det(w) == c and O.norm_res(w) == c, f'N(w) != {c}')
    knownc = [(m, n) for (m, n) in KNOWN[key]
              if sum(a * m ** (4 - i) * n ** i for i, a in enumerate(FORMS[key])) == c]
    check(sorted(knownc) == sorted(orb['known']), f'known solutions with F = c: {knownc} vs printed {orb["known"]}')
    zeros = []
    for (m, n), ex in zip(orb['known'], orb['exps']):
        z = O.mul(O.unit(*ex), w)
        target = (n, -m, 0, 0)
        sgn = 1 if z == target else (-1 if z == tuple(-t for t in target) else 0)
        check(sgn != 0, f'eps^{ex} w = {z} is not +-({n} - {m} alpha)')
        check(O.norm_det(target) == c, f'N({n} - {m} alpha) != c')
        zeros.append(dict(mn=(m, n), exp=ex, z=z))
        print(f'   known solution (m,n) = {(m, n)}: eps1^{ex[0]} eps2^{ex[1]} w = {"+" if sgn > 0 else "-"}({n} - {m} alpha)  [exact]')
    # ---- (1) G_q, Lambda_q, candidates ---------------------------------------------------------------
    data = {}
    for q in orb['used']:
        if (key, q) not in MOD:
            MOD[(key, q)] = ModData(O, q)
        md = MOD[(key, q)]
        cands = md.candidates(O, w)
        data[q] = (md, cands)
        pg, pc = orb['G'][q]
        check(md.order == pg and len(cands) == pc,
              f'q = {q}: |G_q| = {md.order}, candidates = {len(cands)} vs printed {pg}, {pc}')
        print(f'   q = {q:3d} = {md.p}^{md.e}: |G_q| = {md.order:5d} (printed {pg}); Lambda_q = {md.L} '
              f'[index certified]; candidate classes = {len(cands)} (printed {pc})')
    # ---- (2) Proposition 17 at each known zero, each modulus, three bases --------------------------------
    p17 = {}
    for q in orb['used']:
        md = data[q][0]
        for zi, zz in enumerate(zeros):
            bs = three_bases(md.L, rng)
            res = {nm: prop17(O, md.p, md.e, zz['z'], b) for nm, b in bs.items()}
            vals = {(r['lam'], r['nJ'], r['ok']) for r in res.values()}
            oks = {r['ok'] for r in res.values()}
            check(len(oks) == 1, f'Prop. 17 verdict depends on the basis at q = {q}, zero {zz["mn"]}')
            if any(oks):
                check(len(vals) == 1, f'lambda / nu(J) depend on the basis at q = {q}, zero {zz["mn"]}')
            r = res['script']
            p17[(q, zi)] = r['ok']
            print(f'   Prop. 17 at zero {zz["mn"]}, q = {q:3d} (p = {md.p}, e = {md.e}, iota = {int(md.p == 2)}): '
                  + '; '.join(f'{nm}: lam={fmt_nu(res[nm]["lam"])} nuJ={fmt_nu(res[nm]["nJ"])}' for nm in res)
                  + f'  -> bound e-1-iota = {md.e - 1 - (md.p == 2)}: {"HOLDS" if r["ok"] else "fails"}'
                  + f'   [bases: own {bs["own"]}, script {bs["script"]}, unimod {bs["unimod"]}]')
            if r['ok']:
                for nm in ('script', 'unimod'):
                    okm, txt = p17_mechanism(O, md.p, md.e, zz['z'], bs[nm])
                    check(okm, f'mechanism of the proof of Prop. 17 at q = {q}, zero {zz["mn"]}, basis {nm}: {txt}')
                    print(f'      proof mechanism ({nm} basis): {txt}: {okm}')
    # ---- (3) the sieve in the printed order ---------------------------------------------------------------
    sv = run_sieve(orb['used'], data, zeros, p17)
    steps, live, certified, L, sib_uses, used = sv['steps'], sv['live'], sv['certified'], sv['L'], sv['sib'], sv['used']
    CTX.append((orb, data, zeros, p17))
    check(steps == orb['steps'], f'sieve counts {steps} vs printed {orb["steps"]}')
    check(L.index() == orb['index'], f'final lattice index {L.index()} vs printed {orb["index"]}')
    for (q, nl, nc, ns, nr) in steps:
        print(f'   sieve after q = {q:3d}: {nl} candidate classes, {nc} certified, {ns} sibling-empty, {nr} left')
    for (q, v, lst) in sib_uses:
        print(f'   sibling-empty class {v} (step q = {q}): inside the Lambda_qq-class of a certified zero for '
              f'(qq, zero) in {lst}')
    print(f'   final lattice {L}, index {L.index()} (printed {orb["index"]}); live classes left: {len(live)}')
    check(not live, 'live classes remain')
    check(sorted(certified) == list(range(len(zeros))), 'not every known zero certified exactly once')
    for (rep, mn, pq, plam, pnJ, pe) in orb['cert']:
        zi = next(i for i, zz in enumerate(zeros) if zz['mn'] == mn)
        q0, v, Lc, okq = certified[zi]
        md = data[pq][0]
        b = three_bases(md.L, rng)['script']
        r = prop17(O, md.p, md.e, zeros[zi]['z'], b)
        check(q0 == pq and Lc.has((rep[0] - zeros[zi]['exp'][0], rep[1] - zeros[zi]['exp'][1]))
              and (r['lam'], r['nJ'], md.e) == (plam, pnJ, pe),
              f'certified class of {mn}: printed rep {rep}, q {pq}, lam {plam}, nuJ {pnJ}, e {pe}')
        print(f'   printed "class {rep} contains {mn}, Prop. 17 at q = {pq} (lambda = {plam}, nu(J) = {pnJ}, e = {pe})": '
              f'rep == exponent {zeros[zi]["exp"]} mod the lattice of that step: True; recomputed lambda = '
              f'{fmt_nu(r["lam"])}, nu(J) = {fmt_nu(r["nJ"])}; Prop. 17 holds there at q in {okq}')
    # ---- (3') direct check of all classes modulo the final lattice ------------------------------------------
    A, Cc = L.A, L.C
    alive, killed_cand, covered = direct_check(L, used, data, zeros, p17)
    for zz in zeros:                                   # a genuine zero is never a non-candidate
        for q in used:
            check(data[q][0].L.red(zz['exp']) in data[q][1], f'known zero {zz["mn"]} is not a candidate mod {q}')
    check(alive == 0, f'direct check: {alive} classes mod the final lattice are neither excluded nor covered')
    print(f'   direct check over all {A * Cc} classes mod the final lattice: {killed_cand} excluded by a candidate test, '
          f'{covered} inside the Lambda_q-class of a zero certified by Prop. 17, {alive} unaccounted')
    # ---- (4) brute force box (sanity only) -----------------------------------------------------------------
    R = BOX_R
    hits = None
    for P in BOX_PRIMES:
        E1 = np.zeros((2 * R + 1, 4), dtype=np.int64)
        cur = O.pw_mod(O.e1, -R, P)
        for i in range(2 * R + 1):
            E1[i] = cur
            cur = O.mulmod(cur, O.e1, P)
        y = O.mulmod(O.pw_mod(O.e2, -R, P), w, P)
        S = set()
        for j in range(2 * R + 1):
            My = np.array(matmod(O.mat(y), P), dtype=np.int64)
            Z = (E1 @ My.T) % P
            for i in np.nonzero((Z[:, 2] == 0) & (Z[:, 3] == 0))[0]:
                S.add((int(i) - R, j - R))
            y = O.mulmod(y, O.e2, P)
        hits = S if hits is None else hits & S
    exact = sorted(v for v in hits if O.mul(O.unit(*v), w)[2:] == (0, 0))
    check(exact == sorted(orb['exps']), f'box scan zeros {exact} vs known exponents {orb["exps"]}')
    print(f'   box |a|,|b| <= {R}: zeros of both markers mod {BOX_PRIMES[0]} and mod {BOX_PRIMES[1]}: {sorted(hits)}; '
          f'exact zeros {exact} (known exponents {sorted(orb["exps"])})')
    ok = len(FAILS) == nfail0
    SUMMARY.append((key, c, w, len(zeros), ok))
    print(f'   ORBIT {"CONFIRMED" if ok else "NOT CONFIRMED"}: zeros = exponents of the known solutions')

print()
print('=' * 110)
print('NEGATIVE CONTROLS (do the checks above have teeth?)')
nc_ok = True
for orb, data, zeros, p17 in CTX:
    if not zeros:
        continue
    tag = f'class {LABEL[orb["key"]]} c = {orb["c"]} seed {orb["w"]}'
    # (N1) without Proposition 17 the classes of the known zeros cannot be removed
    nop = {k: False for k in p17}
    sv = run_sieve(orb['used'], data, zeros, nop)
    al, _, _ = direct_check(sv['L'], sv['used'], data, zeros, nop)
    n1 = len(sv['live']) >= len(zeros) and al >= len(zeros)
    # (N2) without the sibling rule
    sv2 = run_sieve(orb['used'], data, zeros, p17, use_sib=False)
    # (N3) the mechanism check rejects every (q, zero) at which the hypotheses of Proposition 17 fail
    n3 = True
    nrej = 0
    for q in orb['used']:
        md = data[q][0]
        for zi, zz in enumerate(zeros):
            if not p17[(q, zi)]:
                okm, _ = p17_mechanism(O=FIELDS[orb['key']], p=md.p, e=md.e, z=zz['z'],
                                       basis=three_bases(md.L, rng)['script'])
                n3 = n3 and not okm
                nrej += 1
    nc_ok = nc_ok and n1 and n3
    print(f'   {tag}: without Prop. 17: {len(sv["live"])} live classes after the sieve, {al} unaccounted in the direct '
          f'check (>= {len(zeros)} expected): {n1};  without the sibling rule: {len(sv2["live"])} live classes left;  '
          f'mechanism check rejects all {nrej} (q, zero) pairs where the hypotheses fail: {n3}')
check(nc_ok, 'negative controls')

print()
print('SUMMARY')
for key, c, w, nz, ok in SUMMARY:
    print(f'   class {LABEL[key]} c = {c:2d} seed {w}: {nz} known solution pair(s); independently confirmed: {ok}')
check(all(isprime(P) for P in BOX_PRIMES), 'box primes')
if FAILS:
    print(f'CHECK FAILED: {len(FAILS)} failed check(s)')
else:
    print('ALL 16 ORBITS INDEPENDENTLY CONFIRMED; no discrepancy with 02_exclusion.out.txt')
print(f'[done in {time.time() - T0:.1f} s]')
sys.exit(1 if FAILS else 0)
