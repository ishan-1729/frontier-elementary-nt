"""verify_H_units.py -- independent verification of hypothesis_H.md, Sections 1 and 4 (card triangular_tetrahedral).

Independent verifier (Claude subagent, Opus 5.5), 2026-10-15.  Nothing is imported from cubic.py, core.py, lifting2.py
or any 0x_*.py script, and no code is copied from them: the claims are read from hypothesis_H.md and checked here with
arithmetic written here, by methods chosen to differ from 04_ring.py and 06_units.py where one exists.
Exact arithmetic only: int, fractions.Fraction, sympy.  numpy is used with int64 only, under an explicit overflow
bound.  Real numbers enter only through a rational enclosure of l1, the real root of f(t) = t^3 - 4t + 6, obtained
here from the continued fraction expansion of l1 (Lagrange's method; 04/06 use bisection), and through a centred form
with an exact remainder bound for the values x0 + x1 l1 + x2 l1^2 (04/06 use naive interval arithmetic).

Notation: a triple (a, b, c) is a + b l + c l^2 with l^3 = 4l - 6.  sigma1(u) = a + b l1 + c l1^2; Q(u) is the product
of the two complex embeddings; T2 = sigma1^2 + 2Q; EPS = (-19, -5, 1); E = Q(EPS) = 1/sigma1(EPS).

PART A  Section 1.  Multiplication by reduction mod f; the norm form derived three ways (sympy resultant, own
        Sylvester determinant, Newton's identities on traces) and compared with (1.2); M(u); the adjugate derived by
        Cayley-Hamilton and compared with (1.3); multiplicativity; (1.4) by an exhaustive residue check mod 8 and by
        the coefficient parities quoted in (1.4); the special elements of (1.5) (and the identities of the display of
        Theorem H and of (5.3)); (1.6): uniqueness of the real root (own elementary argument), the enclosure of l1,
        Q derived from z^2 = s z - p, Q = sigma1(u*), multiplicativity, the sum of squares, T2, its minors and
        completed squares, the inverse of the embedding.
PART B  Section 4.  (4.2) values of EPS; (4.3)/(4.3') descent identities and constants rho_lo, rho_hi, A, B, h(A),
        h(B) and the exact threshold K3^(1/3), compared with the values printed in 06_units.out.txt; (4.4) the units
        with T2 <= 114, in the box of (4.4) and in the bounding box of the ellipsoid; THE UNIT THEOREM (4.1) BY A
        DIFFERENT ROUTE: a fundamental domain 1 <= |sigma1| < 1/sigma1(EPS) for multiplication by EPS.  Every integer
        triple with |sigma1| <= 606 and Q <= 1 is enumerated twice (nested ranges from the inverse of the embedding;
        a box from the ellipsoid T2 <= 606^2 + 2, scanned with exact int64 bounds), and the units there are listed;
        Remark (4.6).
PART C  Brute force: every triple with |a|, |b|, |c| <= 60 and N = +-1.
PART D  Negative controls: EPS^2 in place of EPS in the fundamental-domain argument (both orientations of the
        domain), and planted errors that the checks of Parts A and B must catch.

Every check prints PASS or FAIL; the script stops at the first FAIL.  Usage (from the scripts folder):
    python verify_H_units.py > verify_H_units.out.txt 2>&1
"""
import math
import os
import re
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

try:
    sys.stdout.reconfigure(newline='\n')
except (AttributeError, ValueError):
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
T_START = time.time()
NPASS = 0


def check(name, cond):
    """Print PASS or FAIL; stop at the first FAIL."""
    global NPASS
    if cond:
        NPASS += 1
        print('PASS  ' + name, flush=True)
        return
    print('FAIL  ' + name, flush=True)
    print(f'STOPPED at the first FAIL; {NPASS} checks had passed before it.', flush=True)
    sys.exit(1)


def note(s=''):
    print('      ' + s, flush=True)


def part(title):
    print('\n' + '=' * 116 + '\n' + title + '\n' + '=' * 116, flush=True)


def sub(title):
    print('\n--- ' + title, flush=True)


# ======================================================================================================================
part('PART A. Section 1: the ring of integer triples (own arithmetic)')
# ======================================================================================================================
t, L, X, r_ = sp.symbols('t L X r')
a, b, c, d, e, g, h, k, m = sp.symbols('a b c d e g h k m')
a1, b1 = sp.symbols('a1 b1')
FT = t**3 - 4*t + 6
U, V, W = (a, b, c), (d, e, g), (h, k, m)


def reduce_list(cs):
    """Remainder modulo t^3 - 4t + 6 of sum cs[i] t^i, as [r0, r1, r2]: t^j (j >= 3) becomes t^(j-3)(4t - 6),
    from the top degree down.  Works for int, Fraction and sympy coefficients."""
    cs = list(cs)
    while len(cs) < 3:
        cs.append(0)
    for j in range(len(cs) - 1, 2, -1):
        top = cs[j]
        cs[j] = 0
        cs[j - 2] = cs[j - 2] + 4 * top
        cs[j - 3] = cs[j - 3] - 6 * top
    return cs[:3]


def tmul(u, v):
    """Product of two triples: the product of u(t) and v(t), reduced modulo f(t)."""
    prod = [0] * 5
    for i in range(3):
        for j in range(3):
            prod[i + j] = prod[i + j] + u[i] * v[j]
    return tuple(reduce_list(prod))


def smul(u, v):
    return tuple(sp.expand(x) for x in tmul(u, v))


def tadd(u, v):
    return tuple(x + y for x, y in zip(u, v))


def tsub(u, v):
    return tuple(x - y for x, y in zip(u, v))


def tscal(s0, u):
    return tuple(s0 * x for x in u)


def tneg(u):
    return tuple(-x for x in u)


def tpow(u, n):
    out = (1, 0, 0)
    for _ in range(n):
        out = tmul(out, u)
    return out


def same(u, v):
    return all(sp.expand(x - y) == 0 for x, y in zip(u, v))


def as_t(u, var=t):
    return u[0] + u[1] * var + u[2] * var**2


def redL(expr):
    """Remainder of a polynomial in L (coefficients in Q[a, b, ...]) modulo L^3 - 4L + 6, computed by reduce_list."""
    cs = sp.Poly(sp.expand(expr), L).all_coeffs()[::-1]
    r0, r1, r2 = reduce_list(cs)
    return sp.expand(r0 + r1 * L + r2 * L**2)


def sig(u, var=L):
    """sigma1 as a polynomial in L (L stands for l1)."""
    return u[0] + u[1] * var + u[2] * var**2


def poly_fn(expr, gens):
    """An exact integer evaluator for a polynomial with integer coefficients (built from its terms)."""
    terms = []
    for mon, co in sp.Poly(sp.expand(expr), *gens).terms():
        if not co.is_Integer:
            raise ValueError('non-integer coefficient')
        fac = [f'({int(co)})'] + [(f'x{i}**{ex}' if ex > 1 else f'x{i}') for i, ex in enumerate(mon) if ex]
        terms.append('*'.join(fac))
    return eval('lambda x0, x1, x2: ' + (' + '.join(terms) if terms else '0'))


# ---------------------------------------------------------------------------------------------------------------------
sub('A.1 multiplication (1.1)')


def doc_mul(u, v):
    """(1.1) as printed: (a,b,c)(d,e,f) = (ad - 6bf - 6ce, ae + bd + 4bf + 4ce - 6cf, af + be + cd + 4cf)."""
    p0, p1, p2 = u
    q0, q1, q2 = v
    return (p0*q0 - 6*p1*q2 - 6*p2*q1, p0*q1 + p1*q0 + 4*p1*q2 + 4*p2*q1 - 6*p2*q2, p0*q2 + p1*q1 + p2*q0 + 4*p2*q2)


rem_c = sp.Poly(sp.rem(sp.expand(as_t(U) * as_t(V)), FT, t), t).all_coeffs()[::-1]
rem_c = tuple((list(rem_c) + [0, 0, 0])[:3])
check('own reduction (t^j -> t^(j-3)(4t - 6)) of u(t)v(t) agrees with sympy.rem modulo t^3 - 4t + 6',
      same(smul(U, V), rem_c))
check('the product formula printed in (1.1) equals u(t)v(t) reduced modulo f(t)', same(smul(U, V), doc_mul(U, V)))
LAM = (0, 1, 0)
check('l^3 = 4l - 6 and l^4 = 4l^2 - 6l', tpow(LAM, 3) == (-6, 4, 0) and tpow(LAM, 4) == (0, -6, 4))
check('ring laws as polynomial identities: commutative, associative, distributive, unit (1,0,0)',
      same(smul(U, V), smul(V, U)) and same(smul(smul(U, V), W), smul(U, smul(V, W)))
      and same(smul(U, tadd(V, W)), tadd(smul(U, V), smul(U, W))) and same(smul((1, 0, 0), U), U))

# ---------------------------------------------------------------------------------------------------------------------
sub('A.2 the norm form (1.2), derived three ways')
N_res = sp.expand(sp.resultant(FT, as_t(U), t))
SYL = sp.Matrix([[1, 0, -4, 6, 0],
                 [0, 1, 0, -4, 6],
                 [c, b, a, 0, 0],
                 [0, c, b, a, 0],
                 [0, 0, c, b, a]])
N_syl = sp.expand(SYL.det())
e1_, e2_, e3_ = 0, -4, -6                      # f = t^3 - e1 t^2 + e2 t - e3
p1_ = e1_                                      # Newton: power sums of the roots of f
p2_ = e1_ * p1_ - 2 * e2_
check('f = t^3 - e1 t^2 + e2 t - e3 with e1 = 0, e2 = -4, e3 = -6; power sums p1 = 0, p2 = 8 (Newton)',
      sp.expand(t**3 - e1_*t**2 + e2_*t - e3_ - FT) == 0 and (p1_, p2_) == (0, 8))


def trace(x):
    """Tr(x0 + x1 l + x2 l^2) = sum over the three roots r of f of x0 + x1 r + x2 r^2 = 3 x0 + p1 x1 + p2 x2."""
    return 3 * x[0] + p1_ * x[1] + p2_ * x[2]


U2 = smul(U, U)
U3 = smul(U2, U)
P1, P2, P3 = sp.expand(trace(U)), sp.expand(trace(U2)), sp.expand(trace(U3))
N_newton = sp.expand((P1**3 - 3 * P1 * P2 + 2 * P3) / 6)
N_doc = a**3 + 8*a**2*c - 4*a*b**2 + 18*a*b*c + 16*a*c**2 - 6*b**3 + 24*b*c**2 + 36*c**3
note(f'Res_t(t^3 - 4t + 6, a + bt + ct^2) = {N_res}')
check('the resultant (= product of a + b r + c r^2 over the roots r of the monic f) equals the determinant of the '
      'Sylvester matrix built here', sp.expand(N_res - N_syl) == 0)
check('Newton: e3 of the conjugates, (Tr(u)^3 - 3 Tr(u) Tr(u^2) + 2 Tr(u^3))/6, equals the resultant',
      sp.expand(N_res - N_newton) == 0)
check('the cubic form printed in (1.2): N = a^3 + 8a^2c - 4ab^2 + 18abc + 16ac^2 - 6b^3 + 24bc^2 + 36c^3',
      sp.expand(N_res - N_doc) == 0)
cols = [smul(U, (1, 0, 0)), smul(U, (0, 1, 0)), smul(U, (0, 0, 1))]
MU = sp.Matrix(3, 3, lambda i, j: cols[j][i])
M_doc = sp.Matrix([[a, -6*c, -6*b], [b, a + 4*c, 4*b - 6*c], [c, b, a + 4*c]])
colsV = [smul(V, (1, 0, 0)), smul(V, (0, 1, 0)), smul(V, (0, 0, 1))]
MV = sp.Matrix(3, 3, lambda i, j: colsV[j][i])
UV = smul(U, V)
colsUV = [smul(UV, (1, 0, 0)), smul(UV, (0, 1, 0)), smul(UV, (0, 0, 1))]
MUV = sp.Matrix(3, 3, lambda i, j: colsUV[j][i])
check('the matrix with columns u*1, u*l, u*l^2 equals M(u) printed in (1.2); det M(u) = N(u); M(uv) = M(u) M(v)',
      sp.expand(MU - M_doc) == sp.zeros(3, 3) and sp.expand(MU.det() - N_res) == 0
      and sp.expand(MUV - MU * MV) == sp.zeros(3, 3))
NORM = poly_fn(N_res, (a, b, c))


def norm(u):
    return NORM(u[0], u[1], u[2])


def Nsym(u):
    return sp.expand(N_res.xreplace({a: u[0], b: u[1], c: u[2]}))


check('N(uv) = N(u) N(v) (polynomial identity in six variables)', sp.expand(Nsym(UV) - N_res * Nsym(V)) == 0)
check('N(r u) = r^3 N(u) and N(-u) = -N(u)',
      sp.expand(Nsym(tscal(r_, U)) - r_**3 * N_res) == 0 and sp.expand(Nsym(tneg(U)) + N_res) == 0)
check('N(X - theta) = N(X, -6, 0) = X^3 - 144X + 1296 identically, and theta = (0,6,0) satisfies '
      'theta^3 - 144 theta + 1296 = 0',
      sp.expand(Nsym((X, -6, 0)) - (X**3 - 144*X + 1296)) == 0
      and tadd(tadd(tpow((0, 6, 0), 3), tscal(-144, (0, 6, 0))), (1296, 0, 0)) == (0, 0, 0))

# ---------------------------------------------------------------------------------------------------------------------
sub('A.3 the adjugate triple (1.3), derived by Cayley-Hamilton')
s2u = sp.expand((P1**2 - P2) / 2)              # e2 of the conjugates of u
adj_ch = tuple(sp.expand(U2[i] - P1 * U[i] + (s2u if i == 0 else 0)) for i in range(3))
note(f'Tr(u) = {P1},  s2(u) = {s2u}')


def doc_adj(u):
    """(1.3) as printed: u* = (a^2 + 8ac - 4b^2 + 6bc + 16c^2, -ab - 6c^2, b^2 - ac - 4c^2)."""
    x0, x1, x2 = u
    return (x0*x0 + 8*x0*x2 - 4*x1*x1 + 6*x1*x2 + 16*x2*x2, -x0*x1 - 6*x2*x2, x1*x1 - x0*x2 - 4*x2*x2)


check('Cayley-Hamilton: u (u^2 - Tr(u) u + s2(u)) = (N(u), 0, 0) as a polynomial identity',
      same(smul(U, adj_ch), (N_res, 0, 0)))
check('the adjugate printed in (1.3) equals u^2 - Tr(u) u + s2(u), so u u* = (N(u), 0, 0)',
      same(adj_ch, doc_adj(U)) and same(smul(U, doc_adj(U)), (N_res, 0, 0)))
check('(uv)* = u* v* and N(u*) = N(u)^2 (polynomial identities)',
      same(doc_adj(UV), smul(doc_adj(U), doc_adj(V))) and sp.expand(Nsym(doc_adj(U)) - N_res**2) == 0)
check('u* has integer polynomial coordinates, and N(1) = 1',
      all(all(co.is_Integer for co in sp.Poly(x, a, b, c).coeffs()) for x in adj_ch) and norm((1, 0, 0)) == 1)
ADJ = [poly_fn(x, (a, b, c)) for x in adj_ch]


def adj_int(u):
    return (ADJ[0](u[0], u[1], u[2]), ADJ[1](u[0], u[1], u[2]), ADJ[2](u[0], u[1], u[2]))


note('(1.3a): if N(u) = +-1 then u^-1 = N(u) u* is an integer triple; if uv = 1 then N(u)N(v) = N(1) = 1, N(u) = +-1.')
note('(1.3b): a rational triple with N(u) != 0 has the inverse u*/N(u), by u u* = N(u).')
uq = (Fr(1, 2), Fr(-3), Fr(2, 7))
check('example of (1.3b): u = (1/2, -3, 2/7) times u*/N(u) is (1, 0, 0)',
      tmul(uq, tscal(1 / norm(uq), adj_int(uq))) == (1, 0, 0))

# ---------------------------------------------------------------------------------------------------------------------
sub('A.4 the norm vanishes only at 0 (1.4)')
zeros8 = [(x, y, z) for x in range(8) for y in range(8) for z in range(8) if norm((x, y, z)) % 8 == 0]
bad8 = [v for v in zeros8 if any(x % 2 for x in v)]
check(f'exhaustive check mod 8: N(a,b,c) = 0 (mod 8) only if a, b, c are all even ({len(zeros8)} of 512 residue '
      f'triples have N = 0 mod 8, and all of them are even)', not bad8)
note('N has integer coefficients, so N mod 8 depends only on (a,b,c) mod 8: a primitive integer triple has N != 0;')
note('and N(D u) = D^3 N(u) turns a rational zero u != 0 into a primitive integer zero. So N(u) = 0 only for u = 0.')


def all_div(expr, mod, gens):
    return all(co.is_Integer and int(co) % mod == 0 for co in sp.Poly(sp.expand(expr), *gens).coeffs())


N2 = sp.expand(N_res.xreplace({a: 2 * a1}))
PP = sp.expand(N2 / 2)
N4 = sp.expand(N_res.xreplace({a: 2 * a1, b: 2 * b1}))
PP4 = sp.expand(N4 / 4)
check('(1.4) as written: N - a^3 has even coefficients', all_div(N_res - a**3, 2, (a, b, c)))
check('(1.4) as written: N(2a1, b, c) = 2P with P integral and P + 3b^3 with even coefficients',
      all_div(N2, 2, (a1, b, c)) and all_div(PP + 3 * b**3, 2, (a1, b, c)))
check("(1.4) as written: N(2a1, 2b1, c) = 4P' with P' integral and P' - 9c^3 with even coefficients",
      all_div(N4, 4, (a1, b1, c)) and all_div(PP4 - 9 * c**3, 2, (a1, b1, c)))
check('f = t^3 - 4t + 6 is Eisenstein at 2 (-4, 6 even; 6 not divisible by 4)', (-4) % 2 == 0 == 6 % 2 and 6 % 4 != 0)

# ---------------------------------------------------------------------------------------------------------------------
sub('A.5 special elements (1.5)')
EPS = (-19, -5, 1)
EPS_INV = (95, -101, 40)
PI2 = (-4, 1, 1)
P3A = (-1, 1, 0)
P3B = (5, 2, 0)
P3C = (3, -3, 1)
P179 = (107, -114, 45)
THETA = (0, 6, 0)
check('the triples of (1.5) match their polynomials: l^2-5l-19, 40l^2-101l+95, l^2+l-4, l-1, 2l+5, l^2-3l+3, '
      '45l^2-114l+107',
      all(sp.expand(sig(x, L) - pol) == 0 for x, pol in zip(
          (EPS, EPS_INV, PI2, P3A, P3B, P3C, P179),
          (L**2 - 5*L - 19, 40*L**2 - 101*L + 95, L**2 + L - 4, L - 1, 2*L + 5, L**2 - 3*L + 3,
           45*L**2 - 114*L + 107))))
norms_found = [norm(x) for x in (EPS, EPS_INV, PI2, P3A, P3B, P3C, P179)]
check(f'norms of EPS, EPS^-1, PI2, P3A, P3B, P3C, P179: {norms_found} = [1, 1, -2, -3, -3, 3, 179]',
      norms_found == [1, 1, -2, -3, -3, 3, 179])
check('EPS EPS^-1 = 1, and EPS^-1 = EPS* (N(EPS) = 1)', tmul(EPS, EPS_INV) == (1, 0, 0) and adj_int(EPS) == EPS_INV)
check(f'PI2^3 = {tpow(PI2, 3)} = -2 EPS, so 2 = -PI2^3 EPS^-1',
      tpow(PI2, 3) == tscal(-2, EPS) and tmul(tneg(tpow(PI2, 3)), EPS_INV) == (2, 0, 0))
check('P3A P3B P3C = 3', tmul(tmul(P3A, P3B), P3C) == (3, 0, 0))
check('PI2^2 P179 = 3l^2 - 4 = (-4, 0, 3), and N(3l^2 - 4) = 716', tmul(tpow(PI2, 2), P179) == (-4, 0, 3)
      and norm((-4, 0, 3)) == 716)
check('(extra, (5.3)) 36(3l^2 - 4) = EPS^-2 PI2^8 (P3A P3B P3C)^2 P179',
      tmul(tmul(tmul(tpow(EPS_INV, 2), tpow(PI2, 8)), tpow(tmul(tmul(P3A, P3B), P3C), 2)), P179) == (-144, 0, 108))
R1, R2, R3 = (-45, 48, -19), (5, 2, 0), (-8, 8, -3)
check('(extra, display of Theorem H) 9 - theta = EPS P3A P3B r1^2, -15 - theta = -P3A P3C r2^2, '
      '-12 - theta = -EPS P3B P3C r3^2',
      tmul(tmul(tmul(EPS, P3A), P3B), tpow(R1, 2)) == (9, -6, 0)
      and tneg(tmul(tmul(P3A, P3C), tpow(R2, 2))) == (-15, -6, 0)
      and tneg(tmul(tmul(tmul(EPS, P3B), P3C), tpow(R3, 2))) == (-12, -6, 0))

# ---------------------------------------------------------------------------------------------------------------------
sub('A.6 real values (1.6)')


def fval(x):
    return x**3 - 4 * x + 6


s_ = sp.symbols('s')
check('f(-3) = -9 < 0 < 6 = f(-2)', fval(-3) == -9 and fval(-2) == 6)
check('identities for the uniqueness of the real root (own argument): f(t) = t(t - 2)(t + 2) + 6; '
      't(2 - t) = 1 - (t - 1)^2; f(t) - f(s) = (t - s)(t^2 + ts + s^2 - 4)',
      sp.expand(FT - (t * (t - 2) * (t + 2) + 6)) == 0 and sp.expand(t * (2 - t) - (1 - (t - 1)**2)) == 0
      and sp.expand(FT - FT.subs(t, s_) - (t - s_) * (t**2 + t * s_ + s_**2 - 4)) == 0)
note('t >= -2: t(t^2 - 4) >= 0 on [-2, 0] and on [2, oo), and = -(t + 2) t(2 - t) >= -4 on [0, 2]; so f(t) >= 2 > 0.')
note('s < t <= -2: t^2, ts, s^2 >= 4, so f(t) - f(s) >= 8(t - s) > 0: f increases on (-oo, -2] and has exactly one')
note('real root l1 there, in (-3, -2).  (hypothesis_H.md uses the critical values f(+-2/sqrt3) instead; checked next.)')
check("the doc's argument: f(+-2/sqrt3) = 6 -+ 16/(3 sqrt3), and (16/(3 sqrt3))^2 = 256/27 < 36",
      sp.simplify(FT.subs(t, 2 / sp.sqrt(3)) - (6 - 16 / (3 * sp.sqrt(3)))) == 0
      and sp.simplify(FT.subs(t, -2 / sp.sqrt(3)) - (6 + 16 / (3 * sp.sqrt(3)))) == 0 and Fr(256, 27) < 36)
check('f has no rational root (f(r) != 0 for r = +-1, +-2, +-3, +-6): f is irreducible over Q, l1 is irrational and '
      '1, l1, l1^2 are linearly independent over Q', all(fval(x) != 0 for x in (1, -1, 2, -2, 3, -3, 6, -6)))


def peval(cs, x):
    out = 0
    for co in reversed(cs):
        out = out * x + co
    return out


def shift(cs, s0):
    """Coefficients (lowest first) of P(x + s0)."""
    n = len(cs) - 1
    return [sum(cs[j] * math.comb(j, i) * s0**(j - i) for j in range(i, n + 1)) for i in range(n + 1)]


def sgn(x):
    return (x > 0) - (x < 0)


def floor_root(cs):
    """Floor of the unique real root (which is > 1) of the polynomial cs, by sign changes at integers."""
    s1 = sgn(peval(cs, 1))
    if s1 == 0:
        check('continued fraction: no rational root met', False)
    hi_ = 2
    while sgn(peval(cs, hi_)) == s1:
        hi_ *= 2
    lo_ = hi_ // 2
    while hi_ - lo_ > 1:
        mid_ = (lo_ + hi_) // 2
        if sgn(peval(cs, mid_)) == s1:
            lo_ = mid_
        else:
            hi_ = mid_
    return lo_


# Lagrange's method: x = q + 1/y maps the unique real root to the unique real root y > 1 of y^3 P(q + 1/y)
pq = [-3]
cur = shift([6, -4, 0, 1], -3)[::-1]
h_prev, h_cur, k_prev, k_cur = 1, -3, 0, 1
while k_prev * k_cur <= 2**300:
    qn = floor_root(cur)
    pq.append(qn)
    h_prev, h_cur = h_cur, qn * h_cur + h_prev
    k_prev, k_cur = k_cur, qn * k_cur + k_prev
    cur = shift(cur, qn)[::-1]
lo, hi = sorted([Fr(h_prev, k_prev), Fr(h_cur, k_cur)])
check(f'continued fraction of l1 = [{pq[0]}; {", ".join(map(str, pq[1:16]))}, ...] ({len(pq)} partial quotients): '
      f'two consecutive convergents give lo < hi with f(lo) < 0 < f(hi), hi - lo < 2^-300 and -3 < lo < hi < -5/2',
      fval(lo) < 0 < fval(hi) and hi - lo < Fr(1, 2**300) and -3 < lo and hi < Fr(-5, 2))
note(f'so l1 is in (lo, hi) (one real root), l1 = {float(lo):.16f}; 04_ring.out.txt prints -2.525102254814321')
check('consistent with the value printed by 04 (lo and hi round to -2.525102254814321 at 16 digits)',
      abs(lo - Fr('-2.525102254814321')) < Fr(1, 10**15))
KB = 256
SC = 1 << KB
SC2 = SC * SC
MID = math.floor((lo + hi) * SC / 2)
RAD = math.ceil(max(hi * SC - MID, MID - lo * SC)) + 1
MID_SC, MID2, RAD_SC = MID * SC, MID * MID, RAD * SC
check(f'dyadic centred enclosure: [lo, hi] lies in [m - r, m + r] with m = MID/2^256, r = {RAD}/2^256, and '
      f'-3 < m - r < m + r < -5/2', MID - RAD <= lo * SC and hi * SC <= MID + RAD and -3 * SC < MID - RAD
      and 2 * (MID + RAD) < -5 * SC)
note('Centred form: for x = (x0, x1, x2) and t, m in [m - r, m + r] (inside (-3, -2)),')
note('  (x0 + x1 t + x2 t^2) - (x0 + x1 m + x2 m^2) = (t - m)(x1 + x2 (t + m)), of absolute value <= r (|x1| + 6|x2|).')


def enc(x0, x1, x2):
    """sigma1((x0, x1, x2)) * 2^512 lies in [P - ER, P + ER] (exact integers)."""
    return x0 * SC2 + x1 * MID_SC + x2 * MID2, RAD_SC * (abs(x1) + 6 * abs(x2))


def iv(x3):
    P_, ER_ = enc(*x3)
    return Fr(P_ - ER_, SC2), Fr(P_ + ER_, SC2)


def sgn_minus(x3, T=0):
    """Exact sign of sigma1(x3) - T for an integer triple x3 and a rational T.  By the linear independence of
    1, l1, l1^2 the difference is 0 only when x3 = (T, 0, 0); any other undecided case is a FAIL."""
    P_, ER_ = enc(*x3)
    if isinstance(T, int):
        R_ = T * SC2
        if P_ - ER_ > R_:
            return 1
        if P_ + ER_ < R_:
            return -1
    else:
        T = Fr(T)
        lhs, rhs, er = P_ * T.denominator, T.numerator * SC2, ER_ * T.denominator
        if lhs - er > rhs:
            return 1
        if lhs + er < rhs:
            return -1
    if x3[1] == 0 and x3[2] == 0 and x3[0] == T:
        return 0
    check(f'sign of sigma1({x3}) - {T} decided', False)


check('f(t) - (t - L)(t^2 + L t + L^2 - 4) = f(L): z, zbar are the roots of t^2 + l1 t + l1^2 - 4, so s = z + zbar = '
      '-l1, p = z zbar = l1^2 - 4; the discriminant 16 - 3 l1^2 < 0 since l1 < -5/2',
      sp.expand(FT - (t - L) * (t**2 + L * t + L**2 - 4) - (L**3 - 4 * L + 6)) == 0 and 16 - 3 * hi**2 < 0)
sS, pS, zS = sp.symbols('s p z')
check('z^2 = s z - p gives a + bz + cz^2 = (a - cp) + (b + cs) z (the difference is c (z^2 - s z + p))',
      sp.expand(a + b*zS + c*zS**2 - ((a - c*pS) + (b + c*sS)*zS) - c*(zS**2 - sS*zS + pS)) == 0)
Xq, Yq = a - c * pS, b + c * sS
Q_gen = sp.expand(Xq**2 + Xq * Yq * sS + Yq**2 * pS)
check('|X + Y z|^2 = X^2 + XY s + Y^2 p with X = a - cp, Y = b + cs gives the formula of (1.6): '
      'a^2 + s ab + (s^2 - 2p) ac + p b^2 + sp bc + p^2 c^2',
      sp.expand(Q_gen - (a**2 + sS*a*b + (sS**2 - 2*pS)*a*c + pS*b**2 + sS*pS*b*c + pS**2*c**2)) == 0)


def Qd(u, var=L):
    """Q as printed in (1.6): a^2 - L ab + (8 - L^2) ac + (L^2 - 4) b^2 + 6bc + (16 - 6L - 4L^2) c^2."""
    x0, x1, x2 = u
    return (x0**2 - var*x0*x1 + (8 - var**2)*x0*x2 + (var**2 - 4)*x1**2 + 6*x1*x2
            + (16 - 6*var - 4*var**2)*x2**2)


check('with s = -L, p = L^2 - 4 this reduces mod f(L) to the Q printed in (1.6)',
      redL(Q_gen.xreplace({sS: -L, pS: L**2 - 4}) - Qd(U)) == 0)
check('Q(u) = sigma1(u*) (identity in a, b, c, L)', sp.expand(Qd(U) - sig(doc_adj(U))) == 0)
check('u(t)v(t) - (uv)(t) is divisible by f(t), so sigma1(uv) = sigma1(u) sigma1(v)',
      sp.rem(sp.expand(as_t(U) * as_t(V) - as_t(smul(U, V))), FT, t) == 0)
check('mod f(L): sigma1(u) Q(u) = N(u)', redL(sig(U) * Qd(U) - N_res) == 0)
check('mod f(L): Q(uv) = Q(u) Q(v)', redL(Qd(smul(U, V)) - Qd(U) * Qd(V)) == 0)
Qsq = (a - L / 2 * b + (8 - L**2) / 2 * c)**2 + (3 * L**2 - 16) / 4 * (b - L * c)**2
check('mod f(L): Q = (a - (L/2) b + ((8 - L^2)/2) c)^2 + ((3L^2 - 16)/4)(b - L c)^2, and 3 l1^2 - 16 > 11/4 > 0',
      redL(Qsq - Qd(U)) == 0 and 3 * hi**2 - 16 > Fr(11, 4))
check('Q vanishes on the real line through (l1^2 - 4, l1, 1) (mod f(L))', redL(Qd((L**2 - 4, L, 1))) == 0)
note('So Q >= 0, and Q(u) = 0 for a rational u forces b = l1 c, hence c = 0 (l1 irrational), b = 0, a = 0: Q > 0 on')
note('nonzero rational triples (also N = sigma1 Q != 0 by (1.4)).')
T2d = 3*a**2 + 16*a*c + (3*L**2 - 8)*b**2 + 8*L*b*c + (32 - 18*L - 4*L**2)*c**2
check('mod f(L): T2 = sigma1^2 + 2Q equals the form printed in (1.6)', redL(sig(U)**2 + 2 * Qd(U) - T2d) == 0)
G = sp.Matrix([[3, 0, 8], [0, 3*L**2 - 8, 4*L], [8, 4*L, 32 - 18*L - 4*L**2]])
check('T2 has the Gram matrix [[3,0,8],[0,3L^2-8,4L],[8,4L,32-18L-4L^2]]; leading minors 3, 9L^2 - 24, and '
      'det = 716 mod f(L); 9 l1^2 - 24 > 129/4 > 0, so T2 is positive definite',
      sp.expand((sp.Matrix([[a, b, c]]) * G * sp.Matrix([a, b, c]))[0] - T2d) == 0
      and sp.expand(G[:2, :2].det() - (9*L**2 - 24)) == 0 and redL(G.det()) == 716 and 9 * hi**2 - 24 > Fr(129, 4))
d2L = 3 * L**2 - 8
check('completed squares of (1.6): (3L^2 - 8) T2 = 3(3L^2 - 8)(a + 8c/3)^2 + ((3L^2 - 8) b + 4L c)^2 + (716/3) c^2 '
      'mod f(L)', redL(d2L * T2d - (3 * d2L * (a + sp.Rational(8, 3) * c)**2 + (d2L * b + 4 * L * c)**2
                                     + sp.Rational(716, 3) * c**2)) == 0)
cof11 = G[1, 1] * G[2, 2] - G[1, 2]**2
cof22 = G[0, 0] * G[2, 2] - G[0, 2]**2
cof33 = G[0, 0] * G[1, 1] - G[0, 1]**2
check('cofactors of the Gram matrix mod f(L): 64L^2 + 68, 32 - 54L - 12L^2, 9L^2 - 24 (so (G^-1)_ii = cof_ii/716)',
      redL(cof11 - (64*L**2 + 68)) == 0 and redL(cof22 - (32 - 54*L - 12*L**2)) == 0
      and sp.expand(cof33 - (9*L**2 - 24)) == 0)
Y1 = a - L / 2 * b + (8 - L**2) / 2 * c
Y2 = b - L * c
check('inverse of the embedding (identities in a, b, c, L): (3L^2 - 4) c = sigma1 - Y1 - (3L/2) Y2, b = Y2 + L c, '
      'a = Y1 + (L/2) b - ((8 - L^2)/2) c, where Y1 = Re sigma2 and Y2 = Im sigma2 / w, w^2 = (3L^2 - 16)/4, '
      'Q = Y1^2 + w^2 Y2^2',
      sp.expand((3*L**2 - 4) * c - (sig(U) - Y1 - 3 * L / 2 * Y2)) == 0 and sp.expand(b - (Y2 + L * c)) == 0
      and sp.expand(a - (Y1 + L / 2 * b - (8 - L**2) / 2 * c)) == 0
      and redL(Y1**2 + (3 * L**2 - 16) / 4 * Y2**2 - Qd(U)) == 0)
note('values (l1 enclosed as above):   element   N    sigma1                    Q')
ok_vals = True
for nm, x3 in (('EPS', EPS), ('EPS^-1', EPS_INV), ('PI2', PI2), ('P3A', P3A), ('P3B', P3B), ('P3C', P3C),
               ('P179', P179)):
    s_lo, s_hi = iv(x3)
    q_lo, q_hi = iv(adj_int(x3))
    nn = norm(x3)
    prod_lo = min(s_lo * q_lo, s_lo * q_hi, s_hi * q_lo, s_hi * q_hi)
    prod_hi = max(s_lo * q_lo, s_lo * q_hi, s_hi * q_lo, s_hi * q_hi)
    ok_vals = ok_vals and (s_lo > 0 or s_hi < 0) and q_lo > 0 and prod_lo <= nn <= prod_hi
    note(f'                                 {nm:8s} {nn:4d}  {float(s_lo):<24.15g}  {float(q_lo):.15g}')
check('for each special element sigma1 has a definite sign, Q > 0, and sigma1 Q encloses N', ok_vals)

# ======================================================================================================================
part('PART B. Section 4: the units')
# ======================================================================================================================
sub('B.1 values of EPS (4.2)')
se_lo, se_hi = iv(EPS)
E_lo, E_hi = iv(adj_int(EPS))                 # E = Q(EPS) = sigma1(EPS*) = sigma1(EPS^-1)
check('sigma1(EPS) = l1^2 - 5 l1 - 19 lies in (1/606, 1/605)', Fr(1, 606) < se_lo and se_hi < Fr(1, 605))
check('E := Q(EPS) = sigma1(EPS*) = sigma1(EPS^-1) = 40 l1^2 - 101 l1 + 95 lies in (605.08, 605.09)',
      Fr(60508, 100) < E_lo and E_hi < Fr(60509, 100))
check('sigma1(EPS) E = 1: exact (EPS EPS^-1 = 1 and sigma1 is multiplicative), and the interval product contains 1',
      tmul(EPS, EPS_INV) == (1, 0, 0) and se_lo * E_lo <= 1 <= se_hi * E_hi)
note(f'sigma1(EPS) = {float(se_lo):.15g},  E = {float(E_lo):.15g},  enclosure widths < 2^-240')
note('So sigma1(EPS^-1) = E, Q(EPS) = E, Q(EPS^-1) = sigma1(EPS) = 1/E.')

# ---------------------------------------------------------------------------------------------------------------------
sub('B.2 descent identities and constants (4.3), (4.3\')')
xs, qs, Es, ns, ms, ws = sp.symbols('x q E n m w', positive=True)
T2y = xs + 2 * qs
T2yE = xs / Es**2 + 2 * qs * Es               # sigma1(y EPS)^2 + 2 Q(y EPS)
T2yEi = xs * Es**2 + 2 * qs / Es              # sigma1(y EPS^-1)^2 + 2 Q(y EPS^-1)
check('T2(y EPS) - T2(y) = (E - 1)(2q - x(E + 1)/E^2) and T2(y EPS^-1) - T2(y) = (E - 1)(x(E + 1) - 2q/E)',
      sp.cancel(T2yE - T2y - (Es - 1) * (2 * qs - xs * (Es + 1) / Es**2)) == 0
      and sp.cancel(T2yEi - T2y - (Es - 1) * (xs * (Es + 1) - 2 * qs / Es)) == 0)
note('with E > 1, x, q > 0: T2(y EPS) < T2(y) iff x/q > rho_hi = 2E^2/(E+1);  T2(y EPS^-1) < T2(y) iff x/q < rho_lo')
note('= 2/(E(E+1)).')
check('2q - x(E+1)/E^2 = (q (E+1)/E^2)(rho_hi - x/q) and x(E+1) - 2q/E = q(E+1)(x/q - rho_lo) (identities)',
      sp.cancel(2*qs - xs*(Es + 1)/Es**2 - qs*(Es + 1)/Es**2*(2*Es**2/(Es + 1) - xs/qs)) == 0
      and sp.cancel(xs*(Es + 1) - 2*qs/Es - qs*(Es + 1)*(xs/qs - 2/(Es*(Es + 1)))) == 0)
E1s, E2s = sp.symbols('E1 E2', positive=True)
check('monotonicity: rho_hi(E2) - rho_hi(E1) = 2(E2 - E1)(E1 E2 + E1 + E2)/((E1+1)(E2+1)) (rho_hi INCREASES with E), '
      'rho_lo(E1) - rho_lo(E2) = 2(E2 - E1)(E1 + E2 + 1)/(E1 E2 (E1+1)(E2+1)) (rho_lo DECREASES with E)',
      sp.cancel(2*E2s**2/(E2s + 1) - 2*E1s**2/(E1s + 1) - 2*(E2s - E1s)*(E1s*E2s + E1s + E2s)/((E1s + 1)*(E2s + 1)))
      == 0 and sp.cancel(2/(E1s*(E1s + 1)) - 2/(E2s*(E2s + 1))
                         - 2*(E2s - E1s)*(E1s + E2s + 1)/(E1s*E2s*(E1s + 1)*(E2s + 1))) == 0)


def rho_hi_f(x):
    return 2 * x * x / (x + 1)


def rho_lo_f(x):
    return 2 / (x * (x + 1))


rh_lo, rh_hi = rho_hi_f(E_lo), rho_hi_f(E_hi)
rl_lo, rl_hi = rho_lo_f(E_hi), rho_lo_f(E_lo)
Ehb = Fr(60509, 100)
note(f'rho_hi in [{float(rh_lo):.15g}, {float(rh_hi):.15g}],  rho_lo in [{float(rl_lo):.15g}, {float(rl_hi):.15g}]')
check('rho_hi < 1209: rho_hi(E) <= rho_hi(60509/100) = 1208.1833... < 1209 since rho_hi increases and E < 60509/100',
      rh_hi < rho_hi_f(Ehb) < 1209)
check('rho_lo > 5453/10^9: rho_lo(E) >= rho_lo(60509/100) = 5.45347...e-6 > 5453/10^9 since rho_lo DECREASES and '
      'E < 60509/100', rl_lo > rho_lo_f(Ehb) > Fr(5453, 10**9))
note('DISCREPANCY (wording, conclusion unaffected): hypothesis_H.md (4.3) says "Both rho\'s increase with E"; that is')
note('true for rho_hi only. rho_lo = 2/(E(E+1)) decreases with E, and that is what makes E < 60509/100 give the lower')
note('bound rho_lo > 5453e-9. 06_units.py has the correct identity; only the sentence in the text is wrong.')
check('rho_hi/rho_lo = E^3', sp.cancel(2*Es**2/(Es + 1) / (2/(Es*(Es + 1))) - Es**3) == 0)
Aq, Bq = Fr(469, 5000), Fr(2841, 50)
check('A = 469/5000, B = 2841/50: A^3 * 1209 <= 1 <= B^3 * 5453/10^9, so q^3 in (1/1209, 10^9/5453) gives A < q < B',
      Aq**3 * 1209 <= 1 <= Bq**3 * Fr(5453, 10**9))


def hq(x):
    return 1 / x**2 + 2 * x


hA, hB = hq(Aq), hq(Bq)
note(f'h(A) = {hA} = {float(hA):.9f},  h(B) = {hB} = {float(hB):.9f}')
check('h(A) = 62603161709/549902500 < 113.85, h(B) < 113.65, max(h(A), h(B)) = h(A) < 114',
      hA == Fr(62603161709, 549902500) and hA < Fr(11385, 100) and hB < Fr(11365, 100) and max(hA, hB) == hA < 114)
check('for a unit, x q^2 = 1: x/q = q^-3 and T2 = x + 2q = h(q) = q^-2 + 2q; h\'\'(q) = 6/q^4 > 0 (convex)',
      sp.cancel((xs / qs).subs(xs, 1 / qs**2) - qs**-3) == 0 and sp.cancel(T2y.subs(xs, 1 / qs**2) - (qs**-2 + 2 * qs))
      == 0 and sp.cancel(sp.diff(qs**-2 + 2 * qs, qs, 2) - 6 / qs**4) == 0)
check("(4.3') scaled form: with n = m^3, q = m^2 w and x = n^2/q^2: x/q = w^-3 and T2 = m^2 h(w)",
      sp.cancel((ms**6 / (ms**2 * ws)**2) / (ms**2 * ws) - ws**-3) == 0
      and sp.cancel(ms**6 / (ms**2 * ws)**2 + 2 * ms**2 * ws - ms**2 * (ws**-2 + 2 * ws)) == 0)
Rs = sp.symbols('R', positive=True)
K3s = 4 * (Es**2 + Es + 1)**3 / (Es**2 * (Es + 1)**2)
check('exact threshold: with q^3 = R, h(q)^3 = (1 + 2R)^3/R^2; at both ends R = 1/rho_hi = (E+1)/(2E^2) and '
      'R = 1/rho_lo = E(E+1)/2 this is K3 = 4(E^2+E+1)^3/(E^2(E+1)^2), so max h on the admissible q-range is K3^(1/3)',
      sp.cancel(((1 + 2 * qs**3) / qs**2)**3 - ((1 + 2 * Rs)**3 / Rs**2).subs(Rs, qs**3)) == 0
      and sp.cancel(((1 + 2 * Rs)**3 / Rs**2).subs(Rs, (Es + 1) / (2 * Es**2)) - K3s) == 0
      and sp.cancel(((1 + 2 * Rs)**3 / Rs**2).subs(Rs, Es * (Es + 1) / 2) - K3s) == 0)


def K3f(x):
    return 4 * (x * x + x + 1)**3 / (x * x * (x + 1)**2)


K3_lo = 4 * (E_lo**2 + E_lo + 1)**3 / (E_hi**2 * (E_hi + 1)**2)   # numerator and denominator increase with E > 0
K3_hi = 4 * (E_hi**2 + E_hi + 1)**3 / (E_lo**2 * (E_lo + 1)**2)
check(f'K3 in [{float(K3_lo):.10f}, {float(K3_hi):.10f}]: 113.62^3 < K3 < 113.63^3, i.e. 113.62 < C* = K3^(1/3) < '
      f'113.63 (06: "113.62 < C* < 113.63")', Fr(11362, 100)**3 < K3_lo and K3_hi < Fr(11363, 100)**3)
check('K3 <= 113.6238^3, i.e. C0 = K3^(1/3) <= 113.6238 as stated in (2.5) (margin: K3 ~ 1466924.31 vs 1466925.06)',
      K3_hi <= Fr(1136238, 10000)**3)
check('C* < h(A) < 113.85 < 114 (the exact threshold is below the rounded one)', K3_hi < hA**3 and hA < 114)
note('Lemma (4.3) as checked: a unit u with T2(u EPS) >= T2(u) and T2(u EPS^-1) >= T2(u) has rho_lo <= q^-3 <= rho_hi,')
note('so A < q < B and T2(u) = h(q) <= max(h(A), h(B)) = h(A) < 113.85; sharper: T2(u) <= K3^(1/3) < 113.63.')

# comparison with the numbers printed in 06_units.out.txt (read as text; nothing is imported)
vals06 = {}
try:
    with open(os.path.join(HERE, '06_units.out.txt'), encoding='utf-8') as fh:
        txt06 = fh.read()
    pats = {'sigma1(EPS)': (r'sigma1\(EPS\) in \[([-0-9.eE+]+),', (se_lo, se_hi)),
            'E': (r'E = Q\(EPS\) in \[([-0-9.eE+]+),', (E_lo, E_hi)),
            'rho_lo': (r'rho_lo in \[([-0-9.eE+]+),', (rl_lo, rl_hi)),
            'rho_hi': (r'rho_hi in \[([-0-9.eE+]+),', (rh_lo, rh_hi))}
    for key, (pat, enc_) in pats.items():
        mm = re.search(pat, txt06)
        if mm:
            vals06[key] = (Fr(mm.group(1)), enc_)
    mmA = re.search(r'h\(A\) = (\d+)/(\d+)', txt06)
    mmB = re.search(r'h\(B\) = (\d+)/(\d+)', txt06)
except OSError:
    txt06, mmA, mmB = '', None, None
if len(vals06) == 4 and mmA and mmB:
    agree = all(abs(val - (lo_ + hi_) / 2) <= abs(val) * Fr(1, 10**10) for val, (lo_, hi_) in vals06.values())
    check('the values printed in 06_units.out.txt (12 significant digits) for sigma1(EPS), E, rho_lo, rho_hi agree with '
          'these enclosures to 1e-10 (relative), and its h(A), h(B) fractions are equal to these',
          agree and Fr(int(mmA.group(1)), int(mmA.group(2))) == hA and Fr(int(mmB.group(1)), int(mmB.group(2))) == hB)
else:
    note('06_units.out.txt not found or not parsed: comparison with its printed values skipped')

# ---------------------------------------------------------------------------------------------------------------------
sub('B.3 the units of small height (4.4)')


def t2_triple(v):
    """T2(v) = sigma1 of this integer triple: T2 = (3a^2+16ac-8b^2+32c^2) + (8bc-18c^2) l1 + (3b^2-4c^2) l1^2."""
    x0, x1, x2 = v
    return (3*x0*x0 + 16*x0*x2 - 8*x1*x1 + 32*x2*x2, 8*x1*x2 - 18*x2*x2, 3*x1*x1 - 4*x2*x2)


check('T2(a,b,c) = (3a^2 + 16ac - 8b^2 + 32c^2) + (8bc - 18c^2) L + (3b^2 - 4c^2) L^2, the form of (1.6)',
      sp.expand(sig(t2_triple(U)) - T2d) == 0)


def sqrt_up(y):
    """A rational >= sqrt(y), for a rational y >= 0."""
    n_ = math.ceil(Fr(y) * (1 << 128))
    s0 = math.isqrt(n_)
    if s0 * s0 < n_:
        s0 += 1
    return Fr(s0, 1 << 64)


def ellipsoid_box(M):
    """|a|, |b|, |c| bounds for T2 <= M: x_i^2 <= M (G^-1)_ii = M cof_ii/716 (Cauchy-Schwarz for the form T2)."""
    out = []
    for tri in ((68, 0, 64), (32, -54, -12), (-24, 0, 9)):
        out.append(math.floor(sqrt_up(Fr(M) * iv(tri)[1] / 716)))
    return out


def t2_le(v, M):
    return sgn_minus(t2_triple(v), M) <= 0


amax114, bmax114, cmax114 = ellipsoid_box(114)
box_mine = [(x, y, z) for x in range(-amax114, amax114 + 1) for y in range(-bmax114, bmax114 + 1)
            for z in range(-cmax114, cmax114 + 1)]
S114 = sorted(v for v in box_mine if t2_le(v, 114))
check(f'own box from the ellipsoid: T2 <= 114 forces |a| <= {amax114}, |b| <= {bmax114}, |c| <= {cmax114} '
      f'({len(box_mine)} triples, each decided exactly): {len(S114)} have T2 <= 114', (amax114, bmax114, cmax114)
      == (8, 3, 2) and len(S114) == 187)
d2_lo, d2_hi = iv((-8, 0, 3))
lam_abs = Fr(-(MID - RAD), SC)                  # >= |l1|
kabs = 4 * lam_abs / d2_lo                      # >= |k| = 4|l1|/d2
check('the bounds of (4.4) from the completed squares: 9 d3 > 114 (|c| <= 2); d2 (6 - 2|k|)^2 > 114 (|b| <= 5); '
      '3(12 - 16/3)^2 = 400/3 > 114 (|a| <= 11)',
      9 * Fr(716, 3) / d2_hi > 114 and 6 - 2 * kabs > 0 and d2_lo * (6 - 2 * kabs)**2 > 114
      and 3 * (12 - Fr(16, 3))**2 == Fr(400, 3) > 114)
box_doc = [(x, y, z) for x in range(-11, 12) for y in range(-5, 6) for z in range(-2, 3)]
S114_doc = sorted(v for v in box_doc if t2_le(v, 114))
check(f'(4.4): of the {len(box_doc)} triples of the box |a| <= 11, |b| <= 5, |c| <= 2 exactly {len(S114_doc)} have '
      f'T2 <= 114, and they are the same 187 triples as in the own box', len(box_doc) == 1265
      and len(S114_doc) == 187 and S114_doc == S114)
units114 = [v for v in S114 if norm(v) in (1, -1)]
check(f'the units with T2 <= 114 are exactly {units114} (T2(+-1) = 3)', units114 == [(-1, 0, 0), (1, 0, 0)])
note('(4.5) then gives the unit theorem: the least-T2 element v of {u EPS^k} has no smaller neighbour, so T2(v) <= C* <')
note('114 and v = +-1.  The next section proves (4.1) again without the descent.')

# ---------------------------------------------------------------------------------------------------------------------
sub('B.4 the unit theorem (4.1) by a fundamental domain for EPS (different route)')
S64 = 1 << 64
Llo, Lhi = math.floor(lo * S64), math.ceil(hi * S64)              # l1 * 2^64 in [Llo, Lhi]
L2lo, L2hi = math.floor(hi * hi * S64), math.ceil(lo * lo * S64)  # l1^2 * 2^64 in [L2lo, L2hi]
l2min = Fr(L2lo, S64)                                               # <= l1^2
D_lo = 3 * l2min - 4                                                # <= 3 l1^2 - 4
W_lo = 3 * l2min - 16                                               # <= 3 l1^2 - 16 = 4 w^2
Kc2_hi = 4 + 48 / W_lo                                              # >= 1 + 9 l1^2/(3 l1^2 - 16)
check('coarse bounds for the ranges: l1 in [Llo, Lhi]/2^64, l1^2 in [L2lo, L2hi]/2^64, 3 l1^2 - 16 >= W_lo > 0',
      Llo <= lo * S64 and hi * S64 <= Lhi and L2lo <= hi * hi * S64 and lo * lo * S64 <= L2hi and W_lo > 0)
note('Ranges (rigorous, from the inverse of the embedding, A.6): for |sigma1| <= X and Q <= Qm,')
note('  |c| <= (X + sqrt(Qm) sqrt(1 + 9 l1^2/(3 l1^2 - 16)))/(3 l1^2 - 4)   [Cauchy-Schwarz on Y1 + (3 l1/2) Y2],')
note('  |b - l1 c| = |Y2| <= sqrt(Qm)/w,  |a - m(b,c)| = |Y1| <= sqrt(Qm) with m = (l1/2) b - ((8 - l1^2)/2) c,')
note('  |a + b l1 + c l1^2| <= X.')


def cdiv(x, y):
    return -((-x) // y)


def region_candidates(X_, Qm):
    """A superset of the integer triples with |sigma1| <= X_ and Q <= Qm (X_, Qm positive integers)."""
    cmax = math.floor((X_ + sqrt_up(Qm * Kc2_hi)) / D_lo)
    RB = math.ceil(sqrt_up(Fr(4 * Qm) / W_lo) * S64)
    RA = math.ceil(sqrt_up(Qm) * S64)
    XS = X_ * S64
    out = []
    for cc in range(-cmax, cmax + 1):
        lc_lo, lc_hi = min(Llo * cc, Lhi * cc), max(Llo * cc, Lhi * cc)
        l2c_lo, l2c_hi = min(L2lo * cc, L2hi * cc), max(L2lo * cc, L2hi * cc)
        for bb in range(cdiv(lc_lo - RB, S64), (lc_hi + RB) // S64 + 1):
            s_lo_ = min(Llo * bb, Lhi * bb) + l2c_lo
            s_hi_ = max(Llo * bb, Lhi * bb) + l2c_hi
            a_lo = max(cdiv(s_lo_ - 8 * cc * S64 - 2 * RA, 2 * S64), cdiv(-XS - s_hi_, S64))
            a_hi = min((s_hi_ - 8 * cc * S64 + 2 * RA) // (2 * S64), (XS - s_lo_) // S64)
            for aa in range(a_lo, a_hi + 1):
                out.append((aa, bb, cc))
    return out, cmax


def in_region(v, X_, Qm):
    if sgn_minus(v, X_) > 0 or sgn_minus(v, -X_) < 0:
        return False
    return sgn_minus(adj_int(v), Qm) <= 0


def scan_region(X_, Qm):
    cands, cmax = region_candidates(X_, Qm)
    members, units = [], []
    for v in cands:
        if in_region(v, X_, Qm):
            members.append(v)
            if norm(v) in (1, -1):
                units.append(v)
    return cands, cmax, members, units


POW = {0: (1, 0, 0)}
for kk in range(1, 5):
    POW[kk] = tmul(POW[kk - 1], EPS)
    POW[-kk] = tmul(POW[-(kk - 1)], EPS_INV)


def name_unit(v):
    for kk in range(-4, 5):
        if v == POW[kk]:
            return f'EPS^{kk}'
        if v == tneg(POW[kk]):
            return f'-EPS^{kk}'
    return 'NOT +-EPS^k for |k| <= 4'


def abs_sigma_lt(v, w3):
    """|sigma1(v)| < sigma1(w3) (w3 with sigma1 > 0): +1 if true, 0 if equal (exact), -1 if greater."""
    sp_ = sgn_minus(tsub(w3, v), 0)
    sm_ = sgn_minus(tadd(w3, v), 0)
    if sp_ > 0 and sm_ > 0:
        return 1
    if sp_ == 0 or sm_ == 0:
        return 0
    return -1


t_b4 = time.time()
cands606, cmax606, mem606, units606 = scan_region(606, 1)
memb_ge1 = [v for v in mem606 if sgn_minus(v, 1) >= 0 or sgn_minus(v, -1) <= 0]
check(f'enumeration A (nested ranges): |c| <= {cmax606}, {len(cands606)} candidate triples, each decided exactly; '
      f'{len(mem606)} integer triples have |sigma1| <= 606 and Q <= 1, of which {len(memb_ge1)} have |sigma1| >= 1 '
      f'(the other one is 0: a nonzero triple with |sigma1| < 1 and Q <= 1 would have 0 < |N| < 1)',
      len(mem606) == len(memb_ge1) + 1 and (0, 0, 0) in mem606 and (0, 0, 0) not in memb_ge1)
note(f'(the expected number is about 2*606*pi/(sqrt(716)/2) = 284.6; time {time.time() - t_b4:.1f} s)')
# enumeration B: the bounding box of the ellipsoid T2 <= 606^2 + 2, scanned with exact int64 bounds
t_b4 = time.time()
MT = 606**2 + 2
amaxT, bmaxT, cmaxT = ellipsoid_box(MT)
S30 = 1 << 30
L1i = math.floor(lo * S30)
L2i = math.floor(hi * hi * S30)
check(f'enumeration B setup: T2 <= 606^2 + 2 gives |a| <= {amaxT}, |b| <= {bmaxT}, |c| <= {cmaxT}; '
      f'l1 2^30 in [L1, L1 + 1) and l1^2 2^30 in [L2, L2 + 1) with L1 = {L1i}, L2 = {L2i}',
      hi * S30 < L1i + 1 and L1i <= lo * S30 and lo * lo * S30 < L2i + 1 and L2i <= hi * hi * S30)
q0max = amaxT**2 + 8 * amaxT * cmaxT + 4 * bmaxT**2 + 6 * bmaxT * cmaxT + 16 * cmaxT**2
q1max = amaxT * bmaxT + 6 * cmaxT**2
q2max = bmaxT**2 + amaxT * cmaxT + 4 * cmaxT**2
bound64 = q0max * S30 + q1max * abs(L1i) + q2max * abs(L2i) + q1max + q2max
check(f'int64 overflow bound for enumeration B: {bound64} < 2^62', bound64 < 2**62
      and (amaxT + 606) * S30 + bmaxT * abs(L1i) + cmaxT * abs(L2i) < 2**62)
Agrid = np.arange(-amaxT, amaxT + 1, dtype=np.int64)[None, :]
Bgrid = np.arange(-bmaxT, bmaxT + 1, dtype=np.int64)[:, None]
absB = np.abs(Bgrid)
setB = set()
nflagged = 0
for cc in range(-cmaxT, cmaxT + 1):
    c64 = np.int64(cc)
    sgv = Agrid * S30 + Bgrid * L1i + c64 * L2i                    # sigma1 2^30 within |b| + |c| of this
    ok = (np.abs(sgv) - (absB + abs(cc))) <= 606 * S30
    q0 = Agrid * Agrid + 8 * Agrid * c64 - 4 * Bgrid * Bgrid + 6 * Bgrid * c64 + 16 * c64 * c64
    q1 = -Agrid * Bgrid - 6 * c64 * c64
    q2 = Bgrid * Bgrid - Agrid * c64 - 4 * c64 * c64
    qv = q0 * S30 + q1 * L1i + q2 * L2i                             # Q 2^30 within |q1| + |q2| of this
    ok &= (qv - (np.abs(q1) + np.abs(q2))) <= S30
    ii, jj = np.nonzero(ok)
    for i_, j_ in zip(ii.tolist(), jj.tolist()):
        v = (j_ - amaxT, i_ - bmaxT, cc)
        nflagged += 1
        if in_region(v, 606, 1):
            setB.add(v)
nbox = (2 * amaxT + 1) * (2 * bmaxT + 1) * (2 * cmaxT + 1)
check(f'enumeration B: {nbox} triples of the box scanned, {nflagged} kept by the int64 bounds and decided exactly: '
      f'the same {len(setB)} triples as enumeration A', setB == set(mem606))
note(f'(time {time.time() - t_b4:.1f} s)')
uinfo = []
for v in sorted(units606, key=lambda x: (abs(x[0]) + abs(x[1]) + abs(x[2]), x)):
    s_lo, s_hi = iv(v)
    uinfo.append(f'{v} = {name_unit(v)}, sigma1 = {float(s_lo):.10g}')
check('units with |sigma1| <= 606 and Q <= 1: ' + '; '.join(uinfo),
      sorted(units606) == sorted([(1, 0, 0), (-1, 0, 0), EPS_INV, tneg(EPS_INV)]))
dom_units = [v for v in units606 if (sgn_minus(v, 1) >= 0 or sgn_minus(v, -1) <= 0) and abs_sigma_lt(v, EPS_INV) == 1]
bnd_units = [v for v in units606 if abs_sigma_lt(v, EPS_INV) == 0]
check(f'units in the fundamental domain 1 <= |sigma1| < 1/sigma1(EPS) = sigma1(EPS^-1): {sorted(dom_units)}; '
      f'+-EPS^-1 lie exactly on the excluded end (|sigma1(+-EPS^-1)| = 1/sigma1(EPS), exactly, since EPS EPS^-1 = 1)',
      sorted(dom_units) == [(-1, 0, 0), (1, 0, 0)] and sorted(bnd_units) == sorted([EPS_INV, tneg(EPS_INV)]))
note('NOTE on the task statement: the closed region with the rounded bound 606 contains the four units +-1, +-EPS^-1;')
note('only +-1 lie in the half-open domain 1 <= |sigma1| < 1/sigma1(EPS) (= 605.08...), which is what the proof uses.')
note('Proof of (4.1) by this route: let u be a unit and e = sigma1(EPS), 1/606 < e < 1/605. |sigma1(u EPS^k)| =')
note('|sigma1(u)| e^k decreases strictly in k, is <= |sigma1(u)| 605^-k -> 0 (k -> oo) and >= |sigma1(u)| 605^|k| -> oo')
note('(k -> -oo). Let k0 be the largest k with |sigma1(u)| e^k >= 1; then v = u EPS^k0 is a unit with')
note('1 <= |sigma1(v)| < 1/e < 606 and Q(v) = 1/|sigma1(v)| <= 1 (sigma1 Q = N = +-1, Q > 0). By the enumeration')
note('v = +-1, so u = +-EPS^(-k0). k is unique: sigma1(EPS^k) = e^k is not +-1 for k != 0.')
check('uniqueness of k: 0 < sigma1(EPS) < 1, so sigma1(EPS^k) = sigma1(EPS)^k != +-1 for k != 0', 0 < se_lo and se_hi < 1)
check('the four classes of (4.1) are representable: every unit is s EPS^l v^2 with s = +-1, l in {0, 1} '
      '(EPS^k = EPS^l (EPS^k1)^2 for k = 2k1 + l), e.g. EPS^3 = EPS (EPS)^2 and EPS^-3 = EPS (EPS^-2)^2',
      POW[3] == tmul(EPS, tpow(EPS, 2)) and POW[-3] == tmul(EPS, tpow(POW[-2], 2)) and POW[-2] == tpow(EPS_INV, 2))

# ---------------------------------------------------------------------------------------------------------------------
sub('B.5 Remark (4.6)')
check('a square has sigma1 > 0, while sigma1(-1) = -1 < 0 and sigma1(-EPS) < 0 < sigma1(EPS)',
      sgn_minus((-1, 0, 0), 0) < 0 and sgn_minus(tneg(EPS), 0) < 0 < sgn_minus(EPS, 0))
check('a unit v with v^2 = EPS would have sigma1(v)^2 = 1/E and Q(v)^2 = E, so T2(v) = 1/E + 2 sqrt(E) < 114 '
      '(1/E < 1 and 4E < 113^2); by B.3 v = +-1, but (+-1)^2 = 1 != EPS',
      E_lo > 1 and 4 * E_hi < 113**2 and tmul((1, 0, 0), (1, 0, 0)) != EPS)
t2e_lo, t2e_hi = iv(t2_triple(EPS))
t2ei_lo, t2ei_hi = iv(t2_triple(EPS_INV))
note(f'T2(EPS) = {float(t2e_lo):.12g}, T2(EPS^-1) = {float(t2ei_lo):.12g} (06 prints 1210.16196999 and 366123.000052)')
check('T2(EPS) and T2(EPS^-1) agree with the values printed by 06 to 1e-9 (relative), and both exceed 114',
      abs((t2e_lo + t2e_hi) / 2 - Fr('1210.16196999')) < Fr('1210.17') * Fr(1, 10**9)
      and abs((t2ei_lo + t2ei_hi) / 2 - Fr('366123.000052')) < Fr('366124') * Fr(1, 10**9) and t2e_lo > 114)

# ======================================================================================================================
part('PART C. Brute force: every triple with |a|, |b|, |c| <= 60 and N = +-1')
# ======================================================================================================================
t_c = time.time()
ca = sp.Poly(N_res, a).all_coeffs()
check('N = a^3 + k2 a^2 + k1 a + k0 with k2 = 8c, k1 = -4b^2 + 18bc + 16c^2, k0 = -6b^3 + 24bc^2 + 36c^3 '
      '(read off the resultant)', len(ca) == 4 and ca[0] == 1 and sp.expand(ca[1] - 8*c) == 0
      and sp.expand(ca[2] - (-4*b**2 + 18*b*c + 16*c**2)) == 0
      and sp.expand(ca[3] - (-6*b**3 + 24*b*c**2 + 36*c**3)) == 0)
units60 = []
for cc in range(-60, 61):
    for bb in range(-60, 61):
        k2_ = 8 * cc
        k1_ = -4 * bb * bb + 18 * bb * cc + 16 * cc * cc
        k0_ = -6 * bb**3 + 24 * bb * cc * cc + 36 * cc**3
        for aa in range(-60, 61):
            nv = ((aa + k2_) * aa + k1_) * aa + k0_
            if nv == 1 or nv == -1:
                units60.append((aa, bb, cc))
check(f'{121**3} triples scanned exactly; the triples with N = +-1 are: '
      + '; '.join(f'{v} = {name_unit(v)}' for v in sorted(units60)),
      sorted(units60) == sorted([(1, 0, 0), (-1, 0, 0), EPS, tneg(EPS)]) and all(norm(v) in (1, -1) for v in units60))
note(f'(+-EPS^-1 = +-(95, -101, 40) and +-EPS^2 = +-{POW[2]} lie outside the box; time {time.time() - t_c:.1f} s)')

# ======================================================================================================================
part('PART D. Negative controls')
# ======================================================================================================================
sub('D.1 EPS^2 in place of EPS: domain 1 <= |sigma1| < 1/sigma1(EPS^2), rounded bound 606^2')
t_d = time.time()
EPS2, EPS2_INV = POW[2], POW[-2]
s2_lo, s2_hi = iv(EPS2)
check(f'sigma1(EPS^2) = sigma1(EPS)^2 lies in (1/606^2, 1/605^2); EPS^2 = {EPS2}, EPS^-2 = {EPS2_INV}',
      Fr(1, 606**2) < s2_lo and s2_hi < Fr(1, 605**2) and tmul(EPS2, EPS2_INV) == (1, 0, 0))
candsD1, cmaxD1, memD1, unitsD1 = scan_region(606**2, 1)
domD1 = [v for v in unitsD1 if (sgn_minus(v, 1) >= 0 or sgn_minus(v, -1) <= 0) and abs_sigma_lt(v, EPS2_INV) == 1]
bndD1 = [v for v in unitsD1 if abs_sigma_lt(v, EPS2_INV) == 0]
check(f'{len(candsD1)} candidates (|c| <= {cmaxD1}), {len(memD1)} triples with |sigma1| <= 606^2 and Q <= 1; '
      f'units in the domain: ' + ', '.join(f'{v} = {name_unit(v)}' for v in sorted(domD1))
      + f'; on the excluded end: ' + ', '.join(name_unit(v) for v in sorted(bndD1)),
      sorted(domD1) == sorted([(1, 0, 0), (-1, 0, 0), EPS_INV, tneg(EPS_INV)])
      and sorted(bndD1) == sorted([EPS2_INV, tneg(EPS2_INV)]))
check('control caught: with EPS^2 the domain contains the units +-EPS^-1 besides +-1, so the argument cannot conclude '
      '"every unit is +-EPS^(2k)"; it exhibits EPS^-1, whose inverse is the extra unit EPS',
      any(v == EPS_INV for v in domD1))
note(f'(time {time.time() - t_d:.1f} s)')
sub('D.2 EPS^2, mirrored domain sigma1(EPS^2) < |sigma1| <= 1 (so 1 <= Q < 1/sigma1(EPS^2) for units), bound 606^2')
t_d = time.time()
candsD2, cmaxD2, memD2, unitsD2 = scan_region(1, 606**2)
domD2 = []
bndD2 = []
for v in unitsD2:
    vv = v if sgn_minus(v, 0) > 0 else tneg(v)              # |sigma1(v)| = sigma1(vv)
    cmp_ = sgn_minus(tsub(vv, EPS2), 0)                     # sigma1(vv) - sigma1(EPS^2), exact sign
    if cmp_ > 0:
        domD2.append(v)
    elif cmp_ == 0 or vv == EPS2:
        bndD2.append(v)
check(f'{len(candsD2)} candidates (|c| <= {cmaxD2}), {len(memD2)} triples with |sigma1| <= 1 and Q <= 606^2; '
      f'units in the domain: ' + ', '.join(f'{v} = {name_unit(v)}' for v in sorted(domD2))
      + '; on the excluded end: ' + ', '.join(name_unit(v) for v in sorted(bndD2)),
      sorted(domD2) == sorted([(1, 0, 0), (-1, 0, 0), EPS, tneg(EPS)])
      and sorted(bndD2) == sorted([EPS2, tneg(EPS2)]))
check('control caught: the mirrored domain for EPS^2 contains the extra unit EPS itself', EPS in domD2)
note(f'(time {time.time() - t_d:.1f} s)')
sub('D.3 the same enumeration with EPS: mirrored domain sigma1(EPS) < |sigma1| <= 1, bound 606 (no extra unit)')
candsD3, cmaxD3, memD3, unitsD3 = scan_region(1, 606)
domD3 = [v for v in unitsD3 if sgn_minus(tsub(v if sgn_minus(v, 0) > 0 else tneg(v), EPS), 0) > 0]
check(f'{len(memD3)} triples with |sigma1| <= 1 and Q <= 606; units in sigma1(EPS) < |sigma1| <= 1: {sorted(domD3)} '
      f'(all units there: ' + ', '.join(name_unit(v) for v in sorted(unitsD3)) + ')',
      sorted(domD3) == [(-1, 0, 0), (1, 0, 0)])
sub('D.4 planted errors')
N_bad = N_res + a * b * c
check('control caught: the perturbed norm N + abc is not multiplicative',
      sp.expand(N_bad.xreplace({a: UV[0], b: UV[1], c: UV[2]}) - N_bad * N_bad.xreplace({a: d, b: e, c: g})) != 0)
adj_bad = (adj_ch[0], adj_ch[1] + c**2, adj_ch[2])
check('control caught: a perturbed adjugate (second coordinate + c^2) fails u u* = (N(u), 0, 0)',
      not same(smul(U, adj_bad), (N_res, 0, 0)))
rem_bad = sp.Poly(sp.rem(sp.expand(as_t(U) * as_t(V)), t**3 - 4*t - 6, t), t).all_coeffs()[::-1]
rem_bad = tuple((list(rem_bad) + [0, 0, 0])[:3])
check('control caught: the product formula (1.1) does not match the reduction modulo t^3 - 4t - 6 (wrong sign)',
      not same(rem_bad, doc_mul(U, V)))
N0 = sp.expand(sp.resultant(t**3 - 4*t, as_t(U), t))
N0f = poly_fn(N0, (a, b, c))
bad0 = [(x, y, z) for x in range(8) for y in range(8) for z in range(8)
        if N0f(x, y, z) % 8 == 0 and (x % 2 or y % 2 or z % 2)]
check(f'control caught: for t^3 - 4t (not Eisenstein, reducible) the residue test of A.4 finds {len(bad0)} non-even '
      f'residue triples with N = 0 mod 8 (e.g. {bad0[0]}), so it does not certify (1.4) there',
      len(bad0) > 0 and N0f(0, 1, 0) == 0)
S1211 = [v for v in ((x, y, z) for x in range(-30, 31) for y in range(-13, 14) for z in range(-8, 9))
         if norm(v) in (1, -1) and t2_le(v, 1211)]
S1210 = [v for v in S1211 if t2_le(v, 1210)]
amx, bmx, cmx = ellipsoid_box(1211)
check(f'control (detection): the enumeration of B.3 with the bound 1211 (box |a| <= {amx}, |b| <= {bmx}, |c| <= {cmx}) '
      f'finds the units {sorted(S1211)}, i.e. +-1 and +-EPS (T2(EPS) = 1210.16...), and with the bound 1210 only +-1',
      amx <= 30 and bmx <= 13 and cmx <= 8 and sorted(S1211) == sorted([(1, 0, 0), (-1, 0, 0), EPS, tneg(EPS)])
      and sorted(S1210) == [(-1, 0, 0), (1, 0, 0)])

# ======================================================================================================================
part('SUMMARY')
# ======================================================================================================================
print(f'ALL {NPASS} CHECKS PASSED ({time.time() - T_START:.1f} s).  Section 1: the norm form (resultant = Sylvester = '
      f'Newton), the adjugate (Cayley-Hamilton), multiplicativity, (1.4) (mod 8 exhaustive and the quoted parities), the '
      f'identities of (1.5) and the real-embedding facts of (1.6) are confirmed.  Section 4: (4.2) and the descent '
      f'constants are confirmed (rho_lo decreases with E: the sentence "both rho\'s increase" in (4.3) is a wording '
      f'error with no effect); (4.4) is confirmed (187 triples, units +-1); the unit theorem (4.1) is confirmed by a '
      f'fundamental domain for EPS ({len(mem606)} lattice points with |sigma1| <= 606, Q <= 1, found twice; units there '
      f'+-1 and +-EPS^-1, the latter on the excluded end).  Brute force |a|,|b|,|c| <= 60: units +-1, +-EPS.  '
      f'Controls: EPS^2 in place of EPS exposes +-EPS^-1 and +-EPS; planted errors are all caught.', flush=True)
