"""04_ring.py: the ring R = Z[l], l^3 = l - 3, of hypothesis_H.md, Section 1 (exact checks).

1.1 the multiplication of hring.py is the product of polynomials modulo f(t) = t^3 - t + 3, and R is a commutative
    ring with unit;
1.2 the multiplication matrix M(u) (columns u*1, u*l, u*l^2) satisfies M(uv) = M(u)M(v); N(u) = det M(u) is the cubic
    form of hring.py; N(uv) = N(u)N(v), N(r u) = r^3 N(u); the characteristic polynomial of M(u) is monic with integer
    coefficients;
1.3 the adjugate triple: u * adj(u) = (N(u), 0, 0);
1.4 the norm vanishes only at 0 (a descent modulo 2: N(a,b,c) is odd unless a, b, c are all even);
1.5 X - theta, its cofactor, the identity with F'(theta) for F(x) = x^3 - 9x + 81, and N(X - theta) = F(X);
1.6 the special elements EPS, P3A, P3B, P3C, P239 and their identities; ring maps R -> F_p that separate them;
1.7 the real place: the real root l1 of f, sigma1 and Q, the identities sigma1 Q = N, Q >= 0 and T2 = sigma1^2 + 2Q
    positive definite.
"""
import itertools
import time
from fractions import Fraction as Fr

import sympy as sp

from hlib import check, bisect_root, real_forms, iv_eval, Iv
from hring import (mul, matrix, norm, adj, trace, f, x_minus_theta, cofactor, NORM_FORM, F_COEFFS, BASIS_POW,
                   ONE, LAM, THETA, EPS, EPS_INV, P3A, P3B, P3C, P239, power, divides)

T0 = time.time()
a, b, c, d, e, g, h, i_, j = sp.symbols('a b c d e g h i j')
U, V, Wt = (a, b, c), (d, e, g), (h, i_, j)
t = sp.Symbol('t')

print("1.1 multiplication")
fpoly = sp.Poly(t**3 - t + 3, t)


def as_poly(u):
    return u[0] + u[1]*t + u[2]*t**2


check("(a,b,c)(d,e,g) is the product of a + b t + c t^2 and d + e t + g t^2 reduced modulo f(t) = t^3 - t + 3",
      sp.rem(sp.Poly(sp.expand(as_poly(U)*as_poly(V) - as_poly(mul(U, V))), t), fpoly).is_zero)
check("commutative, associative, distributive, with unit (1,0,0)",
      [sp.expand(x) for x in mul(U, V)] == [sp.expand(x) for x in mul(V, U)]
      and [sp.expand(x) for x in mul(mul(U, V), Wt)] == [sp.expand(x) for x in mul(U, mul(V, Wt))]
      and [sp.expand(x) for x in mul(U, tuple(V[k] + Wt[k] for k in range(3)))]
      == [sp.expand(x + y) for x, y in zip(mul(U, V), mul(U, Wt))]
      and mul(U, ONE) == U)

print("1.2 norm")
M = sp.Matrix(matrix(U))
check("matrix(u) has columns u*1, u*l, u*l^2",
      all(list(M[:, k]) == [sp.expand(x) for x in mul(U, (int(k == 0), int(k == 1), int(k == 2)))] for k in range(3)))
check("M(uv) = M(u) M(v)", sp.expand(sp.Matrix(matrix(mul(U, V))) - M*sp.Matrix(matrix(V))) == sp.zeros(3, 3))
check("N(a,b,c) = det M = a^3 + 2a^2c - ab^2 + 9abc + ac^2 - 3b^3 + 3bc^2 + 9c^3 (the form NORM_FORM)",
      sp.expand(M.det() - norm(U)) == 0
      and sp.expand(sum(co*a**m[0]*b**m[1]*c**m[2] for m, co in NORM_FORM.items()) - norm(U)) == 0)
r_ = sp.Symbol('r')
check("N(uv) = N(u)N(v), N(r u) = r^3 N(u), N(-u) = -N(u)",
      sp.expand(norm(mul(U, V)) - norm(U)*norm(V)) == 0 and sp.expand(norm((r_*a, r_*b, r_*c)) - r_**3*norm(U)) == 0
      and sp.expand(norm((-a, -b, -c)) + norm(U)) == 0)
T_ = sp.Symbol('T')
cp = sp.expand(M.charpoly(T_).as_expr())
check(f"characteristic polynomial of M(u): {cp}; monic, with integer coefficients for integer u, and constant term"
      " -N(u)", sp.Poly(cp, T_).LC() == 1 and sp.expand(cp.subs(T_, 0) + norm(U)) == 0
      and sp.expand(sp.Poly(cp, T_).all_coeffs()[1] + trace(U)) == 0)

disc = sp.Matrix(3, 3, lambda k, m: trace(mul(tuple(int(i == k) for i in range(3)), tuple(int(i == m) for i in range(3))))).det()
check(f"the discriminant det(Tr(b_k b_m)) of the basis (1, l, l^2) is {disc} = -239, a prime", disc == -239)

print("1.3 adjugate triple")
check("adj(a,b,c) = (a^2 + 2ac - b^2 + 3bc + c^2, -ab - 3c^2, -ac + b^2 - c^2) and u * adj(u) = (N(u), 0, 0)",
      [sp.expand(x) for x in mul(U, adj(U))] == [sp.expand(norm(U)), 0, 0])

print("1.4 the norm vanishes only at 0")
res = {}
for v in itertools.product((0, 1), repeat=3):
    res[v] = norm(v) % 2
check("modulo 2, N(a,b,c) is odd for each of the 7 nonzero residue triples (a,b,c) mod 2 (it is the norm of the field"
      " F_2[t]/(t^3 + t + 1) = F_8)", all(res[v] == 1 for v in res if v != (0, 0, 0)))
print("    Descent: if N(a,b,c) = 0 with integers not all 0, divide by gcd(a,b,c); then N is still 0 (N(r u) = r^3 N(u)),")
print("    so N is even and a, b, c are all even, against gcd 1. A rational triple scales to an integer one. So N(u) = 0")
print("    only for u = 0.")

print("1.5 X - theta and its cofactor")
X = sp.Symbol('X')
al, be = x_minus_theta(X), cofactor(X)
check("theta = 3l = (0,3,0), X - theta = (X, -3, 0), and N(X - theta) = X^3 - 9X + 81",
      THETA == (0, 3, 0) and sp.expand(norm(al) - (X**3 - 9*X + 81)) == 0)
check("beta := X^2 + theta X + theta^2 - 9 = (X^2 - 9, 3X, 9) and (X - theta) beta = X^3 - 9X + 81",
      [sp.expand(x) for x in mul(al, be)] == [sp.expand(X**3 - 9*X + 81), 0, 0]
      and [sp.expand(x) for x in be] == [sp.expand(y) for y in
                                        (lambda th2: (X**2 + th2[0] - 9, 3*X + th2[1], th2[2]))(mul(THETA, THETA))])
fprime = tuple(sp.expand(x) for x in tuple(be[k] - mul((X, 6, 0), al)[k] for k in range(3)))
check("beta - (X + 2 theta)(X - theta) = (-9, 0, 27) = 3 theta^2 - 9 = F'(theta) with F(x) = x^3 - 9x + 81 (the curve's"
      " cubic, not the f of 1.1), independent of X",
      fprime == (-9, 0, 27) and tuple(3*x - (9 if k == 0 else 0) for k, x in enumerate(mul(THETA, THETA)))
      == (-9, 0, 27))

print("1.6 special elements")
check(f"EPS = l^2 + l - 1 = {EPS} has norm 1, and EPS * EPS_INV = 1 with EPS_INV = l^2 - 2l + 2 = {EPS_INV}",
      norm(EPS) == 1 and mul(EPS, EPS_INV) == ONE)
check(f"P3A = l = {P3A}, P3B = l - 1 = {P3B}, P3C = l + 1 = {P3C} have norm -3, and P3A P3B P3C = l^3 - l = -3",
      norm(P3A) == norm(P3B) == norm(P3C) == -3 and mul(mul(P3A, P3B), P3C) == (-3, 0, 0))
check(f"P239 = 3l^2 - 1 = {P239} has norm 239 (a prime), and F'(theta) = 9 P239 = (P3A P3B P3C)^2 P239",
      norm(P239) == 239 and sp.isprime(239) and (-9, 0, 27) == tuple(9*x for x in P239)
      and mul(power(mul(mul(P3A, P3B), P3C), 2), P239) == (-9, 0, 27))
pairs = [(P3A, P3B), (P3A, P3C), (P3B, P3C)]
check("P3A, P3B, P3C are pairwise not associate: none divides another (their quotients are not integer triples)",
      all(divides(x, y) is None and divides(y, x) is None for x, y in pairs))
SET_S = {'P3A': P3A, 'P3B': P3B, 'P3C': P3C, 'P239': P239}
MAPS = [(3, 0, 'P3A'), (3, 1, 'P3B'), (3, 2, 'P3C'), (239, 124, 'P239')]


def ring_map(u, p, r):
    return (u[0] + u[1]*r + u[2]*r*r) % p


check("a second proof by hand: l -> r is a ring map R -> F_p when f(r) = 0 mod p; the maps l -> 0, 1, -1 mod 3 and"
      " l -> 124 mod 239 each send exactly one of P3A, P3B, P3C, P239 to 0 (the one named), so no one of the four divides"
      " another (if x | y, a map that kills x kills y)",
      all(f(r) % p == 0 for p, r, _ in MAPS)
      and all((ring_map(v, p, r) == 0) == (k == name) for p, r, name in MAPS for k, v in SET_S.items()))

print("1.7 the real place")
lo1, hi1 = bisect_root(f, -2, -1, 200)
check(f"f(-2) = -3 < 0 < 3 = f(-1); bisection gives l1 in [lo, hi] with width {float(hi1 - lo1):.1e} < 2^-199",
      f(Fr(-2)) == -3 and f(Fr(-1)) == 3 and hi1 - lo1 < Fr(1, 2**199) and f(lo1) <= 0 <= f(hi1))
check("f(-3/2) = 9/8 > 0, so l1 lies in (-2, -3/2) and l1^2 > 9/4 > 4/3", f(Fr(-3, 2)) == Fr(9, 8)
      and hi1 < Fr(-3, 2) and lo1**2 > Fr(9, 4))
L = sp.Symbol('L')
sigma1, Q, s_, p_ = real_forms(F_COEFFS, BASIS_POW, L)
fL = sp.Poly(f(L), L)


def zero_mod_f(expr):
    return sp.rem(sp.Poly(sp.expand(expr), L), fL).is_zero


tt = sp.Symbol('tt')
check(f"f(t) = (t - L)(t^2 - s t + p) + f(L) with s = {s_}, p = {p_}; at L = l1 the quadratic has discriminant"
      " s^2 - 4p = 4 - 3 l1^2 < 0 (l1^2 > 4/3), so it has no real root and l1 is the only real root of f",
      sp.expand((tt - L)*(tt**2 - s_*tt + p_) + f(L) - f(tt)) == 0 and lo1**2 > Fr(4, 3) and hi1**2 > Fr(4, 3))
check("sigma1(uv) = sigma1(u) sigma1(v), Q(uv) = Q(u) Q(v) and sigma1(u) Q(u) = N(u) (identities modulo f(L))",
      zero_mod_f(sigma1(mul(U, V)) - sigma1(U)*sigma1(V)) and zero_mod_f(Q(mul(U, V)) - Q(U)*Q(V))
      and zero_mod_f(sigma1(U)*Q(U) - norm(U)))
A_, B_ = sp.symbols('A B')
AU, BU = a - c*p_, b + c*s_
check(f"u = (a, b, c) reduced modulo t^2 - s t + p is A + B t with A = a - c p = a - c({p_}) and B = b + c s = b + c({s_});"
      " Q(u) = A^2 + sAB + pB^2", sp.expand(Q(U) - (AU**2 + s_*AU*BU + p_*BU**2)) == 0)
check("Q = (A + sB/2)^2 + (p - s^2/4) B^2 with p - s^2/4 = 3 l1^2/4 - 1 > 0; so Q >= 0, and Q(u) > 0 when N(u) != 0",
      sp.expand(A_**2 + s_*A_*B_ + p_*B_**2 - (A_ + s_*B_/2)**2 - (p_ - s_**2/4)*B_**2) == 0
      and sp.expand(p_ - s_**2/4 - (sp.Rational(3, 4)*L**2 - 1)) == 0 and Fr(3, 4)*min(lo1**2, hi1**2) - 1 > 0)
T2 = sp.expand(sigma1(U)**2 + 2*Q(U))
gram = [[sp.Rational(1, 2)*sp.diff(T2, U[k], U[m]) for m in range(3)] for k in range(3)]
gram_red = [[sp.rem(sp.Poly(sp.expand(gram[k][m]), L), fL).as_expr() for m in range(3)] for k in range(3)]
print("    T2 = sigma1^2 + 2Q has the Gram matrix (modulo f(L)):", gram_red)
lv = Iv(lo1, hi1)
G_iv = [[iv_eval(gram_red[k][m], L, lv) for m in range(3)] for k in range(3)]
m1 = G_iv[0][0]
m2 = G_iv[0][0]*G_iv[1][1] - G_iv[0][1]*G_iv[1][0]
det3 = (G_iv[0][0]*(G_iv[1][1]*G_iv[2][2] - G_iv[1][2]*G_iv[2][1])
        - G_iv[0][1]*(G_iv[1][0]*G_iv[2][2] - G_iv[1][2]*G_iv[2][0])
        + G_iv[0][2]*(G_iv[1][0]*G_iv[2][1] - G_iv[1][1]*G_iv[2][0]))
check(f"the leading principal minors of the Gram matrix of T2 are positive at l1: {m1}, {m2}, {det3};"
      " so T2 is positive definite", m1.lo > 0 and m2.lo > 0 and det3.lo > 0)
print(f"all checks of 04_ring.py passed ({time.time() - T0:.1f} s)")
