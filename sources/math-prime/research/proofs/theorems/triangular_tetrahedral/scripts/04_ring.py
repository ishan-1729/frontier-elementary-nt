"""04_ring.py: the ring R = Z[l] of integer triples, l^3 = 4l - 6 (hypothesis_H.md, Section 1).

A triple (a, b, c) stands for a + b l + c l^2, where l = theta/6 satisfies l^3 = 4l - 6; the definitions
(mul, matrix, norm, adj, the special elements) are in cubic.py. f(t) = t^3 - 4t + 6.

Checks, by sympy polynomial algebra and exact rational arithmetic (fractions.Fraction):
  1.1 the multiplication rule is polynomial multiplication reduced by f: u(t)v(t) - (uv)(t) = (bf + ce + cf t) f(t);
      ring laws (commutative, associative, distributive, unit (1,0,0)) as polynomial identities;
  1.2 matrix(u) has columns u*1, u*l, u*l^2; matrix(uv) = matrix(u) matrix(v); N(u) = det matrix(u) =
      a^3 + 8a^2c - 4ab^2 + 18abc + 16ac^2 - 6b^3 + 24bc^2 + 36c^3; N(uv) = N(u)N(v), N(r u) = r^3 N(u),
      N(-u) = -N(u); every triple is a root of the monic integer cubic t^3 - Tr(u) t^2 + s2(u) t - N(u);
  1.3 the adjugate triple: u * adj(u) = (N(u), 0, 0); adj(u) is the first column of the adjugate matrix;
      adj(uv) = adj(u) adj(v) and N(adj u) = N(u)^2; units are exactly the elements of norm +-1, with
      u^-1 = N(u) adj(u); a rational triple with N(u) != 0 has the inverse adj(u)/N(u);
  1.4 the norm vanishes only at 0 (2-adic descent; f is Eisenstein at 2): N - a^3 has even coefficients;
      N(2a1,b,c) = 2P with P - b^3 having even coefficients (the b^3 coefficient of P is -3);
      N(2a1,2b1,c) = 4P' with P' - c^3 having even coefficients (the c^3 coefficient of P' is 9);
  1.5 special elements: their norms; EPS*EPS_INV = 1; PI2^3 = -2 EPS, so 2 = -PI2^3 EPS^-1;
      P3A P3B P3C = 3; PI2^2 P179 = 3l^2 - 4; 36(3l^2 - 4) = EPS^-2 PI2^8 (P3A P3B P3C)^2 P179;
      f = t(t-1)(t+1) mod 3, f = t^3 mod 2, f = (t-47)^2 (t-85) mod 179, the reduction maps
      (a,b,c) -> a + b r + c r^2 mod q at these roots, and where each prime element vanishes;
      theta = 6l satisfies theta^3 - 144 theta + 1296 = 0 and N(X - theta) = X^3 - 144X + 1296;
  1.6 real values: f has exactly one real root l1, in (-13/5, -5/2) (factorization with a quadratic factor of
      negative discriminant; cross-checked by the critical values and by Sturm's theorem), enclosed in a rational
      interval of width < 2^-200 by bisection; f(t) = (t - l1)(t^2 + l1 t + l1^2 - 4), so the complex pair has
      sum s = -l1 and product p = l1^2 - 4 = -6/l1; sigma1(u) = a + b l1 + c l1^2;
      Q(u) = |a + bz + cz^2|^2 = a^2 + s ab + (s^2-2p) ac + p b^2 + sp bc + p^2 c^2 (symmetric functions),
      which reduces to a^2 - l1 ab + (8-l1^2) ac + (l1^2-4) b^2 + 6bc + (16-6l1-4l1^2) c^2 = sigma1(adj(u));
      modulo f(L) = 0: sigma1(uv) = sigma1(u)sigma1(v), sigma1(u)Q(u) = N(u), Q(uv) = Q(u)Q(v);
      Q = (a - (l1/2)b + ((8-l1^2)/2)c)^2 + ((3l1^2-16)/4)(b - l1 c)^2 >= 0, with kernel R(l1^2-4, l1, 1):
      Q is positive SEMIdefinite of rank 2 on real triples and positive on nonzero rational triples;
      T2 = sigma1^2 + 2Q = 3a^2 + 16ac + (3l1^2-8)b^2 + 8l1 bc + (32-18l1-4l1^2)c^2 is positive definite
      (leading minors 3, 9l1^2-24, 716); interval values of sigma1, Q, T2 for the special elements.
"""
from fractions import Fraction as Fr
import sympy as sp
from sympy.polys.polyfuncs import symmetrize
from cubic import (mul, add, neg, scal, matrix, norm, adj, inverse_unit, power, check, x_minus_theta,
                   ONE, LAM, EPS, EPS_INV, PI2, P3A, P3B, P3C, P179)

a, b, c, d, e, f, g, h, k, t, r, X, L = sp.symbols('a b c d e f g h k t r X L')
a1, b1 = sp.symbols('a1 b1')
F = t**3 - 4*t + 6
FL = L**3 - 4*L + 6
U, V, W = (a, b, c), (d, e, f), (g, h, k)


def as_poly(u, var=t):
    return u[0] + u[1]*var + u[2]*var**2


def from_poly(p):
    q = sp.Poly(sp.rem(sp.expand(p), F, t), t)
    cs = q.all_coeffs()[::-1] + [0, 0, 0]
    return tuple(sp.expand(x) for x in cs[:3])


def same(u, v):
    return all(sp.expand(x - y) == 0 for x, y in zip(u, v))


def red(expr):
    """Remainder of expr, a polynomial in L with coefficients in Q[a, b, ...], modulo f(L) = L^3 - 4L + 6."""
    return sp.expand(sp.rem(sp.expand(expr), FL, L))


def all_coeffs_divisible(p, m, gens):
    p = sp.Poly(sp.expand(p), *gens)
    return all(co % m == 0 for co in p.coeffs())


print("1.1 multiplication rule and ring laws")
check("mul(u,v) equals u(t) v(t) reduced modulo f(t) = t^3 - 4t + 6", same(mul(U, V), from_poly(as_poly(U)*as_poly(V))))
check("u(t) v(t) - (uv)(t) = (bf + ce + cf t) f(t) exactly (so u -> u(r) respects products at any root r of f)",
      sp.expand(as_poly(U)*as_poly(V) - as_poly(mul(U, V)) - (b*f + c*e + c*f*t)*F) == 0)
check("l^3 = 4l - 6 and l^4 = 4l^2 - 6l", mul(LAM, mul(LAM, LAM)) == (-6, 4, 0) and power(LAM, 4) == (0, -6, 4))
check("commutative", same(mul(U, V), mul(V, U)))
check("associative", same(mul(mul(U, V), W), mul(U, mul(V, W))))
check("distributive", same(mul(U, add(V, W)), add(mul(U, V), mul(U, W))))
check("(1,0,0) is the unit", same(mul(ONE, U), U))

print("1.2 multiplication matrix and the norm")


def M(u):
    return sp.Matrix(3, 3, lambda i, j: matrix(u)[i][j])


cols = [mul(U, (1, 0, 0)), mul(U, (0, 1, 0)), mul(U, (0, 0, 1))]
check("the columns of matrix(u) are u*1, u*l, u*l^2",
      sp.expand(M(U) - sp.Matrix(3, 3, lambda i, j: cols[j][i])) == sp.zeros(3, 3))
print("    matrix(a,b,c) =", M(U).tolist())
check("matrix(uv) = matrix(u) matrix(v)", sp.expand(M(mul(U, V)) - M(U)*M(V)) == sp.zeros(3, 3))
check("N(u) = det matrix(u) = a^3 + 8a^2c - 4ab^2 + 18abc + 16ac^2 - 6b^3 + 24bc^2 + 36c^3",
      sp.expand(M(U).det() - norm(U)) == 0
      and sp.expand(norm(U) - (a**3 + 8*a**2*c - 4*a*b**2 + 18*a*b*c + 16*a*c**2 - 6*b**3 + 24*b*c**2 + 36*c**3)) == 0)
check("N(uv) = N(u) N(v)", sp.expand(norm(mul(U, V)) - norm(U)*norm(V)) == 0)
check("N(r u) = r^3 N(u)", sp.expand(norm(scal(r, U)) - r**3*norm(U)) == 0)
check("N(-u) = -N(u)", sp.expand(norm(neg(U)) + norm(U)) == 0)
tr = 3*a + 8*c
s2 = 3*a**2 + 16*a*c - 4*b**2 + 18*b*c + 16*c**2
check("Tr(u) = trace of matrix(u) = 3a + 8c; s2(u) = sum of the principal 2x2 minors = 3a^2 + 16ac - 4b^2 + 18bc + 16c^2",
      sp.expand(M(U).trace() - tr) == 0
      and sp.expand(M(U).charpoly(t).as_expr() - (t**3 - tr*t**2 + s2*t - norm(U))) == 0)
sq = mul(U, U)
cube = mul(sq, U)
chi = tuple(sp.expand(cube[i] - tr*sq[i] + s2*U[i] - (norm(U) if i == 0 else 0)) for i in range(3))
check("every triple u is a root of t^3 - Tr(u) t^2 + s2(u) t - N(u), a monic cubic with integer coefficients",
      chi == (0, 0, 0))

print("1.3 adjugate triple")
check("u * adj(u) = (N(u), 0, 0)", same(mul(U, adj(U)), (norm(U), 0, 0)))
check("adj(u) is the first column of the adjugate matrix of matrix(u)", same(tuple(M(U).adjugate()[:, 0]), adj(U)))
print("    adj(a,b,c) =", adj(U))
check("adj(uv) = adj(u) adj(v) (polynomial identity)", same(adj(mul(U, V)), mul(adj(U), adj(V))))
check("N(adj u) = N(u)^2", sp.expand(norm(adj(U)) - norm(U)**2) == 0)
check("adj has integer polynomial coefficients, so for N(u) = +-1 the inverse N(u) adj(u) is an integer triple",
      all(all(co.is_integer for co in sp.Poly(x, a, b, c).coeffs()) for x in adj(U)))
check("rational inverse: u * (adj(u)/N(u)) = 1 whenever N(u) != 0 (example u = (1/2, -3, 2/7))",
      mul((Fr(1, 2), Fr(-3), Fr(2, 7)), scal(1/norm((Fr(1, 2), Fr(-3), Fr(2, 7))), adj((Fr(1, 2), Fr(-3), Fr(2, 7)))))
      == (1, 0, 0))

print("1.4 the norm vanishes only at 0 (2-adic descent)")
check("f(t) = t^3 - 4t + 6 is Eisenstein at 2: -4 and 6 are even, 6 is not divisible by 4", (-4) % 2 == 0 and 6 % 2 == 0 and 6 % 4 != 0)
check("N(a,b,c) - a^3 has even coefficients", all_coeffs_divisible(norm(U) - a**3, 2, (a, b, c)))
n1 = sp.expand(norm((2*a1, b, c)))
P = sp.expand(n1/2)
check("N(2a1,b,c) = 2P with P an integer polynomial", all_coeffs_divisible(n1, 2, (a1, b, c)))
print("    P =", P)
check("P - b^3 has even coefficients (P's b^3 coefficient is -3)", all_coeffs_divisible(P - b**3, 2, (a1, b, c)))
n2 = sp.expand(norm((2*a1, 2*b1, c)))
P2 = sp.expand(n2/4)
check("N(2a1,2b1,c) = 4P' with P' an integer polynomial", all_coeffs_divisible(n2, 4, (a1, b1, c)))
print("    P' =", P2)
check("P' - c^3 has even coefficients (the c^3 coefficient of P' is 9)", all_coeffs_divisible(P2 - c**3, 2, (a1, b1, c)))
print("    so a primitive integer triple with N = 0 would have a, then b, then c even: impossible.")

print("1.5 special elements")
expected = {"EPS": (EPS, 1), "EPS_INV": (EPS_INV, 1), "PI2": (PI2, -2), "P3A": (P3A, -3), "P3B": (P3B, -3),
            "P3C": (P3C, 3), "P179": (P179, 179)}
for name, (u, n) in expected.items():
    check(f"N({name}) = N{u} = {n}", norm(u) == n)
check("EPS * EPS_INV = 1, and EPS_INV = adj(EPS) = inverse_unit(EPS)",
      mul(EPS, EPS_INV) == ONE and adj(EPS) == EPS_INV and inverse_unit(EPS) == EPS_INV)
check("PI2^2 = (4,-6,-3) and PI2^3 = (38,10,-2) = -2*EPS", power(PI2, 2) == (4, -6, -3) and power(PI2, 3) == scal(-2, EPS))
check("2 = -PI2^3 EPS^-1", mul(neg(power(PI2, 3)), EPS_INV) == (2, 0, 0))
check("P3A * P3B * P3C = 3", mul(mul(P3A, P3B), P3C) == (3, 0, 0))
check("PI2^2 * P179 = 3l^2 - 4 = (-4,0,3), and N(3l^2 - 4) = 716 = 2^2 * 179",
      mul(power(PI2, 2), P179) == (-4, 0, 3) and norm((-4, 0, 3)) == 716)
check("3l^2 - 4 = f'(l), and disc(f) = -716 = -2^2 * 179", sp.diff(F, t) == 3*t**2 - 4 and sp.discriminant(F, t) == -716)
rhs = mul(mul(mul(power(EPS, -2), power(PI2, 8)), power(mul(mul(P3A, P3B), P3C), 2)), P179)
check("36(3l^2 - 4) = (-144,0,108) = EPS^-2 PI2^8 (P3A P3B P3C)^2 P179", rhs == scal(36, (-4, 0, 3)) == (-144, 0, 108))


def fpoly_mod(q, roots_mult):
    """f(t) - prod (t - r)^m has all coefficients divisible by q."""
    prod = sp.Integer(1)
    for rr, m in roots_mult:
        prod *= (t - rr)**m
    return all_coeffs_divisible(F - prod, q, (t,))


def red_at(u, rr, q):
    return (u[0] + u[1]*rr + u[2]*rr*rr) % q


check("f = t(t - 1)(t + 1) mod 3, so the roots of f mod 3 are 0, 1, 2", fpoly_mod(3, [(0, 1), (1, 1), (-1, 1)]))
check("f = t^3 mod 2, so the only root of f mod 2 is 0", fpoly_mod(2, [(0, 3)]))
check("f = (t - 47)^2 (t - 85) mod 179 (47 is the double root, 85 the simple one)", fpoly_mod(179, [(47, 2), (85, 1)]))
print("    By 1.1, (a,b,c) -> a + b r + c r^2 mod q is a ring homomorphism R -> F_q for each such root r.")
table = {"P3A": [red_at(P3A, rr, 3) for rr in range(3)], "P3B": [red_at(P3B, rr, 3) for rr in range(3)],
         "P3C": [red_at(P3C, rr, 3) for rr in range(3)]}
print("    values mod 3 at r = 0, 1, 2:", table)
check("P3A = l - 1 vanishes mod 3 exactly at r = 1", table["P3A"] == [2, 0, 1])
check("P3B = 2l + 5 vanishes mod 3 exactly at r = 2 = -1", table["P3B"] == [2, 1, 0])
check("P3C = l^2 - 3l + 3 vanishes mod 3 exactly at r = 0", table["P3C"] == [0, 1, 1])
check("PI2 = l^2 + l - 4 vanishes mod 2 at r = 0", red_at(PI2, 0, 2) == 0)
check("P179 = 45l^2 - 114l + 107 vanishes mod 179 at the double root 47 and not at 85 (value 144)",
      red_at(P179, 47, 179) == 0 and red_at(P179, 85, 179) == 144)
THETA = scal(6, LAM)
check("theta = 6l = (0,6,0) satisfies theta^3 - 144 theta + 1296 = 0",
      add(add(power(THETA, 3), scal(-144, THETA)), (1296, 0, 0)) == (0, 0, 0))
check("X - theta = (X,-6,0) and N(X - theta) = X^3 - 144X + 1296 (polynomial identity in X)",
      x_minus_theta(X) == (X, -6, 0) and sp.expand(norm(x_minus_theta(X)) - (X**3 - 144*X + 1296)) == 0)

print("1.6 real values")


def fq(x):
    return x**3 - 4*x + 6


check("f(-3) = -9 < 0 < 6 = f(-2)", fq(-3) == -9 and fq(-2) == 6)
check("f(-13/5) = -147/125 < 0 < 3/8 = f(-5/2), so f has a root l1 in (-13/5, -5/2)",
      fq(Fr(-13, 5)) == Fr(-147, 125) and fq(Fr(-5, 2)) == Fr(3, 8))
lo, hi = Fr(-13, 5), Fr(-5, 2)
steps = 0
while hi - lo >= Fr(1, 2**200):
    mid = (lo + hi)/2
    vm = fq(mid)
    assert vm != 0
    if vm < 0:
        lo = mid
    else:
        hi = mid
    steps += 1
check(f"bisection ({steps} steps): l1 in [lo, hi] with f(lo) < 0 < f(hi), width < 2^-200, l1 = {float(lo):.15f}",
      fq(lo) < 0 < fq(hi) and hi - lo < Fr(1, 2**200))
check("f(t) - (t - L)(t^2 + L t + L^2 - 4) = f(L) (polynomial identity), so f(t) = (t - l1)(t^2 + l1 t + l1^2 - 4)",
      sp.expand(F - (t - L)*(t**2 + L*t + L**2 - 4) - FL) == 0)
check("t^2 + L t + L^2 - 4 = (t + L/2)^2 + (3L^2 - 16)/4 (polynomial identity)",
      sp.expand(t**2 + L*t + L**2 - 4 - ((t + L/2)**2 + (3*L**2 - 16)/4)) == 0)
check("3 l1^2 - 16 > 3(5/2)^2 - 16 = 11/4 > 0, so the quadratic factor has no real root: l1 is the only real root of f",
      3*Fr(5, 2)**2 - 16 == Fr(11, 4) and 3*hi**2 - 16 > 0)
# cross-checks of 'one real root' (not needed after the factorization)
crit = sp.sqrt(3)*2/3
check("cross-check (critical values): f(2/sqrt3) = 6 - 16/(3 sqrt3), f(-2/sqrt3) = 6 + 16/(3 sqrt3), and (16/(3 sqrt3))^2 = 256/27 < 36",
      sp.simplify(F.subs(t, crit) - (6 - 16/(3*sp.sqrt(3)))) == 0 and sp.simplify(F.subs(t, -crit) - (6 + 16/(3*sp.sqrt(3)))) == 0
      and sp.Rational(256, 27) < 36)
seq = sp.sturm(sp.Poly(F, t))


def sign_changes(signs):
    signs = [s_ for s_ in signs if s_ != 0]
    return sum(1 for i in range(len(signs) - 1) if signs[i]*signs[i + 1] < 0)


at_minus = [sp.sign(q.LC())*(-1)**q.degree() for q in seq]
at_plus = [sp.sign(q.LC()) for q in seq]
print("    Sturm sequence:", [q.as_expr() for q in seq])
check("cross-check (Sturm): sign changes at -infinity minus at +infinity = 1 real root",
      sign_changes(at_minus) - sign_changes(at_plus) == 1)
check("f has no rational root (rational root test: f(r) != 0 for r = +-1, +-2, +-3, +-6), so l1 is irrational",
      all(fq(rr) != 0 for rr in (1, -1, 2, -2, 3, -3, 6, -6)))

s_, p_ = -L, L**2 - 4
check("s = -l1 and p = l1^2 - 4: (t - L)(t^2 - s t + p) = f(t) - f(L), i.e. e1 = l1 + s = 0, e2 = l1 s + p = -4, e3 = l1 p = -6 mod f(L)",
      sp.expand((t - L)*(t**2 - s_*t + p_) - (F - FL)) == 0 and red(L + s_) == 0 and red(L*s_ + p_) == -4 and red(L*p_) == -6)
check("p = -6/l1: L p + 6 = f(L)", sp.expand(L*p_ + 6 - FL) == 0)
z1, z2, S, Pp = sp.symbols('z1 z2 S P')
Qsp = a**2 + S*a*b + (S**2 - 2*Pp)*a*c + Pp*b**2 + S*Pp*b*c + Pp**2*c**2
check("(a + b z1 + c z1^2)(a + b z2 + c z2^2) = a^2 + s ab + (s^2-2p) ac + p b^2 + sp bc + p^2 c^2 with s = z1 + z2, p = z1 z2",
      sp.expand((a + b*z1 + c*z1**2)*(a + b*z2 + c*z2**2) - Qsp.subs({S: z1 + z2, Pp: z1*z2})) == 0)


def sigma1(u, var=L):
    return u[0] + u[1]*var + u[2]*var**2


def Qform(u, var=L):
    """Q(a,b,c) = a^2 - L ab + (8 - L^2) ac + (L^2 - 4) b^2 + 6 bc + (16 - 6L - 4L^2) c^2."""
    x, y, z = u
    return (x**2 - var*x*y + (8 - var**2)*x*z + (var**2 - 4)*y**2 + 6*y*z + (16 - 6*var - 4*var**2)*z**2)


QL = Qform(U)
check("with s = -L, p = L^2 - 4 the form reduces mod f(L) to Q = a^2 - L ab + (8-L^2) ac + (L^2-4) b^2 + 6bc + (16-6L-4L^2) c^2",
      red(Qsp.subs({S: s_, Pp: p_}) - QL) == 0)
check("Q(u) = sigma1(adj(u)) exactly (polynomial identity in a, b, c, L)", sp.expand(QL - sigma1(adj(U))) == 0)
check("mod f(L): sigma1(uv) = sigma1(u) sigma1(v)", red(sigma1(mul(U, V)) - sigma1(U)*sigma1(V)) == 0)
check("mod f(L): sigma1(u) Q(u) = N(u)", red(sigma1(U)*QL - norm(U)) == 0)
check("mod f(L): Q(uv) = Q(u) Q(v)", red(Qform(mul(U, V)) - QL*Qform(V)) == 0)
r1, r2, r3 = sp.symbols('r1 r2 r3')
vals = [a + b*x + c*x**2 for x in (r1, r2, r3)]
full = symmetrize(sp.expand(vals[0]*vals[1]*vals[2]), r1, r2, r3, formal=True)
esub = {}
for sname, sdef in full[2]:
    if sp.expand(sdef - (r1 + r2 + r3)) == 0:
        esub[sname] = 0
    elif sp.expand(sdef - (r1*r2 + r1*r3 + r2*r3)) == 0:
        esub[sname] = -4
    elif sp.expand(sdef - r1*r2*r3) == 0:
        esub[sname] = -6
check("derivation: the product of the three values a + b r_i + c r_i^2 is symmetric and equals N when e1 = 0, e2 = -4, e3 = -6",
      full[1] == 0 and len(esub) == len(full[2]) and sp.expand(full[0].subs(esub) - norm(U)) == 0)
e1, e2, e3 = 0, -4, -6
coeffs = {"a^3": (1, 1), "a^2 b: e1": (e1, 0), "a^2 c: e1^2 - 2e2": (e1**2 - 2*e2, 8), "a b^2: e2": (e2, -4),
          "abc: e1 e2 - 3e3": (e1*e2 - 3*e3, 18), "a c^2: e2^2 - 2e1e3": (e2**2 - 2*e1*e3, 16), "b^3: e3": (e3, -6),
          "b^2 c: e3 e1": (e3*e1, 0), "b c^2: e3 e2": (e3*e2, 24), "c^3: e3^2": (e3**2, 36)}
print("    coefficients of N from e1 = 0, e2 = -4, e3 = -6:", {kk: v[0] for kk, v in coeffs.items()})
check("these are the coefficients of N", all(v[0] == v[1] for v in coeffs.values()))
Qsq = (a - L/2*b + (8 - L**2)/2*c)**2 + (3*L**2 - 16)/4*(b - L*c)**2
check("mod f(L): Q = (a - (L/2) b + ((8 - L^2)/2) c)^2 + ((3L^2 - 16)/4) (b - L c)^2", red(Qsq - QL) == 0)
GQ = sp.Matrix([[1, -L/2, (8 - L**2)/2], [-L/2, L**2 - 4, 3], [(8 - L**2)/2, 3, 16 - 6*L - 4*L**2]])
check("Q is the quadratic form of [[1, -L/2, (8-L^2)/2], [-L/2, L^2-4, 3], [(8-L^2)/2, 3, 16-6L-4L^2]]",
      sp.expand((sp.Matrix([[a, b, c]])*GQ*sp.Matrix([a, b, c]))[0] - QL) == 0)
check("its leading minors are 1, (3L^2 - 16)/4 and det = 0 mod f(L): Q is positive SEMIdefinite of rank 2, not definite",
      sp.expand(GQ[:2, :2].det() - (3*L**2 - 16)/4) == 0 and red(GQ.det()) == 0)
check("the kernel: Q(l1^2 - 4, l1, 1) = 0 mod f(L) (an irrational direction, (p, -s, 1))", red(Qform((L**2 - 4, L, 1))) == 0)
print("    So Q(u) >= 0 for all real triples, and Q(u) = 0 only for u in R(l1^2-4, l1, 1). A nonzero rational triple")
print("    is not on that line (l1 is irrational), and also N(u) = sigma1(u) Q(u) != 0 by 1.4: Q(u) > 0 for u in R_Q, u != 0.")
T2L = 3*a**2 + 16*a*c + (3*L**2 - 8)*b**2 + 8*L*b*c + (32 - 18*L - 4*L**2)*c**2
check("mod f(L): T2 := sigma1^2 + 2Q = 3a^2 + 16ac + (3L^2-8) b^2 + 8L bc + (32-18L-4L^2) c^2",
      red(sigma1(U)**2 + 2*QL - T2L) == 0)
GT = sp.Matrix([[3, 0, 8], [0, 3*L**2 - 8, 4*L], [8, 4*L, 32 - 18*L - 4*L**2]])
check("T2 is the quadratic form of [[3,0,8],[0,3L^2-8,4L],[8,4L,32-18L-4L^2]]",
      sp.expand((sp.Matrix([[a, b, c]])*GT*sp.Matrix([a, b, c]))[0] - T2L) == 0)
check("leading minors of T2: 3, 9L^2 - 24, and det = 716 mod f(L) (= |disc f|)",
      sp.expand(GT[:2, :2].det() - (9*L**2 - 24)) == 0 and red(GT.det()) == 716)
check("9 l1^2 - 24 > 9(5/2)^2 - 24 = 129/4 > 0: T2 is positive definite (Sylvester)", 9*Fr(5, 2)**2 - 24 == Fr(129, 4) and 9*hi**2 - 24 > 0)
d2 = 3*L**2 - 8
check("completed squares, mod f(L): (3L^2-8) T2 = 3(3L^2-8)(a + 8c/3)^2 + ((3L^2-8) b + 4L c)^2 + (716/3) c^2,"
      " i.e. T2 = 3(a + 8c/3)^2 + d2 (b + (4L/d2) c)^2 + (716/(3 d2)) c^2 with d2 = 3L^2 - 8",
      red(d2*T2L - (3*d2*(a + sp.Rational(8, 3)*c)**2 + (d2*b + 4*L*c)**2 + sp.Rational(716, 3)*c**2)) == 0)


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

    def contains(self, x):
        return self.lo <= x <= self.hi

    def __str__(self):
        return f"[{float(self.lo):.12g}, {float(self.hi):.12g}]"


LI = Iv(lo, hi)


def sigma1_iv(u):
    return u[0] + u[1]*LI + u[2]*(LI*LI)


def Q_iv(u):
    return sigma1_iv(adj(u))


print("    interval values (l1 in [lo, hi]):")
print(f"    {'element':8s} {'triple':16s} {'N':>4s}  sigma1 in                          Q = sigma1(adj) in                 T2 in")
for name, (u, n) in expected.items():
    s1, qv = sigma1_iv(u), Q_iv(u)
    t2 = s1*s1 + 2*qv
    print(f"    {name:8s} {str(u):16s} {n:4d}  {str(s1):34s} {str(qv):34s} {str(t2)}")
    check(f"{name}: sigma1 and Q have determined signs, Q > 0, and the product interval of sigma1 * Q contains N = {n}",
          (s1.lo > 0 or s1.hi < 0) and qv.lo > 0 and (s1*qv).contains(n))
s1E = sigma1_iv(EPS)
check("sigma1(EPS) in (1/606, 1/605), so Q(EPS) = 1/sigma1(EPS) in (605, 606)",
      Fr(1, 606) < s1E.lo and s1E.hi < Fr(1, 605))
print("all checks of 04_ring.py passed")
