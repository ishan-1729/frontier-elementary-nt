"""07_classes.py: the class of X - theta (hypothesis_H.md, Section 5).

Notation of cubic.py: a triple (a, b, c) stands for a + b*l + c*l^2 with l = theta/6, l^3 = 4l - 6, and R is the
ring of integer triples. For an integral point (X, Y) of E: Y^2 = X^3 - 144X + 1296 put alpha = X - theta = (X, -6, 0).
The lemmas of Section 3 (Bezout, Euclid's lemma, prime elements of prime norm, exact powers, coprime square
allocation, integral closedness) and the unit theorem of Section 4 (every unit is +-EPS^k) are used as proved. This
script checks every identity, norm and congruence that Section 5 uses. All arithmetic is exact (Python integers and
sympy polynomials). mpmath only guides the search for square roots in 5.3 and 5.7, and every result it suggests is
re-verified by exact multiplication.

Checks:
  5.1 for an integer X: X - theta = (X, -6, 0), beta := X^2 + theta X + theta^2 - 144 = (X^2 - 144, 6X, 36),
      alpha*beta = (X^3 - 144X + 1296, 0, 0), N(alpha) = X^3 - 144X + 1296 and
      beta - (X + 2 theta)*alpha = 36(3l^2 - 4) (polynomial identities in X); Y = 0 is impossible, since
      X^3 - 144X + 1296 is never 0 (a parity descent); alpha is never a unit (X^3 - 144X + 1295 has no integer root);
  5.2 2 = -PI2^3 EPS^-1, 3 = P3A P3B P3C and 3l^2 - 4 = PI2^2 P179, hence
      36(3l^2 - 4) = EPS^-2 PI2^8 (P3A P3B P3C)^2 P179; the five S-primes PI2, P3A, P3B, P3C, P179 have norms
      -2, -3, -3, 3, 179 (prime in absolute value) and are pairwise non-associate;
  5.3 the statement "alpha = u PI2^e2 P3A^a P3B^b P3C^c P179^g R0^2 with u a unit and R0 in R", traced step by step
      on the 14 values of X of the 28 known integral points: the exact powers (e, f, w) of each S-prime in alpha,
      beta and Y, with e + f = 2w; alpha1*beta1 = Y1^2 for the S-free parts; an explicit Bezout identity
      lambda*alpha1 + mu*beta1 = 1 (found by integer column reduction of [M(alpha1) | M(beta1)], checked by exact
      multiplication), so alpha1 and beta1 are coprime; alpha1 = u R0^2 with u in {+-1, +-EPS}. Remark: the exact
      power of PI2 in alpha is 0, 3 or 4 as X is odd, 2 mod 4 or 0 mod 4, and 3 is impossible on E;
  5.4 the norm condition N(u)(-2)^e2 (-3)^(a+b) 3^c 179^g N(R0)^2 = Y^2: of the 2^7 vectors
      (sign, EPS, PI2, P3A, P3B, P3C, P179) in {0,1}^7 exactly 8 give a norm that is a positive square, namely those
      with e2 = g = 0, a+b+c even and sign = (-1)^(a+b): the set T = {1, EPS P3A P3B, -P3A P3C, -EPS P3B P3C} and the
      set {EPS, P3A P3B, -EPS P3A P3C, -P3B P3C}, which is EPS*T modulo squares;
  5.5 modulo 8: mu0*S^2 mod 8 depends only on S mod 8; the number of S in (Z/8)^3 for which mu0*S^2 has
      l^2-coordinate = 0 and l-coordinate = -6 (mod 8) is 0 for the four mu0 of EPS*T, and 32, 48, 48, 32 for
      1, EPS P3A P3B, -P3A P3C, -EPS P3B P3C. Remark: modulo 4 the counts are 0 and 8, 16, 16, 8 (so 4 suffices);
      modulo 2 they are 2 and 4 (2 does not); a hand proof modulo 4 (the l-coordinate forces x even, and then a
      table of four residues of (y, z) mod 2) is checked for each of the four mu0 of EPS*T;
  5.6 9 - theta = EPS P3A P3B (-45 + 48l - 19l^2)^2, -15 - theta = -P3A P3C (5 + 2l)^2 and
      -12 - theta = -EPS P3B P3C (-8 + 8l - 3l^2)^2; (9-theta)(-15-theta)(-12-theta) is a square in R.
      The four mu = 1, 9-theta, -15-theta, -12-theta are pairwise distinct modulo K*^2. Ingredients: evaluation at a
      root r of l^3 - 4l + 6 modulo a prime q is a ring homomorphism R -> F_q; the Legendre symbol is multiplicative;
      every triple is a root of its monic integer characteristic polynomial (for integral closedness). The
      characters at (q, r) = (7, 3) and (23, 10) take the four distinct values (+,+), (-,+), (+,-), (-,-). A second
      certificate at the completely split primes, (61, 45) and (31, 19), needs no integral closedness.
      Remark: the eight admissible mu0 of 5.4 are pairwise distinct modulo K*^2 ((31,18), (31,25), (61,45));
  5.7 confirmation, not used by the proof: for each of the 14 values of X, X - theta = mu0 S^2 with an explicit S in R
      and mu0 in T, and (X - theta) T^2 = mu S^2 with mu the matching element of {1, 9-theta, -15-theta,
      -12-theta}; each of the other seven admissible mu0 is excluded by a character; the classes agree with the
      2-descent (targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt, CURVE 1611.a1);
  5.8 confirmation, not used by the proof: no integer X with |X| <= 10^5 other than the 14 values makes
      X^3 - 144X + 1296 a square.
"""
import functools
import itertools
import math
import sympy as sp
import mpmath as mp
from cubic import (mul, add, sub, scal, matrix, norm, adj, power, divides, check, x_minus_theta,
                   ONE, LAM, EPS, EPS_INV, PI2, P3A, P3B, P3C, P179)

X, X1, X2, k = sp.symbols('X X1 X2 k')
THETA = (0, 6, 0)                       # theta = 6l
DELTA = (-144, 0, 108)                  # 36(3l^2 - 4)
KNOWN_X = [-15, -12, -8, 0, 4, 9, 12, 24, 40, 48, 108, 252, 420, 4953]
SPRIMES = [("PI2", PI2), ("P3A", P3A), ("P3B", P3B), ("P3C", P3C), ("P179", P179)]


def same(u, v):
    return all(sp.expand(p - q) == 0 for p, q in zip(u, v))


def prod(factors):
    return functools.reduce(mul, factors, ONE)


def even_coeffs(p, gens, m=2):
    """All coefficients of the polynomial p in gens are divisible by m."""
    return all(co % m == 0 for co in sp.Poly(sp.expand(p), *gens).coeffs())


def exact_power(p, u):
    """(e, v) with u = p^e v and p not dividing v (u nonzero)."""
    e = 0
    while True:
        q = divides(p, u)
        if q is None:
            return e, u
        e, u = e + 1, q


# The admissible representatives of 5.4, by name. T is the set of the four survivors of 5.5.
MU0 = {
    "1": ONE,
    "EPS P3A P3B": prod([EPS, P3A, P3B]),
    "-P3A P3C": scal(-1, prod([P3A, P3C])),
    "-EPS P3B P3C": scal(-1, prod([EPS, P3B, P3C])),
    "EPS": EPS,
    "P3A P3B": prod([P3A, P3B]),
    "-EPS P3A P3C": scal(-1, prod([EPS, P3A, P3C])),
    "-P3B P3C": scal(-1, prod([P3B, P3C])),
}
T_NAMES = ["1", "EPS P3A P3B", "-P3A P3C", "-EPS P3B P3C"]
EPST_NAMES = ["EPS", "P3A P3B", "-EPS P3A P3C", "-P3B P3C"]

# The four classes of (H), with the squares that identify them (5.6).
MU = {"1": (1, 0, 0), "9-theta": (9, -6, 0), "-15-theta": (-15, -6, 0), "-12-theta": (-12, -6, 0)}
IDENT = {"1": ("1", ONE), "9-theta": ("EPS P3A P3B", (-45, 48, -19)),
         "-15-theta": ("-P3A P3C", (5, 2, 0)), "-12-theta": ("-EPS P3B P3C", (-8, 8, -3))}

# ---------------------------------------------------------------- guided square roots (results checked exactly)
mp.mp.dps = 50
_roots = mp.polyroots([1, 0, -4, 6], maxsteps=200, extraprec=200)
L_REAL = mp.re([r for r in _roots if abs(mp.im(r)) < mp.mpf(10)**-30][0])
L_CPX = [r for r in _roots if mp.im(r) > mp.mpf(10)**-30][0]
_VINV = mp.matrix([[1, L_REAL, L_REAL**2], [1, L_CPX, L_CPX**2], [1, mp.conj(L_CPX), mp.conj(L_CPX)**2]])**-1


def sqrt_in_R(g):
    """An integer triple S with S*S = g exactly, found by rounding the guided values; None if the guide finds none."""
    g_real = g[0] + g[1]*L_REAL + g[2]*L_REAL**2
    if g_real <= 0:
        return None
    g_cpx = g[0] + g[1]*L_CPX + g[2]*L_CPX**2
    for sign in (1, -1):
        s_cpx = sign*mp.sqrt(g_cpx)
        coords = _VINV*mp.matrix([mp.sqrt(g_real), s_cpx, mp.conj(s_cpx)])
        S = tuple(int(mp.nint(mp.re(coords[i]))) for i in range(3))
        if mul(S, S) == tuple(g):
            return S
    return None


def bezout(a1, b1):
    """(lam, mu) in R with lam*a1 + mu*b1 = 1, or None if the ideal a1 R + b1 R is not all of R.

    a1 R + b1 R is the Z-span of the six columns of M(a1) and M(b1); integer column operations bring it to a
    triangular basis, whose diagonal gives the index of the span in Z^3."""
    basis = (ONE, LAM, mul(LAM, LAM))
    cols = [list(mul(a1, e)) for e in basis] + [list(mul(b1, e)) for e in basis]
    tags = [[int(i == j) for i in range(6)] for j in range(6)]
    active, piv = list(range(6)), []
    for row in range(3):
        while True:
            nz = [j for j in active if cols[j][row] != 0]
            if not nz:
                return None
            jm = min(nz, key=lambda j: abs(cols[j][row]))
            others = [j for j in nz if j != jm]
            if not others:
                break
            for j in others:
                qt = cols[j][row] // cols[jm][row]
                cols[j] = [x - qt*y for x, y in zip(cols[j], cols[jm])]
                tags[j] = [x - qt*y for x, y in zip(tags[j], tags[jm])]
        piv.append(jm)
        active.remove(jm)
    p0, p1, p2 = (cols[j] for j in piv)
    if abs(p0[0]*p1[1]*p2[2]) != 1:
        return None
    c0 = p0[0]                            # 1/d0 = d0 for d0 = +-1
    c1 = -c0*p0[1]*p1[1]
    c2 = -(c0*p0[2] + c1*p1[2])*p2[2]
    t = [c0*x + c1*y + c2*z for x, y, z in zip(tags[piv[0]], tags[piv[1]], tags[piv[2]])]
    return tuple(t[:3]), tuple(t[3:])


def unit_name(u):
    """'+-EPS^k' if the triple u is +-EPS^k with |k| <= 20 (checked exactly), else ''."""
    for kk in range(-20, 21):
        for s, sg in ((1, ""), (-1, "-")):
            if scal(s, power(EPS, kk)) == tuple(u):
                return f" = {sg}EPS^{kk}"
    return ""


def legendre(x, q):
    """By definition: 0 if q | x, +1 if x is a nonzero square mod q, -1 otherwise."""
    x %= q
    if x == 0:
        return 0
    return 1 if x in {t*t % q for t in range(1, q)} else -1


def ev(u, r, q):
    return (u[0] + u[1]*r + u[2]*r*r) % q


def roots_mod(q):
    return [r for r in range(q) if (r**3 - 4*r + 6) % q == 0]


# =====================================================================================================================
print("5.1 the two factors of X^3 - 144X + 1296")
alpha = x_minus_theta(X)
beta = (X**2 - 144, 6*X, 36)
check("theta = 6l = (0,6,0) satisfies theta^3 - 144 theta + 1296 = 0",
      add(sub(power(THETA, 3), scal(144, THETA)), (1296, 0, 0)) == (0, 0, 0))
check("X - theta = (X, -6, 0)", same(alpha, sub((X, 0, 0), THETA)))
check("beta = X^2 + theta X + theta^2 - 144 = (X^2 - 144, 6X, 36)",
      same(add(add((X**2 - 144, 0, 0), scal(X, THETA)), mul(THETA, THETA)), beta))
check("alpha*beta = (X^3 - 144X + 1296, 0, 0)", same(mul(alpha, beta), (X**3 - 144*X + 1296, 0, 0)))
check("N(alpha) = X^3 - 144X + 1296", sp.expand(norm(alpha) - (X**3 - 144*X + 1296)) == 0)
check("X + 2 theta = (X, 12, 0) and beta - (X + 2 theta)*alpha = (-144, 0, 108) = 36(3l^2 - 4) for every X",
      same(add((X, 0, 0), scal(2, THETA)), (X, 12, 0)) and same(sub(beta, mul((X, 12, 0), alpha)), DELTA)
      and DELTA == scal(36, (-4, 0, 3)))
print("    so for an integral point (X, Y): alpha*beta = Y^2 (as triples), N(alpha) = Y^2, and every common divisor")
print("    of alpha and beta divides 36(3l^2 - 4).")
print("    Y = 0 is impossible: by Section 1, N(alpha) = 0 only for alpha = 0, and alpha has l-coordinate -6.")
print("    Directly, X^3 - 144X + 1296 is never 0 for an integer X:")
check("X = 2k+1: X^3 - 144X + 1296 - 1 has even coefficients in k, so the value is odd",
      even_coeffs((2*k + 1)**3 - 144*(2*k + 1) + 1296 - 1, [k]))
check("X = 2X1: X^3 - 144X + 1296 = 8(X1^3 - 36X1 + 162)",
      sp.expand((2*X1)**3 - 144*2*X1 + 1296 - 8*(X1**3 - 36*X1 + 162)) == 0)
check("X1 = 2k+1: X1^3 - 36X1 + 162 - 1 has even coefficients in k, so it is odd",
      even_coeffs((2*k + 1)**3 - 36*(2*k + 1) + 162 - 1, [k]))
check("X1 = 2X2: X1^3 - 36X1 + 162 = 2(4X2^3 - 36X2 + 81), and 4X2^3 - 36X2 + 81 - 1 has even coefficients",
      sp.expand((2*X2)**3 - 36*2*X2 + 162 - 2*(4*X2**3 - 36*X2 + 81)) == 0
      and even_coeffs(4*X2**3 - 36*X2 + 81 - 1, [X2]))
print("    So Y != 0, and alpha, beta and the triple (Y, 0, 0) are nonzero.")
unit_roots = [d for dd in sp.divisors(1295) for d in (dd, -dd) if d**3 - 144*d + 1295 == 0]
check("alpha is never a unit: N(alpha) = Y^2 >= 0, and N(alpha) = 1 needs an integer root of X^3 - 144X + 1295,"
      " which would divide 1295 = 5*7*37; none of the 16 divisors +-d is a root", unit_roots == [])
print("    (Not needed: if alpha were a unit, 5.3 would hold with all exponents 0 and R0 = 1.)")

# =====================================================================================================================
print("5.2 the common divisor 36(3l^2 - 4) and the five S-primes")
check("EPS*EPS_INV = 1 and N(EPS) = 1", mul(EPS, EPS_INV) == ONE and norm(EPS) == 1)
check("PI2^3 = -2 EPS, so 2 = -PI2^3 EPS^-1", power(PI2, 3) == scal(-2, EPS) and mul(power(PI2, 3), EPS_INV) == (-2, 0, 0))
check("P3A*P3B*P3C = 3", prod([P3A, P3B, P3C]) == (3, 0, 0))
check("PI2^2*P179 = (-4, 0, 3) = 3l^2 - 4", mul(power(PI2, 2), P179) == (-4, 0, 3))
check("36 = EPS^-2 PI2^6 (P3A P3B P3C)^2", prod([power(EPS, -2), power(PI2, 6), power(prod([P3A, P3B, P3C]), 2)]) == (36, 0, 0))
check("36(3l^2 - 4) = EPS^-2 PI2^8 (P3A P3B P3C)^2 P179",
      prod([power(EPS, -2), power(PI2, 8), power(prod([P3A, P3B, P3C]), 2), P179]) == DELTA)
check("N(36(3l^2 - 4)) = 36^3 * 716 = 2^8 * 3^6 * 179 = (-2)^8 ((-3)(-3)3)^2 179",
      norm(DELTA) == 36**3*716 == 2**8*3**6*179 == (-2)**8*((-3)*(-3)*3)**2*179)
NORMS = {"PI2": -2, "P3A": -3, "P3B": -3, "P3C": 3, "P179": 179}
for name, p in SPRIMES:
    check(f"N({name}) = N{p} = {NORMS[name]}, and |N| = {abs(NORMS[name])} is prime, so {name} is a prime element (3.3)",
          norm(p) == NORMS[name] and sp.isprime(abs(NORMS[name])))
for (n1, p1), (n2, p2) in itertools.combinations(SPRIMES, 2):
    if abs(norm(p1)) != abs(norm(p2)):
        check(f"{n1}, {n2} are not associates: |N| = {abs(norm(p1))} and {abs(norm(p2))} differ",
              abs(norm(p1)) != abs(norm(p2)))
    else:
        check(f"{n1}, {n2} are not associates: {n1}/{n2} and {n2}/{n1} are not integral triples",
              divides(p2, p1) is None and divides(p1, p2) is None)
tbl = {name: [ev(p, r, 3) for r in (0, 1, 2)] for name, p in SPRIMES[1:4]}
check("second view: l^3 - 4l + 6 = l(l-1)(l+1) mod 3, and P3A, P3B, P3C vanish mod 3 at the single roots"
      f" l = 1, 2, 0 respectively (values at l = 0,1,2: {tbl}); a unit vanishes at no root",
      even_coeffs(k**3 - 4*k + 6 - k*(k - 1)*(k + 1), [k], 3)
      and [r for r in (0, 1, 2) if tbl["P3A"][r] == 0] == [1] and [r for r in (0, 1, 2) if tbl["P3B"][r] == 0] == [2]
      and [r for r in (0, 1, 2) if tbl["P3C"][r] == 0] == [0])
check("by the way: l = PI2*P3C (so theta = 6l = -EPS^-1 PI2^4 P3A P3B P3C^2)",
      mul(PI2, P3C) == LAM and prod([scal(-1, EPS_INV), power(PI2, 4), P3A, P3B, power(P3C, 2)]) == THETA)

# =====================================================================================================================
print("5.3 alpha = u PI2^e2 P3A^a P3B^b P3C^c P179^g R0^2 on the 28 known integral points")
print("    For each X: exact powers (e, f, w) of each S-prime in alpha, beta, Y (by exact division); the S-free parts")
print("    alpha1, beta1, Y1; an explicit Bezout identity for alpha1, beta1; alpha1 = u R0^2 with u in {+-1, +-EPS}.")
REP = {}
for Xv in KNOWN_X:
    Y2 = Xv**3 - 144*Xv + 1296
    Yv = math.isqrt(Y2)
    al, be, yy = x_minus_theta(Xv), (Xv**2 - 144, 6*Xv, 36), (Yv, 0, 0)
    ok = Yv > 0 and Yv*Yv == Y2 and mul(al, be) == mul(yy, yy)
    a1, b1, y1 = al, be, yy
    exps = {}
    for name, p in SPRIMES:
        e, a1 = exact_power(p, a1)
        f, b1 = exact_power(p, b1)
        w, y1 = exact_power(p, y1)
        exps[name] = (e, f, w)
        ok = ok and e + f == 2*w
    ok = ok and all(divides(p, t) is None for _, p in SPRIMES for t in (a1, b1, y1))
    ok = ok and al == mul(prod([power(p, exps[n][0]) for n, p in SPRIMES]), a1)
    ok = ok and be == mul(prod([power(p, exps[n][1]) for n, p in SPRIMES]), b1)
    ok = ok and yy == mul(prod([power(p, exps[n][2]) for n, p in SPRIMES]), y1)
    ok = ok and mul(a1, b1) == mul(y1, y1)
    bz = bezout(a1, b1)
    ok = ok and bz is not None and add(mul(bz[0], a1), mul(bz[1], b1)) == ONE
    found = None
    for uname, u in (("1", ONE), ("-1", scal(-1, ONE)), ("EPS", EPS), ("-EPS", scal(-1, EPS))):
        R0 = sqrt_in_R(divides(u, a1))
        if R0 is not None:
            found = (uname, u, R0)
            break
    ok = ok and found is not None and mul(found[1], mul(found[2], found[2])) == a1
    e2, a, b, c, g = (exps[n][0] for n, _ in SPRIMES)
    uname, u, R0 = found
    ok = ok and al == prod([u, power(PI2, e2), power(P3A, a), power(P3B, b), power(P3C, c), power(P179, g), R0, R0])
    sigma, kk = (1 if uname in ("1", "EPS") else -1), (1 if "EPS" in uname else 0)
    mu0 = scal(sigma, prod([power(EPS, kk), power(P3A, a % 2), power(P3B, b % 2), power(P3C, c % 2)]))
    mu0_name = [nm for nm, m in MU0.items() if m == mu0]
    parity_ok = e2 % 2 == 0 and g % 2 == 0 and (a + b + c) % 2 == 0 and norm(u)*(-1)**(a + b) == 1
    REP[Xv] = mu0_name[0] if mu0_name else None
    print(f"    X = {Xv:5d}, Y = {Yv:6d}: (e,f,w) = " + ", ".join(f"{n} {exps[n]}" for n, _ in SPRIMES))
    print(f"        N(alpha1) = {norm(a1)}, N(beta1) = {norm(b1)}, gcd of the norms {math.gcd(norm(a1), norm(b1))};"
          f" u = {uname}, R0 = {R0}{unit_name(R0)}")
    check(f"X = {Xv}: e + f = 2w for every S-prime, alpha1*beta1 = Y1^2, lambda*alpha1 + mu*beta1 = 1 (exact),"
          f" alpha = u PI2^{e2} P3A^{a} P3B^{b} P3C^{c} P179^{g} R0^2 (exact)", ok)
    check(f"X = {Xv}: e2, g, a+b+c even and N(u)(-1)^(a+b) = 1, so mu0 = {REP[Xv]} (in T)",
          parity_ok and REP[Xv] in T_NAMES)
print("    remark (not used): the exact power of PI2 in alpha")
a_2 = (2*X1, -6, 0)
q3 = tuple(sp.expand(x/norm(power(PI2, 3))) for x in mul(a_2, adj(power(PI2, 3))))
a_4 = (4*X2, -6, 0)
q4 = tuple(sp.expand(x/norm(power(PI2, 4))) for x in mul(a_4, adj(power(PI2, 4))))
print("    X odd: N(alpha) is odd (5.1), so PI2 does not divide alpha (N(PI2) = -2 would divide N(alpha)).")
check(f"X = 2X1: alpha = PI2^3 * {q3} with N({q3}) = -(X1^3 - 36X1 + 162), odd for X1 odd (5.1): exact power 3",
      all(co.is_integer for x in q3 for co in sp.Poly(x, X1).coeffs()) and same(mul(power(PI2, 3), q3), a_2)
      and sp.expand(norm(q3) + X1**3 - 36*X1 + 162) == 0)
print("    X = 2X1 with X1 odd: N(alpha) = 8(X1^3 - 36X1 + 162) has 2-adic valuation 3, so it is not a square:")
print("    no integral point has X = 2 (mod 4) (the parity argument of 5.4 in this case).")
check(f"X = 4X2: alpha = PI2^4 * {q4} with N({q4}) = 4X2^3 - 36X2 + 81 odd: exact power 4",
      all(co.is_integer for x in q4 for co in sp.Poly(x, X2).coeffs()) and same(mul(power(PI2, 4), q4), a_4)
      and sp.expand(norm(q4) - (4*X2**3 - 36*X2 + 81)) == 0 and even_coeffs(4*X2**3 - 36*X2 + 80, [X2]))
print("    So on E the exact power of PI2 in alpha is 0 (X odd) or 4 (4 | X); in the examples above PI2^8 exactly")
print("    divides beta and PI2^6 exactly divides Y when 4 | X. The argument needs no bound on these exponents.")

# =====================================================================================================================
print("5.4 the norm condition and the eight admissible mu0")
VEC_NAMES = ("sign", "EPS", "PI2", "P3A", "P3B", "P3C", "P179")
admissible, formula_ok, predicted_ok = [], True, True
for vec in itertools.product((0, 1), repeat=7):
    s, kk, e2, a, b, c, g = vec
    mu = scal((-1)**s, prod([power(EPS, kk), power(PI2, e2), power(P3A, a), power(P3B, b), power(P3C, c), power(P179, g)]))
    n = norm(mu)
    formula_ok = formula_ok and n == (-1)**s * (-2)**e2 * (-3)**a * (-3)**b * 3**c * 179**g
    adm = n > 0 and math.isqrt(n)**2 == n
    predicted_ok = predicted_ok and adm == (e2 == 0 and g == 0 and (a + b + c) % 2 == 0 and (-1)**s == (-1)**(a + b))
    if adm:
        admissible.append((vec, mu, n))
check("for all 128 vectors, N(sigma EPS^k PI2^e2 P3A^a P3B^b P3C^c P179^g) = sigma (-2)^e2 (-3)^a (-3)^b 3^c 179^g",
      formula_ok)
check("the norm is a positive square exactly when e2 = g = 0, a+b+c is even and sigma = (-1)^(a+b)", predicted_ok)
check("exactly 8 of the 128 vectors are admissible", len(admissible) == 8)
for vec, mu, n in admissible:
    nm = [x for x, m in MU0.items() if m == mu][0]
    print(f"    {dict(zip(VEC_NAMES, vec))}  mu0 = {nm:13s} = {str(mu):18s} N = {n}"
          f"  ({'T' if nm in T_NAMES else 'EPS*T'})")
check("the 8 admissible mu0 are exactly 1, EPS P3A P3B, -P3A P3C, -EPS P3B P3C, EPS, P3A P3B, -EPS P3A P3C, -P3B P3C",
      sorted(m for _, m, _ in admissible) == sorted(MU0.values()))
check("EPS*T = {EPS, P3A P3B, -EPS P3A P3C, -P3B P3C} modulo squares: EPS*mu0 for mu0 in T is the listed element"
      " times 1 or EPS^2",
      all(mul(EPS, MU0[t]) in (MU0[s], mul(MU0[s], power(EPS, 2))) for t, s in zip(T_NAMES, EPST_NAMES)))
print("    Reduction used in the text: with u = sigma EPS^k (Section 4) and the exponents of 5.3,")
print("    alpha = mu0 S^2, mu0 = sigma EPS^(k mod 2) P3A^(a mod 2) P3B^(b mod 2) P3C^(c mod 2), S in R")
print("    (S = EPS^floor(k/2) PI2^(e2/2) P3A^floor(a/2) P3B^floor(b/2) P3C^floor(c/2) P179^(g/2) R0; EPS^-1 is in R).")

# =====================================================================================================================
print("5.5 the congruence modulo 8")
s0, s1, s2, t0, t1, t2, m0, m1, m2 = sp.symbols('s0 s1 s2 t0 t1 t2 m0 m1 m2')
Sg, Tg, Mg = (s0, s1, s2), (t0, t1, t2), (m0, m1, m2)
shifted = add(Sg, scal(8, Tg))
diff = sub(mul(Mg, mul(shifted, shifted)), mul(Mg, mul(Sg, Sg)))
check("mu*(S + 8T)^2 - mu*S^2 has all coefficients divisible by 8 (indeed 16), for generic integer triples mu, S, T:"
      " mu0*S^2 mod 8 depends only on S mod 8", all(even_coeffs(d, Sg + Tg + Mg, 16) for d in diff))
sq = mul(Sg, Sg)
print(f"    S^2 = {sq}")


def count(mu, q):
    n = 0
    for S in itertools.product(range(q), repeat=3):
        t = mul(mu, mul(S, S))
        if t[2] % q == 0 and (t[1] + 6) % q == 0:
            n += 1
    return n


EXPECTED8 = {"EPS": 0, "P3A P3B": 0, "-EPS P3A P3C": 0, "-P3B P3C": 0,
             "1": 32, "EPS P3A P3B": 48, "-P3A P3C": 48, "-EPS P3B P3C": 32}
for nm in EPST_NAMES + T_NAMES:
    mu = MU0[nm]
    coords = mul(mu, sq)
    print(f"    mu0 = {nm:13s} = {str(mu):18s} l-coordinate of mu0*S^2: {sp.expand(coords[1])}")
    print(f"    {'':38s}l^2-coordinate: {sp.expand(coords[2])}")
    c8 = count(mu, 8)
    check(f"mu0 = {nm}: {c8} residues S mod 8 give l^2-coordinate = 0 and l-coordinate = -6 (mod 8)"
          f" (expected {EXPECTED8[nm]})", c8 == EXPECTED8[nm])
print("    So alpha = (X, -6, 0) = mu0 S^2 with S in R is impossible for the four mu0 of EPS*T, for every X;")
print("    the four mu0 of T survive.")
print("    remark: the same counts modulo 4 and modulo 2")
c4 = {nm: count(MU0[nm], 4) for nm in EPST_NAMES + T_NAMES}
c2 = {nm: count(MU0[nm], 2) for nm in EPST_NAMES + T_NAMES}
check(f"modulo 4: {c4}; the four of EPS*T already have no solution modulo 4",
      all(c4[nm] == 0 for nm in EPST_NAMES) and [c4[nm] for nm in T_NAMES] == [8, 16, 16, 8])
check(f"modulo 2: {c2}; modulo 2 does not exclude them",
      all(c2[nm] == 2 for nm in EPST_NAMES) and all(c2[nm] == 4 for nm in T_NAMES))
xs, ys, zs = sp.symbols('x y z')
S4 = (xs, ys, zs)
shifted4 = add(Sg, scal(4, Tg))
diff4 = sub(mul(Mg, mul(shifted4, shifted4)), mul(Mg, mul(Sg, Sg)))
check("likewise mu*(S + 4T)^2 - mu*S^2 has all coefficients divisible by 8: mu0*S^2 mod 4 depends only on S mod 4",
      all(even_coeffs(d, Sg + Tg + Mg, 8) for d in diff4))
print("    remark: a hand proof modulo 4 for the four mu0 of EPS*T, with S = (x, y, z) (target: (l, l^2) = (2, 0) mod 4)")


def red4(mu, x0, y0, z0):
    t = mul(mu, mul((x0, y0, z0), (x0, y0, z0)))
    return (t[1] % 4, t[2] % 4)


for nm in EPST_NAMES:
    mu = MU0[nm]
    m1 = mu[1]
    step1 = m1 % 2 == 1 and even_coeffs(mul(mu, mul(S4, S4))[1] - m1*xs**2, [xs, ys, zs])
    step2 = all(red4(mu, x0, y0, z0) == red4(mu, 0, y0 % 2, z0 % 2)
                for x0 in (0, 2) for y0 in range(4) for z0 in range(4))
    table = {(y0, z0): red4(mu, 0, y0, z0) for y0 in (0, 1) for z0 in (0, 1)}
    check(f"mu0 = {nm}: l-coordinate = ({m1})x^2 (mod 2) with {m1} odd, so x is even; for x even the"
          f" (l, l^2)-coordinates mod 4 depend only on (y, z) mod 2: {table}; never (2, 0)",
          step1 and step2 and all(v != (2, 0) for v in table.values()))

# =====================================================================================================================
print("5.6 identification with the classes of (H), and their distinctness")
for mname, (m0name, T) in IDENT.items():
    if mname == "1":
        continue
    check(f"{mname} = ({m0name}) * {T}^2, N{T} = {norm(T)}", mul(MU0[m0name], mul(T, T)) == MU[mname])
print("    So X - theta = mu0 S^2 with mu0 in T gives X - theta = mu (S/T)^2 with mu = 1, 9-theta, -15-theta,"
      " -12-theta: this is (H).")
T1, T2, T3 = IDENT["9-theta"][1], IDENT["-15-theta"][1], IDENT["-12-theta"][1]
check("(9-theta)(-15-theta)(-12-theta) = (3 EPS T1 T2 T3)^2, a square in R",
      prod([MU["9-theta"], MU["-15-theta"], MU["-12-theta"]]) == power(prod([(3, 0, 0), EPS, T1, T2, T3]), 2))
check("T is closed modulo squares: (EPS P3A P3B)(-P3A P3C) = (-EPS P3B P3C) P3A^2",
      mul(MU0["EPS P3A P3B"], MU0["-P3A P3C"]) == mul(MU0["-EPS P3B P3C"], mul(P3A, P3A)))
print("    quadratic characters")
a_, b_, c_, d_, e_, f_, r_ = sp.symbols('a b c d e f r')
U, V = (a_, b_, c_), (d_, e_, f_)
prod_r = mul(U, V)
check("u(r)v(r) - (uv)(r) = (r^3 - 4r + 6)(bf + ce + cf r) for triples u = (a,b,c), v = (d,e,f): evaluation at a root"
      " r of l^3 - 4l + 6 mod q is a ring homomorphism R -> F_q",
      sp.expand((a_ + b_*r_ + c_*r_**2)*(d_ + e_*r_ + f_*r_**2) - (prod_r[0] + prod_r[1]*r_ + prod_r[2]*r_**2)
                - (r_**3 - 4*r_ + 6)*(b_*f_ + c_*e_ + c_*f_*r_)) == 0)
Mg_ = sp.Matrix(matrix(U))
tt = sp.symbols('t')
cp = sp.Poly(Mg_.charpoly(tt).as_expr(), tt)
s1c, s2c, s3c = -cp.coeff_monomial(tt**2), cp.coeff_monomial(tt), -cp.coeff_monomial(1)
cube = mul(mul(U, U), U)
ch = tuple(sp.expand(cube[i] - s1c*mul(U, U)[i] + s2c*U[i] - (s3c if i == 0 else 0)) for i in range(3))
check(f"every triple u is a root of t^3 - s1 t^2 + s2 t - N(u) with s1 = {sp.expand(s1c)},"
      f" s2 = {sp.expand(s2c)} (integer polynomials): so eta^2 = u in R makes eta a root of the monic"
      " t^6 - s1 t^4 + s2 t^2 - N(u)",
      ch == (0, 0, 0) and sp.expand(s3c - norm(U)) == 0
      and all(co.is_integer for p in (s1c, s2c) for co in sp.Poly(sp.expand(p), a_, b_, c_).coeffs()))
for q in (7, 23, 31, 61):
    check(f"Legendre symbol mod {q} is multiplicative: leg(xy) = leg(x)leg(y) for all x, y in 1..{q - 1}",
          all(legendre(x*y, q) == legendre(x, q)*legendre(y, q) for x in range(1, q) for y in range(1, q)))
check("roots: l^3 - 4l + 6 = 0 at l = 3 mod 7 (its only root) and l = 10 mod 23 (its only root)",
      roots_mod(7) == [3] and roots_mod(23) == [10])
chi = {nm: (legendre(ev(m, 3, 7), 7), legendre(ev(m, 10, 23), 23)) for nm, m in MU.items()}
check(f"version A: (chi_(7,3), chi_(23,10)) on 1, 9-theta, -15-theta, -12-theta = {list(chi.values())}:"
      " four distinct values, all nonzero",
      list(chi.values()) == [(1, 1), (-1, 1), (1, -1), (-1, -1)])
print("    version A: if mu = mu' xi^2 with xi in K*, then eta = mu' xi has eta^2 = mu mu' in R, so eta is integral")
print("    and lies in R (Section 3); then chi(mu)chi(mu') = leg(eta(r)^2) = +1. Distinct values: distinct classes.")
check("completely split primes: l^3 - 4l + 6 has the three distinct roots 18, 19, 25 mod 31 and 27, 45, 50 mod 61",
      roots_mod(31) == [18, 19, 25] and roots_mod(61) == [27, 45, 50])
check("each of the four mu is nonzero at all three roots mod 31 and mod 61",
      all(ev(m, r, q) != 0 for m in MU.values() for q in (31, 61) for r in roots_mod(q)))
chiB = {nm: (legendre(ev(m, 45, 61), 61), legendre(ev(m, 19, 31), 31)) for nm, m in MU.items()}
check(f"version B: (chi_(61,45), chi_(31,19)) on 1, 9-theta, -15-theta, -12-theta = {list(chiB.values())}:"
      " four distinct values", list(chiB.values()) == [(1, 1), (-1, 1), (1, -1), (-1, -1)])
print("    version B (no integral closedness): if mu = mu' (rho/n)^2 with n >= 1 least, and q | n, then")
print("    mu' rho^2 = 0 mod q gives rho(r_i) = 0 at the three distinct roots, so rho = 0 mod q (a polynomial of degree")
print("    <= 2 over F_q with three roots is 0), against the choice of n; so q does not divide n, and")
print("    n^2 mu(r) = mu'(r) rho(r)^2 gives chi(mu) = chi(mu').")
chi8 = {nm: (legendre(ev(m, 18, 31), 31), legendre(ev(m, 25, 31), 31), legendre(ev(m, 45, 61), 61))
        for nm, m in MU0.items()}
check(f"remark: the 8 admissible mu0 are nonzero at all roots mod 31, 61 and (chi_(31,18), chi_(31,25), chi_(61,45))"
      f" separates them: {chi8}",
      all(ev(m, r, q) != 0 for m in MU0.values() for q in (31, 61) for r in roots_mod(q))
      and len(set(chi8.values())) == 8)

# =====================================================================================================================
print("5.7 confirmation (not used): X - theta = mu0 S^2 at the 14 values of X")
DESCENT = {"1": [-8, 4, 40], "9-theta": [0, 9, 108, 252], "-15-theta": [-15, 12, 48, 4953],
           "-12-theta": [-12, 24, 420]}       # weierstrass_2descent.out.txt, CURVE 1611.a1, "grouped by class"
WITNESS = [(q, r) for q in sp.primerange(5, 200) for r in roots_mod(q)]
MU0_TO_MU = {m0: mname for mname, (m0, _) in IDENT.items()}
for Xv in KNOWN_X:
    al = x_minus_theta(Xv)
    hits = []
    for nm in T_NAMES + EPST_NAMES:
        quo = divides(MU0[nm], al)
        S = sqrt_in_R(quo) if quo is not None else None
        if S is not None and mul(MU0[nm], mul(S, S)) == al:
            hits.append((nm, S))
    if len(hits) != 1:
        check(f"X = {Xv}: exactly one admissible mu0 with X - theta = mu0 S^2 found (found {len(hits)})", False)
    ok = hits[0][0] in T_NAMES
    nm, S = hits[0]
    mname = MU0_TO_MU[nm]
    T = IDENT[mname][1]
    ok = ok and mul(al, mul(T, T)) == mul(MU[mname], mul(S, S)) and Xv in DESCENT[mname] and REP[Xv] == nm
    excluded = []
    for other in MU0:
        if other == nm:
            continue
        w = [(q, r) for q, r in WITNESS if ev(al, r, q) != 0 and ev(MU0[other], r, q) != 0
             and legendre(ev(al, r, q), q) != legendre(ev(MU0[other], r, q), q)]
        excluded.append(bool(w))
    check(f"X = {Xv:5d}: X - theta = ({nm}) * {S}^2; (X - theta) T^2 = ({mname}) S^2; agrees with the descent"
          f" class and with 5.3; the other 7 mu0 excluded by a character", ok and all(excluded))
print("    (A character witness (q, r) with alpha(r), mu0(r) != 0 mod q and leg(alpha(r)) != leg(mu0(r)) shows that")
print("    alpha is not mu0 S^2 with S in R, since S(r) != 0 would give leg(alpha(r)) = leg(mu0(r)).)")

# =====================================================================================================================
print("5.8 confirmation (not used): a direct search")
found = []
for Xv in range(-10**5, 10**5 + 1):
    v = Xv**3 - 144*Xv + 1296
    if v >= 0 and math.isqrt(v)**2 == v:
        found.append(Xv)
check(f"the integers X with |X| <= 10^5 and X^3 - 144X + 1296 a square are exactly {found}", found == KNOWN_X)
print("all checks of 07_classes.py passed")
