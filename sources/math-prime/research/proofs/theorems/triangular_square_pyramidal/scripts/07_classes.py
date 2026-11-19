"""07_classes.py: the class of X - theta for an integral point of E: y^2 = x^3 - 9x + 81 (hypothesis_H.md, Section 5).

5.1 the two factors alpha = X - theta, beta = X^2 + theta X + theta^2 - 9: alpha beta = X^3 - 9X + 81 = Y^2,
    N(alpha) = X^3 - 9X + 81, and beta - (X + 2 theta) alpha = f'(theta) = 9 P239 = (P3A P3B P3C)^2 P239;
5.2 the set S = {P3A, P3B, P3C, P239}: prime elements (prime norms), pairwise not associate;
5.3 the norm condition and the eight candidates mu0 = EPS^l P3A^na P3B^nb P3C^nc with na + nb + nc even;
5.4 the step modulo 4: for the four candidates with l = 1 no S in R gives mu0 S^2 = (X, -3, 0) modulo 4; for the
    other four such S exist;
5.5 the identification with the classes of (H): 3 - theta = P3A P3C P3B^2, -3 - theta = P3A P3B P3C^2,
    -theta = P3B P3C P3A^2; the four classes 1, 3 - theta, -3 - theta, -theta are pairwise distinct modulo squares
    (quadratic characters at roots of f modulo 7 and 13);
5.6 confirmation (not used): the 12 integral X of E (24 points) each have X - theta = mu0 S^2 with S in R explicit,
    in the class that the 2-descent assigns.
"""
import itertools
import time

import sympy as sp

from hring import (mul, norm, adj, power, divides, x_minus_theta, cofactor, ONE, THETA, EPS, EPS_INV, P3A, P3B,
                   P3C, P239)
from hlib import check

T0 = time.time()
X, Y = sp.symbols('X Y')
print("5.1 two factors")
al, be = x_minus_theta(X), cofactor(X)
check("alpha = X - theta = (X, -3, 0), beta = (X^2 - 9, 3X, 9); alpha beta = (X^3 - 9X + 81, 0, 0) and"
      " N(alpha) = X^3 - 9X + 81",
      [sp.expand(t) for t in mul(al, be)] == [sp.expand(X**3 - 9*X + 81), 0, 0]
      and sp.expand(norm(al) - (X**3 - 9*X + 81)) == 0)
rest = tuple(sp.expand(be[k] - mul((X, 6, 0), al)[k]) for k in range(3))
F3 = mul(mul(P3A, P3B), P3C)
check("beta - (X + 2 theta) alpha = (-9, 0, 27) = 9 P239, and 9 = (P3A P3B P3C)^2 (P3A P3B P3C = -3)",
      rest == (-9, 0, 27) and F3 == (-3, 0, 0) and mul(mul(F3, F3), P239) == (-9, 0, 27))

print("5.2 the prime elements")
SS = {'P3A': P3A, 'P3B': P3B, 'P3C': P3C, 'P239': P239}
check("N(P3A) = N(P3B) = N(P3C) = -3 and N(P239) = 239, all prime up to sign",
      [norm(v) for v in SS.values()] == [-3, -3, -3, 239] and sp.isprime(3) and sp.isprime(239))
check("pairwise not associate: for two elements of the same norm, neither divides the other",
      all(divides(SS[u], SS[v]) is None for u in SS for v in SS if u != v and norm(SS[u]) == norm(SS[v])))
print("    (elements of different norms are not associate, since associates have equal |N|)")

print("5.3 the norm condition and the eight candidates")
print("    X - theta = u' P3A^ea P3B^eb P3C^ec P239^e239 R0^2 (hypothesis_H.md 5.4) gives, taking norms,")
print("    Y^2 = N(u') (-3)^(ea+eb+ec) 239^e239 N(R0)^2: so e239 is even, ea + eb + ec is even, and N(u') = 1.")
print("    By Section 4, u' = sigma EPS^k with N(u') = sigma; so sigma = 1. Reducing exponents mod 2:")
cands = []
for l, na, nb, nc in itertools.product((0, 1), repeat=4):
    if (na + nb + nc) % 2:
        continue
    mu0 = mul(mul(mul(power(EPS, l), power(P3A, na)), power(P3B, nb)), power(P3C, nc))
    cands.append(((l, na, nb, nc), mu0))
check(f"{len(cands)} candidates mu0 = EPS^l P3A^na P3B^nb P3C^nc with na + nb + nc even, all of norm 1 or 9",
      len(cands) == 8 and all(norm(m) in (1, 9) for _, m in cands))
for key, m in cands:
    print(f"    (l, na, nb, nc) = {key}: mu0 = {m}, N = {norm(m)}")

print("5.4 the step modulo 4")
print("    S has integer coordinates, so the coordinates of mu0 S^2 modulo 4 depend only on S modulo 4.")
print("    X - theta = (X, -3, 0) needs its l-coordinate == 1 and its l^2-coordinate == 0 (mod 4).")
counts = {}
for key, m in cands:
    hits = 0
    for S in itertools.product(range(4), repeat=3):
        v = mul(m, mul(S, S))
        if (v[1] - 1) % 4 == 0 and v[2] % 4 == 0:
            hits += 1
    counts[key] = hits
check("for the four candidates with l = 1 (the coset EPS T) no residue S mod 4 works; for the four with l = 0 some do",
      all(counts[k] == 0 for k in counts if k[0] == 1) and all(counts[k] > 0 for k in counts if k[0] == 0))
for key in counts:
    print(f"    (l, na, nb, nc) = {key}: residues S mod 4 with the right shape: {counts[key]}")
S_ = sp.symbols('s1 s2 s3')
Z_ = sp.symbols('z1 z2 z3')
lhs = mul(tuple(S_[k] + 2*Z_[k] for k in range(3)), tuple(S_[k] + 2*Z_[k] for k in range(3)))
check("a reason readable by hand: (S + 2Z)^2 - S^2 = 4(SZ + Z^2) for every Z in R, so S^2 mod 4 depends only on S mod 2,"
      " and mu0 S^2 mod 4 (coordinatewise) takes at most 8 values, listed below",
      all(sp.Poly(sp.expand(lhs[k] - mul(S_, S_)[k]), *S_, *Z_).coeffs() and
          all(co % 4 == 0 for co in sp.Poly(sp.expand(lhs[k] - mul(S_, S_)[k]), *S_, *Z_).coeffs()) for k in range(3)))
FORCED = {}
for key, m in cands:
    vals = sorted({tuple(t % 4 for t in mul(m, mul(S, S))) for S in itertools.product(range(2), repeat=3)})
    good = [v for v in vals if v[1] == 1 and v[2] == 0]
    FORCED[key] = sorted({v[0] for v in good})
    print(f"    {key}: mu0 S^2 mod 4 in {vals}; of the shape (*, 1, 0): {good}")
print("    So the step also fixes X mod 4 in each surviving class:",
      "; ".join(f"{k}: X = {' or '.join(str(r) for r in FORCED[k])} mod 4" for k in FORCED if k[0] == 0))

print("5.5 identification with (H) and distinctness")
check("3 - theta = (3, -3, 0) = P3A P3C P3B^2; -3 - theta = (-3, -3, 0) = P3A P3B P3C^2; -theta = (0, -3, 0) ="
      " P3B P3C P3A^2",
      mul(mul(P3A, P3C), mul(P3B, P3B)) == (3, -3, 0) and mul(mul(P3A, P3B), mul(P3C, P3C)) == (-3, -3, 0)
      and mul(mul(P3B, P3C), mul(P3A, P3A)) == (0, -3, 0))
print("    So mu0 in {1, P3A P3C, P3A P3B, P3B P3C} is mu tau^-2 for mu in {1, 3 - theta, -3 - theta, -theta}, with")
print("    tau = 1, P3B, P3C, P3A: X - theta = mu0 S^2 = mu (S/tau)^2, which is (H).")
MU = {'1': (1, 0, 0), '3-th': (3, -3, 0), '-3-th': (-3, -3, 0), '-th': (0, -3, 0)}


def chi(u, p, r):
    """Legendre symbol of u(r) = a + b r + c r^2 modulo p (r a root of f mod p)."""
    v = (u[0] + u[1]*r + u[2]*r*r) % p
    return 0 if v == 0 else (1 if pow(v, (p - 1)//2, p) == 1 else -1)


roots = {p: [r for r in range(p) if (r**3 - r + 3) % p == 0] for p in (7, 13)}
print("    roots of f modulo 7 and 13:", roots)
sig = {k: tuple(chi(v, p, r) for p in (7, 13) for r in roots[p]) for k, v in MU.items()}
print("    quadratic characters (one per root):", sig)
check("the character values are all +-1 and the vectors of the four mu are pairwise different; a character that is"
      " nonzero on mu and mu' takes the same value on both when mu/mu' is a square in K (hypothesis_H.md 5.8), so the"
      " four classes are pairwise distinct", len(set(sig.values())) == 4 and all(0 not in s for s in sig.values()))

print("5.6 confirmation (not used by the proof): the integral points")
XS = [-5, -3, 0, 3, 7, 9, 24, 33, 39, 513, 1099, 5112]
check("these 12 X make X^3 - 9X + 81 a square (24 integral points)",
      all(sp.sqrt(x**3 - 9*x + 81).is_integer for x in XS))


def sqrt_in_R(t):
    """An integer triple S with S^2 = t, if one exists (search guided by the real place, verified exactly)."""
    import numpy as np
    rts = np.roots([1, 0, -1, 3])
    l1 = [r.real for r in rts if abs(r.imag) < 1e-12][0]
    z = [r for r in rts if r.imag > 1e-12][0]
    v1 = t[0] + t[1]*l1 + t[2]*l1**2
    v2 = t[0] + t[1]*z + t[2]*z**2
    if v1 < 0:
        return None
    for s1 in (1, -1):
        for s2 in (1, -1):
            a1, a2 = s1*v1**0.5, s2*complex(v2)**0.5
            Mm = np.array([[1, l1, l1**2], [1, z.real, (z*z).real], [0, z.imag, (z*z).imag]])
            rhs = np.array([a1, a2.real, a2.imag])
            sol = np.linalg.solve(Mm, rhs)
            S = tuple(int(round(x)) for x in sol)
            if mul(S, S) == tuple(t):
                return S
    return None


rep = {}
for x in XS:
    found = []
    for key, m in cands:
        if key[0] == 1:
            continue
        d = divides(m, x_minus_theta(x))          # (X - theta)/mu0, if integral
        if d is not None:
            S = sqrt_in_R(d)
            if S is not None and mul(m, mul(S, S)) == x_minus_theta(x):
                found.append((key, S))
    rep[x] = found
check("for every one of the 12 X exactly one mu0 of the four of (H) has X - theta = mu0 S^2 with S in R (S is unique"
      " up to sign)", all(len(rep[x]) == 1 for x in XS))
check("each X lies in the residue class mod 4 that the step modulo 4 (5.4 above) forces for its mu0",
      all(x % 4 in FORCED[rep[x][0][0]] for x in XS))
# the 2-descent's classes (proof.md section 5): trivial -5, 7, 1099; (1,0): 3, 39; (0,1): -3, 24, 33;
# (1,1): 0, 9, 513, 5112; with mu = 1, 3 - theta, -3 - theta, -theta, that is mu0 = 1, P3A P3C, P3A P3B, P3B P3C
EXPECT = {(0, 0, 0): {-5, 7, 1099}, (1, 0, 1): {3, 39}, (1, 1, 0): {-3, 24, 33}, (0, 1, 1): {0, 9, 513, 5112}}
check("each X is in the class that the 2-descent assigns (proof.md section 5: trivial, (1,0), (0,1), (1,1) <->"
      " mu0 = 1, P3A P3C, P3A P3B, P3B P3C)",
      all(x in EXPECT[rep[x][0][0][1:]] for x in XS))
for x in XS:
    key, S = rep[x][0]
    print(f"    X = {x:5d}: mu0 = P3A^{key[1]} P3B^{key[2]} P3C^{key[3]}, S = {S}")
print(f"all checks of 07_classes.py passed ({time.time() - T0:.1f} s)")
