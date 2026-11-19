"""check_2adic_zero.py -- the markers of the orbit (0,1), c = 36 of 8604.a1 have a common 2-adic zero that is not an
integer zero (proof.md §4.5).  A remark outside the proof: nothing in Theorems A', A or B uses it.

Written by the round-2 global referee (review.md) in its scratch folder; moved here with this docstring, the asserts
and the final lines added by the author.  Independent of the repository code: the field, the basis, the units and the
seed are typed from proof.md (tables 2.1, 2.4 and 3.4).

Field: g(t) = t^4 + 18 t^2 - 72 t + 9, basis omega = 1, b, (3+b^2)/6, (9+9b+b^2+b^3)/12; eps1 = (1,-1,1,0),
eps2 = (1,0,-2,0), seed w = (3,-1,0,0) (coordinates in omega).

Method.
1. Lift the candidate exponents level by level: at level k the classes are those mod 2^k Z x 2^(k-1) Z, which lies in
   Lambda_(2^k) (eps1^(2^k) == eps2^(2^(k-1)) == 1 mod 2^k, checked for k = 2..9; for larger k it follows by squaring,
   since u == 1 mod 2^k with k >= 2 gives u^2 == 1 mod 2^(k+1)).  Keep the classes whose markers are == 0 mod 2^k.
2. Hensel certificate at a lifted point z* of each cluster.  With A = eps1^4 and B = eps2^2 (both == 1 mod 4O, as
   Lambda_4 = 4Z x 2Z), F(i,j) = markers(A^i B^j z*) = markers(exp(i log A + j log B) z*) is a power series in (i,j)
   whose degree-k coefficients have 2-adic valuation >= 2k - v(k!) >= k + 1 (log A == log B == 0 mod 4).  The
   multivariable Hensel lemma (Newton's method) gives a unique zero of F within 2^-v(det DF(0)) of 0 when
   v(F(0)) > 2 v(det DF(0)).  That lemma is p-adic analysis and is used only in this remark.
"""
from fractions import Fraction as Fr
import time
import sympy as sp

T0 = time.time()
b = sp.symbols('b')
g = b**4 + 18*b**2 - 72*b + 9
omega = [sp.Integer(1), b, (3 + b**2)/6, (9 + 9*b + b**2 + b**3)/12]
# change of basis: power basis coords of omega_j
P = sp.Matrix([[sp.Poly(sp.expand(o), b).coeff_monomial(b**k) for o in omega] for k in range(4)])  # columns = omega_j
Pinv = P.inv()


def to_poly(x):
    return sp.expand(sum(x[i]*omega[i] for i in range(4)))


def coords(e):
    r = sp.Poly(sp.rem(sp.expand(e), g, b), b)
    pc = sp.Matrix([r.coeff_monomial(b**k) for k in range(4)])
    c = Pinv * pc
    assert all(ci.is_integer for ci in c), c
    return [int(ci) for ci in c]


def mulmat(x):
    xp = to_poly(x)
    cols = [coords(xp * omega[j]) for j in range(4)]
    return [[cols[j][i] for j in range(4)] for i in range(4)]


def mm(A, B, mod=None):
    C = [[sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
    if mod:
        C = [[c % mod for c in row] for row in C]
    return C


def mv(A, v, mod=None):
    r = [sum(A[i][k]*v[k] for k in range(4)) for i in range(4)]
    return [x % mod for x in r] if mod else r


def mpow(A, e, mod):
    R = [[int(i == j) for j in range(4)] for i in range(4)]
    X = [[a % mod for a in row] for row in A]
    while e:
        if e & 1:
            R = mm(R, X, mod)
        X = mm(X, X, mod)
        e >>= 1
    return R


def v2(x):
    x = int(x)
    if x == 0:
        return 10**9
    k = 0
    while x % 2 == 0:
        x //= 2
        k += 1
    return k


print('check_2adic_zero.py (8604.a1, class (0,1), c = 36): a common 2-adic zero of the markers that is not an integer zero')
e1, e2, w = [1, -1, 1, 0], [1, 0, -2, 0], [3, -1, 0, 0]
M1, M2 = mulmat(e1), mulmat(e2)
assert sp.Matrix(M1).det() == 1 and sp.Matrix(M2).det() == 1 and sp.Matrix(mulmat(w)).det() == 36
print('det M1 = det M2 = 1 (units of norm 1); N(w) = 36')

MOD = 2**200


def elt(a, bb, mod=MOD):
    assert a >= 0 and bb >= 0
    return mv(mm(mpow(M1, a, mod), mpow(M2, bb, mod), mod), w, mod)


# Lambda_{2^k} contains 2^k Z x 2^(k-1) Z: eps1^(2^k) = 1 and eps2^(2^(k-1)) = 1 mod 2^k for k = 2..9
I = [[int(i == j) for j in range(4)] for i in range(4)]
for k in range(2, 10):
    assert mpow(M1, 2**k, 2**k) == I and mpow(M2, 2**(k-1), 2**k) == I
print('Lambda_(2^k) contains 2^k Z x 2^(k-1) Z for k = 2..9: checked (and for all k >= 2 by squaring)')

# lift candidates from mod (64, 32) (markers = 0 mod 2^6), as in verify_exclusion
cands = [(a, bb) for a in range(64) for bb in range(32) if all(x % 64 == 0 for x in elt(a, bb, 64)[2:])]
print('candidate classes mod (64,32):', cands)
K = 120
cl = cands
for k in range(6, K):
    q2 = 2**(k+1)
    cl = [(a + 2**k*i, bb + 2**(k-1)*j) for (a, bb) in cl for i in (0, 1) for j in (0, 1)
          if all(x % q2 == 0 for x in elt(a + 2**k*i, bb + 2**(k-1)*j, q2)[2:])]
clusters = {}
for (a, bb) in cl:
    clusters.setdefault((a % 64, bb % 32), []).append((a, bb))
print(f'after lifting to 2^{K}: classes per cluster (exponents mod (64,32)):', {k2: len(v) for k2, v in clusters.items()})
assert sorted(clusters) == [(0, 0), (46, 6)]
for k2, v in clusters.items():
    a, bb = v[0]
    print('   cluster', k2, 'exponent mod 2^24 =', (a % 2**24, bb % 2**24))

# Hensel certificate at each cluster
A = mpow(M1, 4, MOD)
B = mpow(M2, 2, MOD)
assert all((A[i][j] - I[i][j]) % 4 == 0 and (B[i][j] - I[i][j]) % 4 == 0 for i in range(4) for j in range(4))


def logm(Mx):
    # log(I + X), X = Mx - I == 0 mod 4; the terms X^k/k with 2k - v2(k) > 260 are == 0 mod 2^200 and are dropped
    X = [[(Mx[i][j] - I[i][j]) for j in range(4)] for i in range(4)]
    S = [[Fr(0)]*4 for _ in range(4)]
    Pk = [row[:] for row in I]
    k = 0
    while True:
        k += 1
        Pk = mm(Pk, X)
        Pk = [[x % 2**400 for x in row] for row in Pk]  # 2k <= 400 here, so Pk stays == 0 mod 2^(2k)
        if 2*k - (k.bit_length()) > 260:
            break
        for i in range(4):
            for j in range(4):
                S[i][j] += Fr((-1)**(k+1) * Pk[i][j], k)
    return S


def frac_mod(x, mod):
    # x a Fraction with odd denominator
    assert x.denominator % 2 == 1
    return (x.numerator * pow(x.denominator, -1, mod)) % mod


LA, LB = logm(A), logm(B)
LAm = [[frac_mod(x, MOD) for x in row] for row in LA]
LBm = [[frac_mod(x, MOD) for x in row] for row in LB]
assert all(x % 4 == 0 for row in LAm + LBm for x in row)
print('log A == log B == 0 mod 4 (as 4x4 matrices, mod 2^200)')
ok = {}
for k2, v in clusters.items():
    a, bb = v[0]
    z = elt(a, bb, MOD)
    F0 = [z[2], z[3]]
    c1 = mv(LAm, z, MOD)
    c2 = mv(LBm, z, MOD)
    J = (c1[2]*c2[3] - c1[3]*c2[2]) % MOD
    vF = min(v2(F0[0] % MOD), v2(F0[1] % MOD))
    vJ = v2(J)
    ok[k2] = min(vF, 200) > 2*vJ
    print(f'   cluster {k2}: v(F(0)) >= {min(vF, 200)}, v(det DF(0)) = {vJ}; Hensel condition v(F) > 2 v(det DF): {ok[k2]}')
assert ok[(46, 6)]
print()
print('=> by the multivariable Hensel lemma there is a common 2-adic zero of the markers at an exponent alpha in Z_2^2 with')
print('   alpha == (46, 6) mod (64, 32).  The only integer zero of this orbit is the exponent (0, 0) of (1,3) (Theorem A\'),')
print('   so alpha is not an integer.  alpha lies in the Lambda_4-class (2, 0) + Lambda_4, so its class mod Lambda_(2^k) is a')
print('   candidate for every k, holds no integer zero, and lies in no Lambda_(2^j)-class of (0, 0): no sieve with moduli')
print('   2^k alone can finish this orbit, for any k.')
print(f'[done in {time.time() - T0:.1f} s]')
