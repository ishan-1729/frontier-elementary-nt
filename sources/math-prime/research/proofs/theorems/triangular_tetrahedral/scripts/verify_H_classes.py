#!/usr/bin/env python3
"""verify_H_classes.py: an independent check of Section 5 of hypothesis_H.md (Theorem H, the class of X - theta).

Written 2026-10-15 by Claude (Opus 5.5) as an independent verifier. Run from this folder:
    python verify_H_classes.py > verify_H_classes.out.txt 2>&1

Independence. Nothing is imported from cubic.py, core.py, lifting2.py or any 0x_*.py script, and no code is copied
from them; 07_classes.py and its output were read only to learn WHAT is claimed. Where another method exists, this
script uses it:
  * R = Z[l], l^3 = 4l - 6, is built as Z[t]/(f), f = t^3 - 4t + 6: a product is a product of coefficient lists
    followed by long division by the monic f (no hard-coded product formula);
  * the norm is the determinant (Leibniz expansion) of the matrix of multiplication, cross-checked against the
    resultant Res_t(f(t), u(t)) computed by sympy (no explicit cubic form);
  * division in R is Cramer's rule over Q (fractions.Fraction) followed by an integrality test (no adjugate);
  * square roots in R are found q-adically: the three roots of f modulo a prime q at which f splits completely are
    lifted by Newton's method to q^48, square roots of the values are lifted the same way, all 8 sign patterns are
    interpolated (Lagrange modulo q^48) and reduced to symmetric residues, and a candidate is accepted only if its
    square is checked exactly (07_classes.py uses mpmath complex embeddings instead);
  * non-squares in K are certified by a quadratic character at a root of f modulo a completely split prime other
    than 31 and 61 (the primes of 07); the argument, which needs no integral closedness, is printed in part (d);
  * coprimality (an extra in part (e)) is the gcd of the 3x3 minors of [M(alpha') | M(beta')], i.e. the index of
    alpha'R + beta'R in R = Z^3 (07 uses a Bezout column reduction);
  * the search for integral points sieves by quadratic residues modulo 64, 63, 65 and 11 and then tests with an
    exact integer square root.
No floating-point number is used anywhere: every decision is exact (int, Fraction, sympy). The script prints
PASS/FAIL lines, stops at the first FAIL and ends with a SUMMARY line.

Parts: (0) the ring; (a) 5.2-5.3 identities; (b) 5.5-5.6 the 128 parity vectors; (c) 5.7 the step modulo 4 and 8;
(d) 5.8 identities and characters; (e) Theorem H on every integral point with |X| <= 2*10^6; (f) negative controls.
"""
import itertools
import math
import random
import sys
from fractions import Fraction

import sympy as sp

# =====================================================================================================================
# bookkeeping
N_PASS = 0


def check(label, cond):
    """Print PASS or FAIL; stop at the first FAIL."""
    global N_PASS
    if cond:
        N_PASS += 1
        print('PASS  ' + label, flush=True)
        return
    print('FAIL  ' + label, flush=True)
    print(f'SUMMARY: FAIL (stopped at the first failure, after {N_PASS} passed checks)', flush=True)
    sys.exit(1)


def info(text=''):
    print('    ' + text, flush=True)


# =====================================================================================================================
# the ring R = Z[t]/(f), f = t^3 - 4t + 6; l is the class of t; a triple (a, b, c) is a + b l + c l^2
F_COEFFS = (6, -4, 0, 1)        # f(t) = 6 - 4t + 0t^2 + t^3, lowest degree first, monic
WRONG_F = (-6, -4, 0, 1)        # t^3 - 4t - 6: used only by a planted error in part (f)


def poly_mul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] = out[i + j] + x * y
    return out


def poly_rem(p, fc):
    """Remainder of p (coefficient list, lowest degree first) on division by the monic polynomial fc."""
    n = len(fc) - 1
    if fc[n] != 1:
        raise ValueError('the divisor must be monic')
    p = list(p) + [0] * max(0, n - len(p))
    for d in range(len(p) - 1, n - 1, -1):
        lead = p[d]
        for k in range(n + 1):
            p[d - n + k] = p[d - n + k] - lead * fc[k]
    return tuple(p[:n])


def mul(u, v, fc=F_COEFFS):
    return poly_rem(poly_mul(u, v), fc)


def add(u, v):
    return tuple(x + y for x, y in zip(u, v))


def sub(u, v):
    return tuple(x - y for x, y in zip(u, v))


def neg(u):
    return tuple(-x for x in u)


def smul(k, u):
    return tuple(k * x for x in u)


ONE, LAM = (1, 0, 0), (0, 1, 0)


def prod(*factors, fc=F_COEFFS):
    out = ONE
    for x in factors:
        out = mul(out, x, fc)
    return out


def pw(u, n, fc=F_COEFFS):
    """u^n for an integer n >= 0, by binary powering."""
    out, base = ONE, u
    while n:
        if n & 1:
            out = mul(out, base, fc)
        base = mul(base, base, fc)
        n >>= 1
    return out


def mult_matrix(u, fc=F_COEFFS):
    """Entry (i, j) is coordinate i of u * l^j: the matrix of multiplication by u in the basis 1, l, l^2."""
    cols = [mul(u, e, fc) for e in (ONE, LAM, mul(LAM, LAM, fc))]
    return [[cols[j][i] for j in range(3)] for i in range(3)]


def det3(m):
    """Determinant by the Leibniz expansion over the six permutations."""
    total = 0
    for perm in itertools.permutations(range(3)):
        inversions = sum(1 for i in range(3) for j in range(i + 1, 3) if perm[i] > perm[j])
        term = -1 if inversions % 2 else 1
        for i in range(3):
            term = term * m[i][perm[i]]
        total = total + term
    return total


def norm(u, fc=F_COEFFS):
    return det3(mult_matrix(u, fc))


def res_norm(u):
    """The norm as the resultant Res_t(f(t), u(t)) (f monic), computed by sympy."""
    t = sp.Symbol('t')
    return sp.resultant(t**3 - 4*t + 6, u[0] + u[1]*t + u[2]*t**2, t)


def solve(u, d):
    """The rational triple x with d*x = u (N(d) != 0), by Cramer's rule on M(d) x = u."""
    m = mult_matrix(d)
    D = det3(m)
    if D == 0:
        raise ZeroDivisionError('norm 0')
    out = []
    for j in range(3):
        mj = [[u[i] if k == j else m[i][k] for k in range(3)] for i in range(3)]
        out.append(Fraction(det3(mj), D))
    return tuple(out)


def div_R(u, d):
    """u/d if it lies in R, else None."""
    x = solve(u, d)
    return tuple(int(c) for c in x) if all(c.denominator == 1 for c in x) else None


def exponent(p, u):
    """(e, v) with u = p^e v and p not dividing v in R (u != 0)."""
    e = 0
    while True:
        q = div_R(u, p)
        if q is None:
            return e, u
        e, u = e + 1, q


def minors_gcd(u, v):
    """gcd of the 3x3 minors of the 3x6 matrix [M(u) | M(v)]: the index of uR + vR in R = Z^3 (0 if infinite)."""
    cols = [list(col) for col in zip(*mult_matrix(u))] + [list(col) for col in zip(*mult_matrix(v))]
    g = 0
    for trio in itertools.combinations(range(6), 3):
        g = math.gcd(g, det3([[cols[j][i] for j in trio] for i in range(3)]))
    return g


def same(u, v):
    return all(sp.expand(x - y) == 0 for x, y in zip(u, v))


def coeff_gcd(expr, gens):
    """gcd of the integer coefficients of a polynomial (0 for the zero polynomial)."""
    g = 0
    for co in sp.Poly(sp.expand(expr), *gens).coeffs():
        g = math.gcd(g, int(co))
    return g


# =====================================================================================================================
# quadratic characters and q-adic square roots
_SQUARES = {}


def legendre(x, q):
    """0 if q | x; +1 if x is a nonzero square mod q; -1 otherwise (by definition, from the list of squares)."""
    if q not in _SQUARES:
        _SQUARES[q] = {(s * s) % q for s in range(1, q)}
    x %= q
    if x == 0:
        return 0
    return 1 if x in _SQUARES[q] else -1


def ev(u, r, m):
    return (u[0] + u[1] * r + u[2] * r * r) % m


def roots_mod(q):
    return [r for r in range(q) if (r**3 - 4*r + 6) % q == 0]


SPLIT = []
for _q in sp.primerange(5, 3000):
    _rs = roots_mod(int(_q))
    if len(_rs) == 3:
        SPLIT.append((int(_q), _rs))
SPLIT_IND = [(q, rs) for q, rs in SPLIT if q not in (31, 61)]     # 31 and 61 are the split primes used by 07
HK = 48                                                           # q-adic precision: q^48 > 10^97 for q >= 107
_LIFTED = {}


def lifted_roots(q):
    """The three roots of f in Z/q^HK lifting the roots mod q (Newton's method; the roots are simple)."""
    if q not in _LIFTED:
        out = []
        for r in dict(SPLIT)[q]:
            x, prec = r, 1
            while prec < HK:
                prec = min(2 * prec, HK)
                mod = q ** prec
                x = (x - (x**3 - 4*x + 6) * pow(3*x*x - 4, -1, mod)) % mod
            if (x**3 - 4*x + 6) % q**HK:
                raise RuntimeError('root lifting failed')
            out.append(x)
        _LIFTED[q] = out
    return _LIFTED[q]


def sqrt_mod(g, q):
    """x with x^2 = g mod q^HK, for g a nonzero square mod q (q odd), by Newton's method from a root mod q."""
    g0 = g % q
    x = next(s for s in range(1, q) if (s * s) % q == g0)
    prec = 1
    while prec < HK:
        prec = min(2 * prec, HK)
        mod = q ** prec
        x = (x - (x * x - g) * pow(2 * x, -1, mod)) % mod
    if (x * x - g) % q**HK:
        raise RuntimeError('square-root lifting failed')
    return x


def interpolate(points, values, mod):
    """(s0, s1, s2) with s0 + s1 p + s2 p^2 = value at each of the three points, modulo mod (Lagrange)."""
    c = [0, 0, 0]
    for i in range(3):
        o1, o2 = [points[j] for j in range(3) if j != i]
        w = values[i] * pow((points[i] - o1) * (points[i] - o2), -1, mod) % mod
        c[0] += w * o1 * o2
        c[1] -= w * (o1 + o2)
        c[2] += w
    return [x % mod for x in c]


def symmetric(x, mod):
    x %= mod
    return x - mod if x > mod // 2 else x


def sqrt_in_R(g):
    """('square', S) with S*S == g checked exactly; ('nonsquare', certificate) with g not a square in K*;
    ('unknown', None) if neither is found (every caller treats that as a failure).
    Certificates: ('norm', N(g)) when N(g) is not a positive square; (q, r) when g(r) is a non-residue mod q at a root
    r of f modulo a completely split prime q (valid in K*, see part (d))."""
    g = tuple(g)
    n = norm(g)
    if n == 0:
        raise ValueError('zero')
    if n < 0 or math.isqrt(n) ** 2 != n:
        return 'nonsquare', ('norm', n)
    hensel_done = False
    for q, rs in SPLIT_IND:
        vals = [ev(g, r, q) for r in rs]
        if 0 in vals:
            continue
        for r, v in zip(rs, vals):
            if legendre(v, q) == -1:
                return 'nonsquare', (q, r)
        if not hensel_done:
            hensel_done = True
            mod = q ** HK
            rho = lifted_roots(q)
            xs = [sqrt_mod(ev(g, x, mod), q) for x in rho]
            for signs in itertools.product((1, -1), repeat=3):
                S = tuple(symmetric(c, mod) for c in interpolate(rho, [s * x for s, x in zip(signs, xs)], mod))
                if mul(S, S) == g:
                    return 'square', S
    return 'unknown', None


def classify(alpha, targets):
    """For each (name, mu) in targets decide whether alpha*mu is a square in K (i.e. alpha in mu K*^2)."""
    squares, nonsquares, unknown = [], [], []
    for name, mu in targets:
        verdict, data = sqrt_in_R(mul(alpha, mu))
        if verdict == 'square':
            squares.append((name, data))
        elif verdict == 'nonsquare':
            nonsquares.append((name, data))
        else:
            unknown.append(name)
    return squares, nonsquares, unknown


# =====================================================================================================================
# the elements of hypothesis_H.md (1.5) and Section 5, entered from the text
EPS = (-19, -5, 1)          # l^2 - 5l - 19
PI2 = (-4, 1, 1)            # l^2 + l - 4
P3A = (-1, 1, 0)            # l - 1
P3B = (5, 2, 0)             # 2l + 5
P3C = (3, -3, 1)            # l^2 - 3l + 3
P179 = (107, -114, 45)      # 45l^2 - 114l + 107
THETA = (0, 6, 0)           # theta = 6l
GENS = [('EPS', EPS), ('PI2', PI2), ('P3A', P3A), ('P3B', P3B), ('P3C', P3C), ('P179', P179)]
SPRIMES = GENS[1:]


def signed(sign, *factors):
    v = prod(*factors)
    return v if sign > 0 else neg(v)


LISTED = {      # the eight mu0 of 5.6, built from their names in the text
    '1': ONE,
    'EPS P3A P3B': signed(1, EPS, P3A, P3B),
    '-P3A P3C': signed(-1, P3A, P3C),
    '-EPS P3B P3C': signed(-1, EPS, P3B, P3C),
    'EPS': EPS,
    'P3A P3B': signed(1, P3A, P3B),
    '-EPS P3A P3C': signed(-1, EPS, P3A, P3C),
    '-P3B P3C': signed(-1, P3B, P3C),
}
T_NAMES = ['1', 'EPS P3A P3B', '-P3A P3C', '-EPS P3B P3C']        # the four mu0 of (5.1)
EPST_NAMES = ['EPS', 'P3A P3B', '-EPS P3A P3C', '-P3B P3C']        # the coset EPS*T
MU4 = [('1', ONE), ('9-theta', (9, -6, 0)), ('-15-theta', (-15, -6, 0)), ('-12-theta', (-12, -6, 0))]
MU0_OF = {'1': '1', '9-theta': 'EPS P3A P3B', '-15-theta': '-P3A P3C', '-12-theta': '-EPS P3B P3C'}

# =====================================================================================================================
print('(0) the ring R = Z[t]/(t^3 - 4t + 6), built here independently')
a_, b_, c_, d_, e_, f_, r_ = sp.symbols('a b c d e f r')
U, V = (a_, b_, c_), (d_, e_, f_)
UV = mul(U, V)
_diff = sp.Poly(sp.expand((a_ + b_*r_ + c_*r_**2) * (d_ + e_*r_ + f_*r_**2) - (UV[0] + UV[1]*r_ + UV[2]*r_**2)), r_)
_quo, _rem = sp.div(_diff, sp.Poly(r_**3 - 4*r_ + 6, r_))
_quo_int = all(co.is_integer for co in sp.Poly(_quo.as_expr(), a_, b_, c_, d_, e_, f_, r_).coeffs())
check('u(t)v(t) - (u*v)(t) = f(t) q(t) with q of integer coefficients, for generic u = (a,b,c), v = (d,e,f): the '
      'product of triples is the product of Z[t]/(f); so R is a commutative ring with 1, and for a root r of f modulo a '
      'prime q the evaluation u -> u(r) mod q is a ring homomorphism R -> F_q', _rem.is_zero and _quo_int)
info(f'q(t) = {_quo.as_expr()}')
_MU, _MV, _MUV = (sp.Matrix(mult_matrix(x)) for x in (U, V, UV))
check('M(uv) = M(u) M(v) for generic u, v, so N = det M is multiplicative', (_MU * _MV - _MUV).expand() == sp.zeros(3, 3))
N_TEXT = a_**3 + 8*a_**2*c_ - 4*a_*b_**2 + 18*a_*b_*c_ + 16*a_*c_**2 - 6*b_**3 + 24*b_*c_**2 + 36*c_**3
M_TEXT = sp.Matrix([[a_, -6*c_, -6*b_], [b_, a_ + 4*c_, 4*b_ - 6*c_], [c_, b_, a_ + 4*c_]])
check('M(u) and det M(u) agree with the matrix and the cubic form N printed in hypothesis_H.md (1.2)',
      (_MU - M_TEXT).expand() == sp.zeros(3, 3) and sp.expand(norm(U) - N_TEXT) == 0)
check('det M(u) = Res_t(t^3 - 4t + 6, a + bt + ct^2) (sympy resultant), for generic a, b, c',
      sp.expand(norm(U) - res_norm(U)) == 0)

# =====================================================================================================================
print('(a) 5.2 and 5.3: the two factors alpha, beta and their common divisor')
X = sp.Symbol('X')
check('theta = 6l satisfies theta^3 - 144 theta + 1296 = 0 in R',
      add(sub(pw(THETA, 3), smul(144, THETA)), (1296, 0, 0)) == (0, 0, 0))
alpha_s = sub((X, 0, 0), THETA)
beta_s = add(add((X**2 - 144, 0, 0), smul(X, THETA)), mul(THETA, THETA))
check('alpha = X - theta = (X, -6, 0) and beta = X^2 + theta X + theta^2 - 144 = (X^2 - 144, 6X, 36)',
      same(alpha_s, (X, -6, 0)) and same(beta_s, (X**2 - 144, 6*X, 36)))
check('alpha*beta = (X^3 - 144X + 1296, 0, 0), as polynomials in X',
      same(mul(alpha_s, beta_s), (X**3 - 144*X + 1296, 0, 0)))
_pts = list(range(-6, 8))
check(f'alpha*beta, second method: exact evaluation at the {len(_pts)} integers X = -6..7 (each coordinate has degree'
      ' <= 3 in X, so 4 points determine it)',
      all(mul((x, -6, 0), (x*x - 144, 6*x, 36)) == (x**3 - 144*x + 1296, 0, 0) for x in _pts))
check('N(alpha) = det M(alpha) = X^3 - 144X + 1296', sp.expand(norm(alpha_s) - (X**3 - 144*X + 1296)) == 0)
check('N(alpha), second method: Res_t(t^3 - 4t + 6, X - 6t) = X^3 - 144X + 1296',
      sp.expand(res_norm(alpha_s) - (X**3 - 144*X + 1296)) == 0)
_Xp2t = add((X, 0, 0), smul(2, THETA))
check('X + 2 theta = (X, 12, 0) and beta - (X + 2 theta) alpha = (-144, 0, 108) = 36(3l^2 - 4), identically in X;'
      ' also by evaluation at X = -6..7',
      same(_Xp2t, (X, 12, 0)) and same(sub(beta_s, mul(_Xp2t, alpha_s)), (-144, 0, 108))
      and smul(36, (-4, 0, 3)) == (-144, 0, 108) and add(smul(3, mul(LAM, LAM)), (-4, 0, 0)) == (-4, 0, 3)
      and all(sub((x*x - 144, 6*x, 36), mul((x, 12, 0), (x, -6, 0))) == (-144, 0, 108) for x in _pts))
_int_roots = [d for dd in sp.divisors(1296) for d in (dd, -dd) if d**3 - 144*d + 1296 == 0]
check(f'X^3 - 144X + 1296 has no integer root (rational root test over the {2 * len(sp.divisors(1296))} divisors'
      ' +-d of 1296): so Y != 0 on E, and alpha, beta, (Y,0,0) are nonzero', _int_roots == [])
info('so for an integral point: alpha*beta = (Y^2, 0, 0) and N(alpha) = Y^2; a common divisor of alpha and beta divides'
     ' beta - (X + 2 theta) alpha = 36(3l^2 - 4)')
EPS_INV = div_R(ONE, EPS)
check(f'EPS^-1 by Cramer\'s rule is integral: EPS^-1 = {EPS_INV} = 95 - 101l + 40l^2 (text (1.5)), EPS*EPS^-1 = 1',
      EPS_INV == (95, -101, 40) and mul(EPS, EPS_INV) == ONE)
_expected_norms = {'EPS': 1, 'PI2': -2, 'P3A': -3, 'P3B': -3, 'P3C': 3, 'P179': 179}
for _nm, _g in GENS:
    check(f'N({_nm}) = N{_g} = {_expected_norms[_nm]} by determinant and by resultant',
          norm(_g) == _expected_norms[_nm] and res_norm(_g) == _expected_norms[_nm])
check('PI2^3 = -2 EPS, P3A P3B P3C = 3 and PI2^2 P179 = 3l^2 - 4 = (-4, 0, 3) (text (1.5))',
      pw(PI2, 3) == smul(-2, EPS) and prod(P3A, P3B, P3C) == (3, 0, 0) and mul(pw(PI2, 2), P179) == (-4, 0, 3))
check('36(3l^2 - 4) = EPS^-2 PI2^8 (P3A P3B P3C)^2 P179 (exact product)',
      prod(pw(EPS_INV, 2), pw(PI2, 8), pw(prod(P3A, P3B, P3C), 2), P179) == (-144, 0, 108))
check('N(36(3l^2 - 4)) = 36^3 * 716 = 2^8 3^6 179 = (-2)^8 ((-3)(-3)(3))^2 179',
      norm((-144, 0, 108)) == 36**3 * 716 == 2**8 * 3**6 * 179 == (-2)**8 * ((-3) * (-3) * 3)**2 * 179)
check('|N| of PI2, P3A, P3B, P3C, P179 is a prime number (2, 3, 3, 3, 179): prime elements by (3.3)',
      all(sp.isprime(abs(norm(_g))) for _, _g in SPRIMES))
_assoc_ok = True
for (_n1, _p1), (_n2, _p2) in itertools.combinations(SPRIMES, 2):
    if abs(norm(_p1)) == abs(norm(_p2)):
        _assoc_ok = _assoc_ok and div_R(_p1, _p2) is None and div_R(_p2, _p1) is None
    else:
        _assoc_ok = _assoc_ok and abs(norm(_p1)) != abs(norm(_p2))
check('the five S-primes are pairwise non-associate (different |N|, or for P3A, P3B, P3C: neither quotient is in R,'
      ' by Cramer)', _assoc_ok)

# =====================================================================================================================
print('(b) 5.5-5.6: the 2^7 parity vectors (sign, EPS, PI2, P3A, P3B, P3C, P179) and their norms')
_admissible, _formula_ok, _rule_ok = [], True, True
for _vec in itertools.product((0, 1), repeat=7):
    _s, _k, _e2, _a, _b, _c, _g = _vec
    _rep = ONE
    for (_nm, _gen), _ex in zip(GENS, _vec[1:]):
        if _ex:
            _rep = mul(_rep, _gen)
    if _s:
        _rep = neg(_rep)
    _n = norm(_rep)
    _formula_ok = _formula_ok and _n == (-1)**_s * (-2)**_e2 * (-3)**(_a + _b) * 3**_c * 179**_g
    _adm = _n > 0 and math.isqrt(_n)**2 == _n
    _rule_ok = _rule_ok and _adm == (_e2 == 0 and _g == 0 and (_a + _b + _c) % 2 == 0 and _s == (_a + _b) % 2)
    if _adm:
        _admissible.append((_vec, _rep, _n))
check('for all 128 vectors, det M(rep) = (-1)^s (-2)^e2 (-3)^(a+b) 3^c 179^g (the norm formula of 5.5)', _formula_ok)
check('exactly 8 of the 128 representatives have a norm that is a positive square', len(_admissible) == 8)
check('a representative is admissible exactly when e2 = g = 0, a+b+c is even and sign = (-1)^(a+b)', _rule_ok)
for _vec, _rep, _n in _admissible:
    _names = [k for k, v in LISTED.items() if v == _rep]
    info(f'{_vec}  rep = {str(_rep):17s} N = {_n}  = {(_names or ["?"])[0]}'
         f'  ({"T" if _names and _names[0] in T_NAMES else "EPS*T"})')
check('the 8 admissible representatives are exactly the 8 mu0 listed in 5.6: 1, EPS P3A P3B, -P3A P3C, -EPS P3B P3C,'
      ' EPS, P3A P3B, -EPS P3A P3C, -P3B P3C (compared as triples)',
      sorted(r for _, r, _ in _admissible) == sorted(LISTED.values()) and len(set(LISTED.values())) == 8)
check('EPS * mu0 for mu0 in T is the listed EPS*T element times 1 or EPS^2 (so the second four are EPS*T modulo'
      ' squares; the text\'s "products with EPS" is meant modulo EPS^2)',
      all(mul(EPS, LISTED[t]) in (LISTED[s], mul(LISTED[s], pw(EPS, 2))) for t, s in zip(T_NAMES, EPST_NAMES)))
info('why this is the right list: if X - theta = mu0 S^2 with S in R, then Y^2 = N(X - theta) = N(mu0) N(S)^2, so N(mu0)')
info('is a positive rational square, hence a positive integer square; reducing the exponents of 5.4 mod 2 changes the')
info('representative by a square (of EPS, of S-primes, of R0), which multiplies the norm by a positive square.')

# =====================================================================================================================
print('(c) 5.7: the step modulo 4 (and 8)')
info('target: X - theta = (X, -6, 0), so mu0 S^2 must have l-coordinate = -6 and l^2-coordinate = 0 (mod m)')


def hits(mu0, m):
    count = 0
    for S in itertools.product(range(m), repeat=3):
        v = mul(mu0, mul(S, S))
        if (v[1] + 6) % m == 0 and v[2] % m == 0:
            count += 1
    return count


H2 = {nm: hits(LISTED[nm], 2) for nm in LISTED}
H4 = {nm: hits(LISTED[nm], 4) for nm in LISTED}
H8 = {nm: hits(LISTED[nm], 8) for nm in LISTED}
for _nm in EPST_NAMES + T_NAMES:
    info(f'mu0 = {_nm:13s} {str(LISTED[_nm]):17s} residues S with the target: mod 2: {H2[_nm]:2d}/8,'
         f' mod 4: {H4[_nm]:3d}/64, mod 8: {H8[_nm]:3d}/512')
check('mod 4 (brute force over all 64 residues S): no residue reaches the target for any mu0 of EPS*T',
      all(H4[nm] == 0 for nm in EPST_NAMES))
check('mod 4: each of the four mu0 of T has residues reaching the target', all(H4[nm] > 0 for nm in T_NAMES))
check('mod 8 (brute force over all 512 residues S): counts 0, 0, 0, 0 for EPS*T and 32, 48, 48, 32 for'
      ' 1, EPS P3A P3B, -P3A P3C, -EPS P3B P3C, as stated in 5.7',
      [H8[nm] for nm in EPST_NAMES] == [0, 0, 0, 0] and [H8[nm] for nm in T_NAMES] == [32, 48, 48, 32])
_m = sp.symbols('m0:3')
_s3 = sp.symbols('s0:3')
_t3 = sp.symbols('t0:3')
_g9 = (*_m, *_s3, *_t3)
for _k in (4, 8):
    _Sk = add(_s3, smul(_k, _t3))
    _D = sub(mul(_m, mul(_Sk, _Sk)), mul(_m, mul(_s3, _s3)))
    _gs = [coeff_gcd(x, _g9) for x in _D]
    check(f'symbolic, generic integer triples mu, S, T: every coefficient of every coordinate of mu(S + {_k}T)^2 -'
          f' mu S^2 is divisible by {_k} (coefficient gcds {_gs}); so mu0 S^2 mod {_k} depends only on S mod {_k}'
          f' and the test mod {_k} is legitimate', all(g % _k == 0 for g in _gs))
check('symbolic: (S + 4T)^2 - S^2 = 8ST + 16T^2 in R (the reason: the structure constants are integers)',
      same(sub(mul(add(_s3, smul(4, _t3)), add(_s3, smul(4, _t3))), mul(_s3, _s3)),
           add(smul(8, mul(_s3, _t3)), smul(16, mul(_t3, _t3)))))
info('the hand argument of 5.7 for the coset EPS*T, with S = (x, y, z); target (l, l^2) = (2, 0) mod 4')
_x, _y, _z = sp.symbols('x y z')
_xp, _yp, _zp, _y0, _z0 = sp.symbols('xp yp zp y0 z0')
_CLAIMED_U = {'EPS': -5, 'P3A P3B': 3, '-EPS P3A P3C': 1, '-P3B P3C': 1}
_EXPECTED_VALUES = sorted([(0, 0), (2, 2), (2, 1), (0, 3)])
for _nm in EPST_NAMES:
    _mu0 = LISTED[_nm]
    _sq = [sp.expand(c) for c in mul(_mu0, mul((_x, _y, _z), (_x, _y, _z)))]
    _u = _mu0[1]
    _step1 = (_u == _CLAIMED_U[_nm] and _u % 2 == 1 and coeff_gcd(_sq[1] - _u * _x**2, (_x, _y, _z)) % 2 == 0)

    def _val(x0, y0, z0, mu0=_mu0):
        v = mul(mu0, mul((x0, y0, z0), (x0, y0, z0)))
        return v[1] % 4, v[2] % 4

    _step1b = all(_val(x0, y0, z0)[0] % 2 == 1 for x0 in (1, 3) for y0 in range(4) for z0 in range(4))
    _step2 = all(_val(x0, y0, z0) == _val(0, y0 % 2, z0 % 2) for x0 in (0, 2) for y0 in range(4) for z0 in range(4))
    _sub_full = {_x: 2*_xp, _y: _y0 + 2*_yp, _z: _z0 + 2*_zp}
    _sub_base = {_x: 0, _y: _y0, _z: _z0}
    _step2s = all(coeff_gcd(_sq[i].subs(_sub_full, simultaneous=True) - _sq[i].subs(_sub_base, simultaneous=True),
                            (_xp, _yp, _zp, _y0, _z0)) % 4 == 0 for i in (1, 2))
    _table = {(y0, z0): _val(0, y0, z0) for y0 in (0, 1) for z0 in (0, 1)}
    _step3 = sorted(_table.values()) == _EXPECTED_VALUES and (2, 0) not in _table.values()
    check(f'mu0 = {_nm}: l-coordinate of mu0 S^2 = ({_u}) x^2 mod 2 (all other coefficients even; text: {_CLAIMED_U[_nm]}),'
          f' so x odd gives an odd l-coordinate (also by brute force) and x is even; for x even the pair mod 4 depends'
          f' only on (y, z) mod 2 (brute force, and symbolically with x = 2x\', y = y0 + 2y\', z = z0 + 2z\');'
          f' table {_table}: the values are (0,0), (2,2), (2,1), (0,3) in some order, never (2, 0)',
          _step1 and _step1b and _step2 and _step2s and _step3)
_kk = sp.Symbol('k')
_n42 = sp.expand((4*_kk + 2)**3 - 144*(4*_kk + 2) + 1296)
_odd = sp.expand(_n42 / 8)
_odd_ok = (all(sp.Integer(co).is_integer for co in sp.Poly(_odd, _kk).coeffs())
           and coeff_gcd(_odd - 1, (_kk,)) % 2 == 0)
_pi2_ok = all(exponent(PI2, (x, -6, 0))[0] == 3 for x in range(-1002, 1003) if x % 4 == 2)
check(f'remark 5.9: for X = 4k + 2, X^3 - 144X + 1296 = 8({_odd}) with the bracket odd for every k, so it is not a'
      ' square and no integral point has X = 2 (mod 4); PI2 divides X - theta exactly 3 times for every X = 2 (mod 4)'
      ' with |X| <= 1002 (exact division)', _odd_ok and _pi2_ok)

# =====================================================================================================================
print('(d) 5.8: the identities with the classes of (H), and their distinctness by quadratic characters')
R1, R2, R3 = (-45, 48, -19), (5, 2, 0), (-8, 8, -3)
check('9 - theta = (EPS P3A P3B) * (-45 + 48l - 19l^2)^2 (exact)', mul(LISTED['EPS P3A P3B'], mul(R1, R1)) == (9, -6, 0))
check('-15 - theta = (-P3A P3C) * (5 + 2l)^2 (exact)', mul(LISTED['-P3A P3C'], mul(R2, R2)) == (-15, -6, 0))
check('-12 - theta = (-EPS P3B P3C) * (-8 + 8l - 3l^2)^2 (exact)', mul(LISTED['-EPS P3B P3C'], mul(R3, R3)) == (-12, -6, 0))
info(f'N(r1), N(r2), N(r3) = {norm(R1)}, {norm(R2)}, {norm(R3)} (nonzero), so mu0 = mu r^-2 in K and'
     ' X - theta = mu0 S^2 gives X - theta = mu (S/r)^2')
_v, _w = sqrt_in_R(prod((9, -6, 0), (-15, -6, 0), (-12, -6, 0)))
check(f'(9-theta)(-15-theta)(-12-theta) = W^2 with W = {_w} in R (q-adic root, exact square), and W = +-3 EPS r1 r2 r3,'
      ' so the four mu form a group modulo squares', _v == 'square' and _w in (smul(3, prod(EPS, R1, R2, R3)),
                                                                             neg(smul(3, prod(EPS, R1, R2, R3)))))
info('why a character separates classes (no integral closedness needed). Let q be a prime at which f has three')
info('distinct roots r1, r2, r3 mod q, and phi_i: R -> F_q the evaluation at r_i (a ring homomorphism, part (0)).')
info('Let mu, mu\' in R be nonzero at all three roots, and suppose mu = mu\' xi^2 with xi in K*. Let n >= 1 be least')
info('with rho := n xi in R; then n^2 mu = mu\' rho^2. If q | n, then phi_i(mu\') phi_i(rho)^2 = 0, so phi_i(rho) = 0')
info('for i = 1, 2, 3; a polynomial of degree <= 2 over F_q with three distinct roots is 0, so rho = q rho\' with')
info('rho\' in R and (n/q) xi = rho\', against the choice of n. So q does not divide n, and n^2 phi_i(mu) =')
info('phi_i(mu\') phi_i(rho)^2 with n^2 and phi_i(rho)^2 nonzero squares of F_q: phi_i(mu) is a square of F_q exactly')
info('when phi_i(mu\') is. Different character values therefore mean different classes of K*/K*^2. With mu\' = 1 the')
info('q | n step needs nothing of mu, so a non-residue phi_i(mu) certifies that mu is not a square in K (the')
info('certificates (q, r) of sqrt_in_R). Only the definition of a square mod q is used, not multiplicativity.')
_cands = [(q, r) for q, rs in SPLIT_IND[:40] if all(ev(m, rr, q) != 0 for _, m in MU4 for rr in rs) for r in rs]
_pair = None
for _pr in itertools.combinations(_cands, 2):
    if len({tuple(legendre(ev(m, r, q), q) for q, r in _pr) for _, m in MU4}) == 4:
        _pair = _pr
        break
_vals4 = {nm: tuple(legendre(ev(m, r, q), q) for q, r in _pair) for nm, m in MU4} if _pair else {}
check(f'characters at the split primes chosen here, (q, r) = {_pair} (q not 31, 61): all four mu are nonzero at the'
      f' three roots of each q, and the values {_vals4} are four distinct pairs: 1, 9-theta, -15-theta, -12-theta are'
      ' pairwise distinct modulo K*^2', _pair is not None and len(set(_vals4.values())) == 4
      and all(len(roots_mod(q)) == 3 for q, _ in _pair))
_cands8 = [(q, r) for q, rs in SPLIT_IND[:40] if all(ev(m, rr, q) != 0 for m in LISTED.values() for rr in rs) for r in rs]
_triple = None
for _tr in itertools.combinations(_cands8, 3):
    if len({tuple(legendre(ev(m, r, q), q) for q, r in _tr) for m in LISTED.values()}) == 8:
        _triple = _tr
        break
_vals8 = {nm: tuple(legendre(ev(m, r, q), q) for q, r in _triple) for nm, m in LISTED.items()} if _triple else {}
check(f'the 8 admissible mu0 are pairwise distinct modulo K*^2: characters at {_triple} give 8 distinct sign vectors'
      f' {_vals8}; so X - theta lies in exactly one of the 8 classes', _triple is not None and len(set(_vals8.values())) == 8)
_text_vals = [(legendre(ev(m, 3, 7), 7), legendre(ev(m, 10, 23), 23)) for _, m in MU4]
_mult_ok = all(legendre(x * y, q) == legendre(x, q) * legendre(y, q) for q in (7, 23) for x in range(1, q) for y in range(1, q))
check(f'the text\'s certificate: f has the single roots 3 mod 7 and 10 mod 23, and (chi_(7,3), chi_(23,10)) on'
      f' 1, 9-theta, -15-theta, -12-theta = {_text_vals} = (1,1), (-1,1), (1,-1), (-1,-1) as printed in 5.8; this'
      ' version needs (3.6) (integral closedness) and the multiplicativity of the Legendre symbol (checked by brute force'
      ' mod 7 and 23), since 7 and 23 are not completely split',
      roots_mod(7) == [3] and roots_mod(23) == [10] and _text_vals == [(1, 1), (-1, 1), (1, -1), (-1, -1)] and _mult_ok)

# =====================================================================================================================
print('(e) Theorem H on every integral point with |X| <= 2*10^6: explicit factorization X - theta = mu0 S^2')
BOUND = 2 * 10**6
QRES = {m: {(x * x) % m for x in range(m)} for m in (64, 63, 65, 11)}
FOUND, _survivors, _nonneg = [], 0, 0
for _X in range(-BOUND, BOUND + 1):
    _val3 = _X * _X * _X - 144 * _X + 1296
    if _val3 < 0:
        continue
    _nonneg += 1
    if (_val3 % 64 not in QRES[64] or _val3 % 63 not in QRES[63] or _val3 % 65 not in QRES[65]
            or _val3 % 11 not in QRES[11]):
        continue
    _survivors += 1
    _root = math.isqrt(_val3)
    if _root * _root == _val3:
        FOUND.append(_X)
info(f'{2 * BOUND + 1} integers X; {_nonneg} with X^3 - 144X + 1296 >= 0 (it is < 0 for X <= -16); {_survivors} pass the'
     f' quadratic-residue sieve mod 64, 63, 65, 11; the exact isqrt test leaves')
info(f'X in {FOUND}')
EXPECTED_X = [-15, -12, -8, 0, 4, 9, 12, 24, 40, 48, 108, 252, 420, 4953]
# grouping of the integral points by 2-descent class, from
# research/proofs/targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt lines 548-554 (CURVE 1611.a1);
# used only for a comparison, not for any decision
DESCENT = {'1': [-8, 4, 40], '9-theta': [0, 9, 108, 252], '-15-theta': [-15, 12, 48, 4953], '-12-theta': [-12, 24, 420]}
UNITS = [('1', ONE), ('-1', neg(ONE)), ('EPS', EPS), ('-EPS', neg(EPS))]
ADM8 = list(LISTED.items())
RESULT = {}
for _X in FOUND:
    _Y2 = _X**3 - 144*_X + 1296
    _Y = math.isqrt(_Y2)
    _al, _be, _yt = (_X, -6, 0), (_X * _X - 144, 6 * _X, 36), (_Y, 0, 0)
    _base_ok = _Y > 0 and mul(_al, _be) == mul(_yt, _yt) and norm(_al) == _Y2
    _a1, _b1, _y1, _E = _al, _be, _yt, {}
    for _nm, _p in SPRIMES:
        _e, _a1 = exponent(_p, _a1)
        _f, _b1 = exponent(_p, _b1)
        _w, _y1 = exponent(_p, _y1)
        _E[_nm] = (_e, _f, _w)
    _match_ok = all(e + f == 2 * w for e, f, w in _E.values())
    _free_ok = all(div_R(x, p) is None for _, p in SPRIMES for x in (_a1, _b1, _y1))
    _recon_ok = (_al == prod(*[pw(p, _E[nm][0]) for nm, p in SPRIMES], _a1)
                 and _be == prod(*[pw(p, _E[nm][1]) for nm, p in SPRIMES], _b1)
                 and _yt == prod(*[pw(p, _E[nm][2]) for nm, p in SPRIMES], _y1))
    _sq_ok = mul(_a1, _b1) == mul(_y1, _y1)
    _gm = minors_gcd(_a1, _b1)
    _hits, _certs, _unit_ok = [], [], True
    for _un, _u in UNITS:
        _verdict, _data = sqrt_in_R(div_R(_a1, _u))
        if _verdict == 'square':
            _hits.append((_un, _u, _data))
        elif _verdict == 'nonsquare':
            _certs.append((_un, _data))
        else:
            _unit_ok = False
    info(f'X = {_X}, Y = {_Y}: exponents (in alpha, beta, Y): ' + ', '.join(f'{nm} {_E[nm]}' for nm, _ in SPRIMES))
    check(f'X = {_X}: e + f = 2w for each S-prime, the S-free parts satisfy alpha\'beta\' = Y\'^2, and alpha\', beta\' are'
          f' coprime (gcd of the 3x3 minors of [M(alpha\') | M(beta\')] = {_gm}); alpha\' = u R0^2 for exactly one u in'
          f' {{1, -1, EPS, -EPS}}, the other three certified non-square: {_certs}',
          _base_ok and _match_ok and _free_ok and _recon_ok and _sq_ok and _gm == 1 and _unit_ok and len(_hits) == 1)
    _un, _u, _R0 = _hits[0]
    _e2, _ea, _eb, _ec, _eg = (_E[nm][0] for nm, _ in SPRIMES)
    _fact_ok = prod(_u, pw(PI2, _e2), pw(P3A, _ea), pw(P3B, _eb), pw(P3C, _ec), pw(P179, _eg), _R0, _R0) == _al
    _sigma = norm(_u)
    _parity_ok = _e2 % 2 == 0 and _eg % 2 == 0 and (_ea + _eb + _ec) % 2 == 0 and _sigma * (-1)**(_ea + _eb) == 1
    _mu0 = smul(_sigma, prod(pw(EPS, 1 if 'EPS' in _un else 0), pw(P3A, _ea % 2), pw(P3B, _eb % 2), pw(P3C, _ec % 2)))
    _name0 = ([k for k, v in LISTED.items() if v == _mu0] or [None])[0]
    _S_struct = prod(pw(PI2, _e2 // 2), pw(P3A, _ea // 2), pw(P3B, _eb // 2), pw(P3C, _ec // 2), pw(P179, _eg // 2), _R0)
    _q0 = div_R(_al, _mu0)
    _verdict, _S = sqrt_in_R(_q0) if _q0 is not None else ('unknown', None)
    _S_ok = (_verdict == 'square' and mul(_mu0, mul(_S, _S)) == _al and _S_struct in (_S, neg(_S))
             and mul(_mu0, mul(_S_struct, _S_struct)) == _al)
    check(f'X = {_X}: X - theta = ({_un}) PI2^{_e2} P3A^{_ea} P3B^{_eb} P3C^{_ec} P179^{_eg} R0^2 with R0 = {_R0} (exact);'
          f' e2, g even, a+b+c even, N(u)(-1)^(a+b) = 1; mu0 = {_name0} is one of the four targets; X - theta = mu0 S^2'
          f' with S = {_S} (exact; equal up to sign to the S built from the exponents as in 5.6)',
          _fact_ok and _parity_ok and _name0 in T_NAMES and _S_ok)
    _sq8, _ns8, _unk8 = classify(_al, ADM8)
    _sq4, _ns4, _unk4 = classify(_al, MU4)
    _cls = _sq4[0][0] if len(_sq4) == 1 else None
    _W = _sq4[0][1] if len(_sq4) == 1 else None
    check(f'X = {_X}: among the 8 admissible mu0 only {[n for n, _ in _sq8]} has alpha*mu0 a square (the other 7'
          f' certified: {[(n, c) for n, c in _ns8]}); class {_cls}: (X - theta)({_cls}) = W^2 with W = {_W}; the other'
          f' three classes certified non-square {[(n, c) for n, c in _ns4]}; agrees with mu0 and with the 2-descent'
          f' grouping', not _unk8 and not _unk4 and [n for n, _ in _sq8] == [_name0] and _cls is not None
          and MU0_OF[_cls] == _name0 and _X in DESCENT[_cls])
    RESULT[_X] = {'mu0': _name0, 'class': _cls, 'S': _S}
check(f'the integers X with |X| <= 2*10^6 and X^3 - 144X + 1296 a square are exactly {EXPECTED_X} (14 values, 28'
      ' integral points), and every one has its X - theta in one of the four target classes', FOUND == EXPECTED_X)
info('class of X - theta: ' + '; '.join(f'{cls}: {[x for x in FOUND if RESULT[x]["class"] == cls]}' for cls, _ in MU4))
check('each of the four target classes contains integral points (so none of the four can be dropped)',
      all(any(RESULT[x]['class'] == cls for x in FOUND) for cls, _ in MU4))

# =====================================================================================================================
print('(f) negative controls (planted errors)')
_PLANTED_MU4 = [('1', ONE), ('9-theta', (9, -6, 0)), ('EPS*(-15-theta)', mul(EPS, (-15, -6, 0))),
                ('-12-theta', (-12, -6, 0))]
_PLANTED_T = ['1', 'EPS P3A P3B', '-EPS P3A P3C', '-EPS P3B P3C']
_det_mu, _det_mu0, _cert_example = [], [], None
for _X in FOUND:
    _sq, _ns, _unk = classify((_X, -6, 0), _PLANTED_MU4)
    if _unk:
        check(f'planted list at X = {_X}: every class decided', False)
    if not _sq:
        _det_mu.append(_X)
        if _cert_example is None:
            _cert_example = (_X, _ns)
    if RESULT[_X]['mu0'] not in _PLANTED_T:
        _det_mu0.append(_X)
check(f'planted error 1 (the task\'s control): with -15-theta replaced by EPS*(-15-theta) = {mul(EPS, (-15, -6, 0))}'
      f' in the target list, check (e) finds no target class for X in {_det_mu}, each planted class certified'
      f' non-square (e.g. X = {_cert_example[0] if _cert_example else None}: {_cert_example[1] if _cert_example else None})'
      ': the error is detected exactly at the four points of the class -15-theta',
      _det_mu == [-15, 12, 48, 4953])
check(f'planted error 1, mu0 form: with -P3A P3C replaced by -EPS P3A P3C in T, the mu0 found by the factorization lies'
      f' outside the list for X in {_det_mu0}: detected at the same four points', _det_mu0 == [-15, 12, 48, 4953])
_bw = add(add((X**2 - 144, 0, 0), smul(X, THETA)), mul(THETA, THETA, WRONG_F))
_abw = [sp.expand(c) for c in mul((X, -6, 0), _bw, WRONG_F)]
check(f'planted error 2: in the ring with l^3 = 4l + 6 the product (X - theta)(X^2 + theta X + theta^2 - 144) is {_abw},'
      ' not (X^3 - 144X + 1296, 0, 0): the identity check of (a) depends on the ring',
      not same(_abw, (X**3 - 144*X + 1296, 0, 0)))
check('planted error 3: with r2 = 5 + 2l + l^2 instead of 5 + 2l, the identity -15 - theta = (-P3A P3C) r2^2 fails',
      mul(LISTED['-P3A P3C'], mul((5, 2, 1), (5, 2, 1))) != (-15, -6, 0))
_LISTED_BAD = dict(LISTED)
_LISTED_BAD['-P3A P3C'] = prod(P3A, P3C)
check(f'planted error 4: a sign error in the list of 5.6 (+P3A P3C for -P3A P3C; N = {norm(prod(P3A, P3C))}) makes the'
      ' comparison of (b) fail', sorted(r for _, r, _ in _admissible) != sorted(_LISTED_BAD.values()))
check(f'planted error 5: the step needs modulus 4; modulo 2 every mu0 of EPS*T still has residues'
      f' ({[H2[nm] for nm in EPST_NAMES]} of 8), so a test mod 2 would not exclude the coset',
      all(H2[nm] > 0 for nm in EPST_NAMES))
_rng = random.Random(20261015)
_sq_ok, _ns_ok, _trials = True, True, 0
for _ in range(200):
    _S = tuple(_rng.randint(-10**4, 10**4) for _ in range(3))
    if _S == (0, 0, 0):
        continue
    _trials += 1
    _vv, _TT = sqrt_in_R(mul(_S, _S))
    _sq_ok = _sq_ok and _vv == 'square' and _TT in (_S, neg(_S))
    for _mu in (EPS, neg(ONE), LISTED['P3A P3B'], LISTED['-P3B P3C'], LISTED['-EPS P3A P3C']):
        _vv2, _ = sqrt_in_R(mul(_mu, mul(_S, _S)))
        _ns_ok = _ns_ok and _vv2 == 'nonsquare'
check(f'control of the square-root routine: for {_trials} random S with |coordinates| <= 10^4 it returns +-S from S^2,'
      ' and certifies mu S^2 as a non-square for mu = EPS, -1, P3A P3B, -P3B P3C, -EPS P3A P3C', _sq_ok and _ns_ok)
_known_ns = [('-1', neg(ONE)), ('2', (2, 0, 0)), ('3', (3, 0, 0)), ('PI2', PI2), ('P3A', P3A)] + \
            [(nm, LISTED[nm]) for nm in LISTED if nm != '1'] + MU4[1:]
_ns_res = {nm: sqrt_in_R(m) for nm, m in _known_ns}
check(f'control: every non-trivial class representative is certified non-square in K: '
      f'{ {nm: v[1] for nm, v in _ns_res.items()} }', all(v[0] == 'nonsquare' for v in _ns_res.values()))

# =====================================================================================================================
print(f'SUMMARY: all {N_PASS} checks passed. Section 5 of hypothesis_H.md verified independently: the identities of'
      ' 5.2, 5.3 and 5.8; exactly 8 of the 128 parity vectors have a positive square norm, the 8 mu0 of 5.6; modulo 4 the'
      ' coset EPS*T never reaches (l, l^2) = (2, 0) and the four targets do (mod 8 counts 0,0,0,0 / 32,48,48,32); the'
      ' four classes are pairwise distinct (characters at split primes other than 31, 61); every integral point with'
      f' |X| <= 2*10^6, X in {FOUND}, has X - theta = mu0 S^2 with mu0 in the four targets; the planted errors were'
      ' all detected.')
