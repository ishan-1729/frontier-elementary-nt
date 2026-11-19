"""00_fields.py -- the three quartic Thue fields of 8604.a1 (triangular = square-pyramidal), checked exactly.

Everything below is exact (Python integers, Fractions, sympy).  Floats (mpmath, 60 digits) are used only to print
information (roots, unit logarithms) and to PROPOSE an isomorphism of two fields if no certificate of
non-isomorphism is found; such a proposal would then be verified exactly.

Notation (core.py).  For each class key '10', '01', '11' of the 2-descent of 8604.a1 (y^2 = x^3 - 9x + 81):
  F_cl(m,n)  the CLASS FORM of the 2-descent (core.CLASSES[key]['CLASS_F']); f_cl(t) = F_cl(1,t) is monic, alpha is a
             root of f_cl, and F_cl(m,n) = N(n - m alpha);
  G(m,n)     the FIELD FORM that is solved (core.CLASSES[key]['F']), G(m,n) = a0 m^4 + a1 m^3 n + a2 m^2 n^2
             + a3 m n^3 + n^4; g(t) = G(1,t) is monic, beta is a root of g, and G(m,n) = N(n - m beta).
  Classes (0,1) and (1,1): G = F_cl and beta = alpha.  Class (1,0): beta = (1 + alpha)/2 and F_cl(m,n) = G(2m, n + m).
  An element of the order O is given by its coordinates in the basis omega_0..omega_3 of core.CLASSES[key]['BASIS'].

 1. g is irreducible over Q: it has no integer root (a rational root of a monic integer polynomial is an integer
    dividing g(0)) and no factorization into monic integer quadratics (Gauss's lemma: none over Q either); sympy
    agrees.  G(m,n) = m^4 g(n/m) = Res_t(g(t), n - m t) = det(multiplication by n - m beta) = N(n - m beta) as
    polynomials.  disc(g) is factored.  g has exactly two real roots (own exact Sturm chain; sympy count_roots).
    Galois group S4: the resolvent cubic R(y) = y^3 - b y^2 + (ac - 4d) y - (a^2 d - 4bd + c^2) of
    g = t^4 + a t^3 + b t^2 + c t + d (whose roots are b1b2+b3b4, b1b3+b2b4, b1b4+b2b3: checked symbolically) has
    no rational root, and disc(g) < 0 is not a square.
 2. The basis of core.CLASSES spans an order O (core.Field asserts an integral multiplication table) containing
    Z[beta] (the inverse basis matrix is integral); [O : Z[beta]] = 1/|det B|; disc(O) = disc(g) det(B)^2, equal
    to det(Tr(omega_i omega_j)).
    MAXIMALITY LEMMA.  Let O be an order of K with Z-basis omega_0..omega_3 and O_K the maximal order.  If a prime p
    divides [O_K : O], then O_K/O has an element of order p: y in O_K \\ O with p y in O.  Subtracting an element of
    O, y = (sum c_i omega_i)/p with 0 <= c_i < p, not all 0.  Such a y is integral iff its characteristic polynomial
    lies in Z[t] (Cayley-Hamilton one way; the other way it is a power of the minimal polynomial of an algebraic
    integer), i.e. iff p^k divides the coefficient c_k of t^(4-k) in det(t - M(x)), x = sum c_i omega_i, k = 1..4.
    So if none of the p^4 - 1 classes gives an integral y, p does not divide the index.  Since
    disc(O) = [O_K : O]^2 disc(O_K), only primes with p^2 | disc(O) can divide the index; all of them are tested.
    Cross-checks: the maximal order computed from Z[beta] by the same test with enlargement has the same Hermite
    normal form as the basis of core.CLASSES, and so does the integral basis of sympy's round_two.
    Maximality is NOT used by the proof.  What is used is that O is a ring containing Z[beta] and the units.
 3. MARKER LEMMA.  omega_0 = 1, omega_1 = beta, and the beta^2-, beta^3-power coordinates of x = sum x_i omega_i
    are (x2, x3) D with an invertible 2x2 matrix D.  As 1, beta, beta^2, beta^3 are linearly independent
    (g irreducible), x lies in Q + Q beta iff x2 = x3 = 0, and then x = x0 + x1 beta with x0, x1 in Z.  The rows
    0, 1 of the inverse basis matrix are (1,0,0,0), (0,1,0,0), so n - m beta has coordinates (n, -m, 0, 0).
 4. SUBSTITUTION LEMMA.  (a) F_cl(m,n) = G(2m, n + m) for class (1,0); G = F_cl for (0,1) and (1,1) (sympy).
    (b) Class (1,0): (1 + alpha)/2 has characteristic polynomial g(t), monic in Z[t], so it is an algebraic integer
    in Q + Q alpha outside Z + Z alpha, and the marker lemma does not hold for alpha; it holds for beta = (1+alpha)/2
    (section 3), and Q + Q alpha = Q + Q beta.  Classes (0,1), (1,1): it holds for alpha itself (section 3).
    (c) For c in RHS, core.from_class (m,n) -> (2m, n + m) is a bijection from the coprime solutions of
    F_cl(m,n) = c onto the coprime solutions of G(m',n') = c with m' even, with inverse core.to_class
    (m',n') -> (m'/2, n' - m'/2), compatible with (m,n) -> (-m,-n).  The congruence it needs (F_cl(m,n) == 0 mod 16
    when m, n are both odd) is verified at all residues mod 16.  A coprime solution of G = c with m' odd gives
    F_cl(m', 2n' - m') = 16c, which is not in RHS.  For (0,1), (1,1) both maps are the identity.
 5. N(eps1) = N(eps2) = 1 for core.UNITS (exact determinants and characteristic polynomials); their inverses lie in O.
 6. The known solutions core.KNOWN of G grouped by c; their class-form images (core.to_class); the integral point of
    each class-form solution (core.descent_maps) and, when X == 3 (mod 6), the common value b(b+1)/2 = a(a+1)(2a+1)/6
    with (X, Y) = (6a + 3, 18b + 9); a confirmation box search |m|,|n| <= 1000 for G and for F_cl (not used by any
    proof).
 7. The three quartic fields are pairwise non-isomorphic: different factorization patterns of the g's modulo a prime
    dividing neither discriminant (Dedekind), and the field discriminants of section 2.  Not used by any proof.

Exit status: non-zero iff some CHECK FAILED.
"""
import os
import sys
import time
import itertools
from fractions import Fraction as Fr
from math import gcd, isqrt

sys.dont_write_bytecode = True          # keep the scripts folder free of __pycache__

import numpy as np                      # noqa: E402
import sympy as sp                      # noqa: E402
import mpmath                           # noqa: E402
from sympy.polys.numberfields.basis import round_two      # noqa: E402  (cross-check of section 2 only)

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from core import (CURVE, CURVE_EQ, GA, GB, CLASSES, RHS, KNOWN, UNITS, Field, Fval, Fclass,  # noqa: E402
                  fcoeffs, pmulmod, mat_inv, vec_mat, det_int, descent_maps, to_class, from_class)

T0 = time.time()
mpmath.mp.dps = 60
KEYS = ('10', '01', '11')
LABEL = {k: CLASSES[k]['label'] for k in KEYS}
FAILED = []
t_, m_, n_, y_ = sp.symbols('t m n y')
beta_s = sp.Symbol('beta')

# Tables stated in the write-up (Theorem A' for the field forms G, Theorem A for the class forms F_cl): all integer
# solutions, one of each pair +-(m,n).  They are compared below with core.KNOWN and its images; the proof that they
# are complete is 01_coverage.py + 02_exclusion.py with section 4 of this script.
STATED_G = {'10': {1: [(0, 1)], 4: [(1, 0)], 9: [(1, -1), (2, 1)], 36: []},
            '01': {1: [(0, 1)], 4: [], 9: [(1, 0)], 36: [(1, 3)]},
            '11': {1: [(0, 1), (2, 1)], 4: [(1, 4)], 9: [], 36: [(1, 0)]}}
STATED_F = {'10': {1: [(0, 1)], 4: [], 9: [(1, 0)], 36: []},
            '01': {1: [(0, 1)], 4: [], 9: [(1, 0)], 36: [(1, 3)]},
            '11': {1: [(0, 1), (2, 1)], 4: [(1, 4)], 9: [], 36: [(1, 0)]}}
# Integral points of each non-trivial class, from weierstrass_2descent.out.txt (CURVE 8604.a1, "grouped by class");
# the trivial class (0,0) has X in {-5, 7, 1099} there (a reducible form; not treated here).
DESCENT_X = {'10': [3, 39], '01': [-3, 24, 33], '11': [0, 9, 513, 5112]}
TRIVIAL_X = [-5, 7, 1099]
SEARCH_X = [-5, -3, 0, 3, 7, 9, 24, 33, 39, 513, 1099, 5112]     # same file: search -92 <= x <= 10^6, 24 points


def check(cond, msg):
    if cond:
        print(f'   [ok] {msg}')
    else:
        print(f'   CHECK FAILED: {msg}')
        FAILED.append(msg)
    return cond


def header(s):
    print()
    print('=' * 116)
    print(s)
    print('=' * 116)


# ------------------------------------------------------------------------------------------------------------------
# small exact helpers (from ../../triangular_tetrahedral/scripts/00_fields.py)
# ------------------------------------------------------------------------------------------------------------------
def fexpr(fc, var=t_):
    """sympy expression of the polynomial with coefficients fc (low degree first; integers or Fractions)."""
    return sum(sp.Rational(Fr(c).numerator, Fr(c).denominator) * var ** k for k, c in enumerate(fc))


def form_expr(coef, x, z):
    """Binary quartic sum coef[i] x^(4-i) z^i."""
    return sum(sp.Integer(c) * x ** (4 - i) * z ** i for i, c in enumerate(coef))


def elt(pc, var):
    """sympy expression of an element given by power coordinates (Fractions)."""
    return sum(sp.Rational(Fr(c).numerator, Fr(c).denominator) * var ** k for k, c in enumerate(pc))


def fz(s):
    return dict(sorted(sp.factorint(s).items()))


def sgn(x):
    return (x > 0) - (x < 0)


def p_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def p_rem(a, b):
    a = p_trim([Fr(x) for x in a])
    b = p_trim([Fr(x) for x in b])
    db = len(b) - 1
    while len(a) - 1 >= db and any(a):
        q = a[-1] / b[-1]
        s = len(a) - 1 - db
        for i in range(db + 1):
            a[s + i] -= q * b[i]
        assert a[-1] == 0
        a = p_trim(a[:-1]) if len(a) > 1 else [Fr(0)]
    return a


def sturm_count(fc):
    """Number of distinct real roots of the polynomial fc (low first), by an exact Sturm chain."""
    chain = [[Fr(c) for c in fc], [Fr(k * fc[k]) for k in range(1, len(fc))]]
    while True:
        r = p_rem(chain[-2], chain[-1])
        if not any(r):
            break
        chain.append([-c for c in r])

    def var(signs):
        s = [x for x in signs if x != 0]
        return sum(1 for i in range(len(s) - 1) if s[i] != s[i + 1])
    at_pos = [sgn(p[-1]) for p in chain]
    at_neg = [sgn(p[-1]) * (-1) ** (len(p) - 1) for p in chain]
    last_const = len(chain[-1]) == 1 and chain[-1][0] != 0
    return var(at_neg) - var(at_pos), len(chain), last_const


def divisors_pm(a):
    a = abs(int(a))
    ds = sp.divisors(a)
    return sorted(set(ds + [-d for d in ds]))


def quartic_irreducible_elementary(fc):
    """Monic integer quartic: no integer root and no factorization (t^2+pt+q)(t^2+rt+s) over Z."""
    c0, c1, c2, c3 = (int(c) for c in fc[:4])
    assert c0 != 0
    for r in divisors_pm(c0):
        if sum(c * r ** k for k, c in enumerate(fc)) == 0:
            return False, f'integer root {r}'
    for q in divisors_pm(c0):
        s = c0 // q
        if q != s:
            if (c1 - q * c3) % (s - q):
                continue
            p = (c1 - q * c3) // (s - q)
            r = c3 - p
            if q + s + p * r == c2:
                return False, f'(t^2 + {p} t + {q})(t^2 + {r} t + {s})'
        else:
            if q * c3 != c1:
                continue
            D = c3 * c3 - 4 * (c2 - 2 * q)
            if D >= 0 and isqrt(D) ** 2 == D and (c3 + isqrt(D)) % 2 == 0:
                p = (c3 + isqrt(D)) // 2
                return False, f'(t^2 + {p} t + {q})(t^2 + {c3 - p} t + {s})'
    return True, 'no integer root among the divisors of g(0); no quadratic factorization over Z'


def resolvent(fc):
    """Resolvent cubic (low first) of f = t^4 + a t^3 + b t^2 + c t + d."""
    d, c, b, a = (int(x) for x in fc[:4])
    return [-(a * a * d - 4 * b * d + c * c), a * c - 4 * d, -b, 1]


def charpoly(M):
    """[c_0, ..., c_n] with det(t I - M) = sum c_k t^k, for an integer square matrix (Faddeev-LeVerrier, exact)."""
    n = len(M)
    c = [0] * (n + 1)
    c[n] = 1
    Mk = [[0] * n for _ in range(n)]
    for k in range(1, n + 1):
        AM = [[sum(M[i][l] * Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        Mk = [[AM[i][j] + (c[n - k + 1] if i == j else 0) for j in range(n)] for i in range(n)]
        AMk = [[sum(M[i][l] * Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        tr = sum(AMk[i][i] for i in range(n))
        assert tr % k == 0
        c[n - k] = -tr // k
    return c


def charpoly_q(M):
    """The same for a matrix of rationals (Fractions), exact."""
    n = len(M)
    M = [[Fr(x) for x in row] for row in M]
    c = [Fr(0)] * (n + 1)
    c[n] = Fr(1)
    Mk = [[Fr(0)] * n for _ in range(n)]
    for k in range(1, n + 1):
        AM = [[sum(M[i][l] * Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        Mk = [[AM[i][j] + (c[n - k + 1] if i == j else 0) for j in range(n)] for i in range(n)]
        AMk = [[sum(M[i][l] * Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        c[n - k] = -sum(AMk[i][i] for i in range(n)) / k
    return c


class GenField(Field):
    """core.Field for an arbitrary monic quartic fc (low first) and a basis B (rows = power coordinates)."""

    def __init__(self, fc, B, label=''):
        self.key = None
        self.label = label
        self.fc = [int(c) for c in fc]
        self.B = [[Fr(c) for c in row] for row in B]
        self.Binv = mat_inv(self.B)
        self.mt = [[None] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                co = vec_mat(pmulmod(self.B[i], self.B[j], self.fc), self.Binv)
                assert all(c.denominator == 1 for c in co), 'basis does not span an order'
                self.mt[i][j] = [int(c) for c in co]


def integral_classes(K, p):
    """All nonzero (c_0..c_3) in [0,p)^4 with (sum c_i omega_i)/p integral (characteristic polynomial test)."""
    out = []
    for cs in itertools.product(range(p), repeat=4):
        if not any(cs):
            continue
        cp = charpoly(K.matrix(cs))
        assert cp[0] == K.norm(cs)
        if all(cp[4 - k] % p ** k == 0 for k in range(1, 5)):
            out.append(cs)
    return out


def ext_gcd(a, b):
    """(g, s, t) with s a + t b = g = gcd(a, b) >= 0."""
    s0, s1, t0, t1, r0, r1 = 1, 0, 0, 1, a, b
    while r1:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
        t0, t1 = t1, t0 - q * t1
    if r0 < 0:
        r0, s0, t0 = -r0, -s0, -t0
    assert s0 * a + t0 * b == r0
    return r0, s0, t0


def det_fr(A):
    """Exact determinant of a square matrix of Fractions (Gaussian elimination)."""
    M = [[Fr(c) for c in row] for row in A]
    n, d = len(M), Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return d


def hnf_int(vectors):
    """Lower-triangular Hermite normal form (rows; row k has zeros right of column k, positive diagonal,
    entries left of the diagonal reduced into [0, diagonal of their column)) of a full-rank lattice in Z^4."""
    vecs = [list(map(int, v)) for v in vectors if any(v)]
    basis = [None] * 4
    for col in range(3, -1, -1):
        piv, rest = None, []
        for v in vecs:
            if v[col] == 0:
                rest.append(v)
            elif piv is None:
                piv = v
            else:
                a, b = piv[col], v[col]
                g, s, t = ext_gcd(a, b)
                newp = [s * x + t * y for x, y in zip(piv, v)]
                newv = [(b // g) * x - (a // g) * y for x, y in zip(piv, v)]
                piv = newp
                if any(newv):
                    rest.append(newv)
        assert piv is not None, 'lattice not of full rank'
        if piv[col] < 0:
            piv = [-x for x in piv]
        basis[col] = piv
        vecs = rest
    assert not vecs
    for i in range(4):
        for j in range(i - 1, -1, -1):
            q = basis[i][j] // basis[j][j]
            basis[i] = [x - q * y for x, y in zip(basis[i], basis[j])]
    return basis


def hnf_q(vectors):
    D = 1
    for v in vectors:
        for c in v:
            d = Fr(c).denominator
            D = D * d // gcd(D, d)
    H = hnf_int([[int(Fr(c) * D) for c in v] for v in vectors])
    return [[Fr(c, D) for c in row] for row in H]


def ring_closure(fc, gens):
    """HNF basis of the smallest Z-lattice containing gens (power coordinates, including 1) closed under products."""
    L = hnf_q(gens)
    while True:
        prods = [pmulmod(L[i], L[j], fc) for i in range(4) for j in range(i, 4)]
        L2 = hnf_q(L + prods)
        if L2 == L:
            return L
        L = L2


def maximal_order(fc, log):
    """Maximal order of Q[t]/(fc) from Z[t]/(fc): p-maximality test with enlargement for every p with p^2 | disc."""
    disc = int(sp.discriminant(fexpr(fc), t_))
    primes = [p for p, e in fz(abs(disc)).items() if e >= 2]
    B = [[Fr(int(i == j)) for j in range(4)] for i in range(4)]
    for p in primes:
        while True:
            K = GenField(fc, B)
            found = integral_classes(K, p)
            if not found:
                log.append(f'p={p}: none of the {p ** 4 - 1} classes is integral -> {p}-maximal')
                break
            cs = found[0]
            y = [sum(Fr(cs[i]) * K.B[i][j] for i in range(4)) / p for j in range(4)]
            log.append(f'p={p}: {len(found)} integral classes; adjoin {elt(y, sp.Symbol("t"))}')
            B = ring_closure(fc, B + [y])
    K = GenField(fc, B)
    for p in primes:                          # final pass on the final order
        assert not integral_classes(K, p)
    return B, primes, disc


def disc_of_basis(fc, B):
    """(det B, disc(O) = disc(f) det(B)^2, trace-form determinant) for an order basis B."""
    K = GenField(fc, B)
    detB = det_fr(B)
    disc_f = int(sp.discriminant(fexpr(fc), t_))
    tr1 = [sum(K.matrix(tuple(int(i == j) for j in range(4)))[k][k] for k in range(4)) for i in range(4)]
    TF = [[sum(K.mt[i][j][k] * tr1[k] for k in range(4)) for j in range(4)] for i in range(4)]
    return detB, disc_f * detB * detB, det_int(TF)


def poly_at(hc, P, fc):
    """h(P) in Q[t]/(fc), P given by power coordinates."""
    acc = [Fr(0)] * 4
    for c in reversed(hc):
        acc = pmulmod(acc, P, fc)
        acc[0] += c
    return acc


def mp_roots(fc):
    return mpmath.polyroots(list(reversed([int(c) for c in fc])), maxsteps=800, extraprec=800)


def rationalize(x, maxden=10 ** 12):
    return Fr(mpmath.nstr(x, 55, min_fixed=-10 ** 9, max_fixed=10 ** 9)).limit_denominator(maxden)


def find_embeddings(hc, fc):
    """All roots of h (monic, low first) in Q[t]/(f): candidates from the 24 root bijections (floats), then exact
    verification h(P(beta)) = 0.  Returns the list of power-coordinate vectors P."""
    br, ar = mp_roots(hc), mp_roots(fc)
    V = mpmath.matrix([[a ** j for j in range(4)] for a in ar])
    out = []
    for perm in itertools.permutations(range(4)):
        c = mpmath.lu_solve(V, mpmath.matrix([br[perm[j]] for j in range(4)]))
        if max(abs(mpmath.im(x)) for x in c) > mpmath.mpf(10) ** -35 * (1 + max(abs(x) for x in c)):
            continue
        P = [rationalize(mpmath.re(x)) for x in c]
        if P not in out and not any(poly_at(hc, P, fc)):
            out.append(P)
    return out


def split_pattern(fc, p):
    _, facs = sp.Poly(list(reversed([int(c) for c in fc])), t_, modulus=p).factor_list()
    return tuple(sorted(sp.degree(g, t_) for g, e in facs for _ in range(e)))


def nonisomorphism_certificate(fa, fb, pmax=400):
    da = int(sp.discriminant(fexpr(fa), t_))
    db = int(sp.discriminant(fexpr(fb), t_))
    for p in sp.primerange(2, pmax):
        if (da * db) % p == 0:
            continue
        sa, sb = split_pattern(fa, p), split_pattern(fb, p)
        if sa != sb:
            return p, sa, sb
    return None


def pm_set(pairs):
    return {frozenset({(a, b), (-a, -b)}) for (a, b) in pairs}


# ==================================================================================================================
print(f'00_fields.py -- the quartic Thue fields of {CURVE} ({CURVE_EQ}; triangular = square-pyramidal); exact checks')
print('conventions: G(m,n) = a0 m^4 + a1 m^3 n + a2 m^2 n^2 + a3 m n^3 + n^4 is the FIELD FORM of a class, g(t) = G(1,t),')
print('             beta a root of g; F_cl is the CLASS FORM of the 2-descent, f_cl(t) = F_cl(1,t), alpha a root of f_cl;')
print('             an element is given by its coordinates in the integral basis omega_0..omega_3 of core.CLASSES.')

# ------------------------------------------------------------------------------------------------------------------
header('1. g irreducible; G(m,n) = m^4 g(n/m) = N(n - m beta); disc(g); two real roots; Galois group S4')
x1, x2, x3, x4 = sp.symbols('x1:5')
e1, e2 = x1 + x2 + x3 + x4, x1 * x2 + x1 * x3 + x1 * x4 + x2 * x3 + x2 * x4 + x3 * x4
e3, e4 = x1 * x2 * x3 + x1 * x2 * x4 + x1 * x3 * x4 + x2 * x3 * x4, x1 * x2 * x3 * x4
A_, B_, C_, D_ = -e1, e2, -e3, e4
lhs = (y_ - (x1 * x2 + x3 * x4)) * (y_ - (x1 * x3 + x2 * x4)) * (y_ - (x1 * x4 + x2 * x3))
rhs = y_ ** 3 - B_ * y_ ** 2 + (A_ * C_ - 4 * D_) * y_ - (A_ * A_ * D_ - 4 * B_ * D_ + C_ * C_)
check(sp.expand(lhs - rhs) == 0,
      'resolvent formula: prod over the 3 pairings of (y - (x_i x_j + x_k x_l)) = y^3 - b y^2 + (ac - 4d) y'
      ' - (a^2 d - 4bd + c^2) for g = prod (t - x_i) = t^4 + a t^3 + b t^2 + c t + d (symbolic identity)')
DISC = {}
for key in KEYS:
    Gco = CLASSES[key]['F']
    gc = fcoeffs(key)
    K = Field(key)
    print()
    print(f'CLASS {LABEL[key]}  (mu = {CLASSES[key]["mu"]})   G(m,n) = {sp.expand(form_expr(Gco, m_, n_))}')
    print(f'   class form F_cl(m,n) = {sp.expand(form_expr(CLASSES[key]["CLASS_F"], m_, n_))};'
          f'  core: {CLASSES[key]["SUB"]}')
    print(f'   g(t) = G(1,t) = {fexpr(gc)}')
    ok, why = quartic_irreducible_elementary(gc)
    check(ok, f'g irreducible over Q (elementary: {why})')
    check(sp.Poly(fexpr(gc), t_).is_irreducible, 'g irreducible over Q (sympy)')
    Ge = sp.expand(form_expr(Gco, m_, n_))
    check(sp.expand(sp.cancel(m_ ** 4 * fexpr(gc).subs(t_, n_ / m_))) == Ge, 'm^4 g(n/m) = G(m,n)')
    res = sp.expand(sp.resultant(fexpr(gc), n_ - m_ * t_, t_))
    check(res == Ge, 'Res_t(g(t), n - m t) = G(m,n)   [= prod_i (n - m beta_i) = N(n - m beta), g monic]')
    x = (n_, -m_, 0, 0)
    Msym = sp.Matrix(4, 4, lambda k, j: sum(x[i] * K.mt[i][j][k] for i in range(4)))
    check(sp.expand(Msym.det()) == Ge,
          'det(matrix of multiplication by n - m beta = (n,-m,0,0) in the basis of core.Field) = G(m,n)')
    disc = int(sp.discriminant(fexpr(gc), t_))
    DISC[key] = disc
    print(f'   disc(g) = {disc} = {fz(disc)}')
    nreal, chainlen, lastc = sturm_count(gc)
    check(lastc, f'Sturm chain of length {chainlen} ends in a nonzero constant (g squarefree)')
    check(nreal == 2, f'exactly two real roots (exact Sturm count = {nreal})')
    check(sp.Poly(fexpr(gc), t_).count_roots() == 2, 'exactly two real roots (sympy count_roots)')
    rts = mp_roots(gc)
    rr = sorted([mpmath.re(r) for r in rts if abs(mpmath.im(r)) < 1e-40])
    cr = [r for r in rts if mpmath.im(r) > 1e-40]
    print(f'   roots (floats, for information): real {mpmath.nstr(rr[0], 15)}, {mpmath.nstr(rr[1], 15)};'
          f' complex {mpmath.nstr(cr[0], 15)} and its conjugate')
    R = resolvent(gc)
    print(f'   resolvent cubic R(y) = {fexpr(R, y_)}')
    roots_R = [r for r in divisors_pm(R[0]) if sum(c * r ** k for k, c in enumerate(R)) == 0] if R[0] != 0 else [0]
    check(not roots_R, f'R has no rational root (a rational root is an integer dividing R(0) = {R[0]}: none of'
                       f' the {len(divisors_pm(R[0])) if R[0] else 0} divisors is a root) -> R irreducible over Q')
    check(int(sp.discriminant(fexpr(R, y_), y_)) == disc, 'disc(R) = disc(g)')
    check(disc < 0, 'disc(g) < 0, so disc(g) is not a square in Q')
    print('   => Galois group of g is S4 (g irreducible, R irreducible, disc(g) not a square)')

# ------------------------------------------------------------------------------------------------------------------
header('2. The integral basis: an order containing Z[beta]; index; discriminant; maximality (not used by the proof)')
DISC_O, INDEX = {}, {}
for key in KEYS:
    gc = fcoeffs(key)
    K = Field(key)
    print()
    print(f'CLASS {LABEL[key]}   basis (power coordinates of omega_i):')
    for i, row in enumerate(K.B):
        print(f'      omega_{i} = {elt(row, beta_s)}')
    check(True, 'core.Field built: every omega_i omega_j has integral coordinates (core asserts it), so the'
                ' Z-span O of the basis is a ring')
    print('   multiplication table (coordinates of omega_i omega_j, i <= j):')
    for i in range(4):
        print('      ' + '   '.join(f'w{i}w{j} = {tuple(K.mt[i][j])}' for j in range(i, 4)))
    check(all(c.denominator == 1 for row in K.Binv for c in row),
          'the inverse basis matrix is integral: 1, beta, beta^2, beta^3 lie in O, so Z[beta] <= O')
    detB, discO, trdet = disc_of_basis(gc, K.B)
    index = 1 / abs(detB)
    check(index.denominator == 1, f'[O : Z[beta]] = 1/|det B| = {index}')
    check(discO == trdet, f'disc(O) = disc(g) det(B)^2 = {discO} equals det(Tr(omega_i omega_j)) = {trdet}')
    discO = int(discO)
    DISC_O[key] = discO
    INDEX[key] = int(index)
    fdisc = fz(discO)
    print(f'   disc(O) = {discO} = {fdisc}')
    test_primes = [p for p, e in fdisc.items() if p > 1 and e >= 2]
    other = {p: e for p, e in fdisc.items() if p > 1 and e < 2}
    for p in test_primes:
        found = integral_classes(K, p)
        check(not found, f'p = {p}: none of the {p ** 4 - 1} elements (sum c_i omega_i)/{p}, 0 <= c_i < {p},'
                         f' not all 0, is integral (characteristic polynomial test) -> O is {p}-maximal')
    check(all(e == 1 for e in other.values()), f'every other prime divides disc(O) exactly once: {other}')
    print(f'   => O is the maximal order; disc(O_K) = {discO} = {fdisc}')
    log = []
    Bmax, primes, _ = maximal_order(gc, log)
    for line in log:
        print(f'      from Z[beta]: {line}')
    check(hnf_q(Bmax) == hnf_q(K.B), 'the maximal order computed from Z[beta] (enlargement) has the same HNF as'
                                     ' the basis of core.CLASSES')
    ZK, dK = round_two(sp.Poly(fexpr(gc), t_))
    Mrt = ZK.QQ_matrix.to_Matrix()                       # columns = power coordinates of a basis
    Brt = [[Fr(int(Mrt[i, j].p), int(Mrt[i, j].q)) for i in range(4)] for j in range(4)]
    check(int(dK) == discO and hnf_q(Brt) == hnf_q(K.B), f'sympy round_two (another implementation): disc(O_K) ='
                                                          f' {dK} and an integral basis with the same HNF')
print()
print('Maximality is not used by the proof; it uses only that O is a ring containing Z[beta] (above) and the units of'
      ' section 5.')

# ------------------------------------------------------------------------------------------------------------------
header('3. Marker lemma: x in O lies in Z + Z beta iff x2 = x3 = 0; n - m beta = (n, -m, 0, 0)')
for key in KEYS:
    K = Field(key)
    B = K.B
    print(f'CLASS {LABEL[key]}')
    check(B[0] == [1, 0, 0, 0] and B[1] == [0, 1, 0, 0], 'omega_0 = 1, omega_1 = beta')
    check(all(B[i][j] == 0 for i in (0, 1) for j in (2, 3)),
          'only omega_2, omega_3 involve beta^2, beta^3: the (beta^2, beta^3) power coordinates of'
          ' x = sum x_i omega_i are (x2, x3) D')
    D = [[B[2][2], B[2][3]], [B[3][2], B[3][3]]]
    detD = D[0][0] * D[1][1] - D[0][1] * D[1][0]
    check(detD != 0, f'D = [[{D[0][0]}, {D[0][1]}], [{D[1][0]}, {D[1][1]}]] has det {detD} != 0, so the beta^2,'
                     f' beta^3 coordinates vanish iff x2 = x3 = 0; then x = x0 + x1 beta, x0, x1 in Z')
    check(K.Binv[0] == [1, 0, 0, 0] and K.Binv[1] == [0, 1, 0, 0],
          'rows 0, 1 of the inverse basis matrix are (1,0,0,0), (0,1,0,0): n - m beta has coordinates (n,-m,0,0)')
    ok = all(K.from_power([n, -m, 0, 0]) == [n, -m, 0, 0] for m in range(-3, 4) for n in range(-3, 4))
    check(ok, 'from_power(n - m beta) = (n, -m, 0, 0) on the 49 pairs |m|,|n| <= 3 (spot check of the code)')
    print('   => O cap (Q + Q beta) = Z + Z beta, and the solutions of G(m,n) = c are the elements z of O with N(z) = c'
          ' and x2(z) = x3(z) = 0')

# ------------------------------------------------------------------------------------------------------------------
header('4. Substitution lemma: class form F_cl, field form G, and the bijection of their coprime solutions')
print("core.from_class(key, m, n) sends a solution (m,n) of F_cl to a solution (m',n') of G; core.to_class(key, m', n')")
print("sends it back, or gives None when (m',n') comes from no solution of F_cl.  Right sides: "
      + '; '.join(f'class {LABEL[k]}: {CLASSES[k]["RHS"]}' for k in KEYS) + f' (core.RHS = {RHS}).')
GRID = range(-40, 41)
for key in KEYS:
    Fcl = tuple(CLASSES[key]['CLASS_F'])
    Gco = tuple(CLASSES[key]['F'])
    rhs_k = tuple(CLASSES[key]['RHS'])
    Fe = sp.expand(form_expr(Fcl, m_, n_))
    Ge = sp.expand(form_expr(Gco, m_, n_))
    fcl = [int(c) for c in Fcl]                  # f_cl(t) = F_cl(1,t), low degree first
    gc = fcoeffs(key)
    K = Field(key)
    print()
    print(f'CLASS {LABEL[key]}   F_cl(m,n) = {Fe};   G(m,n) = {Ge}')
    print(f'   core: {CLASSES[key]["SUB"]}')
    check(Fcl[4] == 1, 'F_cl has n^4-coefficient 1: f_cl(t) = F_cl(1,t) is monic')
    check(sp.expand(sp.resultant(fexpr(fcl), n_ - m_ * t_, t_)) == Fe,
          'Res_t(f_cl(t), n - m t) = F_cl(m,n): F_cl(m,n) = N(n - m alpha) (sympy)')
    if key == '10':
        print('   (a) the identity')
        check(sp.expand(form_expr(Gco, 2 * m_, n_ + m_)) == Fe,
              '(a) F_cl(m,n) = G(2m, n + m) (sympy expand: identity of forms)')
        check(sp.expand(form_expr(Fcl, m_ / 2, n_ - m_ / 2)) == Ge,
              '    equivalently G(m,n) = F_cl(m/2, n - m/2) (identity over Q)')
        check(sp.expand(fexpr(fcl).subs(t_, 2 * t_ - 1) - 16 * fexpr(gc)) == 0,
              '    f_cl(2t - 1) = 16 g(t): if g(beta) = 0 then alpha = 2 beta - 1 is a root of f_cl, i.e.'
              ' beta = (1 + alpha)/2')
        dfcl = int(sp.discriminant(fexpr(fcl), t_))
        check(dfcl == 2 ** 12 * DISC[key], f'    disc(f_cl) = {dfcl} = {fz(dfcl)} = 2^12 disc(g)'
                                           ' (the roots are alpha_i = 2 beta_i - 1)')
        xs = (n_ + m_, -2 * m_, 0, 0)
        Msym = sp.Matrix(4, 4, lambda k, j: sum(xs[i] * K.mt[i][j][k] for i in range(4)))
        check(sp.expand(Msym.det()) == Fe,
              '    n - m alpha = (n + m) - 2m beta has omega-coordinates (n + m, -2m, 0, 0), and the determinant of'
              ' multiplication by it is F_cl(m,n): F_cl(m,n) = N((n + m) - 2m beta) = G(2m, n + m)')
        print('   (b) why the substitution is needed')
        h = [Fr(1, 2), Fr(1, 2), Fr(0), Fr(0)]                   # (1 + alpha)/2 in the power basis of alpha
        cols = [pmulmod(h, [Fr(int(i == j)) for i in range(4)], fcl) for j in range(4)]   # h * alpha^j
        Mh = [[cols[j][k] for j in range(4)] for k in range(4)]
        cp = charpoly_q(Mh)
        check(cp == [Fr(c) for c in gc],
              f'(b) the characteristic polynomial of (1 + alpha)/2 (exact, from its multiplication matrix in the basis'
              f' 1, alpha, alpha^2, alpha^3) is {fexpr(cp)} = g(t): monic with integer coefficients')
        xx = sp.Symbol('x')
        check(sp.expand(sp.resultant(fexpr(fcl), xx - (1 + t_) / 2, t_) - fexpr(gc, xx)) == 0,
              '    the same by sympy: Res_t(f_cl(t), x - (1 + t)/2) = prod_i (x - (1 + alpha_i)/2) = g(x)')
        print('    => (1 + alpha)/2 is an algebraic integer (a root of the monic integer polynomial g) in Q + Q alpha that')
        print('       is not in Z + Z alpha.  So O_K cap (Q + Q alpha) != Z + Z alpha, and the marker lemma does not hold')
        print('       for alpha: an element of O_K with vanishing alpha^2-, alpha^3-coordinates need not have the form')
        print('       n - m alpha with m, n integers, and no Z-basis of O_K has the shape 1, alpha, omega_2, omega_3 of')
        print('       section 3 (that shape forces O_K cap (Q + Q alpha) = Z + Z alpha).  Since Q + Q alpha = Q + Q beta,')
        print('       section 3 (omega_0 = 1, omega_1 = beta) gives O cap (Q + Q alpha) = Z + Z beta, which contains')
        print('       Z + Z alpha = {x0 + x1 beta : x1 even} with index 2 (a + b alpha = (a - b) + 2b beta).  So the field')
        print('       form G, with beta, is the one that is solved.')
        print('   (c) the bijection of coprime solutions, in steps')
        print("      Step 1. By (a), (m,n) -> (m',n') = (2m, n + m) [core.from_class] gives G(m',n') = F_cl(m,n), and for m'")
        print("              even (m',n') -> (m,n) = (m'/2, n' - m'/2) [core.to_class] gives F_cl(m,n) = G(m',n').  The two")
        print("              linear maps are inverse to each other ((2m)/2 = m, (n + m) - m = n; 2(m'/2) = m',")
        print("              (n' - m'/2) + m'/2 = n') and commute with (m,n) -> (-m,-n).")
        ok1 = all(to_class(key, *from_class(key, m, n)) == (m, n) for m in GRID for n in GRID)
        ok2 = all(from_class(key, *to_class(key, a, b)) == (a, b) for a in GRID for b in GRID if a % 2 == 0)
        ok3 = all((to_class(key, a, b) is None) == (a % 2 == 1) for a in GRID for b in GRID)
        ok4 = all(Fclass(key, m, n) == Fval(key, *from_class(key, m, n)) for m in GRID for n in GRID)
        ok5 = all(from_class(key, -m, -n) == tuple(-v for v in from_class(key, m, n)) for m in GRID for n in GRID)
        check(ok1 and ok2 and ok3 and ok4 and ok5,
              "step 1 on the code (|m|,|n| <= 40): to_class(from_class(m,n)) = (m,n); from_class(to_class(m',n')) ="
              " (m',n') for m' even; to_class(m',n') is None iff m' is odd; F_cl(m,n) = G(from_class(m,n));"
              " from_class(-m,-n) = -from_class(m,n)")
        print('      Step 2. If m and n are both odd, then F_cl(m,n) == 0 (mod 16): m + n = 2k and, by (a),')
        print('              F_cl(m,n) = G(2m, 2k) = 16 G(m,k).')
        check(sp.expand(form_expr(Gco, 2 * m_, 2 * n_) - 16 * Ge) == 0, 'step 2: G(2m, 2k) = 16 G(m,k) (sympy)')
        res16 = [(a, b) for a in range(1, 16, 2) for b in range(1, 16, 2)]
        check(all(Fclass(key, a, b) % 16 == 0 for a, b in res16),
              f'step 2 by residues: F_cl(m,n) == 0 (mod 16) at all {len(res16)} residue pairs (m,n) mod 16 with m, n'
              ' odd (F_cl has integer coefficients, so F_cl(m,n) mod 16 depends only on m, n mod 16)')
        print('      Step 3. No c in RHS is divisible by 16.  So a coprime solution (m,n) of F_cl(m,n) = c, c in RHS, has m')
        print('              and n of opposite parity: not both even (coprime), not both odd (step 2).')
        check(all(c % 16 != 0 for c in rhs_k), f'step 3: RHS = {rhs_k} is {[c % 16 for c in rhs_k]} mod 16; none is 0')
        print('      Step 4. Then n + m is odd, and gcd(2m, n + m) = 1: an odd prime dividing 2m and n + m divides m and')
        print('              n = (n + m) - m, which are coprime, and 2 does not divide n + m.  So (m\',n\') = (2m, n + m) is a')
        print("              coprime solution of G(m',n') = c with m' even.")
        check(all((a + b) % 2 == 1 for a, b in ((0, 1), (1, 0))),
              'step 4 by residues: (m, n) == (0, 1) or (1, 0) mod 2 gives n + m == 1 mod 2')
        ok = all(gcd(*from_class(key, m, n)) == (2 if (m % 2 and n % 2) else 1)
                 for m in GRID for n in GRID if gcd(m, n) == 1)
        check(ok, 'steps 2-4 on the code: for coprime (m,n), |m|,|n| <= 40, gcd(from_class(m,n)) is 1, or 2 exactly when'
                  ' m and n are both odd')
        print("      Step 5. Conversely, let (m',n') be a coprime solution of G(m',n') = c with m' even.  Then (m,n) =")
        print("              (m'/2, n' - m'/2) is integral; a common divisor of m and n divides m' = 2m and n' = n + m, so")
        print("              gcd(m,n) = 1; and F_cl(m,n) = G(m',n') = c by step 1.")
        ok = all(gcd(*to_class(key, a, b)) == 1 for a in GRID for b in GRID if a % 2 == 0 and gcd(a, b) == 1)
        check(ok, "step 5 on the code: for coprime (m',n') with m' even, |m'|,|n'| <= 40, to_class(m',n') is coprime")
        print('      Step 6. For every c in RHS, from_class is therefore a bijection from the coprime solutions of F_cl(m,n) = c')
        print("              onto the coprime solutions of G(m',n') = c with m' even, with inverse to_class.  It maps -(m,n)")
        print("              to -(m',n'), so it is one-to-one on pairs +-(m,n).  A coprime solution of G = c with m' odd")
        print('              corresponds to no solution of F_cl = c (to_class gives None); it gives no integral point.')
        print("              (Such a solution gives the coprime pair (m', 2n' - m') of odd numbers with F_cl(m', 2n' - m') =")
        print("              G(2m', 2n') = 16 G(m',n') = 16c, by (a): a solution of F_cl = 16c, and 16c is not in RHS by step 3.)")
        ok = all(Fclass(key, a, 2 * b - a) == 16 * Fval(key, a, b) and gcd(a, 2 * b - a) == gcd(a, b)
                 for a in GRID for b in GRID if a % 2)
        check(ok, "step 6, remark, on |m'|,|n'| <= 40 with m' odd: F_cl(m', 2n' - m') = 16 G(m',n') and"
                  " gcd(m', 2n' - m') = gcd(m',n')")
        print('      In the language of section 3: the solutions of F_cl(m,n) = c are the elements z of O with N(z) = c,')
        print('      markers x2 = x3 = 0 and x1 even (z = n - m alpha = (n + m) - 2m beta), with the coprimality of steps 3-5.')
    else:
        check(Fcl == Gco, '(a) G = F_cl (the same coefficients): the substitution is the identity, beta = alpha')
        print('   (b) no substitution is needed: by section 3 (omega_0 = 1, omega_1 = alpha) and section 2 (O = O_K),')
        print('       O_K cap (Q + Q alpha) = Z + Z alpha, so the marker lemma holds for the class form itself.')
        ok = all(from_class(key, m, n) == (m, n) and to_class(key, m, n) == (m, n) for m in GRID for n in GRID)
        check(ok, '(c) core.from_class and core.to_class are the identity (checked on |m|,|n| <= 40): the coprime'
                  ' solutions of F_cl(m,n) = c are those of G(m,n) = c, for every c in RHS, with no condition')

# ------------------------------------------------------------------------------------------------------------------
header('5. Units: N(eps1) = N(eps2) = 1 for core.UNITS')
for key in KEYS:
    K = Field(key)
    print(f'CLASS {LABEL[key]}')
    for nm, u in zip(('eps1', 'eps2'), UNITS[key]):
        nu_ = K.norm(u)
        cp = charpoly(K.matrix(u))
        check(nu_ == 1 and cp[0] == 1, f'N({nm}) = {nu_}   {nm} = {u} = {elt(K.to_power(u), beta_s)}'
                                        f'   power coordinates {[str(c) for c in K.to_power(u)]}'
                                        f'   charpoly {fexpr(cp)}')
        inv = K.pow(u, -1)
        check(K.mul(u, inv) == (1, 0, 0, 0), f'{nm}^(-1) = {inv} (the adjugate, integral coordinates) lies in O and'
                                              f' {nm} {nm}^(-1) = 1: {nm} is a unit of O')
    ar = mp_roots(fcoeffs(key))
    rr = sorted([mpmath.re(r) for r in ar if abs(mpmath.im(r)) < 1e-40])

    def logs(u, K=K, rr=rr):
        pc = K.to_power(u)
        return [mpmath.log(abs(sum(mpmath.mpf(c.numerator) / c.denominator * r ** k for k, c in enumerate(pc))))
                for r in rr]
    L1, L2 = logs(UNITS[key][0]), logs(UNITS[key][1])
    print(f'   (information, floats) log|sigma_1|, log|sigma_2| of eps1: {mpmath.nstr(L1[0], 10)}, '
          f'{mpmath.nstr(L1[1], 10)}; of eps2: {mpmath.nstr(L2[0], 10)}, {mpmath.nstr(L2[1], 10)};'
          f' |det| = {mpmath.nstr(abs(L1[0] * L2[1] - L1[1] * L2[0]), 12)}')
print('Independence of eps1, eps2 follows from the cone certificate of 01_coverage.py; fundamentality is not used.')

# ------------------------------------------------------------------------------------------------------------------
header('6. Known solutions: G-solutions, class-form images, integral points; confirmation box search |m|,|n| <= 1000')
print('A class-form solution (m,n) gives an integral point: (u,v,w) = conic map, g = gcd(u,v,w), F_cl(m,n) = g^2,')
print(f'X = X4(m,n)/g^2, Y = |Y6(m,n)|/g^3 on {CURVE_EQ} (core.descent_maps).  If (X, Y) = (6a + 3, 18b + 9), then')
print('b(b+1)/2 = a(a+1)(2a+1)/6: a triangular number equal to a square-pyramidal number, by the identity below.')
a_s, b_s = sp.symbols('a b')
ident = sp.expand((18 * b_s + 9) ** 2 - ((6 * a_s + 3) ** 3 + GA * (6 * a_s + 3) + GB)
                  - 648 * (b_s * (b_s + 1) / 2 - a_s * (a_s + 1) * (2 * a_s + 1) / 6))
check(GA == -9 and GB == 81 and ident == 0,
      'identity: (18b + 9)^2 - ((6a + 3)^3 - 9(6a + 3) + 81) = 648 (b(b+1)/2 - a(a+1)(2a+1)/6)')
POINTS = {}
VALUES = set()
CLASS_SOLS = {}
for key in KEYS:
    rhs_k = tuple(CLASSES[key]['RHS'])
    print(f'CLASS {LABEL[key]}   G = {CLASSES[key]["F"]}, F_cl = {CLASSES[key]["CLASS_F"]}')
    uvw, X4, Y6 = descent_maps(key)
    byc_G = {c: [] for c in rhs_k}
    byc_F = {c: [] for c in rhs_k}
    for (a, b) in KNOWN[key]:
        c = Fval(key, a, b)
        check(gcd(a, b) == 1 and c in rhs_k, f"G-solution (m',n') = {(a, b)}: coprime, G = {c}")
        if c not in byc_G:
            continue
        byc_G[c].append((a, b))
        img = to_class(key, a, b)
        if img is None:
            check(key == '10' and a % 2 == 1,
                  f"   to_class gives None (m' = {a} odd): a solution of G = {c} that gives no solution of F_cl = {c}"
                  ' and no integral point (section 4(c), step 6)')
            continue
        m, n = img
        cF = Fclass(key, m, n)
        check(gcd(m, n) == 1 and cF == c and from_class(key, m, n) == (a, b),
              f"   class-form image (m,n) = to_class(m',n') = {img}: coprime, F_cl = {cF}, from_class(m,n) = {(a, b)}")
        byc_F[c].append(img)
        u, v, w = uvw(m, n)
        g = gcd(gcd(u, v), w)
        ok = g * g == cF and X4(m, n) % (g * g) == 0 and Y6(m, n) % g ** 3 == 0
        X, Y = X4(m, n) // (g * g), abs(Y6(m, n)) // g ** 3
        ok = ok and Y * Y == X ** 3 + GA * X + GB
        txt = f'   (m,n) = {img}: (u,v,w) = {(u, v, w)}, g = {g}, F_cl = g^2, (X,Y) = ({X}, +-{Y}) on the curve'
        if X % 6 == 3:
            aa, bb = (X - 3) // 6, (Y - 9) // 18
            ok = ok and Y % 18 == 9 and 3 * bb * (bb + 1) == aa * (aa + 1) * (2 * aa + 1)
            val = bb * (bb + 1) // 2
            txt += f'; X = 6a + 3, Y = 18b + 9 with a = {aa}, b = {bb}: b(b+1)/2 = a(a+1)(2a+1)/6 = {val}'
            VALUES.add(val)
        else:
            txt += f'; X == {X % 6} (mod 6), not 3: no common value'
        POINTS[(key, m, n)] = (X, Y)
        check(ok, txt)
    CLASS_SOLS[key] = byc_F
    print("   G-solutions by c (Theorem A'):     " + '; '.join(f'c = {c}: {byc_G[c]}' for c in rhs_k))
    print('   F_cl-solutions by c (Theorem A):   ' + '; '.join(f'c = {c}: {byc_F[c]}' for c in rhs_k))
    check(byc_G == STATED_G[key], "core.KNOWN grouped by c is the stated table of Theorem A' for this class")
    check(byc_F == STATED_F[key], 'their class-form images grouped by c are the stated table of Theorem A for this class')
    xs = sorted(X for (k, _, _), (X, _) in POINTS.items() if k == key)
    check(xs == DESCENT_X[key], f'the X of these points, {xs}, are the integral points of class {LABEL[key]} listed in'
                                ' weierstrass_2descent.out.txt (CURVE 8604.a1)')
check(all(c % d ** 4 for c in RHS for d in range(2, isqrt(isqrt(max(RHS))) + 1)),
      f'no c in RHS = {RHS} is divisible by d^4 for d >= 2 (16 and 81 divide none): every integer solution of G = c'
      ' or F_cl = c, c in RHS, is coprime, so "all solutions" and "all coprime solutions" are the same lists')
# box search
R = 1000
print(f'box search: every coprime (m,n), |m|,|n| <= {R}, with G(m,n) in RHS, resp. F_cl(m,n) in RHS (numpy int64, then'
      ' exact recheck)')
nvec = np.arange(-R, R + 1, dtype=np.int64)
n2 = nvec * nvec


def box(co, evalfn, rhs_k):
    hits = []
    for m in range(-R, R + 1):
        mm = np.int64(m)
        Fv = (co[0] * mm ** 4 + co[1] * mm ** 3 * nvec + co[2] * mm * mm * n2 + co[3] * mm * n2 * nvec
              + co[4] * n2 * n2)
        idx = np.nonzero(np.isin(Fv, list(rhs_k)))[0]
        for i in idx:
            n = int(nvec[i])
            if gcd(m, n) == 1 and evalfn(m, n) in rhs_k:
                hits.append((m, n))
    return hits


for key in KEYS:
    rhs_k = tuple(CLASSES[key]['RHS'])
    images = [to_class(key, a, b) for (a, b) in KNOWN[key] if to_class(key, a, b) is not None]
    for name, co, evalfn, expected in (('G', CLASSES[key]['F'], lambda m, n, k=key: Fval(k, m, n), KNOWN[key]),
                                       ('F_cl', CLASSES[key]['CLASS_F'], lambda m, n, k=key: Fclass(k, m, n), images)):
        bound = sum(abs(c) for c in co) * R ** 4
        check(bound < 2 ** 62, f'class {LABEL[key]}, {name}: |{name}| <= {sum(abs(c) for c in co)} * {R}^4 = {bound}'
                               ' < 2^62 (no int64 overflow)')
        hits = box(co, evalfn, rhs_k)
        found = pm_set(hits)
        want = pm_set(expected)
        check(found == want, f'class {LABEL[key]}, {name}: {len(hits)} solutions = {len(found)} pairs +-(m,n);'
                             f' equal to {"core.KNOWN" if name == "G" else "the images to_class(core.KNOWN)"} up to'
                             f' sign: {found == want}')

# ------------------------------------------------------------------------------------------------------------------
header('7. The three quartic fields are pairwise non-isomorphic (not used by any proof)')
print('Certificate: a prime p dividing neither disc(g_1) nor disc(g_2) at which g_1 and g_2 factor with different degree')
print('patterns.  For such p (which does not divide [O_K : Z[beta]]) the pattern is the splitting type of p in the field')
print('(Dedekind), an invariant of the field.  The field discriminants disc(O_K) of section 2 are a second invariant.')
ISOS = []
for k1, k2 in itertools.combinations(KEYS, 2):
    cert = nonisomorphism_certificate(fcoeffs(k1), fcoeffs(k2))
    if cert is not None:
        check(True, f'class {LABEL[k1]} vs {LABEL[k2]}: at p = {cert[0]} (not dividing {DISC[k1]} * {DISC[k2]}) g'
                    f' factors with degrees {cert[1]} and {cert[2]}: not isomorphic')
    else:
        P = find_embeddings(fcoeffs(k2), fcoeffs(k1))
        if P:
            P = P[0]
            ok = not any(poly_at(fcoeffs(k2), P, fcoeffs(k1)))
            check(ok, f'class {LABEL[k1]} vs {LABEL[k2]}: ISOMORPHIC; g_{LABEL[k2]}(P(beta)) = 0 exactly modulo'
                      f' g_{LABEL[k1]} with P(beta) = {elt(P, beta_s)}')
            ISOS.append((k1, k2, P))
        else:
            check(False, f'class {LABEL[k1]} vs {LABEL[k2]}: no certificate up to p = 400 and no isomorphism found')
    same = DISC_O[k1] == DISC_O[k2]
    print(f'      field discriminants {DISC_O[k1]} and {DISC_O[k2]}: '
          + ('equal' if same else 'different (a second certificate of non-isomorphism)'))

# ------------------------------------------------------------------------------------------------------------------
header('SUMMARY')
for key in KEYS:
    print(f'class {LABEL[key]} (mu = {CLASSES[key]["mu"]}): G = {sp.expand(form_expr(CLASSES[key]["F"], m_, n_))},'
          f' g = {fexpr(fcoeffs(key))}; disc(g) = {DISC[key]} = {fz(DISC[key])};'
          f' [O : Z[beta]] = {INDEX[key]}; disc(O_K) = {DISC_O[key]} = {fz(DISC_O[key])}; {CLASSES[key]["SUB"]}')
print("Theorem A' (field forms), one of each pair +-(m',n').  Completeness is proved by 01_coverage.py and"
      ' 02_exclusion.py; this script checks the listed solutions and confirms them by the box search:')
for key in KEYS:
    print(f'   class {LABEL[key]}: ' + '; '.join(f'c = {c}: {STATED_G[key][c] or "none"}' for c in CLASSES[key]['RHS']))
print('Theorem A (class forms; from Theorem A\' by section 4(c)), one of each pair +-(m,n):')
for key in KEYS:
    print(f'   class {LABEL[key]}: ' + '; '.join(f'c = {c}: {STATED_F[key][c] or "none"}' for c in CLASSES[key]['RHS']))
allX = sorted({X for (X, _) in POINTS.values()})
print(f'integral points from the three non-trivial classes: X in {allX} ({2 * len(allX)} points, both signs of Y)')
ok = True
for X in TRIVIAL_X:
    Y2 = X ** 3 + GA * X + GB
    ok = ok and Y2 >= 0 and isqrt(Y2) ** 2 == Y2 and X % 6 != 3
check(ok, f'the trivial class (weierstrass_2descent.out.txt) has X in {TRIVIAL_X}: each X^3 - 9X + 81 is a square, and'
          ' no such X is 3 mod 6 (no common value from that class)')
allX_B = sorted(set(allX) | set(TRIVIAL_X))
check(allX_B == SEARCH_X, f'with the trivial class: X in {allX_B} ({2 * len(allX_B)} points), the integral points found'
                          ' by the search -92 <= x <= 10^6 of weierstrass_2descent.out.txt (equal to the LMFDB list);'
                          ' under hypothesis (H) there are no others (Theorem B)')
check(sorted(VALUES) == [0, 1, 55, 91, 208335], f'common values b(b+1)/2 = a(a+1)(2a+1)/6 from these points:'
                                                 f' {sorted(VALUES)}')
print(f'fields: pairwise non-isomorphic: {not ISOS}')
print(f'[done in {time.time() - T0:.1f} s]')
if FAILED:
    print(f'{len(FAILED)} CHECK FAILED')
    sys.exit(1)
print('all checks passed')
