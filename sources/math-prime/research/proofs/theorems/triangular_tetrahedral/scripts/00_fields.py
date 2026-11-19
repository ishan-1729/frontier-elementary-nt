"""00_fields.py -- the three quartic Thue fields of 1611.a1 (triangular = tetrahedral), checked exactly, and the
comparison with E. T. Avanesov, Acta Arith. 12 (1966/67) 409-420.

Everything below is exact (Python integers, Fractions, sympy).  Floats (mpmath, 60 digits) are used only to PROPOSE
candidates (an isomorphism, a substitution, a pair of unit exponents), which are then verified exactly.

For each class key '10', '01', '11' of core.CLASSES, with F(m,n) = a0 m^4 + a1 m^3 n + a2 m^2 n^2 + a3 m n^3 + n^4
and the monic f(t) = F(1,t):

 1. f is irreducible over Q: it has no integer root (a rational root of a monic integer polynomial is an integer
    dividing f(0)) and no factorization into monic integer quadratics (Gauss's lemma: none over Q either); sympy
    agrees.  F(m,n) = m^4 f(n/m) = Res_t(f(t), n - m t) = det(multiplication by n - m alpha) = N(n - m alpha) as
    polynomials.  disc(f) is factored.  f has exactly two real roots (own exact Sturm chain; sympy count_roots).
    Galois group S4: the resolvent cubic R(y) = y^3 - b y^2 + (ac - 4d) y - (a^2 d - 4bd + c^2) of
    f = t^4 + a t^3 + b t^2 + c t + d (whose roots are a1a2+a3a4, a1a3+a2a4, a1a4+a2a3: checked symbolically) has
    no rational root, and disc(f) < 0 is not a square.
 2. The basis of core.CLASSES spans an order O (core.Field asserts an integral multiplication table) containing
    Z[alpha] (the inverse basis matrix is integral); [O : Z[alpha]] = 1/|det B|; disc(O) = disc(f) det(B)^2, equal
    to det(Tr(omega_i omega_j)).
    MAXIMALITY LEMMA.  Let O be an order of K with Z-basis omega_0..omega_3 and O_K the maximal order.  If a prime p
    divides [O_K : O], then O_K/O has an element of order p: y in O_K \\ O with p y in O.  Subtracting an element of
    O, y = (sum c_i omega_i)/p with 0 <= c_i < p, not all 0.  Such a y is integral iff its characteristic polynomial
    lies in Z[t] (Cayley-Hamilton one way; the other way it is a power of the minimal polynomial of an algebraic
    integer), i.e. iff p^k divides the coefficient c_k of t^(4-k) in det(t - M(x)), x = sum c_i omega_i, k = 1..4.
    So if none of the p^4 - 1 classes gives an integral y, p does not divide the index.  Since
    disc(O) = [O_K : O]^2 disc(O_K), only primes with p^2 | disc(O) can divide the index; all of them are tested.
    Cross-check: the maximal order computed from Z[alpha] by the same test with enlargement has the same Hermite
    normal form as the basis of core.CLASSES.
 3. MARKER LEMMA.  omega_0 = 1, omega_1 = alpha, and the alpha^2-, alpha^3-power coordinates of x = sum x_i omega_i
    are (x2, x3) D with an invertible 2x2 matrix D.  As 1, alpha, alpha^2, alpha^3 are linearly independent
    (f irreducible), x lies in Q + Q alpha iff x2 = x3 = 0, and then x = x0 + x1 alpha with x0, x1 in Z.  The rows
    0, 1 of the inverse basis matrix are (1,0,0,0), (0,1,0,0), so n - m alpha has coordinates (n, -m, 0, 0).
 4. N(eps1) = N(eps2) = 1 for core.UNITS (exact determinants).
 5. The known solutions core.KNOWN grouped by the right side c; the integral point and the original pair of each
    (core.descent_maps); and a confirmation box search |m|,|n| <= 1000 (not used by any proof).
 6. Avanesov: (a) explicit isomorphisms of his three quartic fields with ours, and certificates of
    non-isomorphism (different factorization patterns mod a prime p not dividing either discriminant, Dedekind);
    (b) the maximal orders and discriminants of his fields; (c) GL2(Q) substitutions F(T(N,M)) = lambda L_i(N,M)
    for his six forms, his solutions mapped to core.KNOWN, and the congruences that make L_i = 1 equivalent to
    F = c; (d) his units: identities, norms, transport into our field, comparison with core.UNITS;
    (e) further statements of the paper (Skolem data, the cubic field, lemma solution lists by box search);
    (f) which equations F = c can contain points of his model 6Y^2 = X^3 - 4X + 6.

Exit status: non-zero iff some CHECK FAILED (a check of our own data).  A statement of Avanesov's paper that is
tested is reported as 'AVANESOV ...: TRUE/FALSE'; FALSE there is a finding about the paper, not a failure.
"""
import os
import sys
import time
import itertools
from fractions import Fraction as Fr
from math import gcd, isqrt

import numpy as np
import sympy as sp
import mpmath

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from core import (CLASSES, RHS, KNOWN, UNITS, Field, Fval, fcoeffs, pmulmod, mat_inv, vec_mat, det_int,  # noqa: E402
                  descent_maps)

T0 = time.time()
mpmath.mp.dps = 60
KEYS = ('10', '01', '11')
LABEL = {k: CLASSES[k]['label'] for k in KEYS}
FAILED = []
CLAIMS = []
t_, m_, n_, N_, M_, y_ = sp.symbols('t m n N M y')


def check(cond, msg):
    if cond:
        print(f'   [ok] {msg}')
    else:
        print(f'   CHECK FAILED: {msg}')
        FAILED.append(msg)
    return cond


def claim(label, cond, detail=''):
    print(f'   AVANESOV {label}: {"TRUE" if cond else "FALSE"}' + (f'   [{detail}]' if detail else ''))
    CLAIMS.append((label, bool(cond)))


def header(s):
    print()
    print('=' * 116)
    print(s)
    print('=' * 116)


# ------------------------------------------------------------------------------------------------------------------
# small exact helpers
# ------------------------------------------------------------------------------------------------------------------
def fexpr(fc, var=t_):
    """sympy expression of the polynomial with coefficients fc (low degree first)."""
    return sum(sp.Integer(c) * var ** k for k, c in enumerate(fc))


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
    return True, 'no integer root among the divisors of f(0); no quadratic factorization over Z'


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
    verification h(P(alpha)) = 0.  Returns the list of power-coordinate vectors P."""
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


def mobius_relations(Fco, Gco):
    """All primitive integer T = [[a,b],[c,d]] and rational lam with F(aN+bM, cN+dM) = lam G(N,M) (verified
    exactly).  Candidates: Mobius maps sending the roots n/m of F to the roots N/M of G (all 24 bijections)."""
    ar = mp_roots(list(Fco))                    # F(1,t) = sum Fco[i] t^i (low first): the roots n/m
    xr = mpmath.polyroots(list(Gco), maxsteps=800, extraprec=800)   # G(x,1), high first
    out = []
    Fe = form_expr(Fco, m_, n_)
    Ge = form_expr(Gco, N_, M_)
    for perm in itertools.permutations(range(4)):
        rows = [[ar[j], 1, -ar[j] * xr[perm[j]], -xr[perm[j]]] for j in range(4)]
        det4 = mpmath.det(mpmath.matrix(rows))
        scale = 1
        for r in rows:
            scale *= max(abs(x) for x in r)
        if abs(det4) > mpmath.mpf(10) ** -40 * scale:
            continue
        v = [(-1) ** k * mpmath.det(mpmath.matrix([[rows[i][c] for c in range(4) if c != k] for i in range(3)]))
             for k in range(4)]
        kk = max(range(4), key=lambda c: abs(v[c]))
        v = [x / v[kk] for x in v]
        if max(abs(mpmath.im(x)) for x in v) > 1e-35:
            continue
        q = [rationalize(mpmath.re(x), 10 ** 6) for x in v]
        D = 1
        for x in q:
            D = D * x.denominator // gcd(D, x.denominator)
        ints = [int(x * D) for x in q]
        g = 0
        for x in ints:
            g = gcd(g, x)
        A, B, C, Dd = (x // g for x in ints)
        T = [[-C, A], [Dd, -B]]                  # (m, n) = T (N, M)
        if T[0][0] < 0 or (T[0][0] == 0 and T[0][1] < 0):
            T = [[-x for x in row] for row in T]
        lam = Fr(Fval_co(Fco, T[0][0], T[1][0]), Gco[0])
        sub = Fe.subs({m_: T[0][0] * N_ + T[0][1] * M_, n_: T[1][0] * N_ + T[1][1] * M_}, simultaneous=True)
        if sp.expand(sub - sp.Rational(lam.numerator, lam.denominator) * Ge) == 0 and (T, lam) not in out:
            out.append((T, lam))
    return out


def Fval_co(co, m, n):
    return sum(c * m ** (4 - i) * n ** i for i, c in enumerate(co))


def mn_str(T):
    def lin(a, b):
        s = []
        for c, v in ((a, 'N'), (b, 'M')):
            if c == 0:
                continue
            s.append(f'{"+" if c > 0 else "-"} {"" if abs(c) == 1 else abs(c)}{v}')
        r = ' '.join(s)
        return r[2:] if r.startswith('+ ') else '-' + r[2:]
    return f'(m,n) = ({lin(T[0][0], T[0][1])}, {lin(T[1][0], T[1][1])})'


# ==================================================================================================================
print('00_fields.py -- the quartic Thue fields of 1611.a1 (C(n,2) = C(m,3)) and Avanesov\'s fields; exact checks')
print('conventions: F(m,n) = a0 m^4 + a1 m^3 n + a2 m^2 n^2 + a3 m n^3 + n^4, f(t) = F(1,t), alpha a root of f;')
print('             an element is given by its coordinates in the integral basis omega_0..omega_3 of core.CLASSES.')

# ------------------------------------------------------------------------------------------------------------------
header('1. f irreducible; F(m,n) = m^4 f(n/m) = N(n - m alpha); disc(f); two real roots; Galois group S4')
x1, x2, x3, x4 = sp.symbols('x1:5')
e1, e2 = x1 + x2 + x3 + x4, x1 * x2 + x1 * x3 + x1 * x4 + x2 * x3 + x2 * x4 + x3 * x4
e3, e4 = x1 * x2 * x3 + x1 * x2 * x4 + x1 * x3 * x4 + x2 * x3 * x4, x1 * x2 * x3 * x4
A_, B_, C_, D_ = -e1, e2, -e3, e4
lhs = (y_ - (x1 * x2 + x3 * x4)) * (y_ - (x1 * x3 + x2 * x4)) * (y_ - (x1 * x4 + x2 * x3))
rhs = y_ ** 3 - B_ * y_ ** 2 + (A_ * C_ - 4 * D_) * y_ - (A_ * A_ * D_ - 4 * B_ * D_ + C_ * C_)
check(sp.expand(lhs - rhs) == 0,
      'resolvent formula: prod over the 3 pairings of (y - (x_i x_j + x_k x_l)) = y^3 - b y^2 + (ac - 4d) y'
      ' - (a^2 d - 4bd + c^2) for f = prod (t - x_i) = t^4 + a t^3 + b t^2 + c t + d (symbolic identity)')
DISC = {}
for key in KEYS:
    F = CLASSES[key]['F']
    fc = fcoeffs(key)
    K = Field(key)
    print()
    print(f'CLASS {LABEL[key]}  (mu = {CLASSES[key]["mu"]})   F(m,n) = {sp.expand(form_expr(F, m_, n_))}')
    print(f'   f(t) = F(1,t) = {fexpr(fc)}')
    ok, why = quartic_irreducible_elementary(fc)
    check(ok, f'f irreducible over Q (elementary: {why})')
    check(sp.Poly(fexpr(fc), t_).is_irreducible, 'f irreducible over Q (sympy)')
    Fe = sp.expand(form_expr(F, m_, n_))
    check(sp.expand(sp.cancel(m_ ** 4 * fexpr(fc).subs(t_, n_ / m_))) == Fe, 'm^4 f(n/m) = F(m,n)')
    res = sp.expand(sp.resultant(fexpr(fc), n_ - m_ * t_, t_))
    check(res == Fe, 'Res_t(f(t), n - m t) = F(m,n)   [= prod_i (n - m alpha_i) = N(n - m alpha), f monic]')
    x = (n_, -m_, 0, 0)
    Msym = sp.Matrix(4, 4, lambda k, j: sum(x[i] * K.mt[i][j][k] for i in range(4)))
    check(sp.expand(Msym.det()) == Fe,
          'det(matrix of multiplication by n - m alpha = (n,-m,0,0) in the basis of core.Field) = F(m,n)')
    disc = int(sp.discriminant(fexpr(fc), t_))
    DISC[key] = disc
    print(f'   disc(f) = {disc} = {fz(disc)}')
    nreal, chainlen, lastc = sturm_count(fc)
    check(lastc, f'Sturm chain of length {chainlen} ends in a nonzero constant (f squarefree)')
    check(nreal == 2, f'exactly two real roots (exact Sturm count = {nreal})')
    check(sp.Poly(fexpr(fc), t_).count_roots() == 2, 'exactly two real roots (sympy count_roots)')
    rts = mp_roots(fc)
    rr = sorted([mpmath.re(r) for r in rts if abs(mpmath.im(r)) < 1e-40])
    cr = [r for r in rts if mpmath.im(r) > 1e-40]
    print(f'   roots (floats, for information): real {mpmath.nstr(rr[0], 15)}, {mpmath.nstr(rr[1], 15)};'
          f' complex {mpmath.nstr(cr[0], 15)} and its conjugate')
    R = resolvent(fc)
    print(f'   resolvent cubic R(y) = {fexpr(R, y_)}')
    roots_R = [r for r in divisors_pm(R[0]) if sum(c * r ** k for k, c in enumerate(R)) == 0] if R[0] != 0 else [0]
    check(not roots_R, f'R has no rational root (a rational root is an integer dividing R(0) = {R[0]}: none of'
                       f' the {len(divisors_pm(R[0]))} divisors is a root) -> R irreducible over Q')
    check(int(sp.discriminant(fexpr(R, y_), y_)) == disc, 'disc(R) = disc(f)')
    check(disc < 0, 'disc(f) < 0, so disc(f) is not a square in Q')
    print('   => Galois group of f is S4 (f irreducible, R irreducible, disc(f) not a square)')

# ------------------------------------------------------------------------------------------------------------------
header('2. The integral basis: an order containing Z[alpha]; index; discriminant; maximality')
DISC_O = {}
for key in KEYS:
    fc = fcoeffs(key)
    K = Field(key)
    print()
    print(f'CLASS {LABEL[key]}   basis (power coordinates of omega_i):')
    for i, row in enumerate(K.B):
        print(f'      omega_{i} = {elt(row, sp.Symbol("alpha"))}')
    check(True, 'core.Field built: every omega_i omega_j has integral coordinates (core asserts it), so the'
                ' Z-span O of the basis is a ring')
    print('   multiplication table (coordinates of omega_i omega_j, i <= j):')
    for i in range(4):
        print('      ' + '   '.join(f'w{i}w{j} = {tuple(K.mt[i][j])}' for j in range(i, 4)))
    check(all(c.denominator == 1 for row in K.Binv for c in row),
          'the inverse basis matrix is integral: 1, alpha, alpha^2, alpha^3 lie in O, so Z[alpha] <= O')
    detB, discO, trdet = disc_of_basis(fc, K.B)
    index = 1 / abs(detB)
    check(index.denominator == 1, f'[O : Z[alpha]] = 1/|det B| = {index}')
    check(discO == trdet, f'disc(O) = disc(f) det(B)^2 = {discO} equals det(Tr(omega_i omega_j)) = {trdet}')
    discO = int(discO)
    DISC_O[key] = discO
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
    Bmax, primes, _ = maximal_order(fc, log)
    for line in log:
        print(f'      from Z[alpha]: {line}')
    check(hnf_q(Bmax) == hnf_q(K.B), 'the maximal order computed from Z[alpha] (enlargement) has the same HNF as'
                                     ' the basis of core.CLASSES')

# ------------------------------------------------------------------------------------------------------------------
header('3. Marker lemma: x in O lies in Z + Z alpha iff x2 = x3 = 0; n - m alpha = (n, -m, 0, 0)')
for key in KEYS:
    K = Field(key)
    B = K.B
    print(f'CLASS {LABEL[key]}')
    check(B[0] == [1, 0, 0, 0] and B[1] == [0, 1, 0, 0], 'omega_0 = 1, omega_1 = alpha')
    check(all(B[i][j] == 0 for i in (0, 1) for j in (2, 3)),
          'only omega_2, omega_3 involve alpha^2, alpha^3: the (alpha^2, alpha^3) power coordinates of'
          ' x = sum x_i omega_i are (x2, x3) D')
    D = [[B[2][2], B[2][3]], [B[3][2], B[3][3]]]
    detD = D[0][0] * D[1][1] - D[0][1] * D[1][0]
    check(detD != 0, f'D = [[{D[0][0]}, {D[0][1]}], [{D[1][0]}, {D[1][1]}]] has det {detD} != 0, so the alpha^2,'
                     f' alpha^3 coordinates vanish iff x2 = x3 = 0; then x = x0 + x1 alpha, x0, x1 in Z')
    check(K.Binv[0] == [1, 0, 0, 0] and K.Binv[1] == [0, 1, 0, 0],
          'rows 0, 1 of the inverse basis matrix are (1,0,0,0), (0,1,0,0): n - m alpha has coordinates (n,-m,0,0)')
    ok = all(K.from_power([n, -m, 0, 0]) == [n, -m, 0, 0] for m in range(-3, 4) for n in range(-3, 4))
    check(ok, 'from_power(n - m alpha) = (n, -m, 0, 0) on the 49 pairs |m|,|n| <= 3 (spot check of the code)')

# ------------------------------------------------------------------------------------------------------------------
header('4. Units: N(eps1) = N(eps2) = 1 for core.UNITS')
alpha_s = sp.Symbol('alpha')
for key in KEYS:
    K = Field(key)
    print(f'CLASS {LABEL[key]}')
    for nm, u in zip(('eps1', 'eps2'), UNITS[key]):
        nu_ = K.norm(u)
        cp = charpoly(K.matrix(u))
        check(nu_ == 1 and cp[0] == 1, f'N({nm}) = {nu_}   {nm} = {u} = {elt(K.to_power(u), alpha_s)}'
                                        f'   power coordinates {[str(c) for c in K.to_power(u)]}'
                                        f'   charpoly {fexpr(cp)}')
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

# ------------------------------------------------------------------------------------------------------------------
header('5. Known solutions; their integral points; confirmation box search |m|,|n| <= 1000')
print('integral point of (m,n): (u,v,w) = conic map, g = gcd(u,v,w), X = X4(m,n)/g^2, Y = |Y6(m,n)|/g^3 on')
print('Y^2 = X^3 - 144X + 1296 (core.descent_maps); Avanesov\'s model 6Y\'^2 = X\'^3 - 4X\' + 6 is X = 6X\', Y = 36Y\','
      ' and X\' = 2x + 2, Y\' = 2y + 1 with 3y(y+1) = x(x+1)(x+2), common value C(x+2,3) = C(y+1,2).')
POINTS = {}
for key in KEYS:
    print(f'CLASS {LABEL[key]}')
    uvw, X4, Y6 = descent_maps(key)
    byc = {c: [] for c in RHS}
    for (m, n) in KNOWN[key]:
        c = Fval(key, m, n)
        check(gcd(m, n) == 1 and c in RHS, f'(m,n) = {(m, n)}: coprime, F = {c}')
        byc[c].append((m, n))
        u, v, w = uvw(m, n)
        g = gcd(gcd(u, v), w)
        ok = (g * g == c and X4(m, n) % (g * g) == 0 and Y6(m, n) % g ** 3 == 0)
        X, Y = X4(m, n) // (g * g), abs(Y6(m, n)) // g ** 3
        ok = ok and Y * Y == X ** 3 - 144 * X + 1296
        txt = f'(m,n) = {(m, n)}: g = {g}, F = g^2, (X,Y) = ({X}, {Y}) on the monic model'
        if X % 6 == 0 and Y % 36 == 0:
            Xa, Ya = X // 6, Y // 36
            xx, yy = (Xa - 2) // 2, (Ya - 1) // 2
            val = xx * (xx + 1) * (xx + 2) // 6
            ok = ok and (Xa - 2) % 2 == 0 and (Ya - 1) % 2 == 0 and 3 * yy * (yy + 1) == xx * (xx + 1) * (xx + 2)
            txt += f"; Avanesov X' = {Xa}, Y' = {Ya}; x = {xx}, y = {yy}, value {val}"
        else:
            txt += "; not on Avanesov's model (6 does not divide X or 36 does not divide Y)"
        POINTS[(key, m, n)] = (X, Y)
        check(ok, txt)
    print('   grouped by c: ' + '; '.join(f'c = {c}: {byc[c]}' for c in RHS))
# box search
R = 1000
print(f'box search: every coprime (m,n), |m|,|n| <= {R}, with F(m,n) in {RHS} (numpy int64, then exact recheck)')
for key in KEYS:
    a = CLASSES[key]['F']
    bound = sum(abs(c) for c in a) * R ** 4
    check(bound < 2 ** 62, f'class {LABEL[key]}: |F| <= {sum(abs(c) for c in a)} * {R}^4 = {bound} < 2^62'
                           ' (no int64 overflow)')
    nvec = np.arange(-R, R + 1, dtype=np.int64)
    n2 = nvec * nvec
    hits = []
    for m in range(-R, R + 1):
        mm = np.int64(m)
        Fv = (a[0] * mm ** 4 + a[1] * mm ** 3 * nvec + a[2] * mm * mm * n2 + a[3] * mm * n2 * nvec + a[4] * n2 * n2)
        idx = np.nonzero(np.isin(Fv, RHS))[0]
        for i in idx:
            n = int(nvec[i])
            if gcd(m, n) == 1 and Fval(key, m, n) in RHS:
                hits.append((m, n))
    found = {frozenset({(m, n), (-m, -n)}) for (m, n) in hits}
    known = {frozenset({(m, n), (-m, -n)}) for (m, n) in KNOWN[key]}
    check(found == known, f'class {LABEL[key]}: {len(hits)} solutions = {len(found)} pairs +-(m,n); equal to'
                          f' core.KNOWN up to sign: {found == known}')

# ------------------------------------------------------------------------------------------------------------------
header("6. Avanesov, Acta Arith. 12 (1966/67) 409-420: fields, discriminants, forms, units")
print('His quartic fields, labelled H1, H2, H3 here (in the paper H1 is K1(eta), and K2(lambda) is the CUBIC field):')
print('H1: eta^4 - 36 eta + 36 = 0 (his (3)); H2: eta^4 - 12 eta^3 + 18 eta^2 + 9 = 0, H3: eta^4 - 18 eta^2 - 36 eta + 9 = 0 (his Remark, p. 413).')
print('His forms (his variables (n, m) written here (N, M)); his lemmas 2-7 state that L_i(N,M) = 1 has exactly the')
print('listed solutions, up to sign:')
HIS_F = {'H1': [36, -36, 0, 0, 1], 'H2': [9, 0, 18, -12, 1], 'H3': [9, -36, -18, 0, 1]}
HIS_FORMS = {   # coefficients of N^4, N^3 M, N^2 M^2, N M^3, M^4; his solution list (N, M)
    'L2': ((1, 0, 0, 36, 36), [(1, 0), (1, -1)]),
    'L3': ((9, 0, 0, 12, 4), [(1, -1)]),
    'L4': ((1, -12, 18, 0, 9), [(1, 0), (2, 1)]),
    'L5': ((1, 0, 18, 36, 9), [(1, 0)]),
    'L6': ((1, -12, -18, 0, 9), [(1, 0)]),
    'L7': ((1, 0, -18, -36, 9), [(1, 0)]),
}
for nm, (co, sols) in HIS_FORMS.items():
    print(f'   {nm}(N,M) = {sp.expand(form_expr(co, N_, M_))} = 1;  solutions {sols}')
    check(all(Fval_co(co, a, b) == 1 for a, b in sols), f'{nm} = 1 at each listed solution (substitution)')

print()
print('(a) which of our fields each of his polynomials defines')
for hk, hc in HIS_F.items():
    print(f'   {hk}: h(t) = {fexpr(hc)}; irreducible {sp.Poly(fexpr(hc), t_).is_irreducible};'
          f' real roots {sp.Poly(fexpr(hc), t_).count_roots()}; disc(h) = {fz(int(sp.discriminant(fexpr(hc), t_)))}')
claim('p. 410: (3) has two real roots and a pair of complex conjugate roots',
      sp.Poly(fexpr(HIS_F['H1']), t_).count_roots() == 2)
ISO = {}
for hk, hc in HIS_F.items():
    for key in KEYS:
        fc = fcoeffs(key)
        P = find_embeddings(hc, fc)
        if P:
            check(len(P) == 1, f'{hk} -> class {LABEL[key]}: exactly one root of h in Q(alpha) found')
            P = P[0]
            K = Field(key)
            co = K.from_power(P)
            check(not any(poly_at(hc, P, fc)), f'{hk} = class {LABEL[key]}: h(P(alpha)) = 0 exactly modulo f, with'
                                                f' P(alpha) = {elt(P, alpha_s)}  (coordinates {[int(c) for c in co]})')
            # inverse: alpha as a polynomial in beta = P(alpha)
            W = [[Fr(1), Fr(0), Fr(0), Fr(0)]]
            for k in range(1, 4):
                W.append(pmulmod(W[-1], P, fc))
            Q = vec_mat([Fr(0), Fr(1), Fr(0), Fr(0)], mat_inv(W))
            check(not any(poly_at(fc, Q, hc)), f'inverse: alpha = Q(eta) with Q(eta) = {elt(Q, sp.Symbol("eta"))};'
                                                f' f(Q(eta)) = 0 exactly modulo h')
            ISO[hk] = (key, P, Q)
        else:
            cert = nonisomorphism_certificate(hc, fc)
            check(cert is not None, f'{hk} is not class {LABEL[key]}: at p = {cert[0]} (not dividing either'
                                    f' discriminant) h and f factor with degrees {cert[1]} and {cert[2]}'
                                    if cert else f'{hk} vs class {LABEL[key]}: no certificate found')
check(sorted(v[0] for v in ISO.values()) == sorted(KEYS), 'his three fields are our three fields, one each: '
      + ', '.join(f'{hk} = class {LABEL[v[0]]}' for hk, v in ISO.items()))
# the same isomorphisms as Mobius maps (they come from the substitutions of (c)): (c alpha + d) eta = a alpha + b
for hk, (num, den, txt) in {'H1': ([6, 0], [0, 1], 'eta = 6/alpha'),
                            'H2': ([3, 0], [1, 1], 'eta = 3/(alpha + 1)'),
                            'H3': ([3, 3], [1, -1], 'eta = 3(1 + alpha)/(1 - alpha)')}.items():
    if hk in ISO:
        key, P, _ = ISO[hk]
        lhs_ = pmulmod([Fr(den[0]), Fr(den[1]), Fr(0), Fr(0)], P, fcoeffs(key))
        check(lhs_ == [Fr(num[0]), Fr(num[1]), Fr(0), Fr(0)], f'{hk} = class {LABEL[key]}: equivalently {txt}'
                                                              ' (exact, modulo f)')
print('   our three fields are pairwise non-isomorphic:')
for k1, k2 in itertools.combinations(KEYS, 2):
    cert = nonisomorphism_certificate(fcoeffs(k1), fcoeffs(k2))
    check(cert is not None, f'class {LABEL[k1]} vs {LABEL[k2]}: at p = {cert[0]} the factorization degrees are'
                            f' {cert[1]} and {cert[2]}' if cert else f'class {LABEL[k1]} vs {LABEL[k2]}: no certificate')
distinct = True
for h1k, h2k in itertools.combinations(HIS_F, 2):
    cert = nonisomorphism_certificate(HIS_F[h1k], HIS_F[h2k])
    distinct = distinct and cert is not None
    print(f'      {h1k} vs {h2k}: p = {cert[0]}, degrees {cert[1]} vs {cert[2]}' if cert else f'      {h1k} vs {h2k}: ?')
claim('Remark (p. 413), first half: the three fields are distinct', distinct)

print()
print('(b) maximal orders of his fields (p-maximality test with enlargement, from Z[eta])')
HIS_DISC = {}
for hk, hc in HIS_F.items():
    log = []
    Bmax, primes, dh = maximal_order(hc, log)
    for line in log:
        print(f'   {hk}: {line}')
    detB, discO, trdet = disc_of_basis(hc, Bmax)
    check(discO == trdet, f'{hk}: disc = disc(h) det(B)^2 = det(trace form) = {discO}')
    HIS_DISC[hk] = int(discO)
    print(f'   {hk}: maximal order basis ' + ', '.join(str(elt(r, sp.Symbol('eta'))) for r in Bmax)
          + f';  [O : Z[eta]] = {1 / abs(detB)};  disc = {int(discO)} = {fz(int(discO))}')
    check(int(discO) == DISC_O[ISO[hk][0]], f'{hk}: equals disc(O) of the isomorphic class {LABEL[ISO[hk][0]]}'
                                            f' ({DISC_O[ISO[hk][0]]})')
K1basis = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, Fr(1, 6), 0], [0, 0, 0, Fr(1, 6)]]
log = []
B1, _, _ = maximal_order(HIS_F['H1'], log)
claim('p. 410: [1, eta, eta^2/6, eta^3/6] is an integral basis of K1(eta) = H1', hnf_q(B1) == hnf_q(K1basis))
target = -2 ** 4 * 3 ** 2 * 179
claim('p. 410: D(eta) = -2^4 3^2 179 for K1(eta) = H1', HIS_DISC['H1'] == target, f'computed {HIS_DISC["H1"]}')
claim('Remark (p. 413), second half: all three fields have discriminant -2^4 3^2 179 = -25776',
      all(v == target for v in HIS_DISC.values()), ', '.join(f'{k}: {v} = {fz(v)}' for k, v in HIS_DISC.items()))
B2 = [[1, 0, 0, 0], [0, 1, 0, 0], [Fr(1, 2), 0, Fr(1, 6), 0], [Fr(1, 4), Fr(1, 4), Fr(1, 12), Fr(1, 12)]]
K2o = GenField(HIS_F['H2'], B2, 'H2')                      # asserts an integral multiplication table
d2 = disc_of_basis(HIS_F['H2'], B2)
y2 = (0, 0, 0, 1)
check(d2[1] == d2[2] == -6444,
      'certificate for H2: O = Z + Z eta + Z (eta^2 + 3)/6 + Z (eta^3 + eta^2 + 3 eta + 3)/12 is an order (integral'
      f' multiplication table) with disc(O) = {d2[1]}; as disc(O) = [O_K : O]^2 disc(O_K), disc(H2) divides -6444,'
      f' so it is not -25776; the last generator has characteristic polynomial {fexpr(charpoly(K2o.matrix(y2)))}')
claim('p. 410: the cubic resolvent of (3) is z^3 - 144 z - 1296', resolvent(HIS_F['H1']) == [-1296, -144, 0, 1],
      f'computed {fexpr(resolvent(HIS_F["H1"]), sp.Symbol("z"))}')

print()
print('(c) his six forms as GL2(Q)-transforms of ours: F(T(N,M)) = lambda L_i(N,M), T primitive integral')
RELS = {}
for nm, (co, sols) in HIS_FORMS.items():
    found = []
    for key in KEYS:
        for (T, lam) in mobius_relations(CLASSES[key]['F'], co):
            found.append((key, T, lam))
    check(len(found) == 1, f'{nm}: exactly one relation found (over the three classes, all 24 root bijections)')
    if not found:
        continue
    key, T, lam = found[0]
    detT = T[0][0] * T[1][1] - T[0][1] * T[1][0]
    hk_of = [h for h, v in ISO.items() if v[0] == key]
    inv = [[sp.Rational(T[1][1], detT), sp.Rational(-T[0][1], detT)], [sp.Rational(-T[1][0], detT),
                                                                         sp.Rational(T[0][0], detT)]]
    print(f'   {nm}: class {LABEL[key]} (his field {hk_of[0] if hk_of else "?"}):  F_{LABEL[key]}(T(N,M)) = {lam} *'
          f' {nm}(N,M)   with {mn_str(T)}, det T = {detT}   (identity of forms: sympy expand, exact);'
          f'  inverse (N,M) = ({inv[0][0] * m_ + inv[0][1] * n_}, {inv[1][0] * m_ + inv[1][1] * n_})')
    cs = set()
    for (Nn, Mm) in sols:
        mm, nn = T[0][0] * Nn + T[0][1] * Mm, T[1][0] * Nn + T[1][1] * Mm
        g = gcd(mm, nn)
        mm, nn = mm // g, nn // g
        c = Fval(key, mm, nn)
        cs.add((c, g))
        inK = (mm, nn) in KNOWN[key] or (-mm, -nn) in KNOWN[key]
        check(inK and c * g ** 4 == lam, f'his solution (N,M) = {(Nn, Mm)} -> T(N,M) = {g} * {(mm, nn)};'
                                          f' F{(mm, nn)} = {c} = lambda/g^4; in core.KNOWN (up to sign): {inK}')
    # converse: our known solutions of the class that are his solutions
    adj = [[T[1][1], -T[0][1]], [-T[1][0], T[0][0]]]
    back = []
    for (mm, nn) in KNOWN[key]:
        Nn, Mm = adj[0][0] * mm + adj[0][1] * nn, adj[1][0] * mm + adj[1][1] * nn
        g = gcd(Nn, Mm)
        Nn, Mm = Nn // g, Mm // g
        if Fval_co(co, Nn, Mm) == 1:
            back.append(((mm, nn), (Nn, Mm)))
    print(f'      core.KNOWN of class {LABEL[key]} giving solutions of {nm} = 1, (m,n) -> (N,M) = primitive part of'
          f' adj(T)(m,n): ' + ', '.join(f'{a} -> {b}' for a, b in back))
    check(len(back) == len(sols), f'{nm}: the known solutions of the class account for exactly his {len(sols)}'
                                  f' solution(s)')
    # forced congruence: F = c (coprime) implies det T | g adj(T)(m,n), i.e. (m,n) = T(N,M)/g with L_i(N,M) = 1
    check(len(cs) == 1, f'{nm}: all his solutions give the same right side c and gcd g: {sorted(cs)}')
    if len(cs) == 1:
        (c, g), = cs
        RELS[nm] = (key, T, lam, c)
        forced_at = None
        for Q in sorted({abs(detT), abs(detT) * 2, abs(detT) * 3, abs(detT) ** 2, abs(detT) ** 3}):
            if Q > 1000:
                continue
            primes_Q = [p for p in sp.primefactors(Q)]
            okQ = True
            for a in range(Q):
                for b in range(Q):
                    if any(a % p == 0 and b % p == 0 for p in primes_Q):
                        continue
                    if (Fval(key, a, b) - c) % Q:
                        continue
                    if (g * (adj[0][0] * a + adj[0][1] * b)) % detT or (g * (adj[1][0] * a + adj[1][1] * b)) % detT:
                        okQ = False
                        break
                if not okQ:
                    break
            if okQ:
                forced_at = Q
                break
        # and g is forced for the solutions of L_i = 1: g = gcd(T(N,M), det T) depends on (N,M) mod det T
        gs = set()
        Qd = abs(detT)
        for a in range(max(Qd, 1)):
            for b in range(max(Qd, 1)):
                if Qd > 1 and (Fval_co(co, a, b) - 1) % Qd:
                    continue
                mm, nn = T[0][0] * a + T[0][1] * b, T[1][0] * a + T[1][1] * b
                gs.add(gcd(gcd(mm, nn), Qd) if Qd > 1 else 1)
        check(forced_at is not None and gs == {g},
              f'{nm} = 1  <=>  F_{LABEL[key]}(m,n) = {c} with (m,n) coprime: every residue (m,n) mod {forced_at}'
              f' with F == {c}, primitive, satisfies det T | {g} adj(T)(m,n); and L_i = 1 forces'
              f' gcd(T(N,M)) = {g} (residues mod {Qd})')
print('   our equations F = c not among his six:')
cov = {(v[0], v[3]) for v in RELS.values()}
for key in KEYS:
    for c in RHS:
        if (key, c) not in cov:
            ks = [mn for mn in KNOWN[key] if Fval(key, *mn) == c]
            print(f'      class {LABEL[key]}, c = {c}: known solutions {ks}'
                  + (f' (integral points (X,Y) = {[POINTS[(key,) + mn] for mn in ks]})' if ks else ''))

print()
print('(d) his units of K1(eta) = H1 (proof of Lemma 2)')
eta_s = sp.Symbol('eta')
hc1 = HIS_F['H1']
K1 = GenField(hc1, K1basis, 'H1')
pw = lambda pc: [Fr(c) for c in pc]
ae1 = K1.from_power(pw([-6, 0, Fr(1, 6), Fr(1, 6)]))           # eta^3/6 + eta^2/6 - 6
ae2 = K1.from_power(pw([-33, 22, Fr(23, 3), Fr(8, 3)]))        # (8/3)eta^3 + (23/3)eta^2 + 22 eta - 33
check(all(c.denominator == 1 for c in ae1 + ae2), f'eps1, eps2 lie in O_H1: coordinates {[int(c) for c in ae1]},'
                                                   f' {[int(c) for c in ae2]} in [1, eta, eta^2/6, eta^3/6]')
ae1, ae2 = tuple(int(c) for c in ae1), tuple(int(c) for c in ae2)
claim('N(eps1) = N(eps2) = +-1 (p. 411)', abs(K1.norm(ae1)) == 1 and abs(K1.norm(ae2)) == 1,
      f'N(eps1) = {K1.norm(ae1)}, N(eps2) = {K1.norm(ae2)}')
inv1 = K1.pow(ae1, -1)
inv2 = K1.pow(ae2, -1)
print(f'   exact inverses: eps1^(-1) = {elt(K1.to_power(inv1), eta_s)};  eps2^(-1) = {elt(K1.to_power(inv2), eta_s)}')


def ppow(pc, e, fcp):
    r = [Fr(1), Fr(0), Fr(0), Fr(0)]
    for _ in range(e):
        r = pmulmod(r, pc, fcp)
    return r


mid1 = [-c for c in ppow(pw([0, 0, Fr(-1, 6), 0]), 2, hc1)]           # -(-(1/6) eta^2)^2
lhs1 = K1.to_power(K1.mul((-1, 0, 0, 0), K1.pow(ae1, -2)))              # -eps1^(-2)
claim('(5): -eps1^(-2) = -(-(1/6)eta^2)^2 = 1 - eta', lhs1 == mid1 == pw([1, -1, 0, 0]),
      f'-eps1^(-2) = {elt(lhs1, eta_s)}, -(-(1/6)eta^2)^2 = {elt(mid1, eta_s)}')
mid2 = [-c for c in ppow(pw([3, -2, Fr(1, 3), 0]), 6, hc1)]           # -((1/3)eta^2 - 2 eta + 3)^6
lhs2 = K1.to_power(K1.mul((-1, 0, 0, 0), K1.pow(ae2, -6)))              # -eps2^(-6)
claim('(5): -eps2^(-6) = -((1/3)eta^2 - 2eta + 3)^6 = -24156 eta^3 + 22494 eta^2 + 207144 eta - 211085',
      lhs2 == mid2 == pw([-211085, 207144, 22494, -24156]),
      f'-eps2^(-6) = {elt(lhs2, eta_s)}, -((1/3)eta^2 - 2eta + 3)^6 = {elt(mid2, eta_s)}')
key1, P1, Q1 = ISO['H1']
K = Field(key1)
fc = fcoeffs(key1)


def transport(pc_eta):
    """sum c_k eta^k -> sum c_k P(alpha)^k in power coordinates of class key1."""
    acc = [Fr(0)] * 4
    for c in reversed(pc_eta):
        acc = pmulmod(acc, P1, fc)
        acc[0] += Fr(c)
    return acc


te = []
for nm, u in (('eps1', ae1), ('eps2', ae2)):
    co = K.from_power(transport(K1.to_power(u)))
    check(all(c.denominator == 1 for c in co), f'his {nm} transported into class {LABEL[key1]} (eta -> P(alpha)):'
                                               f' coordinates {[int(c) for c in co]}')
    te.append(tuple(int(c) for c in co))
    check(K.norm(te[-1]) == K1.norm(u), f'norm preserved: {K.norm(te[-1])}')
ar = mp_roots(fc)
rr = sorted([mpmath.re(r) for r in ar if abs(mpmath.im(r)) < 1e-40])


def loglist(u):
    pc = K.to_power(u)
    return [mpmath.log(abs(sum(mpmath.mpf(c.numerator) / c.denominator * r ** k for k, c in enumerate(pc))))
            for r in rr]


def express(target, basis):
    """Integers (a,b), sign s with target = s * basis0^a basis1^b, verified exactly (or None)."""
    L = mpmath.matrix([[loglist(basis[0])[0], loglist(basis[1])[0]], [loglist(basis[0])[1], loglist(basis[1])[1]]])
    ab = mpmath.lu_solve(L, mpmath.matrix(loglist(target)))
    a, b = int(mpmath.nint(ab[0])), int(mpmath.nint(ab[1]))
    g = K.mul(K.pow(basis[0], a), K.pow(basis[1], b))
    if g == tuple(target):
        return a, b, 1
    if g == tuple(-x for x in target):
        return a, b, -1
    return None


ours = UNITS[key1]
Amat, Bmat = [], []
for nm, u in zip(('eps1', 'eps2'), te):
    r = express(u, ours)
    check(r is not None, f'his {nm} = {"" if r and r[2] == 1 else "-"}eps1^{r[0] if r else "?"} eps2^'
                         f'{r[1] if r else "?"} in our units of class {LABEL[key1]} (exact)')
    Amat.append(r[:2] if r else (0, 0))
for nm, u in zip(('eps1', 'eps2'), ours):
    r = express(u, te)
    check(r is not None, f'our {nm} = {"" if r and r[2] == 1 else "-"}e1^{r[0] if r else "?"} e2^'
                         f'{r[1] if r else "?"} in his transported units e1, e2 (exact)')
    Bmat.append(r[:2] if r else (0, 0))
prod = [[sum(Amat[i][k] * Bmat[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
detA = Amat[0][0] * Amat[1][1] - Amat[0][1] * Amat[1][0]
check(prod == [[1, 0], [0, 1]] and abs(detA) == 1,
      f'exponent matrices {Amat} (his in ours) and {Bmat} (ours in his) are inverse, det = {detA}: his units and'
      f' core.UNITS generate the same group modulo +-1 in class {LABEL[key1]}')

print()
print('(e) further statements of the paper, checked')
# (7), (8): Skolem data in Z[eta] (power coordinates)
b1 = pw([1, -1, 0, 0])                                        # bar eps1 = 1 - eta
b2 = pw([-211085, 207144, 22494, -24156])                     # bar eps2
p9 = [Fr(1), Fr(0), Fr(0), Fr(0)]
for _ in range(9):
    p9 = pmulmod(p9, b1, hc1)
A9 = [(c - (1 if k == 0 else 0)) / 3 for k, c in enumerate(p9)]
claim('(7): bar eps1^9 = 1 + 3A with A in Z[eta], A == 2 eta^3 (mod 3)',
      all(c.denominator == 1 for c in A9) and [int(c) % 3 for c in A9] == [0, 0, 0, 2],
      f'A mod 3 = {[int(c) % 3 for c in A9]}')
Bq = [(c - (1 if k == 0 else 0)) / 3 for k, c in enumerate(b2)]
claim('(7): bar eps2 = 1 + 3B with B in Z[eta], B == eta^2 (mod 3)',
      all(c.denominator == 1 for c in Bq) and [int(c) % 3 for c in Bq] == [0, 0, 1, 0],
      f'B mod 3 = {[int(c) % 3 for c in Bq]}')
table8 = {0: [1, 0, 0, 0], 1: [1, -1, 0, 0], 2: [1, 1, 1, 0], 3: [1, 0, 0, -1], 4: [1, -1, 0, -1],
          5: [1, 1, 1, -1], 6: [1, 0, 0, 1], 7: [1, -1, 0, 1], 8: [1, 1, 1, 1]}
ok8 = True
cur = [Fr(1), Fr(0), Fr(0), Fr(0)]
two_term = []
for r in range(9):
    ok8 = ok8 and [int(c) % 3 for c in cur] == [c % 3 for c in table8[r]]
    if [int(c) % 3 for c in cur][2:] == [0, 0]:
        two_term.append(r)
    cur = pmulmod(cur, b1, hc1)
claim('(8): the table of V(r, beta) = bar eps1^r bar eps2^beta mod 3, r = 0..8', ok8)
claim('p. 412: only r = 0, 1 give eta^2-, eta^3-coefficients == 0 mod 3', two_term == [0, 1], f'r = {two_term}')
# lemma solution lists by box search (confirmation only)
Rb = 300
for nm, (co, sols) in HIS_FORMS.items():
    hits = set()
    for a in range(-Rb, Rb + 1):
        bv = np.arange(-Rb, Rb + 1, dtype=np.int64)
        aa = np.int64(a)
        v = co[0] * aa ** 4 + co[1] * aa ** 3 * bv + co[2] * aa * aa * bv * bv + co[3] * aa * bv ** 3 + co[4] * bv ** 4
        for i in np.nonzero(v == 1)[0]:
            b = int(bv[i])
            if Fval_co(co, a, b) == 1:
                hits.add(frozenset({(a, b), (-a, -b)}))
    want = {frozenset({s, (-s[0], -s[1])}) for s in sols}
    claim(f'Lemma {nm[1]}: the solutions of {nm} = 1 are {sols} (box search |N|,|M| <= {Rb}, confirmation only)',
          hits == want, f'found {sorted(min(h) for h in hits)}')
# the cubic field
lam_s = sp.Symbol('lambda')
fl = lam_s ** 3 - 4 * lam_s + 6
red = lambda e: sp.expand(sp.rem(sp.expand(e), fl, lam_s))
th_s = sp.Symbol('theta')
check(sp.expand((6 * lam_s) ** 3 - 144 * (6 * lam_s) + 1296 - 216 * fl) == 0,
      'theta = 6 lambda: theta^3 - 144 theta + 1296 = 216 (lambda^3 - 4 lambda + 6) (our cubic = his (17))')
claim('p. 413: D(lambda) = -2^2 179', int(sp.discriminant(fl, lam_s)) == -716, f'disc = {sp.discriminant(fl, lam_s)}')
a_, b_, c_ = sp.symbols('a b c')
claim('(18) norm formula', sp.expand(sp.resultant(fl, a_ + b_ * lam_s + c_ * lam_s ** 2, lam_s)
                                     - (a_ ** 3 - 6 * b_ ** 3 + 36 * c_ ** 3 + 8 * a_ ** 2 * c_ - 4 * a_ * b_ ** 2
                                        + 16 * a_ * c_ ** 2 + 24 * b_ * c_ ** 2 + 18 * a_ * b_ * c_)) == 0)
claim('(19) square formula', red((a_ + b_ * lam_s + c_ * lam_s ** 2) ** 2 - (a_ ** 2 - 12 * b_ * c_ + lam_s * (
    2 * a_ * b_ + 8 * b_ * c_ - 6 * c_ ** 2) + lam_s ** 2 * (b_ ** 2 + 2 * a_ * c_ + 4 * c_ ** 2))) == 0)
eps0 = 40 * lam_s ** 2 - 101 * lam_s + 95
epsb = lam_s ** 2 - 5 * lam_s - 19
claim('p. 413-414: eps0 = 40 lambda^2 - 101 lambda + 95 is a unit and 1/eps0 = lambda^2 - 5 lambda - 19',
      red(eps0 * epsb) == 1 and abs(sp.resultant(fl, eps0, lam_s)) == 1, f'N(eps0) = {sp.resultant(fl, eps0, lam_s)}')
claim('(21): 2 = -(lambda^2 + lambda - 4)^3 eps^(-1), eps = 1/eps0', red(-(lam_s ** 2 + lam_s - 4) ** 3 * eps0) == 2)
claim('(22): 3 = (lambda - 1)(2 lambda + 5)(3 lambda^2 + 4 lambda - 9) eps^(-1) = (lambda-1)(2 lambda+5)'
      '(lambda^2 - 3 lambda + 3)',
      red((lam_s - 1) * (2 * lam_s + 5) * (3 * lam_s ** 2 + 4 * lam_s - 9) * eps0) == 3
      and red((lam_s - 1) * (2 * lam_s + 5) * (lam_s ** 2 - 3 * lam_s + 3)) == 3)
claim('p. 414: 3 lambda^2 - 4 = (lambda^2 + lambda - 4)^2 (45 lambda^2 - 114 lambda + 107), N(45 lambda^2 - ...) = 179',
      red((lam_s ** 2 + lam_s - 4) ** 2 * (45 * lam_s ** 2 - 114 * lam_s + 107)) == red(3 * lam_s ** 2 - 4)
      and sp.resultant(fl, 45 * lam_s ** 2 - 114 * lam_s + 107, lam_s) == 179)
avx = sorted(XY[0] // 6 for XY in POINTS.values() if XY[0] % 6 == 0 and XY[1] % 36 == 0)
claim('p. 419: the integral solutions of (16) are X\' = -2, 0, 2, 4, 8, 18, 42, 70; the known solutions core.KNOWN give'
      ' exactly these points of his model', avx == [-2, 0, 2, 4, 8, 18, 42, 70], f'from core.KNOWN: X\' in {avx}')

print()
print("(f) which equations F = c can contain points of Avanesov's model (6 | X), by residues mod 6c")
print('    an integral point needs F(m,n) = c = g^2 with g | u, v, w and X = X4(m,n)/g^2; 6 | X needs X4 == 0 mod 6c')
can = {}
for key in KEYS:
    uvw, X4, Y6 = descent_maps(key)
    for c in RHS:
        g = isqrt(c)
        Q = 6 * c
        poss = False
        for a in range(Q):
            for b in range(Q):
                if (a % 2 == 0 and b % 2 == 0) or (a % 3 == 0 and b % 3 == 0):
                    continue
                if (Fval(key, a, b) - c) % Q:
                    continue
                u, v, w = uvw(a, b)
                if u % g or v % g or w % g:
                    continue
                if X4(a, b) % Q == 0:
                    poss = True
                    break
            if poss:
                break
        can[(key, c)] = poss
    print(f'   class {LABEL[key]}: ' + ', '.join(f'c = {c}: {"possible" if can[(key, c)] else "impossible"}'
                                              for c in RHS))
possible = {kc for kc, v in can.items() if v}
check(cov <= possible, 'consistency: the six (class, c) of his equations L2..L7 pass the residue test: '
      + ', '.join(f'({LABEL[k]}, {c})' for k, c in sorted(cov)))
if possible == cov:
    print('   => exactly the six (class, c) of his equations L2..L7 can contain points of his model; the other six'
          ' equations F = c contain none (residues mod 6c).')
else:
    print('   => residue test inconclusive for: ' + ', '.join(f'({LABEL[k]}, {c})' for k, c in sorted(possible - cov)))

# ------------------------------------------------------------------------------------------------------------------
header('SUMMARY')
for key in KEYS:
    print(f'class {LABEL[key]}: f = {fexpr(fcoeffs(key))}; disc(f) = {fz(DISC[key])}; disc(O_K) = {DISC_O[key]} ='
          f' {fz(DISC_O[key])}')
for hk, (key, P, Q) in ISO.items():
    print(f'{hk} (h = {fexpr(HIS_F[hk])}) = class {LABEL[key]} via eta = {elt(P, alpha_s)}; disc = {HIS_DISC[hk]}')
for nm, (key, T, lam, c) in RELS.items():
    print(f'{nm}: F_{LABEL[key]}(T(N,M)) = {lam} {nm}(N,M), {mn_str(T)};  {nm} = 1 <=> F_{LABEL[key]}(m,n) = {c}'
          f' (m,n coprime)')
if possible == cov:
    print('for points of Avanesov\'s model only these six (class, c) can occur (6(f)); a solution of one of the other'
          ' six equations F = c gives no integral point with 6 | X.')
print('his units eps1, eps2 of K1(eta) = H1 and core.UNITS[\'10\'] generate the same group modulo +-1 (exponent matrix'
      f' {Amat}); whether either pair is fundamental is not checked here.')
print(f'{sum(1 for _, v in CLAIMS if v)} of {len(CLAIMS)} tested statements of Avanesov\'s paper are true.')
for lab, v in CLAIMS:
    if not v:
        print(f'   FALSE: {lab}')
print(f'[done in {time.time() - T0:.1f} s]')
if FAILED:
    print(f'{len(FAILED)} CHECK FAILED')
    sys.exit(1)
print('all checks passed')
