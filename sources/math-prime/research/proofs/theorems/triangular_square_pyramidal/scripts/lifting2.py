"""Helpers for the exclusion step: kernel lattices, coset systems, and the hypotheses of the two-index lemma.

The lemma is research/proofs/methods/lifting.binomial.multi/proof.md.  Notation there: T_l = I + q K_l
(q = p^e), seed v, markers l_1, l_2, f_k(m,n) = l_k T_1^m T_2^n v.  Here:
  * the markers are the coordinates x2 and x3 in the integral basis (rows (0,0,1,0) and (0,0,0,1));
  * T_j = matrix of multiplication by a unit eta_j = +- eps1^a eps2^b with eta_j == 1 (mod q);
  * a "coset system" is the orbit {T_1^m T_2^n v : (m,n) in Z^2} of a seed v.

Orbit decomposition.  Let Gamma = <-1, eps1, eps2>, w a seed of norm c, and q = p^e (q >= 3).
Lambda_q := {(a,b) in Z^2 : eps1^a eps2^b == 1 mod q}, with Hermite basis (h1, 0), (-a1, b0)
(h1 = order of eps1 mod q; b0 = least b > 0 with eps2^b in <eps1> mod q, eps2^b0 == eps1^a1).
  * mode 'eps':    the family <eps1,eps2> w splits into the |Z^2/Lambda_q| coset systems
                   {eps^(a,b) w : (a,b) in (a0,b0) + Lambda_q}, with eta_1 = eps1^h1, eta_2 = eps1^-a1 eps2^b0.
                   Zeros in this family correspond one-to-one to solution pairs +-(m,n).
  * mode 'signed': if eps^(a*,b*) == -1 (mod q) for some (a*,b*), put Lambda'' = Lambda_q + Z (a*,b*) and use
                   signed units (-1)^s eps^lambda == 1 (mod q).  The coset systems then cover Gamma w = +-<eps> w,
                   and each solution pair gives two zeros, in different coset systems (q does not divide 2).
Every element of a coset system is == its seed (mod q), so a coset system whose seed has a marker not == 0
(mod q) has no zero.  The others are the candidates.

Precision: seeds are kept mod p^N and the matrices K_j = (T_j - I)/q mod p^Nk, with N = e + PREC_EXTRA and
Nk = N - e = PREC_EXTRA.  nu() returns the cap Nk for a residue == 0 mod p^Nk, so a printed value equal to the cap
means "at least Nk".  Such values never certify anything: Proposition 17 is accepted only when nu(J) < Nk and
lambda <= e - 1 - iota < Nk.
"""
from math import gcd
from core import Field, UNITS

PREC_EXTRA = 40


def nu(x, p, cap):
    """p-adic valuation of the integer x, capped at cap (x == 0 gives cap)."""
    if x == 0:
        return cap
    v = 0
    while x % p == 0:
        x //= p
        v += 1
        if v >= cap:
            return cap
    return v


def mat_mul_mod(A, B, M):
    return [[sum(A[i][k] * B[k][j] for k in range(4)) % M for j in range(4)] for i in range(4)]


def mat_vec_mod(A, v, M):
    return [sum(A[i][k] * v[k] for k in range(4)) % M for i in range(4)]


def det4_mod(cols, M):
    """Determinant (mod M) of the 4x4 matrix with the given columns."""
    A = [[cols[j][i] % M for j in range(4)] for i in range(4)]
    # exact integer determinant by cofactor expansion, then reduce
    def det(B):
        n = len(B)
        if n == 1:
            return B[0][0]
        s = 0
        for j in range(n):
            if B[0][j]:
                minor = [row[:j] + row[j + 1:] for row in B[1:]]
                s += (-1) ** j * B[0][j] * det(minor)
        return s
    return det(A) % M


class Kernel:
    """Lambda_q for <eps1, eps2> acting on O/q, plus the optional sign twist."""

    def __init__(self, K, q, units=None):
        self.K, self.q = K, q
        e1, e2 = units if units is not None else UNITS[K.key]
        self.e1, self.e2 = e1, e2
        one = (1, 0, 0, 0)
        # powers of eps1 mod q
        pw = [one]
        x = tuple(v % q for v in e1)
        while x != one:
            pw.append(x)
            x = K.mulmod(x, e1, q)
            assert len(pw) < 10 ** 7
        self.h1 = len(pw)
        index = {g: a for a, g in enumerate(pw)}
        # least b0 > 0 with eps2^b0 in <eps1> mod q
        y = tuple(v % q for v in e2)
        b = 1
        while y not in index:
            y = K.mulmod(y, e2, q)
            b += 1
            assert b < 10 ** 7
        self.b0, self.a1 = b, index[y]          # eps2^b0 == eps1^a1 (mod q)
        self.order = self.h1 * self.b0          # |G_q|
        self._pw1 = pw
        # sign twist: is -1 in G_q ?
        minus = tuple((-1) % q if i == 0 else 0 for i in range(4))
        self.star = None
        z = one
        for bb in range(self.b0):
            for aa in range(self.h1):
                if K.mulmod(pw[aa], z, q) == minus:
                    self.star = (aa, bb)
                    break
            if self.star:
                break
            z = K.mulmod(z, e2, q)

    def reduce(self, a, b):
        """Canonical representative of (a,b) mod Lambda_q: 0 <= a < h1, 0 <= b < b0."""
        k = b // self.b0
        b2 = b - k * self.b0
        a2 = (a + k * self.a1) % self.h1
        return (a2, b2)

    def reps(self):
        return [(a, b) for b in range(self.b0) for a in range(self.h1)]

    def residue(self, a, b, w):
        K, q = self.K, self.q
        g = K.mulmod(K.powmod(self.e1, a % self.h1, q), K.powmod(self.e2, b, q), q) if b >= 0 else None
        if g is None:
            a, b = self.reduce(a, b)
            g = K.mulmod(K.powmod(self.e1, a, q), K.powmod(self.e2, b, q), q)
        return K.mulmod(g, w, q)

    def candidates(self, w):
        """Representatives (a,b) of Z^2/Lambda_q whose residue eps^(a,b) w has both markers == 0 (mod q)."""
        K, q = self.K, self.q
        out = []
        z2 = tuple(v % q for v in w)
        for b in range(self.b0):
            for a in range(self.h1):
                r = K.mulmod(self._pw1[a], z2, q)
                if r[2] == 0 and r[3] == 0:
                    out.append((a, b))
            z2 = K.mulmod(z2, self.e2, q)
        return out

    def eta_exponents(self, mode='eps'):
        """Exponent vectors (s, a, b) of eta_1, eta_2 (unit (-1)^s eps1^a eps2^b == 1 mod q) for the mode."""
        if mode == 'eps':
            return [(0, self.h1, 0), (0, -self.a1, self.b0)]
        assert mode == 'signed' and self.star is not None
        # Lambda'' = Lambda_q + Z(a*,b*): Hermite basis of the lattice spanned by (h1,0), (-a1,b0), (a*,b*)
        gens = [(0, self.h1, 0), (0, -self.a1, self.b0), (1,) + self.star]
        return _hermite2_signed(gens, self)


def _hermite2_signed(gens, ker):
    """Basis of the rank-2 lattice spanned by the (a,b) parts of gens, carrying the signs (mod 2) along."""
    # integer row reduction on (a,b) with sign bookkeeping; rows (s,a,b)
    rows = [list(g) for g in gens]
    # make the b-column have a single nonzero entry (gcd), Euclid on column 2
    def euclid(rows, col):
        while True:
            nz = [r for r in rows if r[col] != 0]
            if len(nz) <= 1:
                return rows
            nz.sort(key=lambda r: abs(r[col]))
            piv = nz[0]
            for r in nz[1:]:
                k = r[col] // piv[col]
                r[0] = (r[0] - k * piv[0]) % 2
                r[1] -= k * piv[1]
                r[2] -= k * piv[2]
    rows = euclid(rows, 2)
    rb = [r for r in rows if r[2] != 0]
    ra = [r for r in rows if r[2] == 0]
    ra = euclid(ra, 1)
    ra = [r for r in ra if r[1] != 0]
    # rows with a = b = 0 must have sign 0, otherwise -1 == 1 (mod q)
    assert all(r[0] % 2 == 0 for r in rows if r[1] == 0 and r[2] == 0)
    assert len(ra) == 1 and len(rb) == 1
    return [tuple(ra[0]), tuple(rb[0])]


def unit_from_exponents(K, e1, e2, s, a, b):
    u = K.mul(K.pow(e1, a), K.pow(e2, b))
    return tuple(-c for c in u) if s % 2 else u


class CosetSystem:
    """The system (T_1, T_2, v) of one coset, reduced mod p^N.  K_j = (T_j - I)/q mod p^(N-e)."""

    def __init__(self, K, p, e, seed, etas, N=None):
        self.K, self.p, self.e = K, p, e
        self.q = p ** e
        self.N = N if N is not None else e + PREC_EXTRA
        PN = p ** self.N
        self.PN = PN
        self.Kmats = []
        for eta in etas:
            T = K.matrix(eta)
            for i in range(4):
                for j in range(4):
                    d = T[i][j] - (1 if i == j else 0)
                    assert d % self.q == 0, 'eta is not == 1 mod q'
            Km = [[((T[i][j] - (1 if i == j else 0)) // self.q) % (p ** (self.N - e)) for j in range(4)] for i in range(4)]
            self.Kmats.append(Km)
        self.v = [c % PN for c in seed]
        self.Nk = self.N - e          # precision of K-matrices and of everything derived from them

    # --- Proposition 17 (Jacobian criterion) at the seed, which must be a common zero ------------------
    def prop17(self):
        p, e, Nk = self.p, self.e, self.Nk
        M = p ** Nk
        assert self.v[2] % M == 0 and self.v[3] % M == 0, 'seed is not a common zero (mod p^N)'
        K1v = mat_vec_mod(self.Kmats[0], self.v, M)
        K2v = mat_vec_mod(self.Kmats[1], self.v, M)
        x = (K1v[2], K1v[3])            # l_k K_1 w
        y = (K2v[2], K2v[3])            # l_k K_2 w
        lam = min(nu(t, p, Nk) for t in (x[0], x[1], y[0], y[1]))
        J = (x[0] * y[1] - x[1] * y[0]) % M
        nJ = nu(J, p, Nk)
        iota = 1 if p == 2 else 0
        delta = nJ - lam
        ok = lam <= e - 1 - iota and delta <= e - 1 - iota and nJ < Nk
        return dict(lam=lam, nuJ=nJ, delta=delta, ok=ok, bound=1 if ok else None)

    # --- Theorem 18 (frame criterion) ---------------------------------------------------------------------
    def frame_form(self):
        """(Q20, Q11, Q02) mod p^Nk with Q_v(c) = c1^2 Q20 + 2 c1 c2 Q11 + c2^2 Q02 (Proposition 12(c))."""
        p, Nk = self.p, self.Nk
        M = p ** Nk
        K1, K2 = self.Kmats
        v = [t % M for t in self.v]
        K1v, K2v = mat_vec_mod(K1, v, M), mat_vec_mod(K2, v, M)
        K11v = mat_vec_mod(K1, K1v, M)
        K12v = mat_vec_mod(K1, K2v, M)
        K22v = mat_vec_mod(K2, K2v, M)
        Q20 = det4_mod([v, K1v, K2v, K11v], M)
        Q11 = det4_mod([v, K1v, K2v, K12v], M)
        Q02 = det4_mod([v, K1v, K2v, K22v], M)
        self.Q = (Q20, Q11, Q02)
        return self.Q

    def frame_depth(self, smax=12, cap=10 ** 5):
        """Least s_F <= smax with nu(Q_v(c)) <= s_F for every primitive c in Z^2, or None.

        A primitive c is, up to a p-adic unit factor (which does not change nu(Q(c))), either (1, t) or (u, 1)
        with p | u.  So nu(Q(c)) >= j for some primitive c iff g1(t) = Q(1,t) or g2(u) = Q(u,1) (p | u) has a
        root mod p^j.  The root sets mod p^j are lifted level by level (t -> t + p^(j-1) s), exactly; s_F is one
        less than the first level at which both are empty."""
        p, Nk = self.p, self.Nk
        A, B, C = self.frame_form()
        g1 = lambda t: A + 2 * B * t + C * t * t
        g2 = lambda u: A * u * u + 2 * B * u + C
        R1 = [t for t in range(p) if g1(t) % p == 0]
        R2 = [0] if g2(0) % p == 0 else []
        j = 1
        while R1 or R2:
            if j > smax or j + 1 > Nk:
                return None
            m = p ** (j + 1)
            R1 = [t + p ** j * s for t in R1 for s in range(p) if g1(t + p ** j * s) % m == 0]
            R2 = [u + p ** j * s for u in R2 for s in range(p) if g2(u + p ** j * s) % m == 0]
            if len(R1) + len(R2) > cap:
                return None
            j += 1
        return j - 1

    def thm18(self, smax=6):
        sF = self.frame_depth(smax)
        iota2 = 1 if self.p in (2, 3) else 0
        ok = sF is not None and self.e >= sF + 1 + iota2
        return dict(sF=sF, ok=ok, bound=2 if ok else None)

    def markers_mod(self, M=None):
        M = M or self.q
        return (self.v[2] % M, self.v[3] % M)


def coset_seed(K, w, e1, e2, a, b, s=0):
    """Exact element (-1)^s eps1^a eps2^b w."""
    g = unit_from_exponents(K, e1, e2, s, a, b)
    return K.mul(g, w)


def locate(K, z, w, e1, e2, Em=None, search=60):
    """Find (s, a, b) with z = (-1)^s eps1^a eps2^b w exactly, or None.  Uses real logarithms only to propose
    (a, b); the identity itself is checked exactly."""
    import math
    from core import Embeddings
    q = K.divides_exactly(w, z)
    if q is None or abs(K.norm(q)) != 1:
        return None
    Em = Em or Embeddings(K)
    f = lambda x: (float(Em.sig1(x)), float(Em.sig2(x)))
    l1 = [math.log(abs(t)) for t in f(e1)]
    l2 = [math.log(abs(t)) for t in f(e2)]
    lq = [math.log(abs(t)) for t in f(q)]
    det = l1[0] * l2[1] - l1[1] * l2[0]
    a = (lq[0] * l2[1] - lq[1] * l2[0]) / det
    b = (l1[0] * lq[1] - l1[1] * lq[0]) / det
    for da in (0, -1, 1):
        for db in (0, -1, 1):
            A, B = round(a) + da, round(b) + db
            u = K.mul(K.pow(e1, A), K.pow(e2, B))
            if u == q:
                return (0, A, B)
            if u == tuple(-t for t in q):
                return (1, A, B)
    return None
