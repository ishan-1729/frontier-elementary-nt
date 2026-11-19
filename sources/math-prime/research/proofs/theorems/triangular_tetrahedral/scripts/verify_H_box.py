"""verify_H_box.py: an independent check of hypothesis_H.md, Section 6 (unique factorization of R = Z[l] by the box
principle; 08_ufd_box.py is the author's script for it).

Independence. Nothing is imported from cubic.py, core.py, lifting2.py or any 0x_*.py script, and no code is copied from
them. The arithmetic of R is rebuilt here from the definition l^3 = 4l - 6 (product of polynomials, then reduction of
l^4 and l^3); the norm is the determinant of the multiplication matrix built from that product; division is Cramer's
rule; lattice equality is decided by a Hermite normal form written here; f mod p is factored by exhaustive search
(small p) or by Rabin's gcd method (large p). The text's claims are typed in from hypothesis_H.md as data and compared.
The only file of the author that is read is 08_ufd_box.out.txt, as data, for the cross-check (e7).

Exact arithmetic only: int, Fraction, sympy. numpy is used only with int64 arrays (exact integer arithmetic, values
far below 2^63). No float enters any decision; decimals are printed for reading only.

Parts:
 (a) 6.3, the constant. F(x) = N(x0 + x1(1 + l) + x2(l^2 + l - 3)); the text claims |F| <= C0 = 535/24 on the box
     |x0| <= 2, |x1| <= 1, |x2| <= 1/2 (by Bernstein coefficients). Here, by two other methods:
     (a1) exact critical-point analysis: interior (Euler's identity), the 6 faces (resultant, Sturm isolation, exact
          interval enclosures), the 12 edges (exact quadratic surds), the 8 vertices (exact);
     (a2) branch and bound with the exact Taylor expansion of the cubic at box centres, in Fractions;
     (a3) control: the exact maximum over an integer grid lies below the maximum found.
 (b) 6.3, the lemma, for every prime p < 10^5: B_i = floor(t_i p^(1/3)), B_i^3 <= t_i^3 p < (B_i + 1)^3,
     prod(B_i + 1) > p; and the pigeonhole run for every root r of f mod p, giving beta != 0 in P(p, t - r) with
     |x_i| <= B_i and 1 <= |N(beta)|/p <= 22.
 (c) 6.2 for every prime p <= 60, by all p^3 residues: p | N(beta) iff beta lies in some P(p,g), with f mod p
     factored here by exhaustive search over monic polynomials of degree 1 and 2. Also p = 179 (numpy).
 (d) 6.4, the table: membership, norm, and P(p,g) = pi R, by equal Hermite normal forms and by the index argument.
 (e) 6.5, the descent, for every prime 23 <= p < 2000 with a root mod p (run on to p < 3000 for the text's counts).
 (f) the S-primes PI2, P3A, P3B, P3C, P179 and Q179 = (-23, 11, 8); 179 = (unit) P179^2 Q179.
Lines marked CONTROL plant an error and PASS when the check catches it.
"""
import sys
sys.dont_write_bytecode = True

import heapq
import itertools
import math
import os
import re
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

T_START = time.time()
N_CHECKS = 0


def check(name, cond):
    """Print PASS or FAIL; stop at the first FAIL."""
    global N_CHECKS
    N_CHECKS += 1
    print(("PASS  " if cond else "FAIL  ") + name, flush=True)
    if not cond:
        print(f"SUMMARY: FAIL at check {N_CHECKS} ({time.time() - T_START:.1f} s): {name}")
        sys.exit(1)


def control(name, caught):
    """A planted error: PASS when it is caught."""
    check("CONTROL " + name + (" -> caught" if caught else " -> NOT caught"), bool(caught))


def dec(q, n=12):
    """Decimal string of a Fraction rounded to n places, for reading only."""
    q = Fr(q)
    s = "-" if q < 0 else ""
    r = round(abs(q)*10**n)
    return f"{s}{r // 10**n}.{r % 10**n:0{n}d}"


# ====================================================================================================================
# The ring R = Z[l], l^3 = 4l - 6, rebuilt from the definition.
# ====================================================================================================================
def rmul(u, v):
    """Product of two triples: multiply as polynomials in l, then replace l^k (k = 4, 3) by 4 l^(k-2) - 6 l^(k-3)."""
    w = [0]*5
    for i in range(3):
        for j in range(3):
            w[i + j] += u[i]*v[j]
    for k in (4, 3):
        ck = w[k]
        w[k] = 0
        w[k - 2] += 4*ck
        w[k - 3] -= 6*ck
    return (w[0], w[1], w[2])


E3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
LAM = (0, 1, 0)


def mult_matrix(u):
    """3x3 matrix (rows = coordinates) whose columns are u*1, u*l, u*l^2."""
    cols = [rmul(u, e) for e in E3]
    return [[cols[j][i] for j in range(3)] for i in range(3)]


def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
            + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]))


A_, B_, C_ = sp.symbols('a b c')
NORM_EXPR = sp.expand(sp.Matrix(mult_matrix((A_, B_, C_))).det())
_norm_fast = eval("lambda a, b, c: " + str(NORM_EXPR))      # exact on ints: built from the expression just derived


def N(u):
    return _norm_fast(u[0], u[1], u[2])


def rdiv(beta, pi):
    """beta / pi in R, or None: solve M(pi) y = beta by Cramer's rule and require y integral."""
    M = mult_matrix(pi)
    d = det3(M)
    y = []
    for i in range(3):
        Mi = [row[:] for row in M]
        for r in range(3):
            Mi[r][i] = beta[r]
        di = det3(Mi)
        if di % d:
            return None
        y.append(di//d)
    y = tuple(y)
    if rmul(pi, y) != tuple(beta):
        check(f"internal: Cramer quotient {beta}/{pi} multiplies back", False)
    return y


def is_unit(u):
    return N(u) in (1, -1)


def associates(u, v):
    q = rdiv(u, v)
    return q is not None and is_unit(q)


def rpow(u, k):
    r = (1, 0, 0)
    for _ in range(k):
        r = rmul(r, u)
    return r


def tstr(u):
    """a + b l + c l^2 as text."""
    out = []
    for co, mon in ((u[2], "l^2"), (u[1], "l"), (u[0], "")):
        if co == 0:
            continue
        mag = abs(co)
        body = (str(mag) if (mag != 1 or mon == "") else "") + mon
        out.append(("-" if co < 0 else "+") + body)
    if not out:
        return "0"
    s = " ".join(out).replace("+", "+ ").replace("-", "- ")
    return s[2:] if s.startswith("+ ") else "-" + s[2:]


# ---------------------------------------------------------------- polynomials over F_p (lists, constant term first)
def ptrim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def pdivmod(a, g, p):
    """Quotient and remainder of a by g over F_p (leading coefficient of g invertible mod p)."""
    a = ptrim([x % p for x in a])
    g = ptrim([x % p for x in g])
    dg = len(g) - 1
    inv = pow(g[-1], -1, p)
    q = [0]*max(len(a) - dg, 1)
    while a and len(a) - 1 >= dg:
        co = a[-1]*inv % p
        sh = len(a) - 1 - dg
        q[sh] = co
        for i in range(dg + 1):
            a[sh + i] = (a[sh + i] - co*g[i]) % p
        a = ptrim(a)
    return ptrim(q), a


def prem(a, g, p):
    return pdivmod(a, g, p)[1]


def pmul(u, v, p):
    if not u or not v:
        return []
    w = [0]*(len(u) + len(v) - 1)
    for i, x in enumerate(u):
        for j, y in enumerate(v):
            w[i + j] = (w[i + j] + x*y) % p
    return ptrim(w)


def psub(u, v, p):
    n = max(len(u), len(v))
    return ptrim([((u[i] if i < len(u) else 0) - (v[i] if i < len(v) else 0)) % p for i in range(n)])


def pgcd(u, v, p):
    """Monic gcd over F_p."""
    u, v = ptrim([x % p for x in u]), ptrim([x % p for x in v])
    while v:
        u, v = v, prem(u, v, p)
    if not u:
        return []
    inv = pow(u[-1], -1, p)
    return [x*inv % p for x in u]


def ppowmod(base, e, m, p):
    r, b = [1], prem(base, m, p)
    while e:
        if e & 1:
            r = prem(pmul(r, b, p), m, p)
        b = prem(pmul(b, b, p), m, p)
        e >>= 1
    return r


def peval(g, x, p):
    s = 0
    for co in reversed(g):
        s = (s*x + co) % p
    return s


def f_mod(p):
    return ptrim([6 % p, (-4) % p, 0, 1])


def gstr(g, p=None):
    """A monic polynomial over F_p as text; a linear factor t - r is shown with its root r (0 <= r < p)."""
    d = len(g) - 1
    if d == 1 and p is not None:
        r = (-g[0]) % p
        return "t" if r == 0 else f"t - {r}"
    parts = []
    for i in range(d, -1, -1):
        co = g[i]
        if co == 0:
            continue
        mon = "" if i == 0 else ("t" if i == 1 else f"t^{i}")
        parts.append((str(co) if (co != 1 or i == 0) else "") + mon)
    return " + ".join(parts)


def factor_exhaustive(p):
    """Monic irreducible factors of f mod p with multiplicities, by trying every monic polynomial of degree 1 and 2
    (dividing out each divisor found); a remaining cubic has no linear factor and is therefore irreducible."""
    rest = f_mod(p)
    facs = []
    for d in (1, 2):
        for low in itertools.product(range(p), repeat=d):
            g = list(low) + [1]
            e = 0
            while len(rest) - 1 >= d:
                q, r = pdivmod(rest, g, p)
                if r:
                    break
                rest, e = q, e + 1
            if e:
                facs.append((tuple(g), e))
    if len(rest) - 1 >= 1:
        facs.append((tuple(rest), 1))
    return facs


def factor_ok(p, facs):
    """The factors multiply back to f mod p, are monic, and those of degree 2 or 3 have no root in F_p."""
    prod = [1]
    for g, e in facs:
        for _ in range(e):
            prod = pmul(prod, list(g), p)
    irr = all(len(g) - 1 == 1 or all(peval(g, x, p) for x in range(p)) for g, e in facs)
    return prod == f_mod(p) and irr and all(g[-1] == 1 for g, e in facs)


def roots_brute(p):
    return [r for r in range(p) if (r**3 - 4*r + 6) % p == 0]


def roots_rabin(p):
    """Distinct roots of f mod p: g = gcd(f, t^p - t) is the product of the distinct linear factors; split g by
    gcd(g, (t + delta)^((p-1)/2) - 1) for delta = 0, 1, 2, ... Each root is re-checked by evaluation."""
    if p < 5:
        return roots_brute(p)
    f = f_mod(p)
    g = pgcd(f, psub(ppowmod([0, 1], p, f, p), [0, 1], p), p)
    out = []

    def split(h):
        d = len(h) - 1
        if d == 0:
            return
        if d == 1:
            out.append((-h[0]) % p)
            return
        for delta in range(p):
            w = pgcd(h, psub(ppowmod([delta, 1], (p - 1)//2, h, p), [1], p), p)
            if 0 < len(w) - 1 < d:
                split(w)
                split(pdivmod(h, w, p)[0])
                return
        check(f"internal: Rabin splitting at p = {p}", False)
    split(g)
    out = sorted(out)
    if len(out) != len(g) - 1 or len(set(out)) != len(out) or any(peval(f, r, p) for r in out):
        check(f"internal: Rabin roots at p = {p}", False)
    return out


def primes_below(n):
    sieve = bytearray([1])*n
    sieve[0:2] = b"\x00\x00"
    for i in range(2, math.isqrt(n - 1) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(range(i*i, n, i)))
    return [i for i in range(n) if sieve[i]]


# ---------------------------------------------------------------- lattices: P(p,g), pi R, Hermite normal form
def in_P(beta, p, g):
    """beta in P(p,g): g divides beta(t) in F_p[t] (remainder by long division)."""
    return not prem(list(beta), list(g), p)


def hnf(rows):
    """Row Hermite normal form of the Z-span of the given vectors of Z^3 (rank 3 required, else None): upper
    triangular, positive diagonal, each entry above a pivot reduced into [0, pivot)."""
    rest = [list(r) for r in rows if any(r)]
    piv = []
    for j in range(3):
        while True:
            nz = [r for r in rest if r[j] != 0]
            if len(nz) <= 1:
                break
            m = min(nz, key=lambda r: abs(r[j]))
            for r in nz:
                if r is not m:
                    q = r[j]//m[j]
                    for k in range(3):
                        r[k] -= q*m[k]
            rest = [r for r in rest if any(r)]
        nz = [r for r in rest if r[j] != 0]
        if not nz:
            return None
        pr = nz[0]
        rest = [r for r in rest if r is not pr]
        piv.append([-x for x in pr] if pr[j] < 0 else pr)
    for j in range(3):
        for i in range(j):
            q = piv[i][j]//piv[j][j]
            piv[i] = [piv[i][k] - q*piv[j][k] for k in range(3)]
    return tuple(tuple(r) for r in piv)


def nullspace_mod_p(A, p):
    """Basis of {x in F_p^3 : A x = 0} (A a list of rows), by row reduction mod p."""
    A = [[x % p for x in row] for row in A]
    pivots, r = [], 0
    for c in range(3):
        pr = next((i for i in range(r, len(A)) if A[i][c]), None)
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        inv = pow(A[r][c], -1, p)
        A[r] = [x*inv % p for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c]:
                co = A[i][c]
                A[i] = [(x - co*y) % p for x, y in zip(A[i], A[r])]
        pivots.append(c)
        r += 1
    basis = []
    for fc in (c for c in range(3) if c not in pivots):
        x = [0, 0, 0]
        x[fc] = 1
        for i, pc in enumerate(pivots):
            x[pc] = (-A[i][fc]) % p
        basis.append(tuple(x))
    return basis


def lattice_P_hnf(p, g):
    """HNF of P(p,g) = ker(Z^3 -> F_p[t]/(g)): generated by p e_i and lifts of an F_p-basis of the null space of the
    map beta -> beta(t) mod g (whose matrix has the columns t^i mod g)."""
    d = len(g) - 1
    cols = []
    for i in range(3):
        r = prem([0]*i + [1], list(g), p)
        cols.append(r + [0]*(d - len(r)))
    A = [[cols[i][row] for i in range(3)] for row in range(d)]
    gens = [tuple(p if i == j else 0 for j in range(3)) for i in range(3)] + nullspace_mod_p(A, p)
    H = hnf(gens)
    if H is None or not all(in_P(h, p, g) for h in H):
        check(f"internal: HNF of P({p}, {gstr(g, p)})", False)
    return H


def principal_hnf(pi):
    """HNF of pi R = Z-span of pi*1, pi*l, pi*l^2."""
    return hnf([rmul(pi, e) for e in E3])


def hnf_index(H):
    return H[0][0]*H[1][1]*H[2][2]


def generates(pi, p, g):
    """P(p,g) = pi R, decided two ways: (i) equal HNFs; (ii) the index argument: pi*1, pi*l, pi*l^2 lie in P(p,g)
    (so pi R is inside P(p,g)) and [Z^3 : pi R] = |N(pi)| equals [Z^3 : P(p,g)] (from the HNF of P)."""
    HP = lattice_P_hnf(p, g)
    Hpi = principal_hnf(pi)
    same_hnf = Hpi is not None and Hpi == HP
    contained = all(in_P(rmul(pi, e), p, g) for e in E3)
    index_eq = abs(N(pi)) == hnf_index(HP) == p**(len(g) - 1) and abs(det3(mult_matrix(pi))) == abs(N(pi))
    return same_hnf and contained and index_eq, (same_hnf, contained, index_eq)


# ---------------------------------------------------------------- univariate polynomials over Q, Sturm, intervals
def q_trim(a):
    a = [Fr(x) for x in a]
    while a and a[-1] == 0:
        a.pop()
    return a


def q_eval(a, x):
    s = Fr(0)
    for co in reversed(a):
        s = s*x + co
    return s


def q_deriv(a):
    return q_trim([i*a[i] for i in range(1, len(a))])


def q_divmod(a, b):
    a, b = q_trim(a), q_trim(b)
    q = [Fr(0)]*max(len(a) - len(b) + 1, 1)
    while a and len(a) >= len(b):
        co = a[-1]/b[-1]
        sh = len(a) - len(b)
        q[sh] = co
        for i in range(len(b)):
            a[sh + i] -= co*b[i]
        a = q_trim(a)
    return q_trim(q), a


def q_gcd(a, b):
    a, b = q_trim(a), q_trim(b)
    while b:
        a, b = b, q_divmod(a, b)[1]
    return [x/a[-1] for x in a] if a else a


def q_from_sympy(expr, var):
    co = sp.Poly(expr, var).all_coeffs()[::-1]
    return q_trim([Fr(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in co])


def sturm_chain(a):
    seq = [q_trim(a), q_deriv(a)]
    while seq[-1]:
        r = q_divmod(seq[-2], seq[-1])[1]
        if not r:
            break
        seq.append([-x for x in r])
    return seq


def sign_var(seq, x):
    sg = [s for s in ((q_eval(p, x) > 0) - (q_eval(p, x) < 0) for p in seq) if s]
    return sum(1 for i in range(len(sg) - 1) if sg[i] != sg[i + 1])


def real_roots(a, width=Fr(1, 2**110)):
    """All real roots of the squarefree part of a, as rational intervals [lo, hi] of width < width, each holding
    exactly one root (lo == hi for an exact rational root found on the way). Sturm's theorem plus bisection."""
    a = q_trim(a)
    sf = q_divmod(a, q_gcd(a, q_deriv(a)))[0]
    seq = sturm_chain(sf)
    bound = 1 + max(abs(x/sf[-1]) for x in sf[:-1]) if len(sf) > 1 else Fr(1)
    lo0, hi0 = -bound - 1, bound + 1
    out = []

    def split_point(lo, hi):
        for k in (Fr(1, 2), Fr(1, 3), Fr(2, 3), Fr(1, 5), Fr(4, 5), Fr(2, 7)):
            m = lo + (hi - lo)*k
            if q_eval(sf, m) != 0:
                return m, False
        return lo + (hi - lo)/2, True

    def rec(lo, hi):          # lo, hi not roots
        n = sign_var(seq, lo) - sign_var(seq, hi)
        if n == 0:
            return
        if n == 1:
            out.append(refine(lo, hi))
            return
        m, exact = split_point(lo, hi)
        if exact:
            check("internal: Sturm split point", False)
        rec(lo, m)
        rec(m, hi)

    def refine(lo, hi):
        slo = q_eval(sf, lo) > 0
        while hi - lo >= width:
            m = (lo + hi)/2
            v = q_eval(sf, m)
            if v == 0:
                return (m, m)
            if (v > 0) == slo:
                lo = m
            else:
                hi = m
        return (lo, hi)
    rec(lo0, hi0)
    return sorted(out), sf


def i_mul(A, B):
    ps = (A[0]*B[0], A[0]*B[1], A[1]*B[0], A[1]*B[1])
    return (min(ps), max(ps))


def i_pow(A, n):
    if n == 0:
        return (Fr(1), Fr(1))
    lo, hi = A[0]**n, A[1]**n
    if n % 2 == 0 and A[0] <= 0 <= A[1]:
        return (Fr(0), max(lo, hi))
    return (min(lo, hi), max(lo, hi))


def i_terms(terms, U, V):
    """Enclosure of sum c u^i v^j over the box U x V (exact rational interval arithmetic)."""
    lo = hi = Fr(0)
    for i, j, c in terms:
        T = i_mul(i_pow(U, i), i_pow(V, j))
        lo += c*T[0] if c >= 0 else c*T[1]
        hi += c*T[1] if c >= 0 else c*T[0]
    return (lo, hi)


class Surd:
    """a + b sqrt(D), with a, b, D rational and D > 0 (D not a rational square when b != 0)."""

    def __init__(self, a, b=0, D=0):
        self.a, self.b, self.D = Fr(a), Fr(b), Fr(D)

    def _same(self, o):
        if not isinstance(o, Surd):
            o = Surd(o, 0, self.D)
        if self.b and o.b and self.D != o.D:
            check("internal: surds over different fields", False)
        return o, (self.D if self.b else o.D)

    def __add__(self, o):
        o, D = self._same(o)
        return Surd(self.a + o.a, self.b + o.b, D)

    def __sub__(self, o):
        o, D = self._same(o)
        return Surd(self.a - o.a, self.b - o.b, D)

    def __neg__(self):
        return Surd(-self.a, -self.b, self.D)

    def __mul__(self, o):
        o, D = self._same(o)
        return Surd(self.a*o.a + self.b*o.b*D, self.a*o.b + self.b*o.a, D)

    def sign(self):
        a, b, D = self.a, self.b, self.D
        sa, sb = (a > 0) - (a < 0), (b > 0) - (b < 0)
        if sb == 0 or D == 0:
            return sa
        if sa == 0 or sa == sb:
            return sb if sa == 0 else sa
        d = a*a - b*b*D
        return sa if d > 0 else (-sa if d < 0 else 0)

    def encl(self, bits=200):
        if self.b == 0 or self.D == 0:
            return (self.a, self.a)
        n, d = self.D.numerator, self.D.denominator
        K = 1 << bits
        r = math.isqrt(n*d*K*K)
        s_lo, s_hi = Fr(r, d*K), Fr(r + 1, d*K)
        if self.b > 0:
            return (self.a + self.b*s_lo, self.a + self.b*s_hi)
        return (self.a + self.b*s_hi, self.a + self.b*s_lo)

    def __str__(self):
        if self.b == 0:
            return str(self.a)
        return f"{self.a} + ({self.b})*sqrt({self.D})"


def is_rat_square(q):
    q = Fr(q)
    if q < 0:
        return False
    n, d = q.numerator, q.denominator
    return math.isqrt(n)**2 == n and math.isqrt(d)**2 == d


def rat_sqrt(q):
    return Fr(math.isqrt(q.numerator), math.isqrt(q.denominator))


# ====================================================================================================================
print("(0) the ring, rebuilt from l^3 = 4l - 6")
U3 = sp.symbols('u0 u1 u2')
V3 = sp.symbols('v0 v1 v2')
W3 = sp.symbols('w0 w1 w2')
check("rmul is commutative and associative, l^3 = 4l - 6 in it (symbolic)",
      all(sp.expand(x - y) == 0 for x, y in zip(rmul(U3, V3), rmul(V3, U3)))
      and all(sp.expand(x - y) == 0 for x, y in zip(rmul(rmul(U3, V3), W3), rmul(U3, rmul(V3, W3))))
      and rmul(rmul(LAM, LAM), LAM) == (-6, 4, 0))
d_, e_, f_ = V3
TEXT_MUL = (U3[0]*d_ - 6*U3[1]*f_ - 6*U3[2]*e_, U3[0]*e_ + U3[1]*d_ + 4*U3[1]*f_ + 4*U3[2]*e_ - 6*U3[2]*f_,
            U3[0]*f_ + U3[1]*e_ + U3[2]*d_ + 4*U3[2]*f_)
check("the product formula of hypothesis_H.md 1.1 (lines 60) equals rmul",
      all(sp.expand(x - y) == 0 for x, y in zip(rmul(U3, V3), TEXT_MUL)))
TEXT_NORM = A_**3 + 8*A_**2*C_ - 4*A_*B_**2 + 18*A_*B_*C_ + 16*A_*C_**2 - 6*B_**3 + 24*B_*C_**2 + 36*C_**3
print("    N(a,b,c) := det of the matrix with columns u*1, u*l, u*l^2 =", NORM_EXPR)
check("this determinant is the norm form of hypothesis_H.md 1.2 (line 73)", sp.expand(NORM_EXPR - TEXT_NORM) == 0)
check("N(uv) = N(u) N(v) (symbolic)",
      sp.expand(sp.Matrix(mult_matrix(rmul(U3, V3))).det() - sp.Matrix(mult_matrix(U3)).det()
                * sp.Matrix(mult_matrix(V3)).det()) == 0)
check("Res_t(t^3 - 4t + 6, a + bt + ct^2) = N(a,b,c) (sympy resultant, an independent route to the norm)",
      sp.expand(sp.resultant(sp.Symbol('t')**3 - 4*sp.Symbol('t') + 6,
                             A_ + B_*sp.Symbol('t') + C_*sp.Symbol('t')**2, sp.Symbol('t')) - NORM_EXPR) == 0)

# ====================================================================================================================
print()
print("(a) the constant of 6.3: max |F| on the box |x0| <= 2, |x1| <= 1, |x2| <= 1/2")
X = sp.symbols('x0 x1 x2')
BASIS = ((1, 0, 0), (1, 1, 0), (-3, 1, 1))        # 1, 1 + l, l^2 + l - 3 (hypothesis_H.md 6.3, line 496)
TW = (Fr(2), Fr(1), Fr(1, 2))                      # the weights t_i (line 496)
check("b = (1, 1 + l, l^2 + l - 3) is a Z-basis of R: the matrix of coordinates has determinant 1",
      det3([[BASIS[j][i] for j in range(3)] for i in range(3)]) == 1)
check("t0 t1 t2 = 2 * 1 * 1/2 = 1", TW[0]*TW[1]*TW[2] == 1)
ELEM = tuple(sum(X[i]*BASIS[i][j] for i in range(3)) for j in range(3))
F_EXPR = sp.expand(NORM_EXPR.subs({A_: ELEM[0], B_: ELEM[1], C_: ELEM[2]}, simultaneous=True))
print("    F(x) = N(x0 b0 + x1 b1 + x2 b2) =", F_EXPR)
FP = sp.Poly(F_EXPR, *X)
F_TERMS = [(m, int(c)) for m, c in zip(FP.monoms(), FP.coeffs())]


def F_int(x):
    return sum(c*x[0]**m[0]*x[1]**m[1]*x[2]**m[2] for m, c in F_TERMS)


def F_rat(x):
    return sum(c*Fr(x[0])**m[0]*Fr(x[1])**m[1]*Fr(x[2])**m[2] for m, c in F_TERMS)


check("F is an integer cubic form, F(-x) = -F(x), and F(x) = N(sum x_i b_i) at 240 integer points by rmul/det",
      FP.is_homogeneous and FP.total_degree() == 3 and all(isinstance(c, int) for _, c in F_TERMS)
      and sp.expand(F_EXPR.subs({v: -v for v in X}, simultaneous=True) + F_EXPR) == 0
      and all(F_int(x) == det3(mult_matrix(tuple(sum(x[i]*BASIS[i][j] for i in range(3)) for j in range(3))))
              for x in itertools.product(range(-3, 3), range(-2, 3), range(-3, 5))))
GRAD = [sp.diff(F_EXPR, v) for v in X]
check("Euler: x0 F_x0 + x1 F_x1 + x2 F_x2 = 3F; so grad F = 0 forces F = 0, and the open box holds no extremum"
      " of |F| other than the value 0", sp.expand(sum(X[i]*GRAD[i] for i in range(3)) - 3*F_EXPR) == 0)
print("    Stratification: the maximum of |F| on the closed box is attained at some x*; on the open face (vertex, edge,")
print("    2-face, interior) containing x*, the partial derivatives in the free coordinates vanish. So max |F| is the")
print("    largest |F| over: the 8 vertices, the critical points inside the 12 edges, those inside the 6 faces, and the")
print("    interior critical points (value 0 by Euler). F is odd, so max F = -min F and max |F| = max F.")
CANDS = []       # (lo, hi, label, exact value or None)
for sg in itertools.product((1, -1), repeat=3):
    pt = tuple(s*TW[i] for i, s in enumerate(sg))
    v = F_rat(pt)
    CANDS.append((v, v, f"vertex {tuple(str(c) for c in pt)}", Surd(v)))
print("    vertices: F =", ", ".join(str(c[3]) for c in CANDS))

# ---- edges
S = sp.Symbol('s')
n_edge_cp = 0
EDGE_LINES = []
for i in range(3):
    others = [j for j in range(3) if j != i]
    for sg in itertools.product((1, -1), repeat=2):
        sub = {X[j]: sg[n]*sp.Rational(TW[j].numerator, TW[j].denominator) for n, j in enumerate(others)}
        phi = q_from_sympy(sp.expand(F_EXPR.subs(sub).subs(X[i], S)), S)
        phi = phi + [Fr(0)]*(4 - len(phi))
        c0, c1, c2, c3 = phi
        Aq, Bq, Cq = 3*c3, 2*c2, c1                       # phi'(s) = A s^2 + B s + C
        Disc = Bq*Bq - 4*Aq*Cq
        roots = []
        if Aq == 0:
            check("internal: the edge cubic has degree 3", False)
        if Disc == 0:
            roots = [Surd(-Bq/(2*Aq))]
        elif Disc > 0 and is_rat_square(Disc):
            sq = rat_sqrt(Disc)
            roots = [Surd((-Bq + sq)/(2*Aq)), Surd((-Bq - sq)/(2*Aq))]
        elif Disc > 0:
            roots = [Surd(-Bq/(2*Aq), 1/(2*Aq), Disc), Surd(-Bq/(2*Aq), -1/(2*Aq), Disc)]
        dphi = [Cq, Bq, Aq]
        rem_lin = q_divmod(phi, dphi)[1]
        rem_lin = rem_lin + [Fr(0)]*(2 - len(rem_lin))
        t = TW[i]
        for r in roots:
            if not ((r - t).sign() < 0 and (r + t).sign() > 0):
                continue
            val = Surd(c0) + r*Surd(c1) + r*r*Surd(c2) + r*r*r*Surd(c3)
            val2 = Surd(rem_lin[0]) + r*Surd(rem_lin[1])
            if (val - val2).sign() != 0 or (r*r*Surd(Aq) + r*Surd(Bq) + Surd(Cq)).sign() != 0:
                check("internal: edge critical value (direct = remainder of phi mod phi')", False)
            lo, hi = val.encl()
            fixed = ", ".join(f"x{j} = {sg[n]*TW[j]}" for n, j in enumerate(others))
            lab = f"edge {fixed}: x{i} = {r}"
            CANDS.append((lo, hi, lab, val))
            n_edge_cp += 1
            EDGE_LINES.append(f"    {lab} ~ {dec(r.encl()[0], 6)}; F = {val} ~ {dec(lo, 10)}")
check(f"the 12 edges: phi'(s) = 0 solved exactly (quadratic surds); {n_edge_cp} critical points inside the open"
      f" edges, each value computed exactly in Q(sqrt(D)) two ways (direct, and phi mod phi')", n_edge_cp > 0)
for ln in EDGE_LINES:
    print(ln)

# ---- faces
n_face_cp = 0
FACE_LINES = []
for k in range(3):
    for sg in (1, -1):
        free = [i for i in range(3) if i != k]
        Uv, Vv = X[free[0]], X[free[1]]
        tu, tv = TW[free[0]], TW[free[1]]
        G = sp.expand(F_EXPR.subs(X[k], sg*sp.Rational(TW[k].numerator, TW[k].denominator)))
        Gu, Gv = sp.diff(G, Uv), sp.diff(G, Vv)
        a = [sp.expand(sp.Poly(Gu, Uv).as_expr().coeff(Uv, i)) for i in range(3)]
        b = [sp.expand(sp.Poly(Gv, Uv).as_expr().coeff(Uv, i)) for i in range(3)]
        ok_shape = (sp.degree(Gu, Uv) == 2 and sp.degree(Gv, Uv) <= 2 and a[2].is_Rational and a[2] != 0
                    and sp.expand(sum(a[i]*Uv**i for i in range(3)) - Gu) == 0
                    and sp.expand(sum(b[i]*Uv**i for i in range(3)) - Gv) == 0)
        R = sp.expand((a[2]*b[0] - a[0]*b[2])**2 - (a[2]*b[1] - a[1]*b[2])*(a[1]*b[0] - a[0]*b[1]))
        L1 = sp.expand(b[2]*a[1] - a[2]*b[1])
        L0 = sp.expand(b[2]*a[0] - a[2]*b[0])
        Rq, L1q, L0q = q_from_sympy(R, Vv), q_from_sympy(L1, Vv), q_from_sympy(L0, Vv)
        coprime = len(q_gcd(Rq, L1q)) == 1
        roots, Rsf = real_roots(Rq)
        name = f"x{k} = {sg*TW[k]}"
        n_sympy = sp.Poly(R, Vv).sqf_part().count_roots()
        check(f"face {name}: G_u quadratic in u = x{free[0]} with constant leading coefficient {a[2]}; R(v) ="
              f" Res_u(G_u, G_v) (two-quadratic formula) equals sympy's resultant; gcd(R, L1) = 1 with"
              f" L1 = b2 a1 - a2 b1; Sturm finds {len(roots)} distinct real roots of R (sympy count_roots: {n_sympy})",
              ok_shape and sp.expand(R - sp.resultant(Gu, Gv, Uv)) == 0 and Rq != [] and coprime
              and len(roots) == len(set(roots)) and len(roots) == n_sympy)
        print(f"    face {name}: R(x{free[1]}) = {sp.factor(R)}")
        GT = [(m[0], m[1], Fr(int(sp.Rational(c).p), int(sp.Rational(c).q)))
              for m, c in zip(sp.Poly(G, Uv, Vv).monoms(), sp.Poly(G, Uv, Vv).coeffs())]
        GuT = [(m[0], m[1], Fr(int(sp.Rational(c).p), int(sp.Rational(c).q)))
               for m, c in zip(sp.Poly(Gu, Uv, Vv).monoms(), sp.Poly(Gu, Uv, Vv).coeffs())]
        GvT = [(m[0], m[1], Fr(int(sp.Rational(c).p), int(sp.Rational(c).q)))
               for m, c in zip(sp.Poly(Gv, Uv, Vv).monoms(), sp.Poly(Gv, Uv, Vv).coeffs())]
        L0T = [(0, i, c) for i, c in enumerate(L0q)]
        L1T = [(0, i, c) for i, c in enumerate(L1q)]
        for (vlo, vhi) in roots:
            # a common root u of G_u(., v*) and G_v(., v*) solves L1(v*) u + L0(v*) = 0 (from b2 G_u - a2 G_v)
            while True:
                Vi = (vlo, vhi)
                l1 = i_terms(L1T, (Fr(0), Fr(0)), Vi)
                l0 = i_terms(L0T, (Fr(0), Fr(0)), Vi)
                if not l1[0] <= 0 <= l1[1]:          # L1(v*) != 0 (gcd(R, L1) = 1); refine until its enclosure is
                    inv = (1/l1[1], 1/l1[0])          # free of 0
                    Ui = i_mul((-l0[1], -l0[0]), inv)
                    v_in = vlo >= -tv and vhi <= tv
                    v_out = vhi < -tv or vlo > tv
                    u_in = Ui[0] >= -tu and Ui[1] <= tu
                    u_out = Ui[1] < -tu or Ui[0] > tu
                    if (v_in or v_out) and (u_in or u_out or v_out):
                        break
                if vhi - vlo < Fr(1, 2**400):
                    check(f"internal: face {name} point undecided (on the boundary?)", False)
                # refine the root interval of R's squarefree part by bisection
                slo = q_eval(Rsf, vlo) > 0
                for _ in range(64):
                    m = (vlo + vhi)/2
                    vm = q_eval(Rsf, m)
                    if vm == 0:
                        vlo = vhi = m
                        break
                    if (vm > 0) == slo:
                        vlo = m
                    else:
                        vhi = m
            inside = v_in and u_in
            gu, gv = i_terms(GuT, Ui, Vi), i_terms(GvT, Ui, Vi)
            Fi = i_terms(GT, Ui, Vi)
            crit = gu[0] <= 0 <= gu[1] and gv[0] <= 0 <= gv[1]
            if not crit:
                check(f"internal: face {name}: G_u, G_v enclose 0 at the point found", False)
            pos = f"x{free[0]} ~ {dec(Ui[0], 6)}, x{free[1]} ~ {dec(Vi[0], 6)}"
            if inside:
                n_face_cp += 1
                CANDS.append((Fi[0], Fi[1], f"face {name}: {pos}", None))
                FACE_LINES.append(f"    face {name}: critical point {pos}, inside; F in [{dec(Fi[0], 10)},"
                                  f" {dec(Fi[1], 10)}]")
            else:
                FACE_LINES.append(f"    face {name}: critical point {pos}, outside the face (not a candidate)")
check(f"the 6 faces: every real critical point located exactly (Sturm isolation to width < 2^-110, coordinates and"
      f" value by exact interval arithmetic, G_u and G_v enclose 0); {n_face_cp} lie in the closed faces",
      n_face_cp > 0)
for ln in FACE_LINES:
    print(ln)

# ---- the maximum
best = max(CANDS, key=lambda c: c[0])
worst = min(CANDS, key=lambda c: c[1])
others_max = [c for c in CANDS if c is not best]
others_min = [c for c in CANDS if c is not worst]
M_EXACT = best[3]
check(f"max F over all {len(CANDS)} candidates is attained only at [{best[2]}]: its lower bound exceeds every other"
      f" candidate's upper bound", M_EXACT is not None and all(c[1] < best[0] for c in others_max))
check(f"min F is attained only at [{worst[2]}] and equals -max F exactly",
      worst[3] is not None and all(c[0] > worst[1] for c in others_min) and (worst[3] + M_EXACT).sign() == 0)
M_LO, M_HI = M_EXACT.encl(240)
C0_TEXT = Fr(535, 24)
print(f"    TRUE MAXIMUM of |F| on the box: M = {M_EXACT} = {dec(M_LO, 15)}..., attained at +-({best[2]})")
check(f"M < 535/24 = {dec(C0_TEXT, 4)} exactly (535/24 - M is a positive surd); 535/24 - M ~ {dec(C0_TEXT - M_HI, 4)};"
      f" M < 20, so |N(beta)| <= M p would even give k <= 19",
      (Surd(C0_TEXT) - M_EXACT).sign() > 0 and (Surd(20) - M_EXACT).sign() > 0)

# ---- (a2) branch and bound with the exact Taylor expansion
W = sp.symbols('w0 w1 w2')
HESS = [[sp.diff(F_EXPR, X[i], X[j]) for j in range(3)] for i in range(3)]
taylor = (F_EXPR + sum(GRAD[i]*W[i] for i in range(3)) + sp.Rational(1, 2)*sum(HESS[i][j]*W[i]*W[j]
                                                                                 for i in range(3) for j in range(3))
          + F_EXPR.subs({X[i]: W[i] for i in range(3)}, simultaneous=True))
check("exact Taylor identity of the cubic form: F(v + w) = F(v) + grad F(v).w + w^T H(v) w / 2 + F(w)",
      sp.expand(F_EXPR.subs({X[i]: X[i] + W[i] for i in range(3)}, simultaneous=True) - taylor) == 0)


def terms3(expr):
    P = sp.Poly(expr, *X)
    return [(m, int(c)) for m, c in zip(P.monoms(), P.coeffs())]


GRAD_T = [terms3(g) for g in GRAD]
HESS_T = [[terms3(HESS[i][j]) for j in range(3)] for i in range(3)]


def ev3(tt, v):
    return sum(c*v[0]**m[0]*v[1]**m[1]*v[2]**m[2] for m, c in tt)


def bb_bound(target, max_nodes=200000):
    """Branch and bound: True if |F| <= target is proved on every sub-box (bound |F(v)| + sum |g_i| h_i +
    sum |H_ij| h_i h_j / 2 + sum |n_m| h^m at the centre v, half-widths h); False with a point if some box centre
    has |F| > target (the bound is then false). Boxes are taken largest |F(centre)| first; whether a box is split
    depends only on the box, so a proof uses the same boxes in any order."""
    heap, cnt = [], itertools.count()

    def push(lo, hi, dpt):
        v = tuple((lo[i] + hi[i])/2 for i in range(3))
        Fv = ev3(F_TERMS, v)
        heapq.heappush(heap, (-abs(Fv), next(cnt), lo, hi, dpt, v, Fv))
    push((Fr(-2), Fr(-1), Fr(-1, 2)), (Fr(2), Fr(1), Fr(1, 2)), 0)
    nodes, maxd = 0, 0
    while heap:
        _, _, lo, hi, dpt, v, Fv = heapq.heappop(heap)
        nodes += 1
        if nodes > max_nodes:
            return None, nodes, maxd, None
        h = tuple((hi[i] - lo[i])/2 for i in range(3))
        if abs(Fv) > target:
            return False, nodes, dpt, v
        bnd = (abs(Fv) + sum(abs(ev3(GRAD_T[i], v))*h[i] for i in range(3))
               + sum(abs(ev3(HESS_T[i][j], v))*h[i]*h[j] for i in range(3) for j in range(3))/2
               + sum(abs(c)*h[0]**m[0]*h[1]**m[1]*h[2]**m[2] for m, c in F_TERMS))
        if bnd <= target:
            continue
        maxd = max(maxd, dpt)
        rel = [(hi[i] - lo[i])/TW[i] for i in range(3)]
        ax = max(range(3), key=lambda q: rel[q])
        mid = (lo[ax] + hi[ax])/2
        hi1 = list(hi)
        hi1[ax] = mid
        lo2 = list(lo)
        lo2[ax] = mid
        push(lo, tuple(hi1), dpt + 1)
        push(tuple(lo2), hi, dpt + 1)
    return True, nodes, maxd, None


ok1, n1, d1, _ = bb_bound(C0_TEXT)
check(f"(a2) branch and bound proves |F| <= 535/24 on the box ({n1} boxes, depth {d1})", ok1 is True)
tight = M_HI + Fr(1, 10**12)
ok2, n2, d2, _ = bb_bound(tight)
check(f"(a2) branch and bound proves |F| <= M + 10^-12 on the box ({n2} boxes, depth {d2}): the maximum M of (a1)"
      f" is confirmed from above by a second method", ok2 is True)
ok3, n3, d3, pt3 = bb_bound(M_LO - Fr(1, 10**6))
control(f"branch and bound asked to prove the false bound |F| <= M - 10^-6 finds a box centre with larger |F|"
        f" (after {n3} boxes, at x ~ ({', '.join(dec(c, 6) for c in pt3) if pt3 else '-'}))", ok3 is False)

# ---- (a3) exact integer grid, as a control
NG = 60
ii = np.arange(-NG, NG + 1, dtype=np.int64)
I0, I1, I2 = np.meshgrid(4*ii, 2*ii, ii, indexing='ij')      # x = (4i, 2j, k)/(2 NG): the grid on the box
GV = np.zeros(I0.shape, dtype=np.int64)
for m, c in F_TERMS:
    GV += c*I0**m[0]*I1**m[1]*I2**m[2]
gmax = int(GV.max())
gmin = int(GV.min())
scale = (2*NG)**3
check(f"(a3) control: on the grid x = (2i/{NG}, j/{NG}, k/{2*NG}), |i|,|j|,|k| <= {NG} (int64, exact), max F ="
      f" {gmax}/{scale} ~ {dec(Fr(gmax, scale), 6)} <= M and min F = -max F; the gap M - grid max ~"
      f" {dec(M_LO - Fr(gmax, scale), 6)}",
      (M_EXACT - Surd(Fr(gmax, scale))).sign() >= 0 and gmin == -gmax and abs(GV).max() < 2**40)
control("the grid maximum exceeds the false value M - 10^-2 (so the grid would expose an underestimate of M)",
        Fr(gmax, scale) > M_HI - Fr(1, 100))

# ====================================================================================================================
print()
print("(b) the box principle 6.3 for every prime p < 10^5")


def icbrt(n):
    """Largest integer B >= 0 with B^3 <= n, by bisection on integers."""
    lo, hi = 0, 1
    while hi**3 <= n:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi)//2
        if mid**3 <= n:
            lo = mid
        else:
            hi = mid
    return lo


def box_sizes(p):
    """B_i = floor(t_i p^(1/3)) for t = (2, 1, 1/2): the largest B with B^3 <= t_i^3 p, i.e. B^3 <= 8p, p, p/8."""
    return (icbrt(8*p), icbrt(p), icbrt(p//8))


T3 = tuple(w**3 for w in TW)
VAL_B = lambda r, p: (1 % p, (1 + r) % p, (r*r + r - 3) % p)     # b_i(r) mod p: 1, 1 + r, r^2 + r - 3


def pigeonhole(p, r, B, vals=VAL_B):
    """Two points c != d of the box 0 <= c_i <= B_i with the same residue sum c_i b_i(r) mod p.
    The box is cut into rows (c1, c2) fixed, c0 = 0..B0; a row occupies the B0 + 1 consecutive residues s, ..., s + B0
    (s = c1 b1(r) + c2 b2(r)). Sort the row starts cyclically: if every cyclic gap were > B0, the gaps would add up to
    at least (number of rows)(B0 + 1) = prod(B_i + 1) > p, but they add up to p. So some gap is <= B0, and the point
    (gap, c1, c2) of the earlier row has the residue of the point (0, d1, d2) of the next one."""
    v = vals(r, p)
    rows = sorted(((c1*v[1] + c2*v[2]) % p, c1, c2) for c1 in range(B[1] + 1) for c2 in range(B[2] + 1))
    n = len(rows)
    for i in range(n):
        s1, c1, c2 = rows[i]
        s2, d1, d2 = rows[(i + 1) % n]
        gap = (s2 - s1) % p
        if gap <= B[0]:
            c, d = (gap, c1, c2), (0, d1, d2)
            if c != d and (sum(ci*vi for ci, vi in zip(c, v)) - sum(di*vi for di, vi in zip(d, v))) % p == 0:
                return tuple(ci - di for ci, di in zip(c, d))
    return None


def beta_of(x):
    return tuple(sum(x[i]*BASIS[i][j] for i in range(3)) for j in range(3))


PR5 = primes_below(10**5)
check(f"own sieve: {len(PR5)} primes below 10^5 (sympy primepi: {sp.primepi(10**5 - 1)})",
      len(PR5) == sp.primepi(10**5 - 1))
ok_ineq = all(Fr(Bi**3) <= t3*p < Fr((Bi + 1)**3) for p in PR5 for Bi, t3 in zip(box_sizes(p), T3))
ok_prod = all((box_sizes(p)[0] + 1)*(box_sizes(p)[1] + 1)*(box_sizes(p)[2] + 1) > p for p in PR5)
check("for every prime p < 10^5: B_i^3 <= t_i^3 p < (B_i + 1)^3 (so B_i = floor(t_i p^(1/3)), exactly)", ok_ineq)
check("for every prime p < 10^5: prod(B_i + 1) > p (the lemma's count; in general (B_i + 1) > t_i p^(1/3) and"
      " t0 t1 t2 = 1)", ok_prod)
small_box_fail = [p for p in PR5 if box_sizes(p)[0]*box_sizes(p)[1]*box_sizes(p)[2] <= p]
control(f"with the boxes shrunk to B_i - 1, the count prod(B_i) > p fails for {len(small_box_fail)} of the"
        f" {len(PR5)} primes", len(small_box_fail) > 0)
print("    The bound: |x_i| <= B_i <= t_i p^(1/3), so x / p^(1/3) lies in the box of (a) and |N(beta)| = |F(x)| <=")
print("    M p < (535/24) p < 23p by homogeneity; p | N(beta) (beta in P(p, t - r)) and beta != 0 give")
print("    |N(beta)| = kp with 1 <= k <= 22 (with the true M of (a): k <= 19).")
T_B = time.time()
n_roots, ks, ROOTS5 = 0, {}, {}
bad_b = []
for p in PR5:
    rts = roots_rabin(p)
    ROOTS5[p] = rts
    B = box_sizes(p)
    for r in rts:
        n_roots += 1
        x = pigeonhole(p, r, B)
        if x is None:
            bad_b.append((p, r, "no collision"))
            continue
        beta = beta_of(x)
        nb = N(beta)
        ok = (beta != (0, 0, 0) and all(abs(x[i]) <= B[i] for i in range(3)) and in_P(beta, p, (-r % p, 1))
              and nb == F_int(x) and nb % p == 0 and 1 <= abs(nb)//p <= 22 and Fr(abs(nb)) <= C0_TEXT*p)
        if not ok:
            bad_b.append((p, r, x, nb))
        ks[abs(nb)//p] = ks.get(abs(nb)//p, 0) + 1
check(f"roots of f mod p by Rabin's method agree with brute force for every prime p < 3000",
      all(ROOTS5[p] == roots_brute(p) for p in PR5 if p < 3000))
check(f"for each of the {n_roots} roots r of f mod p (p < 10^5 prime) the pigeonhole gives x != 0 with |x_i| <= B_i,"
      f" beta = sum x_i b_i in P(p, t - r) (remainder of beta(t) by t - r is 0), N(beta) = F(x), p | N(beta),"
      f" 1 <= |N(beta)|/p <= 22 ({time.time() - T_B:.1f} s)", not bad_b)
print(f"    k = |N(beta)|/p over all {n_roots} roots: {dict(sorted(ks.items()))}; largest k = {max(ks)}")
p_c, r_c = 99991, None
for p in reversed(PR5):
    if ROOTS5[p]:
        p_c, r_c = p, ROOTS5[p][0]
        break
xw = pigeonhole(p_c, r_c, box_sizes(p_c), vals=lambda r, p: (1 % p, (1 + r) % p, (r*r + r + 3) % p))
control(f"pigeonhole run at p = {p_c}, r = {r_c} with a wrong residue for b2 (r^2 + r + 3): the beta found is not in"
        f" P(p, t - r)", xw is not None and not in_P(beta_of(xw), p_c, (-r_c % p_c, 1)))

# ====================================================================================================================
print()
print("(c) lemma 6.2 for every prime p <= 60, by all p^3 residues; f mod p factored by exhaustive search")
FACT = {}
T_C = time.time()
for p in primes_below(61):
    facs = factor_exhaustive(p)
    FACT[p] = facs
    if not factor_ok(p, facs):
        check(f"p = {p}: factorization of f", False)
    n_div = n_in = n_bad = 0
    for beta in itertools.product(range(p), repeat=3):
        div = N(beta) % p == 0
        inside = any(in_P(beta, p, g) for g, e in facs)
        n_div += div
        n_in += inside
        n_bad += div != inside
    units = 1
    for g, e in facs:
        dg = len(g) - 1
        units *= (p**dg - 1)*p**(dg*(e - 1))
    fs = " * ".join(f"({gstr(g, p)})" + (f"^{e}" if e > 1 else "") for g, e in facs)
    check(f"p = {p:2d}: f = {fs} mod {p}; of the {p**3} residues, {n_div} have p | N(beta), and exactly these lie in a"
          f" P(p,g); count = p^3 - #(R/pR)* = {p**3 - units}", n_bad == 0 and n_div == p**3 - units)
print(f"    ({time.time() - T_C:.1f} s)")
facs179 = factor_exhaustive(179)
FACT[179] = facs179
check(f"p = 179: f = {' * '.join(f'({gstr(g, 179)})' + (f'^{e}' if e > 1 else '') for g, e in facs179)} mod 179"
      f" (exhaustive search), as hypothesis_H.md 6.2 and 6.4 say: (t - 47)^2 (t - 85)",
      factor_ok(179, facs179) and sorted(facs179) == sorted([((132, 1), 2), ((94, 1), 1)]))
ar = np.arange(179, dtype=np.int64)
Bm, Cm = np.meshgrid(ar, ar, indexing='ij')
ok179 = True
for a in range(179):
    Nv = (a**3 + 8*a*a*Cm - 4*a*Bm**2 + 18*a*Bm*Cm + 16*a*Cm**2 - 6*Bm**3 + 24*Bm*Cm**2 + 36*Cm**3) % 179
    inP = ((a + 47*Bm + 47*47*Cm) % 179 == 0) | ((a + 85*Bm + 85*85*Cm) % 179 == 0)
    ok179 &= bool(np.all((Nv == 0) == inP))
check("p = 179: for all 179^3 residues (numpy int64, exact), 179 | N(beta) iff beta(47) = 0 or beta(85) = 0 mod 179,"
      " i.e. beta in P(179, t - 47) or P(179, t - 85)", ok179)
check("6.2: (36t + 81) f + (-12t^2 - 27t + 32) f' = 358 = 2 * 179 (sympy), so a repeated root forces p in {2, 179};"
      " and among primes p < 10^5 exactly p = 2 and p = 179 have a root r with f(r) = f'(r) = 0 mod p",
      sp.expand((36*sp.Symbol('t') + 81)*(sp.Symbol('t')**3 - 4*sp.Symbol('t') + 6)
                + (-12*sp.Symbol('t')**2 - 27*sp.Symbol('t') + 32)*(3*sp.Symbol('t')**2 - 4) - 358) == 0
      and [p for p in PR5 if any((3*r*r - 4) % p == 0 for r in ROOTS5[p])] == [2, 179])
n_wrong = sum(1 for beta in itertools.product(range(7), repeat=3)
              if (N(beta) % 7 == 0) != any(in_P(beta, 7, g) for g, e in [((4, 1), 1), ((6, 3, 1), 1)]))
control(f"p = 7 with a wrong factor list [t - 3, t^2 + 3t + 6] (the true quadratic is t^2 + 3t + 5): the 'iff' fails"
        f" for {n_wrong} residues", n_wrong > 0)

# ====================================================================================================================
print()
print("(d) the table 6.4 (hypothesis_H.md lines 504-513), typed in from the text")
PI2, P3A, P3B, P3C = (-4, 1, 1), (-1, 1, 0), (5, 2, 0), (3, -3, 1)          # 1.5, lines 101-104
P179, Q179 = (107, -114, 45), (-23, 11, 8)                                   # 1.5 line 105; 6.4 line 513; task (f)
EPS, EPS_INV = (-19, -5, 1), (95, -101, 40)                                  # 1.5, line 100
TEXT_TABLE = [   # (p, g as coefficients mod p with constant term first, generator, name)
    (2, (0, 1), PI2, "PI2"),
    (3, (0, 1), P3C, "P3C"), (3, (2, 1), P3A, "P3A"), (3, (1, 1), P3B, "P3B"),
    (5, (1, 1, 0, 1), (5, 0, 0), "5"),
    (7, (4, 1), (-1, 2, 1), "pi_7 = l^2 + 2l - 1"), (7, (5, 3, 1), (5, -4, 1), "7/pi_7 = l^2 - 4l + 5"),
    (11, (5, 1), (-1, -3, -1), "pi_11 = -l^2 - 3l - 1"), (11, (10, 6, 1), (7, -9, 4), "11/pi_11 = 4l^2 - 9l + 7"),
    (13, (6, 9, 0, 1), (13, 0, 0), "13"), (17, (6, 13, 0, 1), (17, 0, 0), "17"), (19, (6, 15, 0, 1), (19, 0, 0), "19"),
    (179, (132, 1), P179, "P179"), (179, (94, 1), Q179, "Q179 = 8l^2 + 11l - 23"),
]
TEXT_PATTERN = {2: [((0, 1), 3)], 3: [((0, 1), 1), ((1, 1), 1), ((2, 1), 1)], 5: [((1, 1, 0, 1), 1)],
                7: [((4, 1), 1), ((5, 3, 1), 1)], 11: [((5, 1), 1), ((10, 6, 1), 1)],
                13: [((6, 9, 0, 1), 1)], 17: [((6, 13, 0, 1), 1)], 19: [((6, 15, 0, 1), 1)],
                179: [((94, 1), 1), ((132, 1), 2)]}
for p in (2, 3, 5, 7, 11, 13, 17, 19, 179):
    check(f"p = {p}: the factorization of f mod p in the text's table equals the one found in (c):"
          f" {' * '.join(f'({gstr(g, p)})' + (f'^{e}' if e > 1 else '') for g, e in sorted(FACT[p]))}",
          sorted(FACT[p]) == sorted(TEXT_PATTERN[p]))
TABLE = {}
for p, g, pi, name in TEXT_TABLE:
    ok, parts = generates(pi, p, g)
    H = lattice_P_hnf(p, g)
    check(f"P({p}, {gstr(g, p)}) = pi R for pi = {name} = {pi}: pi in P {in_P(pi, p, g)}, N(pi) = {N(pi)} ="
          f" +-{p}^{len(g) - 1}, HNF(P) = {H} equals HNF(pi R), and [Z^3 : P] = {hnf_index(H)} = |N(pi)| with"
          f" pi*1, pi*l, pi*l^2 in P (index argument)",
          ok and in_P(pi, p, g) and abs(N(pi)) == p**(len(g) - 1))
    TABLE[(p, g)] = pi
check("the table covers every prime p <= 22 and every monic irreducible factor g of f mod p, and p = 179",
      set(TABLE) == {(p, g) for p in (2, 3, 5, 7, 11, 13, 17, 19, 179) for g, e in FACT[p]})
check("7/pi_7 and 11/pi_11 as listed: pi_7 (l^2 - 4l + 5) = 7 and pi_11 (4l^2 - 9l + 7) = 11",
      rmul((-1, 2, 1), (5, -4, 1)) == (7, 0, 0) and rmul((-1, -3, -1), (7, -9, 4)) == (11, 0, 0))
control("pi = 2 PI2 (in P(2,t), norm -16) is rejected for P(2, t)", not generates(rmul((2, 0, 0), PI2), 2, (0, 1))[0])
control("pi_7 + 7 = (6, 2, 1) (in P(7, t - 3)) is rejected: N = "
        f"{N((6, 2, 1))}", not generates((6, 2, 1), 7, (4, 1))[0])
control("P3A is rejected for P(3, t) (not a member)", not generates(P3A, 3, (0, 1))[0])
control("pi_7 (norm 7) is rejected for P(7, t^2 + 3t + 5)", not generates((-1, 2, 1), 7, (5, 3, 1))[0])

# ====================================================================================================================
print()
print("(e) the descent 6.5 for every prime 23 <= p < 2000 with a root of f mod p (continued to p < 3000)")
DESC_USED = {}
REC = {}
ks_e = {}
counts = {"inert": 0, "split": 0, "lin*quad": 0, "double": 0}
lat = {2000: 0, 3000: 0}
box_built = {2000: 0, 3000: 0}
bad_e = []
T_E = time.time()


def least_prime_factor(k):
    q = 2
    while k % q:
        q += 1
    return q


def descend(p, r):
    """6.5 for one root: the box principle, then division by table generators until |N| = p."""
    g = ((-r) % p, 1)
    x = pigeonhole(p, r, box_sizes(p))
    beta = beta_of(x)
    k0 = abs(N(beta))//p
    if not (in_P(beta, p, g) and N(beta) % p == 0 and 1 <= k0 <= 22):
        return None, k0, []
    steps = []
    while abs(N(beta)) != p:
        k = abs(N(beta))//p
        q = least_prime_factor(k)
        hits = [(gq, TABLE[(q, gq)]) for gq, e in FACT[q] if in_P(beta, q, gq)]
        if not hits or q > 22:
            return None, k0, steps
        gq, piq = hits[0]
        nb = rdiv(beta, piq)
        if nb is None or not in_P(nb, p, g) or abs(N(nb))*abs(N(piq)) != abs(N(beta)):
            return None, k0, steps
        beta = nb
        steps.append((q, gq))
    return beta, k0, steps


for p in primes_below(3000):
    if p < 23:
        continue
    rts = ROOTS5[p]
    f = f_mod(p)
    if len(rts) == 3:
        kind = "split"
    elif len(rts) == 0:
        kind = "inert"
    elif len(rts) == 1:
        kind = "lin*quad"
    else:
        kind = "double"
    counts[kind] += 1
    gens = {}
    for r in rts:
        gen, k0, steps = descend(p, r)
        if gen is None:
            bad_e.append((p, r, "descent"))
            continue
        ok, parts = generates(gen, p, ((-r) % p, 1))
        if not ok:
            bad_e.append((p, r, "generates", parts))
        gens[r] = gen
        REC[(p, r)] = gen
        ks_e[k0] = ks_e.get(k0, 0) + 1
        for st in steps:
            DESC_USED[st] = DESC_USED.get(st, 0) + 1
    nl = 0
    if kind == "split":
        nl = 3
    elif kind == "double":
        prod = pmul(pmul([(-rts[0]) % p, 1], [(-rts[0]) % p, 1], p), [(-rts[1]) % p, 1], p)
        prod2 = pmul(pmul([(-rts[1]) % p, 1], [(-rts[1]) % p, 1], p), [(-rts[0]) % p, 1], p)
        if prod != f and prod2 != f:
            bad_e.append((p, "double root pattern"))
        nl = 2
    elif kind == "lin*quad":
        r = rts[0]
        gq, rem = pdivmod(f, [(-r) % p, 1], p)
        if rem or any(peval(gq, z, p) == 0 for z in range(p)):
            bad_e.append((p, "quadratic factor"))
        if r in gens:
            piq = rdiv((p, 0, 0), gens[r])
            okq = piq is not None and generates(piq, p, tuple(gq))[0]
            if not okq:
                bad_e.append((p, "P(p, quadratic) = (p/pi_r) R"))
        nl = 2
    else:
        okf = generates((p, 0, 0), p, tuple(f))[0]
        if not okf:
            bad_e.append((p, "P(p, f) = pR"))
        nl = 1
    for lim in (2000, 3000):
        if p < lim:
            lat[lim] += nl
            box_built[lim] += len(rts)
check(f"every root r of f mod p, 23 <= p < 2000 ({box_built[2000]} degree-one lattices): the box principle gives beta"
      f" in P(p, t - r) with k <= 22; dividing by the table's generators (each step: exact quotient in R, still in"
      f" P(p, t - r), |N| divided by q^deg g) ends at |N| = p, and that element generates P(p, t - r): equal HNFs and"
      f" the index argument", not [b for b in bad_e if b[0] < 2000])
check(f"continued to 23 <= p < 3000: {box_built[3000]} degree-one lattices, all generated; the quadratic factor"
      f" of each partially split p is generated by p/pi_r, each inert p by p; prime lattices in all:"
      f" {lat[3000]} (the text, line 526, and 08_ufd_box.out.txt: 781, of which 428 by the box principle)",
      not bad_e and lat[3000] == 781 and box_built[3000] == 428)
print(f"    primes 23 <= p < 3000 by pattern: {counts}; prime lattices for p < 2000: {lat[2000]}"
      f" ({box_built[2000]} by the box principle); k after the box principle: {dict(sorted(ks_e.items()))};"
      f" {time.time() - T_E:.1f} s")
SIX = {(2, (0, 1)), (3, (0, 1)), (3, (1, 1)), (3, (2, 1)), (7, (4, 1)), (11, (5, 1))}
print("    lattices used by the descent (q, g): " + ", ".join(f"P({q}, {gstr(g, q)}) x{n}" for (q, g), n in
                                                            sorted(DESC_USED.items())))
check("the descent used only the six degree-one lattices of norm <= 22 named in 6.5 (above 2, 3 (three), 7, 11)",
      set(DESC_USED) <= SIX)
check("p = 179 through the descent: the generators built for the roots 47 and 85 are associates of P179 and Q179",
      associates(REC[(179, 47)], P179) and associates(REC[(179, 85)], Q179))
gen23 = REC[(23, ROOTS5[23][0])]
wrong = rdiv(rmul(gen23, (-1, 2, 1)), P3A)
control("dividing an element of P(7, t - 3) * P(23, ...) by P3A (wrong table lattice) has no quotient in R",
        wrong is None)

# ---- (e7) cross-check with the generators printed by 08_ufd_box.py for p < 200 (read as data)
OUT08 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "08_ufd_box.out.txt")


def parse_l(s):
    s = s.replace(" ", "").replace("-", "+-")
    a = b = c = 0
    for part in (x for x in s.split("+") if x):
        if part.endswith("l^2"):
            co = part[:-3]
            c += 1 if co == "" else (-1 if co == "-" else int(co))
        elif part.endswith("l"):
            co = part[:-1]
            b += 1 if co == "" else (-1 if co == "-" else int(co))
        else:
            a += int(part)
    return (a, b, c)


rx = re.compile(r"p=(\d+)\s+r=(\d+)\s+box=\((\d+), (\d+), (\d+)\)\s+pts=\d+\s+x=\((-?\d+), (-?\d+), (-?\d+)\)\s+"
                r"k=(\d+)\s+divided by: .*; pi = (.*), N = (-?\d+)")
rows08 = []
with open(OUT08, encoding="utf-8") as fh:
    for line in fh:
        m = rx.search(line)
        if m:
            rows08.append(m.groups())
mine_lt200 = sorted((p, r) for p in primes_below(200) if p >= 23 for r in ROOTS5[p])
ok_x = True
ok_assoc = True
for (ps, rs, b0, b1, b2, y0, y1, y2, kk, pis, ns) in rows08:
    p, r = int(ps), int(rs)
    their = parse_l(pis)
    xx = (int(y0), int(y1), int(y2))
    Bt = (int(b0), int(b1), int(b2))
    ok_x &= (Bt == box_sizes(p) and all(abs(xx[i]) <= Bt[i] for i in range(3))
             and abs(N(beta_of(xx))) == int(kk)*p and in_P(beta_of(xx), p, ((-r) % p, 1)))
    ok_assoc &= N(their) == int(ns) and associates(REC[(p, r)], their)
check(f"08_ufd_box.out.txt lists {len(rows08)} roots for 23 <= p < 200, the same (p, r) as found here; its boxes"
      f" equal B_i here, its x give |N(beta)| = kp with beta in P(p, t - r); its generators have the printed norms",
      sorted((int(g[0]), int(g[1])) for g in rows08) == mine_lt200 and ok_x)
check("each generator printed by 08 for p < 200 is an associate of the one built here (quotient in R of norm +-1)",
      ok_assoc)

# ====================================================================================================================
print()
print("(f) the S-primes")
for name, pi, p, g in (("PI2", PI2, 2, (0, 1)), ("P3C", P3C, 3, (0, 1)), ("P3A", P3A, 3, (2, 1)),
                       ("P3B", P3B, 3, (1, 1)), ("P179", P179, 179, (132, 1)), ("Q179", Q179, 179, (94, 1))):
    H = lattice_P_hnf(p, g)
    ok, parts = generates(pi, p, g)
    check(f"{name} = {pi} = {tstr(pi)}, N = {N(pi)}: generates P({p}, {gstr(g, p)}) (HNF {H})", ok)
check("these are all the prime lattices above 2, 3, 179: f = t^3 mod 2, t(t - 1)(t - 2) mod 3, (t - 47)^2 (t - 85)"
      " mod 179 (from (c))",
      sorted(FACT[2]) == [((0, 1), 3)] and sorted(FACT[3]) == [((0, 1), 1), ((1, 1), 1), ((2, 1), 1)]
      and sorted(FACT[179]) == [((94, 1), 1), ((132, 1), 2)])
check("EPS EPS_INV = 1, N(EPS) = 1; PI2^3 = -2 EPS (so 2 = unit * PI2^3); P3A P3B P3C = 3; PI2^2 P179 = 3l^2 - 4",
      rmul(EPS, EPS_INV) == (1, 0, 0) and N(EPS) == 1 and rpow(PI2, 3) == tuple(-2*x for x in EPS)
      and rmul(rmul(P3A, P3B), P3C) == (3, 0, 0) and rmul(rmul(PI2, PI2), P179) == (-4, 0, 3))
prod179 = rmul(rmul(P179, P179), Q179)
u179 = rdiv((179, 0, 0), prod179)
ident = None
for kexp in range(0, 41):
    for sgn in (1, -1):
        for base, nm in ((EPS, "EPS"), (EPS_INV, "EPS^-1")):
            if u179 is not None and tuple(sgn*x for x in rpow(base, kexp)) == u179:
                ident = f"{'-' if sgn < 0 else ''}{nm}^{kexp}"
check(f"P179^2 * Q179 = {prod179}; 179 / (P179^2 Q179) = {u179} has norm {N(u179) if u179 else None}, a unit,"
      f" = {ident}: 179 = (unit) P179^2 Q179", u179 is not None and is_unit(u179) and ident is not None)
check("N(P179)^2 N(Q179) = 179^3 = N(179), and P179, Q179 are not associates (different lattices)",
      N(P179)**2*N(Q179) == 179**3 and not associates(P179, Q179))
control("P179 does not generate P(179, t - 85)", not generates(P179, 179, (94, 1))[0])
control("Q179 does not generate P(179, t - 47)", not generates(Q179, 179, (132, 1))[0])
control("P3A does not generate P(3, t - 2)", not generates(P3A, 3, (1, 1))[0])

print()
print(f"SUMMARY: all {N_CHECKS} checks PASS (controls included); max |F| on the box = {M_EXACT} ~ {dec(M_LO, 12)}"
      f" < 535/24 ~ 22.2917; lemma 6.3 holds for all {len(PR5)} primes < 10^5 ({n_roots} roots, largest k = {max(ks)});"
      f" lemma 6.2 for p <= 60 and p = 179; table 6.4 and the S-primes confirmed; descent confirmed for"
      f" 23 <= p < 3000 ({time.time() - T_START:.1f} s)")
