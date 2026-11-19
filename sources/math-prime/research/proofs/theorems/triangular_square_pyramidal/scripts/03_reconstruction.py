#!/usr/bin/env python3
"""03_reconstruction.py -- from the Thue solutions back to the integral points of 8604.a1, and on to
n(n+1)/2 = m(m+1)(2m+1)/6 (numbers both triangular and square-pyramidal) and to the model 3Y^2 = X^3 - X + 3.

Curve 8604.a1 (LMFDB): y^2 = x^3 - 9x + 81.  It is already monic: E: Y^2 = X^3 - 9X + 81 with X = x, Y = y.
K = Q(th), th^3 - 9 th + 81 = 0.

(i)   THE REDUCTION.  Every item of the 2-descent that the reduction uses is read from
      research/proofs/targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt (section CURVE 8604.a1) and
      re-checked here by exact arithmetic: the maximal order O_K = Z + Z th/3 + Z th^2/9, the four classes and their
      distinctness, the line identity behind the homomorphism, the lattices Xi_mu = {xi in K : mu xi^2 in O_K}, the
      forms A0, A1, A2, the conic, the completeness of its parametrisation, H = -A1(X(m,n)) = c F_cl, the possible
      values of g, k = 1, X4 = A0(X(m,n)), Y6 = sqrt(N mu) N(xi), and Y6^2 = X4^3 - 9 X4 H^2 + 81 H^3.  The chain of
      reasoning and the facts it uses are printed; the one non-elementary input is rank E(Q) = 2 (LMFDB, mwrank).
(ii)  THE TRIVIAL CLASS.  F0(m,n) = n (m^3 - m n^2 + 3 n^3) = 3, solved completely and elementarily.
(iii) THE INTEGRAL POINTS, assuming Theorem A' (proved separately by 01_coverage.py and 02_exclusion.py): the coprime
      solutions of the FIELD forms G(m,n) in {1,4,9,36} are +-core.KNOWN.  Through core.to_class they give all coprime
      solutions of the CLASS forms F_cl(m,n) in {1,4,9,36} (Theorem A), and these give exactly 24 points.
(iv)  TRANSFER to n(n+1)/2 = m(m+1)(2m+1)/6 (via (x,y) = (6m+3, 18n+9)) and to the model 3Y^2 = X^3 - X + 3.

Exact arithmetic only (int, fractions.Fraction, sympy); no floating point anywhere.
A failed check prints CHECK FAILED and the script exits with status 1.
Run in this folder:  python 03_reconstruction.py > 03_reconstruction.out.txt
"""
import ast
import os
import re
import sys
import time
from fractions import Fraction as Fr
from math import gcd, isqrt

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import core  # noqa: E402
from core import CLASSES, KNOWN, RHS, Fval, Fclass, to_class, from_class, descent_maps  # noqa: E402

T0 = time.time()
TARGET_DIR = os.path.normpath(os.path.join(HERE, '..', '..', '..', 'targets', 'rank2_signatures_2026-10'))
DESC = os.path.join(TARGET_DIR, 'weierstrass_2descent.out.txt')
LMFDB_MD = os.path.join(TARGET_DIR, 'lmfdb_data.md')
GA, GB = -9, 81                            # g(t) = t^3 + GA t + GB;  E: Y^2 = g(X)
WD = (1, 3, 9)                             # O_K = Z w0 + Z w1 + Z w2 with w_k = th^k / WD[k] (re-checked below)
LW = 9                                     # lcm of WD

m, n, u, v, w, th = sp.symbols('m n u v w th')
NCHECK = [0]


def out(*a):
    print(*a)
    sys.stdout.flush()


def check(cond, msg):
    NCHECK[0] += 1
    if not cond:
        out('CHECK FAILED:', msg)
        sys.exit(1)


def fstr(e):
    return str(sp.expand(e)).replace('**', '^')


def onE(X, Y):
    return Y * Y == X ** 3 + GA * X + GB


def gcd3(a, b, c):
    return gcd(gcd(a, b), c)


def pm_rep(a, b):
    """the representative of +-(a,b) with b > 0, or b = 0 and a > 0"""
    return (a, b) if (b > 0 or (b == 0 and a > 0)) else (-a, -b)


# ----------------------------------------------------------------------------------------------------------
# reading the descent output (text is data: only the characters of polynomials and integer lists are accepted)
# ----------------------------------------------------------------------------------------------------------
_EXPR_CHARS = re.compile(r'^[0-9mnuvwth \+\-\*/\^\(\)]+$')
_LIST_CHARS = re.compile(r'^[\[\]\(\)0-9, \-]*$')


def pexpr(s):
    s = s.strip()
    check(_EXPR_CHARS.match(s) is not None, 'unexpected text in an expression read from a file: %r' % s)
    return sp.expand(sp.sympify(s.replace('^', '**'), locals=dict(m=m, n=n, u=u, v=v, w=w, th=th)))


def pK(s):
    """an element of K written as a polynomial in th of degree <= 2 -> [c0, c1, c2] (Fractions)"""
    P = sp.Poly(pexpr(s), th)
    check(P.degree() <= 2, 'element of K of degree <= 2 in th')
    cs = [sp.Rational(P.coeff_monomial(th ** k)) for k in range(3)]
    return [Fr(int(c.p), int(c.q)) for c in cs]


def plist(s):
    s = s.strip()
    check(_LIST_CHARS.match(s) is not None, 'unexpected text in a list read from a file: %r' % s)
    return ast.literal_eval(s)


def grab(pat, text, what):
    mt = re.search(pat, text, re.M)
    check(mt is not None, 'descent output: line for %s found' % what)
    return mt


# ----------------------------------------------------------------------------------------------------------
# arithmetic in K = Q[th]/(th^3 - 9 th + 81); an element is [c0, c1, c2] = c0 + c1 th + c2 th^2
# (entries may be int, Fraction or sympy expressions)
# ----------------------------------------------------------------------------------------------------------
def kmul(p, q):
    r = [0] * 5
    for i in range(3):
        for j in range(3):
            r[i + j] = r[i + j] + p[i] * q[j]
    # th^3 = -GA th - GB,  th^4 = -GA th^2 - GB th
    return [r[0] - GB * r[3], r[1] - GA * r[3] - GB * r[4], r[2] - GA * r[4]]


def kmat(al):
    cols = [kmul(al, e) for e in ([1, 0, 0], [0, 1, 0], [0, 0, 1])]
    return [[cols[j][i] for j in range(3)] for i in range(3)]


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def knorm(al):
    return det3(kmat(al))


def charpoly(al):
    """(e1, e2, e3) with characteristic polynomial T^3 - e1 T^2 + e2 T - e3"""
    M = kmat(al)
    e1 = M[0][0] + M[1][1] + M[2][2]
    e2 = (M[0][0] * M[1][1] - M[0][1] * M[1][0] + M[0][0] * M[2][2] - M[0][2] * M[2][0]
          + M[1][1] * M[2][2] - M[1][2] * M[2][1])
    return e1, e2, det3(M)


def conj_prod(al):
    """N(al)/al (= product of the other two conjugates): first column of the adjugate of kmat(al)"""
    M = kmat(al)
    return [M[1][1] * M[2][2] - M[1][2] * M[2][1],
            -(M[1][0] * M[2][2] - M[1][2] * M[2][0]),
            M[1][0] * M[2][1] - M[1][1] * M[2][0]]


def ok_coords(al):
    """coordinates in the basis w0 = 1, w1 = th/3, w2 = th^2/9 of O_K"""
    return [WD[k] * Fr(al[k]) for k in range(3)]


def in_OK(al):
    return all(c.denominator == 1 for c in ok_coords(al))


def to_sym(al):
    return [sp.Rational(Fr(c).numerator, Fr(c).denominator) for c in al]


def kstr(al):
    terms = []
    for k, c in enumerate(al):
        c = Fr(c)
        if c == 0:
            continue
        mon = ['', 'th', 'th^2'][k]
        if mon == '':
            terms.append(str(c))
        elif c == 1:
            terms.append(mon)
        elif c == -1:
            terms.append('-' + mon)
        else:
            terms.append('%s*%s' % (c, mon))
    return ' + '.join(terms).replace('+ -', '- ') if terms else '0'


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def nonsquare_witness(al):
    """(p, r): p prime >= 5 (odd, not dividing 3), g(r) = 0 mod p, al(r) a non-residue mod p (al in Z[th])"""
    for p in sp.primerange(5, 3000):
        for r in range(p):
            if (r ** 3 + GA * r + GB) % p == 0:
                val = (int(al[0]) + int(al[1]) * r + int(al[2]) * r * r) % p
                if legendre(val, p) == -1:
                    return p, r
    return None


def ec_add(P, Q):
    """group law on Y^2 = X^3 + GA X + GB (exact rationals; None = O)"""
    if P is None:
        return Q
    if Q is None:
        return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and y1 == -y2:
        return None
    if x1 == x2:
        lam = (3 * x1 * x1 + GA) / (2 * y1)
    else:
        lam = (y2 - y1) / (x2 - x1)
    x3 = lam * lam - x1 - x2
    return (x3, lam * (x1 - x3) - y1)


def int_roots(cf):
    """all integer roots of the monic integer polynomial cf[0] + cf[1] t + ... + t^d (rational-root test):
    if cf[0] != 0 every integer root divides cf[0]; if cf[0] = 0, t = 0 is a root and the others are the roots of
    the quotient by t.  Returns (sorted roots, number of candidates tried)."""
    if len(cf) == 1:
        return [], 0
    if cf[0] == 0:
        rr, nt = int_roots(cf[1:])
        return sorted(set([0] + rr)), nt + 1
    cand = [d for d0 in sp.divisors(abs(cf[0])) for d in (int(d0), -int(d0))]
    roots = sorted(d for d in cand if sum(c * d ** k for k, c in enumerate(cf)) == 0)
    return roots, len(cand)


# ==========================================================================================================
out('03_reconstruction.py -- from the Thue solutions back to the integral points of 8604.a1 (exact arithmetic)')
out('python %s, sympy %s' % (sys.version.split()[0], sp.__version__))
out('2-descent data read from research/proofs/targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt')
out('(section CURVE 8604.a1); every item used is re-checked below (a failing check stops the script, exit status 1).')
out()
with open(DESC, encoding='utf-8') as fh:
    TXT = fh.read()
i0, i1 = TXT.index('CURVE 8604.a1'), TXT.index('CURVE 10404.f2')
SEC = TXT[i0:i1]

# ==========================================================================================================
out('=' * 110)
out('(i) THE REDUCTION: every integral point of E gives coprime (m,n) with H(m,n) = g^2')
out('=' * 110)
out('--- (i.a) the curve, the field K and its maximal order')
mt = grab(r'^CURVE 8604\.a1\s+LMFDB a-invariants \[(.*?)\]', SEC, 'the a-invariants')
check(plist('[' + mt.group(1) + ']') == [0, 0, 0, -9, 81], 'a-invariants [0, 0, 0, -9, 81]')
mt = grab(r'monic model  E: Y\^2 = X\^3 \+ \((-?\d+)\) X \+ \((-?\d+)\),  map X = x, Y = y;', SEC, 'the monic model')
check((int(mt.group(1)), int(mt.group(2))) == (GA, GB), 'monic model Y^2 = X^3 - 9 X + 81')
check((core.GA, core.GB) == (GA, GB), 'core.py has the same curve')
out('E: Y^2 = X^3 - 9 X + 81 is the LMFDB model [a1,a2,a3,a4,a6] = [0,0,0,-9,81] itself (a1 = a2 = a3 = 0): X = x, '
    'Y = y  (file and core.py agree)')
divs = sp.divisors(GB)
check(all(r ** 3 + GA * r + GB != 0 for r in divs + [-d for d in divs]), 'g has no rational root')
out('g(t) = t^3 - 9 t + 81 has no rational root (rational-root test over the divisors of 81): g is irreducible, '
    'K = Q(th) is a cubic field, E(Q)[2] = 0')
mt = grab(r'O_K = Z-span of  w0 = (.+),  w1 = (.+),  w2 = (.+)$', SEC, 'the basis of O_K')
WB = [pK(mt.group(k)) for k in (1, 2, 3)]
check(WB == [[1, 0, 0], [0, Fr(1, WD[1]), 0], [0, 0, Fr(1, WD[2])]], 'O_K basis read from the file')
check(all(in_OK(kmul(WB[i], WB[j])) for i in range(3) for j in range(3)), 'O_K closed under multiplication')
e1, e2, e3 = charpoly(WB[1])
check((e1, e2, e3) == (0, -1, -3), 'lambda = th/3 is a root of T^3 - T + 3')
check(kmul(WB[1], WB[1]) == WB[2], 'w2 = w1^2')
out('file: O_K = Z w0 + Z w1 + Z w2 with w0 = 1, w1 = th/3, w2 = th^2/9.  Re-checked: the structure constants '
    'w_i w_j are integral (an order);')
out('   lambda := w1 = th/3 has characteristic polynomial T^3 - T + 3, and w2 = lambda^2, so O_K = Z[lambda]: the ring '
    'of the cubic lambda^3 - lambda + 3 = 0')
out('   (the cubic of the model 3Y^2 = X^3 - X + 3 of part (iv), whose X is x/3)')
tr = lambda al: charpoly(al)[0]                                      # noqa: E731
dO = det3([[tr(kmul(WB[i], WB[j])) for j in range(3)] for i in range(3)])
dZ = det3([[tr(kmul([int(i == 0), int(i == 1), int(i == 2)], [int(j == 0), int(j == 1), int(j == 2)]))
            for j in range(3)] for i in range(3)])
check(dZ == -4 * GA ** 3 - 27 * GB ** 2 == -174231 == -(3 ** 6) * 239, 'disc Z[th] = -3^6 * 239')
check(dO == -239 and Fr(dZ, dO) == 27 ** 2, 'disc O_K = -239, index 27')
check(int(grab(r'disc\(O_K\) = (-?\d+)', SEC, 'disc O_K').group(1)) == dO, 'disc O_K agrees with the file')
check(int(grab(r'\[O_K : Z\[th\]\] = (\d+)', SEC, 'index').group(1)) == 27, 'index 27 agrees with the file')
sqdiv = [k for k in range(1, isqrt(239) + 1) if 239 % (k * k) == 0]
check(sqdiv == [1] and sp.isprime(239), '239 is prime, hence squarefree')
out('disc(Z[th]) = -174231 = -3^6*239, disc(O_K) = -239 (trace form, exact), [O_K : Z[th]] = 27.  Maximality (N1): '
    'disc(O_K) = [R : O_K]^2 disc(R)')
out('   for the ring of integers R, and -239 is squarefree (239 is prime), so [R : O_K] = 1: O_K is the ring of '
    'integers of K.')
check(int(grab(r'\(D1 = (\d+)\)', SEC, 'D1').group(1)) == 3, 'D1 = 3 in the file')
out('the th-coordinate of a w0 + b w1 + c w2 is b/3: th-coordinates of elements of O_K lie in (1/3)Z  (D1 = 3)')

# ----------------------------------------------------------------------------------------------------------
out()
out('--- (i.b) the four classes of phi(P) = [X(P) - th] in K*/K*^2')
reps = {}
for mt in re.finditer(r'^class \((\d), (\d)\): mu = X0 - th from the integral point \((-?\d+),(-?\d+)\) = '
                      r'(-?\d+)\*G1 \+ (-?\d+)\*G2;  N\(mu\) = (\d+);  mu\(r\) is a non-residue mod p = (\d+) '
                      r'at the root r = (\d+)', SEC, re.M):
    cl = (int(mt.group(1)), int(mt.group(2)))
    X0, Y0 = int(mt.group(3)), int(mt.group(4))
    reps[cl] = dict(P=(X0, Y0), ij=(int(mt.group(5)), int(mt.group(6))), N=int(mt.group(7)),
                    p=int(mt.group(8)), r=int(mt.group(9)), mu=[X0, -1, 0])
check(sorted(reps) == [(0, 1), (1, 0), (1, 1)], 'three non-trivial classes in the file')
mtg = grab(r'LMFDB generators \(original model\) \[\((-?\d+), (-?\d+)\), \((-?\d+), (-?\d+)\)\]', SEC, 'generators')
G1 = (Fr(int(mtg.group(1))), Fr(int(mtg.group(2))))
G2 = (Fr(int(mtg.group(3))), Fr(int(mtg.group(4))))
check(G1 == (3, 9) and G2 == (-3, 9) and onE(*G1) and onE(*G2), 'LMFDB generators G1 = (3,9), G2 = (-3,9) on E')
for cl in sorted(reps):
    R_ = reps[cl]
    check(onE(*R_['P']), 'representative point on E')
    # the decomposition i G1 + j G2 of the representative, by the exact group law
    i_, j_ = R_['ij']
    S_ = None
    for _ in range(abs(i_)):
        S_ = ec_add(S_, G1 if i_ > 0 else (G1[0], -G1[1]))
    for _ in range(abs(j_)):
        S_ = ec_add(S_, G2 if j_ > 0 else (G2[0], -G2[1]))
    check(S_ == (Fr(R_['P'][0]), Fr(R_['P'][1])), 'representative = i G1 + j G2 (exact group law)')
    check((i_ % 2, j_ % 2) == cl, 'class = (i mod 2, j mod 2)')
    check(knorm(R_['mu']) == R_['N'], 'N(mu) agrees with the file')
    p, r = R_['p'], R_['r']
    check((r ** 3 + GA * r + GB) % p == 0 and legendre(R_['mu'][0] - r, p) == -1, 'file witness re-checked')
    out('class %s: mu = %s from the point (%d,%d) on E (= %d*G1 + %d*G2, re-checked by the group law); N(mu) = %d; '
        'file witness (p,r) = (%d,%d) re-checked' % (cl, kstr(R_['mu']), R_['P'][0], R_['P'][1], i_, j_, R_['N'],
                                                      p, r))
MU = {(0, 0): [1, 0, 0], (1, 0): reps[(1, 0)]['mu'], (0, 1): reps[(0, 1)]['mu'], (1, 1): reps[(1, 1)]['mu']}
check(MU[(1, 0)] == [3, -1, 0] and MU[(0, 1)] == [-3, -1, 0] and MU[(1, 1)] == [0, -1, 0], 'mu = 3-th, -3-th, -th')
names = {(1, 0): 'mu10', (0, 1): 'mu01', (1, 1): 'mu11'}
tests = [((1, 0),), ((0, 1),), ((1, 1),), ((1, 0), (0, 1)), ((1, 0), (1, 1)), ((0, 1), (1, 1))]
for tcl in tests:
    al = MU[tcl[0]] if len(tcl) == 1 else kmul(MU[tcl[0]], MU[tcl[1]])
    wit = nonsquare_witness(al)
    check(wit is not None, 'non-square witness for %s' % (tcl,))
    out('   %-12s = %-28s is not a square in K: at the root r = %d of g mod p = %d its value is a non-residue'
        % ('*'.join(names[c] for c in tcl), kstr(al), wit[1], wit[0]))
out('   (if al in O_K were xi^2 with xi in K, xi would be integral (N2), so in O_K, which lies in Z[th][1/3]; the '
    'ring map Z[th][1/3] -> F_p, th -> r (p not dividing 3) would make al(r) a square mod p.)')
out('=> the classes [1], [mu10], [mu01], [mu11] are four DISTINCT elements of phi(E(Q)) (phi(O) = [1], and each mu '
    'is X0 - th of a point of E).')
# the line identity behind the homomorphism property, on a sample
sample = [(Fr(3), Fr(9)), (Fr(-3), Fr(9)), (Fr(0), Fr(9)), (Fr(-5), Fr(1)), (Fr(7), Fr(-19)), (Fr(9), Fr(27)),
          (Fr(1, 4), Fr(71, 8)), (Fr(57, 16), Fr(-621, 64))]
check(all(onE(*P_) for P_ in sample), 'sample points on E')
nline = 0
for P1 in sample:
    for P2 in sample:
        S = ec_add(P1, P2)
        if S is None:
            continue
        lam = (3 * P1[0] ** 2 + GA) / (2 * P1[1]) if P1 == P2 else (P2[1] - P1[1]) / (P2[0] - P1[0])
        nu = P1[1] - lam * P1[0]
        check(onE(*S), 'P1 + P2 on E')
        lhs = kmul(kmul([P1[0], -1, 0], [P2[0], -1, 0]), [S[0], -1, 0])
        check(lhs == kmul([nu, lam, 0], [nu, lam, 0]), 'line identity')
        nline += 1
out('phi is a homomorphism (2-descent): for P1 + P2 + P3 = O on Y = lam X + nu, g(X) - (lam X + nu)^2 = '
    '(X - x1)(X - x2)(X - x3), and X = th gives')
out('   (x1 - th)(x2 - th)(x3 - th) = (lam th + nu)^2  [identity re-checked exactly on %d pairs of points]' % nline)
out('COMPLETENESS OF THE CLASS LIST (the non-elementary input): |phi(E(Q))| <= |E(Q)/2E(Q)| = 2^rank * |E(Q)[2]| '
    '= 2^rank (E(Q)[2] = 0 above;')
out('   E(Q) finitely generated by Mordell-Weil).  With rank E(Q) = 2 (LMFDB 8604.a1, computed by mwrank; only '
    'rank <= 2 is used) this is 4,')
out('   so phi(E(Q)) = {[1], [3 - th], [-3 - th], [-th]}: every point P != O has X(P) - th = mu xi^2 with xi in K* for '
    'one of these four mu.')

# ----------------------------------------------------------------------------------------------------------
cls_pos = [(mt.start(), (int(mt.group(1)), int(mt.group(2))))
           for mt in re.finditer(r'^CLASS \((\d), (\d)\)', SEC, re.M)]
iextra = SEC.index('EXTRA for 8604.a1')
BLOCKS = {}
for k_, (pos, cl) in enumerate(cls_pos):
    BLOCKS[cl] = SEC[pos:(cls_pos[k_ + 1][0] if k_ + 1 < len(cls_pos) else iextra)]
check(sorted(BLOCKS) == [(0, 0), (0, 1), (1, 0), (1, 1)], 'four class blocks in the file')

KEYS = [('00', (0, 0)), ('10', (1, 0)), ('01', (0, 1)), ('11', (1, 1))]
RED = {}
for key, cl in KEYS:
    Bk = BLOCKS[cl]
    out()
    out('--- (i.c) class %s  [key %r]' % (cl, key))
    mu = pK(grab(r'^CLASS .*\)\s+mu = (.+)$', Bk, 'mu').group(1))
    check(mu == MU[cl], 'mu of the block = mu of the class list')
    if key != '00':
        check(pK(CLASSES[key]['mu'].replace('theta', 'th')) == mu, 'core.CLASSES mu = file mu')
    check(in_OK(mu) and all(Fr(c).denominator == 1 for c in mu), 'mu in Z[th]')
    Nmu = Fr(knorm(mu))
    check(Nmu.denominator == 1, 'N(mu) integral')
    Nmu = int(Nmu)
    check(Nmu == int(grab(r'^N\(mu\) = (\d+);', Bk, 'N(mu)').group(1)), 'N(mu) agrees with the file')
    s = isqrt(Nmu)
    check(s * s == Nmu, 'N(mu) is a square')
    check(s == int(grab(r'O_K <= Xi_mu <= \(1/(\d+)\) O_K', Bk, 'the bound for Xi_mu').group(1)), 'bound s')
    cp = conj_prod(mu)
    check(kmul(mu, cp) == [Nmu, 0, 0] and in_OK(cp), 'N(mu)/mu in O_K')
    out('mu = %s, N(mu) = %d = %d^2.  N(mu)/mu = %s lies in O_K, so mu xi^2 in O_K gives (%d xi)^2 = (mu xi^2)'
        '(N(mu)/mu) in O_K,' % (kstr(mu), Nmu, s, kstr(cp), s))
    out('   hence %d xi is integral (N2): Xi_mu := {xi in K : mu xi^2 in O_K} lies in (1/%d) O_K' % (s, s))
    # enumeration of (1/s) O_K / O_K:  xi = (c0 w0 + c1 w1 + c2 w2)/s = y/(LW s) with y_k = (LW/WD[k]) c_k
    S = []
    Dk = [LW * LW * s * s // WD[k] for k in range(3)]      # mu xi^2 = z/(LW s)^2 is in O_K iff Dk | z_k
    mui = [int(c) for c in mu]
    for c0 in range(s):
        for c1 in range(s):
            for c2 in range(s):
                y = [LW // WD[0] * c0, LW // WD[1] * c1, LW // WD[2] * c2]
                z = kmul(mui, kmul(y, y))
                if all(z[k] % Dk[k] == 0 for k in range(3)):
                    S.append((c0, c1, c2))
    mt = grab(r'\[Xi_mu : O_K\] = (\d+);\s+Z-basis b1 = (.+?), b2 = (.+?), b3 = (.+)$', Bk, 'the basis of Xi_mu')
    idx = int(mt.group(1))
    bas = [pK(mt.group(k)) for k in (2, 3, 4)]
    check(len(S) == idx, '[Xi_mu : O_K] = number of classes found')
    Bw = [ok_coords(bj) for bj in bas]                                 # rows: b_j in w-coordinates
    check(Fr(1) / abs(det3(Bw)) == idx, '[L : O_K] = idx for L = Z b1 + Z b2 + Z b3')
    BwM = sp.Matrix([[sp.Rational(c.numerator, c.denominator) for c in row] for row in Bw])
    BwI = BwM.inv()

    def in_L(wc):
        co = sp.Matrix([[sp.Rational(Fr(c).numerator, Fr(c).denominator) for c in wc]]) * BwI
        return all(x.is_integer for x in co)
    check(all(in_L([int(i == j) for j in range(3)]) for i in range(3)), 'O_K inside L')
    check(all((s * c).denominator == 1 for row in Bw for c in row), 'L inside (1/s) O_K')
    check(all(in_L([Fr(c, s) for c in sc]) for sc in S), 'every class found lies in L')
    check(all(in_OK(kmul(mu, kmul(bj, bj))) for bj in bas), 'mu b_j^2 in O_K')
    check(bas[0] == [1, 0, 0], 'b1 = 1')
    out('enumeration of the %d classes of (1/%d)O_K/O_K: exactly %d satisfy mu xi^2 in O_K (a property of the class '
        'mod O_K, by N2)' % (s ** 3, s, len(S)))
    out('   file basis b1 = %s, b2 = %s, b3 = %s: L = Z b1 + Z b2 + Z b3 contains O_K with [L : O_K] = %d, lies in '
        '(1/%d) O_K' % (kstr(bas[0]), kstr(bas[1]), kstr(bas[2]), idx, s))
    out('   and contains every class found; so L/O_K = the %d classes and Xi_mu = L  [exact]' % idx)
    # the forms A0, A1, A2
    bs = [to_sym(bj) for bj in bas]
    xi_s = [u * bs[0][k] + v * bs[1][k] + w * bs[2][k] for k in range(3)]
    A0, A1, A2 = [sp.expand(c) for c in kmul(to_sym(mu), kmul(xi_s, xi_s))]
    for nm, Ak, pat in (('A2', A2, r'A2 \(conic\)\s+= (.+)$'), ('A1', A1, r'A1 \(th-coeff\)\s+= (.+)$'),
                        ('A0', A0, r'A0 \(1-coeff\)\s+= (.+)$')):
        check(sp.expand(Ak - pexpr(grab(pat, Bk, nm).group(1))) == 0, '%s agrees with the file' % nm)
    check(all(sp.Rational(c).is_integer for c in sp.Poly(A0, u, v, w).coeffs()), 'A0 has integer coefficients')
    out('mu (u b1 + v b2 + w b3)^2 = A0 + A1 th + A2 th^2 (recomputed, = file):')
    out('   A2 = %s;   A1 = %s;   A0 = %s' % (fstr(A2), fstr(A1), fstr(A0)))
    out('an integral point (X,Y) with X - th = mu xi^2 has xi = u b1 + v b2 + w b3, (u,v,w) in Z^3, with '
        'A2 = 0, A1 = -1, A0 = X')
    # conic
    mt = grab(r'^conic C_mu: (.+) = 0 \(A2 divided by its content (\S+)\); Gram determinant (-?\d+) != 0', Bk,
              'the conic')
    Q = pexpr(mt.group(1))
    cont = sp.Rational(mt.group(2))
    check(sp.expand(A2 - cont * Q) == 0, 'A2 = content * Q')
    Qco = sp.Poly(Q, u, v, w).coeffs()
    check(all(sp.Rational(c).is_integer for c in Qco) and sp.gcd_list([int(c) for c in Qco]) == 1, 'Q primitive')
    VV = [u, v, w]
    Gm = sp.Matrix(3, 3, lambda i, j: sp.Poly(Q, u, v, w).coeff_monomial(VV[i] ** 2) if i == j else
                   sp.Rational(sp.Poly(Q, u, v, w).coeff_monomial(VV[i] * VV[j]), 2))
    check(Gm.det() == int(mt.group(3)) and Gm.det() != 0, 'Gram determinant')
    P0 = plist('[' + grab(r'conic point P0 = \[(.*?)\] from', Bk, 'P0').group(1) + ']')
    check(P0 == [1, 0, 0] and Q.subs({u: 1, v: 0, w: 0}) == 0, 'P0 = (1,0,0) (xi = 1) on the conic')
    mt = grab(r'basis completion \[P0 e f\] = (\[\[.*?\]\]) \(det (-?\d+)\)', Bk, 'basis completion')
    check(plist(mt.group(1)) == [[1, 0, 0], [0, 1, 0], [0, 0, 1]] and int(mt.group(2)) == 1, '[P0 e f] = I')
    out('conic C: Q = %s = 0 (A2 = %s Q), Gram determinant %s != 0 (smooth); P0 = (1,0,0) from xi = 1'
        % (fstr(Q), cont, Gm.det()))
    # parametrisation
    Xu = pexpr(grab(r'^\s+u = (.+)$', Bk, 'u').group(1))
    Xv = pexpr(grab(r'^\s+v = (.+)$', Bk, 'v').group(1))
    Xw = pexpr(grab(r'^\s+w = (.+)$', Bk, 'w').group(1))
    Xp = [Xu, Xv, Xw]
    if key != '00':
        uvw_f, X4_f, Y6_f = descent_maps(key)
        check([sp.expand(e) for e in uvw_f(m, n)] == Xp, 'core.descent_maps (u,v,w) = file')
    subX = {u: Xu, v: Xv, w: Xw}
    check(sp.expand(Q.subs(subX, simultaneous=True)) == 0, 'X(m,n) lies on the conic')
    Tl = plist('[' + grab(r'GL2\(Z\) change \(m,n\) -> \[(.*?)\] \(cosmetic', Bk, 'T').group(1) + ']')
    Tm = sp.Matrix([[Tl[0], Tl[1]], [Tl[2], Tl[3]]])
    check(abs(Tm.det()) == 1, 'T in GL2(Z)')
    mt = grab(r'with l = (.+), q = (.+)$', Bk, 'l and q')
    lpost, qpost = pexpr(mt.group(1)), pexpr(mt.group(2))
    mp_, np_ = Tm[0, 0] * m + Tm[0, 1] * n, Tm[1, 0] * m + Tm[1, 1] * n
    check(Xp == [qpost, sp.expand(-lpost * mp_), sp.expand(-lpost * np_)], 'X(m,n) = (q, -l m\', -l n\')')
    BP = sp.expand(Q.subs({u: u + 1}, simultaneous=True) - Q)              # B(P0, P), linear
    Ti = Tm.inv()
    sub_back = {m: Ti[0, 0] * v + Ti[0, 1] * w, n: Ti[1, 0] * v + Ti[1, 1] * w}
    lhs = [sp.expand(e.subs(sub_back, simultaneous=True)) for e in Xp]
    rhs = [sp.expand(Q * [1, 0, 0][i] - BP * VV[i]) for i in range(3)]
    check(lhs == rhs, 'inverse identity X(T^-1 (v,w)) = Q(P) P0 - B(P0,P) P')
    # tangent case: B(P0,P) = 0 and Q(P) = 0 force P ~ P0
    lv, lw = BP.coeff(v), BP.coeff(w)
    check(BP.coeff(u) == 0 and (lv, lw) != (0, 0), 'B(P0,.) vanishes at P0, is non-zero')
    Tvec = [0, lw, -lv]
    QT = Q.subs({u: Tvec[0], v: Tvec[1], w: Tvec[2]}, simultaneous=True)
    check(QT != 0, 'Q(T0) != 0 on the tangent plane')
    al_, be_ = sp.Poly(lpost, m, n).coeff_monomial(m), sp.Poly(lpost, m, n).coeff_monomial(n)
    m0, n0 = (int(be_), -int(al_))
    g0 = gcd(m0, n0)
    m0, n0 = m0 // g0, n0 // g0
    Xt = [e.subs({m: m0, n: n0}) for e in Xp]
    check(Xt[1] == 0 and Xt[2] == 0 and Xt[0] != 0, 'X(m0,n0) ~ P0')
    out('parametrisation (file%s): u = %s, v = %s, w = %s' % ('' if key == '00' else ' = core.descent_maps',
                                                               fstr(Xu), fstr(Xv), fstr(Xw)))
    out('   = (q, -l m\', -l n\') with l = %s, q = %s, (m\',n\') = T(m,n), T = %s (det %s)'
        % (fstr(lpost), fstr(qpost), Tl, Tm.det()))
    out('   identity in Z[u,v,w]:  X(T^-1(v,w)) = Q(u,v,w) P0 - B(P0,P) P   with B(P0,P) = %s  [exact];' % fstr(BP))
    out('   so a point P of C with B(P0,P) != 0 is proportional to X(T^-1(v,w)); if B(P0,P) = 0 then P ~ P0 (on the '
        'tangent plane Q(a P0 + b T0) = b^2 Q(T0),')
    out('   T0 = %s, Q(T0) = %s != 0), and P0 ~ X(%d,%d).  Every rational point of C is ~ X(m,n) with (m,n) coprime.'
        % (Tvec, QT, m0, n0))
    # g, k
    Pl = sp.Poly(lpost, m, n)
    check(Pl.total_degree() == 1 and len(Pl.terms()) == 1, 'l = c_l * (m or n)')
    (ex, c_l) = Pl.terms()[0]
    xv, yv = (m, n) if ex == (1, 0) else (n, m)
    qx0 = sp.expand(qpost.subs(xv, 0))
    check(qx0 == yv ** 2, 'q = %s^2 mod %s' % (yv, xv))
    gvals = sorted(sp.divisors(abs(int(c_l))))
    gfile = plist('[' + grab(r'possible values of g \(exact, per prime power of R, by CRT\): \[(.*)\]', Bk,
                             'g values').group(1) + ']')
    check(gvals == gfile, 'possible g agree with the file')
    Rfile = int(grab(r'divides R = \|q\(beta,-alpha\)\| = (\d+)', Bk, 'R').group(1))
    Rval = abs(int(qpost.subs({m: be_, n: -al_}, simultaneous=True)))
    check(Rval == Rfile and Rval % abs(int(c_l)) == 0,
          'R = |q(beta,-alpha)| for l = alpha m + beta n agrees with the file')
    mt = grab(r'so k\^2 \| D1 = (\d+): k in \[(.*)\]', Bk, 'k')
    check(int(mt.group(1)) == 3 and plist('[' + mt.group(2) + ']') == [1], 'k = 1 in the file')
    check([k for k in range(1, 3) if 3 % (k * k) == 0] == [1], 'k^2 | 3 => k = 1')
    out('g := gcd(u,v,w) for coprime (m,n): T is unimodular, so gcd(m\',n\') = 1 and gcd(v,w) = |l| = %d|%s|; '
        'q = %s^2 mod %s gives gcd(q,%s) = 1;'
        % (abs(int(c_l)), xv, yv, xv, xv))
    out('   so g = gcd(q, %d) divides %d: g in %s (file, by CRT from R = |q(beta,-alpha)| = %d for '
        'l = alpha m + beta n: %s)' % (abs(int(c_l)), abs(int(c_l)), gvals, Rval, gfile))
    out('k = 1: if (u,v,w) = k (u\',v\',w\'), then (X - th)/k^2 = mu (xi/k)^2 is in O_K, its th-coordinate -1/k^2 '
        'is in (1/3)Z, so k^2 | 3 and k = 1')
    out('   => (u,v,w) is primitive and (u,v,w) = +-X(m,n)/g for coprime (m,n)')
    # H = c F
    cH = sp.Rational(grab(r'with c = (\S+) and primitive', Bk, 'c').group(1))
    Fp = pexpr(grab(r'^\s+F\(m,n\) = (.+)$', Bk, 'F').group(1))
    if key != '00':
        check(sp.expand(Fclass(key, m, n) - Fp) == 0, 'core class form F_cl = file F')
    H = sp.expand(-A1.subs(subX, simultaneous=True))
    check(sp.expand(H - cH * Fp) == 0, 'H = -A1(X(m,n)) = c F')
    rhs_vals = sorted(set(sp.Rational(g * g) / cH for g in gvals))
    rhs_int = [int(r) for r in rhs_vals if r.is_integer]
    rfile = plist('[' + grab(r'so F\(m,n\) in \[(.*)\]', Bk, 'F values').group(1) + ']')
    check(rhs_int == rfile, 'right-hand sides agree with the file')
    out('A1(u,v,w) = -1 and (u,v,w) = +-X(m,n)/g give H(m,n) := -A1(X(m,n)) = g^2, with H = c F, c = %s,' % cH)
    out('   F(m,n) = %s;  so F(m,n) = g^2/c with g in %s: the integral values are F(m,n) in %s' % (fstr(Fp), gvals,
                                                                                               rhs_int))
    GU = None
    if key != '00':
        check(cH == 1 and rhs_int == list(RHS), 'c = 1 and F_cl in {1,4,9,36}')
        check(Xu == qpost and abs(int(c_l)) == 6, 'u = q and |c_l| = 6, so g = gcd(u, 6)')
        Dq = sp.Poly(sp.expand(Xu ** 2 - Fp), m, n)
        cont_D = int(sp.gcd_list([int(c) for c in Dq.coeffs()]))
        check(cont_D % 12 == 0, 'u^2 = F mod 12 identically')
        for c in RHS:
            check(all(gcd(r, 6) == isqrt(c) for r in range(6) if (r * r - c) % 6 == 0), 'u^2 = c mod 6 => g = sqrt c')
            check(any((r * r - c) % 6 == 0 for r in range(6)), 'c is a square mod 6')
        GU = cont_D
        out('   u = q here, so g = gcd(u, 6).  u^2 - F = %d * (%s) identically (content %d): u^2 = F (mod %d), in '
            'particular mod 6;' % (cont_D, fstr(sp.expand((Xu ** 2 - Fp) / cont_D)), cont_D, cont_D))
        out('   residues mod 6: r^2 = c (mod 6) forces gcd(r,6) = sqrt(c) for each c in {1,4,9,36} (r^2 = 1, 4, 3, 0 '
            'mod 6 <=> gcd(r,6) = 1, 2, 3, 6);')
        out('   so EVERY coprime solution of F(m,n) = c in {1,4,9,36} has g = gcd(u,6) = sqrt(c), i.e. H = g^2 holds '
            'automatically')
        # the field form G actually solved by 01/02
        if key == '10':
            check(sp.expand(Fclass(key, m, n) - Fval(key, 2 * m, n + m)) == 0, 'F_cl(m,n) = G(2m, n+m)')
            check(sp.expand(Fclass(key, m, 2 * n - m) - 16 * Fval(key, m, n)) == 0, 'F_cl(m, 2n-m) = 16 G(m,n)')
            out('   field form (the Thue equation of 01/02): G(m,n) = %s, beta = (1 + alpha)/2;'
                % fstr(Fval(key, m, n)))
            out('      identities  F_cl(m,n) = G(2m, n+m)  and  F_cl(m, 2n-m) = 16 G(m,n)  in Z[m,n]  [exact]')
        else:
            check(tuple(CLASSES[key]['F']) == tuple(CLASSES[key]['CLASS_F']), 'G = F_cl')
            check(sp.expand(Fclass(key, m, n) - Fval(key, m, n)) == 0, 'G = F_cl (polynomials)')
            out('   field form (the Thue equation of 01/02): G = F_cl, beta = alpha')
    else:
        check(cH == Fr(4, 3) and gvals == [1, 2] and rhs_int == [3], 'trivial class: c = 4/3, g | 2, F0 = 3')
        check(all(((mm ** 2 + nn ** 2) % 2 == 0) == ((mm - nn) % 2 == 0)
                  for mm in range(2) for nn in range(2)), 'u even iff m = n mod 2')
        out('   here H = (4/3) F0 = g^2 with g | 2: g = 1 would give F0 = 3/4, so F0(m,n) = 3 and g = 2, i.e. '
            'u = m^2 + n^2 even, i.e. m = n (mod 2)')
    # X4, Y6
    X4 = sp.expand(A0.subs(subX, simultaneous=True))
    check(sp.expand(X4 - pexpr(grab(r'X4 = (.+)$', Bk, 'X4').group(1))) == 0, 'X4 = file')
    xiP = [sp.expand(e.subs(subX, simultaneous=True)) for e in xi_s]
    Y6 = sp.expand(s * det3(kmat(xiP)))
    check(sp.expand(Y6 - pexpr(grab(r'Y6 = (.+)$', Bk, 'Y6').group(1))) == 0, 'Y6 = file')
    if key != '00':
        check(sp.expand(X4 - X4_f(m, n)) == 0 and sp.expand(Y6 - Y6_f(m, n)) == 0, 'core X4, Y6 = recomputed')
    check(all(sp.Rational(c).is_integer for c in sp.Poly(X4, m, n).coeffs() + sp.Poly(Y6, m, n).coeffs()),
          'X4, Y6 have integer coefficients')
    check(sp.expand(Y6 ** 2 - (X4 ** 3 + GA * X4 * H ** 2 + GB * H ** 3)) == 0, 'Y6^2 = X4^3 - 9 X4 H^2 + 81 H^3')
    if key != '00':
        check(sp.expand(Y6 ** 2 - (X4 ** 3 + GA * X4 * Fp ** 2 + GB * Fp ** 3)) == 0, 'identity with F')
    out('X = A0(u,v,w) = X4(m,n)/g^2 with X4 = A0(X(m,n)) = %s' % fstr(X4))
    out('Y^2 = N(X - th) = N(mu) N(xi)^2, Y = +-%d N(xi) = +-Y6(m,n)/g^3 with Y6 = %d N(X(m,n).b) = %s'
        % (s, s, fstr(Y6)))
    out('IDENTITY  Y6^2 = X4^3 - 9 X4 %s^2 + 81 %s^3  holds in %s[m,n]  [exact]%s'
        % ('F' if key != '00' else 'H', 'F' if key != '00' else 'H', 'Z' if key != '00' else 'Q',
           '' if key != '00' else '   (H = (4/3) F0; with F0 itself in place of H it does not hold)'))
    if key == '00':
        check(sp.expand(Y6 ** 2 - (X4 ** 3 + GA * X4 * Fp ** 2 + GB * Fp ** 3)) != 0, 'not with F0')
    out('   so whenever H(m,n) = g^2:  (X,Y) = (X4/g^2, +-Y6/g^3) satisfies Y^2 = X^3 - 9X + 81; X is an '
        'integer (A0 has integer coefficients, X(m,n)/g is integral)')
    out('   and Y is then an integer (Y^2 is).  Conversely X - th = mu xi^2 with xi = (X(m,n)/g).b, so the point '
        'lies in the class [mu].')
    # the file's table of known points of the class
    ftab = [(int(a), int(b), (int(c), int(d)), int(e), int(f)) for a, b, c, d, e, f in
            re.findall(r'^\s+(-?\d+)\s+(-?\d+)\s+\((-?\d+), (-?\d+)\)\s+(\d+)\s+(\d+)\s*$', Bk, re.M)]
    check(len(ftab) > 0, 'file table of the class')
    RED[key] = dict(cl=cl, mu=mu, s=s, Nmu=Nmu, idx=idx, bas=bas, Q=Q, cont=cont, Xp=Xp, l=lpost, q=qpost, T=Tl,
                    X4=X4, Y6=Y6, H=H, F=Fp, cH=cH, gvals=gvals, ftab=ftab, GU=GU)

out()
out('--- (i.d) SUMMARY OF THE REDUCTION (chain of reasoning, facts used)')
for line in [
    '1. [class list; NON-ELEMENTARY INPUT] phi: E(Q) -> K*/K*^2, P -> [X(P) - th], is a homomorphism (line identity) '
    'with 2E(Q) in its kernel.',
    '   E(Q)[2] = 0 (g irreducible).  rank E(Q) = 2 (LMFDB, mwrank) and Mordell-Weil give |E(Q)/2E(Q)| = 4.  The four '
    'classes [1], [3-th],',
    '   [-3-th], [-th] are distinct (Legendre witnesses) and lie in phi(E(Q)); hence they are all of it.  Only '
    'rank <= 2 is used, through the class list (H) it implies.',
    '2. [lattice; N1, N2] O_K = Z[th/3] is the ring of integers (index-discriminant formula: disc O_K = -239 is '
    'squarefree).  For an integral point',
    '   in class [mu], X - th lies in O_K, so xi in Xi_mu = Z b1 + Z b2 + Z b3 (enumeration of (1/s)O_K/O_K, which '
    'uses integral closure), xi = u b1 + v b2 + w b3.',
    '3. [conic; elementary] A2(u,v,w) = 0, A1 = -1, A0 = X.  The conic A2 = 0 is smooth with the point P0 = (1,0,0); '
    'by the explicit inverse identity',
    '   every rational point is ~ X(m,n) for coprime integers (m,n) (unique up to sign).',
    '4. [primitivity; elementary] k^2 | 3 (th-coordinates of O_K in (1/3)Z) forces (u,v,w) primitive: (u,v,w) = '
    '+-X(m,n)/g, g = gcd of X(m,n).',
    '5. [equation; elementary] A1 = -1 gives H(m,n) = g^2.  Non-trivial classes: H = F_cl and g | 6, so F_cl(m,n) = '
    'g^2 in {1,4,9,36};',
    '   trivial class: H = (4/3) F0 and g | 2, so F0(m,n) = 3 with g = 2.  Then X = X4/g^2 and Y = +-Y6/g^3.',
    '6. [converse; elementary] each coprime (m,n) with H = g^2 gives an integral point of the class; in the three '
    'non-trivial classes u^2 = F_cl (mod 12),',
    '   so EVERY coprime solution of F_cl = c in {1,4,9,36} has g = sqrt(c) and gives an integral point.  The map '
    '+-(m,n) -> {(X,+-Y)} is a bijection',
    '   (xi is fixed up to sign by X, and P^1(Q) -> C(Q) is bijective by the inverse identity).',
    '7. [field forms; elementary] the Thue equations solved in 01/02 are those of the field forms G: G = F_cl for '
    '(0,1) and (1,1);',
    '   for (1,0), F_cl(m,n) = G(2m, n+m) (part iii.a turns Theorem A\' for G into Theorem A for F_cl).',
    'FACTS USED:',
    '   external, non-elementary: rank E(Q) = 2 for 8604.a1 (LMFDB; only rank <= 2 is used), with the Mordell-Weil '
    'theorem.  LMFDB\'s reliability page',
    '      (see lmfdb_data.md): the analytic rank 2 is rigorous, but for analytic rank >= 2 the Mordell-Weil rank '
    'comes from mwrank and is not called proven.',
    '   standard algebraic number theory (N1, N2 of the descent script): the index-discriminant formula; integral '
    'elements form a ring, alpha is',
    '      integral iff its characteristic polynomial lies in Z[T], alpha^2 integral => alpha integral.  Used in '
    'steps 1 (Legendre witnesses) and 2.',
    '   classical: the 2-descent map phi is a homomorphism (line identity).',
    '   elementary: rational-root test, Legendre symbols, lines through a rational point of a smooth conic, gcd '
    'arguments, polynomial identities.',
]:
    out(line)

# ==========================================================================================================
out()
out('=' * 110)
out('(ii) THE TRIVIAL CLASS (0,0): F0(m,n) = n (m^3 - m n^2 + 3 n^3) = 3, solved completely')
out('=' * 110)
R0 = RED['00']
Cub = m ** 3 - m * n ** 2 + 3 * n ** 3
check(sp.expand(R0['F'] - n * Cub) == 0, 'F0 = n * cubic (file form)')
out('F0 = %s = n * (%s)  (= the file\'s F of class (0, 0))' % (fstr(R0['F']), fstr(Cub)))
out('n divides F0(m,n) = 3, so n in {+-1, +-3}; for each n, m is an integer root of the monic cubic')
out('   m^3 + 0 m^2 - n^2 m + (3n^3 - 3/n) = 0.  Rational-root test: if the constant term is non-zero, m divides it; '
    'if it is 0,')
out('   then m = 0 or m is a root of m^2 - n^2 (constant term -n^2 != 0, so m divides n^2):')
SOL0 = []
for nn in [1, -1, 3, -3]:
    check(3 % nn == 0, 'n | 3')
    c0 = 3 * nn ** 3 - 3 // nn
    cf = [c0, -nn ** 2, 0, 1]
    roots, ntried = int_roots(cf)
    for r in roots:
        check(int((n * Cub).subs({m: r, n: nn})) == 3, 'F0 = 3')
        SOL0.append((r, nn))
    # every root found by the test is a root, and every root is found (the test lists all candidates)
    check(all(sum(c * r ** k for k, c in enumerate(cf)) == 0 for r in roots), 'roots')
    out('   n = %2d: m^3 + (%d) m + (%d) = 0;  %s;  %d candidates tried;  integer roots m = %s'
        % (nn, -nn ** 2, c0, 'constant term 0: m (m^2 - %d)' % (nn ** 2) if c0 == 0 else 'constant term != 0',
           ntried, roots))
SOL0.sort()
file0 = plist(grab(r'coprime solutions of F\(m,n\) in \[3\]: (\[.*\])$', SEC, 'trivial class solutions').group(1))
check(SOL0 == sorted(file0), 'solutions agree with the file')
check(all(gcd(a, b) == 1 for a, b in SOL0), 'all solutions coprime')
check(SOL0 == [(-5, 3), (-1, -1), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 1), (5, -3)], 'the 8 solutions')
out('all integer solutions (automatically coprime: d^4 | 3 forces d = 1): %s  (= file)' % SOL0)
out('map with u = m^2 + n^2, v = -2mn, w = -2n^2 (file), g = gcd(u,v,w); an integral point needs H = (4/3)*3 = 4 '
    '= g^2, i.e. g = 2:')
check([sp.expand(e) for e in R0['Xp']] == [sp.expand(m ** 2 + n ** 2), sp.expand(-2 * m * n), sp.expand(-2 * n ** 2)],
      'trivial-class parametrisation')
out('   %-10s %-16s %-3s %-22s %s' % ('(m,n)', '(u,v,w)', 'g', 'X = X4/4', 'Y = +-Y6/8'))
PTS = {}          # X -> dict(Y=|Y|, src=[...])
NOT_INTEGRAL = []
for (a, b) in SOL0:
    uu, vv, ww = [int(e.subs({m: a, n: b})) for e in R0['Xp']]
    g = gcd3(uu, vv, ww)
    X4v, Y6v = int(R0['X4'].subs({m: a, n: b})), int(R0['Y6'].subs({m: a, n: b}))
    Hv = int(R0['H'].subs({m: a, n: b}))
    check(Hv == 4, 'H = 4')
    check(g == (2 if (a - b) % 2 == 0 else 1), 'g = 2 iff m = n (mod 2)')
    Xq, Yq = Fr(X4v, 4), Fr(Y6v, 8)
    check(onE(Xq, Yq), 'point on E')
    if g == 2:
        check(Xq.denominator == 1 and Yq.denominator == 1, 'integral point')
        X, Y = int(Xq), abs(int(Yq))
        xi = [Fr(uu, g) * R0['bas'][0][k] + Fr(vv, g) * R0['bas'][1][k] + Fr(ww, g) * R0['bas'][2][k]
              for k in range(3)]
        check(kmul(xi, xi) == [X, -1, 0], 'X - th = xi^2 (trivial class)')
        check(knorm(xi) in (Y, -Y), 'Y = +-N(xi)')
        PTS.setdefault(X, dict(Y=Y, src=[]))
        check(PTS[X]['Y'] == Y, 'same |Y|')
        PTS[X]['src'].append(('00', 3, (a, b)))
        hits = [r for r in R0['ftab'] if r[0] == X]
        check(len(hits) == 2 and all(r[2] in ((a, b), (-a, -b)) and r[3] == g and r[4] == 3 for r in hits),
              'file table of class (0, 0) agrees for X = %d' % X)
        note = 'integral point (%d, +-%d); X - th = xi^2 checked; = file table' % (X, Y)
    else:
        check(Xq.denominator != 1, 'not integral')
        NOT_INTEGRAL.append(((a, b), g, Xq, Yq))
        note = 'g = %d != 2: H = 4 != g^2, NOT from an integral point; rational point (%s, +-%s) [non-integral]' % (
            g, Xq, abs(Yq))
    out('   %-10s %-16s %-3d %-22s %s' % ((a, b), (uu, vv, ww), g, Xq, note))
check(sorted(PTS) == [-5, 7, 1099], 'trivial class gives X in {-5, 7, 1099}')
check([s_[0] for s_ in NOT_INTEGRAL] == [(0, -1), (0, 1)], 'only +-(0,1) fails')
out('=> trivial class: X in {-5, 7, 1099} (from +-(1,1), +-(1,-1), +-(5,-3)); the solution +-(0,1) of F0 = 3 has g = 1 '
    'and gives no integral point.')
mt = grab(r'class \(0, 0\) \(reducible F = .*\n.*\n\s+integral points with phi\(P\) in this class: \(X,\|Y\|\) in '
          r'(\[.*?\])', SEC, 'trivial class points')
check(sorted((X, PTS[X]['Y']) for X in PTS) == plist(mt.group(1)), 'trivial class points agree with the file')

# ==========================================================================================================
out()
out('=' * 110)
out('(iii) ALL INTEGRAL POINTS, assuming Theorem A\': the coprime solutions of the field forms G(m,n) in {1,4,9,36} '
    'are +-core.KNOWN')
out('=' * 110)
out('--- (iii.a) Theorem A\' (field forms G) => Theorem A (class forms F_cl)')
check(all(c % 16 != 0 for c in RHS), '16 divides no c in {1,4,9,36}')
out('classes (0,1), (1,1): G = F_cl (checked in i.c), so the two theorems coincide.')
out('class (1,0): F_cl(m,n) = G(2m, n+m) (i.c).  (a) A coprime solution (m,n) of F_cl = c gives the solution (2m, n+m) '
    'of G = c;')
out('   gcd(2m, n+m) divides 2 (gcd(m, n+m) = 1), and if it were 2, then G(2m, n+m) = 16 G(m, (n+m)/2) would be '
    'divisible by 16,')
out('   which no c in {1,4,9,36} is.  So (2m, n+m) = from_class(m,n) is a coprime solution of G = c, hence +- an '
    'element of KNOWN,')
out('   and (m,n) = to_class(2m, n+m).  (b) Conversely, a solution (m,n) of G = c with m even gives the coprime '
    'solution to_class(m,n) =')
out('   (m/2, n - m/2) of F_cl = c (gcd(m/2, n - m/2) = gcd(m/2, n) divides gcd(m,n) = 1); one with m odd gives NO '
    'solution of F_cl = c:')
out('   F_cl(m, 2n-m) = 16 G(m,n) = 16c with gcd(m, 2n-m) = gcd(m, 2n) = 1, so the point (m : 2n-m) of P^1(Q) has the '
    'value 16c at its')
out('   primitive representative, and 16c is not in {1,4,9,36} (not g^2 with g | 6: no integral point).')
CLASS_SOL = {}
SKIPPED = []
for key in ('10', '01', '11'):
    CLASS_SOL[key] = []
    for (a, b) in KNOWN[key]:
        check(gcd(a, b) == 1, 'KNOWN coprime')
        cG = Fval(key, a, b)
        check(cG in CLASSES[key]['RHS'] and cG in RHS, 'G(m,n) in {1,4,9,36}')
        mc = to_class(key, a, b)
        if mc is None:
            check(key == '10' and a % 2 == 1, 'to_class None only for class (1,0) with m odd')
            m2, n2 = a, 2 * b - a
            check(gcd(m2, n2) == 1 and Fclass(key, m2, n2) == 16 * cG, 'F_cl(m, 2n-m) = 16 c, coprime')
            uu, vv, ww = descent_maps(key)[0](m2, n2)
            g = gcd3(uu, vv, ww)
            h = 4 * isqrt(cG)
            check(h * h == 16 * cG and g != h, 'H = 16c = h^2 with h != g')
            Xq, Yq = Fr(descent_maps(key)[1](m2, n2), h * h), Fr(descent_maps(key)[2](m2, n2), h ** 3)
            check(onE(Xq, Yq) and Xq.denominator != 1, 'a rational, non-integral point of E')
            SKIPPED.append((key, (a, b), cG, (m2, n2), g, Xq, Yq))
            out('   class %s: G-solution +-%s, G = %d, m odd: no class-form solution (F_cl%s = %d = 16*%d; g = %d; it '
                'gives the non-integral rational point (%s, +-%s))'
                % (key, (a, b), cG, (m2, n2), 16 * cG, cG, g, Xq, abs(Yq)))
            continue
        a2, b2 = mc
        check(gcd(a2, b2) == 1, 'to_class coprime')
        check(Fclass(key, a2, b2) == cG, 'F_cl(to_class) = G')
        check(from_class(key, a2, b2) == (a, b), 'from_class(to_class) = id')
        check(to_class(key, -a, -b) == (-a2, -b2), 'to_class commutes with -1')
        CLASS_SOL[key].append(((a2, b2), cG, (a, b)))
        out('   class %s: G-solution +-%s, G = %d  ->  class-form solution +-%s, F_cl = %d'
            % (key, (a, b), cG, (a2, b2), cG))
for key in ('10', '01', '11'):
    for (mc, c, _) in CLASS_SOL[key]:
        check(gcd(*from_class(key, *mc)) == 1, 'from_class coprime')
out('THEOREM A (from Theorem A\'): the coprime solutions of F_cl(m,n) = c, c in {1,4,9,36}, up to sign:')
for key in ('10', '01', '11'):
    out('   class %s: %s' % (CLASSES[key]['label'],
                             ', '.join('%s -> %d' % (mc, c) for (mc, c, _) in sorted(CLASS_SOL[key], key=lambda t: (
                                 t[1], t[0])))))
check(sorted((mc, c) for mc, c, _ in CLASS_SOL['10']) == [((0, 1), 1), ((1, 0), 9)], 'Theorem A, class (1,0)')
check(sorted((mc, c) for mc, c, _ in CLASS_SOL['01']) == [((0, 1), 1), ((1, 0), 9), ((1, 3), 36)],
      'Theorem A, class (0,1)')
check(sorted((mc, c) for mc, c, _ in CLASS_SOL['11']) == [((0, 1), 1), ((1, 0), 36), ((1, 4), 4), ((2, 1), 1)],
      'Theorem A, class (1,1)')

out()
out('--- (iii.b) the points of the non-trivial classes')
out('%-6s %-4s %-12s %-10s %-20s %-3s %-6s %-9s %s' % ('class', 'c', 'G: +-(m,n)', 'F_cl:+-(m,n)', '(u,v,w)', 'g', 'X',
                                                       '|Y|', 'checks'))
for key in ('10', '01', '11'):
    Rk = RED[key]
    uvw_f, X4_f, Y6_f = descent_maps(key)
    for ((a, b), cG, (aG, bG)) in CLASS_SOL[key]:
        rows = []
        for sg in (1, -1):
            mm, nn = sg * a, sg * b
            check(gcd(mm, nn) == 1, 'coprime')
            c = Fclass(key, mm, nn)
            check(c in RHS and c == cG, 'F_cl in {1,4,9,36}')
            uu, vv, ww = uvw_f(mm, nn)
            g = gcd3(uu, vv, ww)
            check(c == g * g, 'F_cl = g^2')
            check(g == gcd(uu, 6), 'g = gcd(u, 6)')
            X4v, Y6v = X4_f(mm, nn), Y6_f(mm, nn)
            check(X4v % (g * g) == 0 and Y6v % g ** 3 == 0, 'divisibility by g^2, g^3')
            X, Y = X4v // (g * g), Y6v // g ** 3
            check(onE(X, Y) and onE(X, -Y), 'on E')
            xi = [Fr(uu, g) * Rk['bas'][0][k] + Fr(vv, g) * Rk['bas'][1][k] + Fr(ww, g) * Rk['bas'][2][k]
                  for k in range(3)]
            check(kmul(Rk['mu'], kmul(xi, xi)) == [X, -1, 0], 'X - th = mu xi^2')
            check(Rk['s'] * knorm(xi) == Y, 'Y = s N(xi)')
            rows.append((c, (uu, vv, ww), g, X, abs(Y)))
        check(rows[0] == rows[1], '(m,n) and (-m,-n) give the same (u,v,w), g, X, |Y|')
        c, uvw_, g, X, Y = rows[0]
        check(X not in PTS, 'X not produced twice')
        PTS[X] = dict(Y=Y, src=[(key, c, (a, b))])
        out('%-6s %-4d %-12s %-10s %-20s %-3d %-6d %-9d F_cl = g^2; X - th = mu xi^2; Y = %d N(xi); (-m,-n) the same'
            % (key, c, (aG, bG), (a, b), uvw_, g, X, Y, Rk['s']))
        # agreement with the file's table of the class
        hits = [r for r in RED[key]['ftab'] if r[0] == X]
        check(len(hits) == 2 and all(r[2] in ((a, b), (-a, -b)) and r[3] == g and r[4] == c for r in hits),
              'file table of the class agrees for X = %d' % X)
ALLX = sorted(PTS)
EXPECT = [-5, -3, 0, 3, 7, 9, 24, 33, 39, 513, 1099, 5112]
check(ALLX == EXPECT, 'the 12 X values')
POINTS24 = sorted([(X, PTS[X]['Y']) for X in ALLX] + [(X, -PTS[X]['Y']) for X in ALLX])
check(len(POINTS24) == 24 and all(PTS[X]['Y'] > 0 for X in ALLX), '24 points')
out()
out('trivial class (0,0) (part ii): X = -5, 7, 1099 from +-(1,1), +-(1,-1), +-(5,-3) with c = F0 = 3, g = 2.')
# agreement with the file's table of integral points (and their classes from the Mordell-Weil decomposition)
seg = SEC[SEC.index('integral points (monic coordinates), decomposition'):SEC.index('grouped by class (X values')]
ftab = [(int(a), int(b), (int(c), int(d)), (int(e), int(f))) for a, b, c, d, e, f in
        re.findall(r'^\s+(-?\d+)\s+(-?\d+)\s+\((-?\d+), (-?\d+)\)\s+\((\d), (\d)\)\s*$', seg, re.M)]
check(sorted((a, b) for a, b, _, _ in ftab) == POINTS24, 'the 24 points = the file\'s list')
keycl = {'00': (0, 0), '10': (1, 0), '01': (0, 1), '11': (1, 1)}
check(all(keycl[PTS[a]['src'][0][0]] == cl for a, b, ij, cl in ftab), 'classes agree with the MW decomposition')
check(all((ij[0] % 2, ij[1] % 2) == cl for a, b, ij, cl in ftab), 'class = (i mod 2, j mod 2) in the file table')
IJ = {(a, b): ij for a, b, ij, cl in ftab}
check(all(IJ[(X, -Y)] == (-IJ[(X, Y)][0], -IJ[(X, Y)][1]) for X, Y in POINTS24 if Y > 0), '(X,-Y) = -(X,Y)')
# the decomposition, re-checked by the exact group law
for (X, Y) in POINTS24:
    i_, j_ = IJ[(X, Y)]
    S_ = None
    for _ in range(abs(i_)):
        S_ = ec_add(S_, G1 if i_ > 0 else (G1[0], -G1[1]))
    for _ in range(abs(j_)):
        S_ = ec_add(S_, G2 if j_ > 0 else (G2[0], -G2[1]))
    check(S_ == (Fr(X), Fr(Y)), 'file decomposition i G1 + j G2 re-checked by the group law')
grp = {}
for mt in re.finditer(r'^\s+class \((\d), (\d)\), phi = \[(.*?)\]: X in (\[.*\])$', SEC, re.M):
    grp[(int(mt.group(1)), int(mt.group(2)))] = plist(mt.group(4))
check(sorted(grp) == [(0, 0), (0, 1), (1, 0), (1, 1)], 'grouped-by-class lines')
check(all(grp[keycl[k]] == sorted(X for X in ALLX if PTS[X]['src'][0][0] == k) for k in keycl), 'grouping = file')
out('X values: %s; each with both signs of Y: %d integral points (X, Y):' % (ALLX, len(POINTS24)))
out('   ((m,n) and (-m,-n) give the same (u,v,w), X4, Y6 (forms of even degree); the two signs of Y come from +-xi)')
for X in ALLX:
    src = PTS[X]['src'][0]
    out('   X = %-5d Y = +-%-7d  from class %s, c = %-2d (m,n) = +-%-9s  file: (X,+|Y|) = %s G1 + %s G2 (re-checked)'
        % (X, PTS[X]['Y'], src[0], src[1], src[2], IJ[(X, PTS[X]['Y'])][0], IJ[(X, PTS[X]['Y'])][1]))
out('   per class: ' + '; '.join('%s: X in %s' % (keycl[k], sorted(X for X in ALLX if PTS[X]['src'][0][0] == k))
                                 for k in ('00', '10', '01', '11')))
out('agreement with the file: the same 24 points, and the class of each = (i mod 2, j mod 2) of its Mordell-Weil '
    'decomposition i G1 + j G2: True')
out()
out('Solutions that do NOT give integral points:')
out('   class forms, non-trivial classes: none (every listed solution has g = sqrt(c); by u^2 = F_cl (mod 12) this '
    'holds for every coprime solution).')
for (key, (a, b), cG, (m2, n2), g, Xq, Yq) in SKIPPED:
    out('   field form of class %s: G%s = %d with m odd: no class-form solution (F_cl%s = %d); rational point '
        '(%s, +-%s), not integral' % (key, (a, b), cG, (m2, n2), 16 * cG, Xq, abs(Yq)))
for ((a, b), g, Xq, Yq) in NOT_INTEGRAL:
    out('   trivial class: (m,n) = %s, F0 = 3 but g = %d (H = 4 != 1): rational point (%s, +-%s), not integral'
        % ((a, b), g, Xq, abs(Yq)))
# consistency checks (not proofs): Thue box and integral-point search
out()
out('consistency checks (finite searches, not part of the argument):')
BOX = 300
for key in ('10', '01', '11'):
    foundG, foundC = set(), set()
    for mm in range(-BOX, BOX + 1):
        for nn in range(0, BOX + 1):
            if (nn == 0 and mm <= 0) or gcd(mm, nn) != 1:
                continue
            if Fval(key, mm, nn) in RHS:
                foundG.add((mm, nn))
            if Fclass(key, mm, nn) in RHS:
                foundC.add((mm, nn))
    knownG = set(pm_rep(a, b) for (a, b) in KNOWN[key])
    knownC = set(pm_rep(*mc) for (mc, c, _) in CLASS_SOL[key])
    check(foundG == knownG, 'box search G = KNOWN for %s' % key)
    check(foundC == knownC, 'box search F_cl = Theorem A list for %s' % key)
    out('   class %s, coprime (m,n), |m|,|n| <= %d: G(m,n) in {1,4,9,36} exactly +-%s (= core.KNOWN); F_cl(m,n) in '
        '{1,4,9,36} exactly +-%s'
        % (key, BOX, sorted(foundG), sorted(foundC)))
XMAX = 10 ** 6
check((-6) ** 3 + GA * (-6) + GB < 0, 'g(-6) < 0')
srch = []
for X in range(-5, XMAX + 1):
    D = X ** 3 + GA * X + GB
    if D >= 0:
        r = isqrt(D)
        if r * r == D:
            srch.append(X)
check(srch == EXPECT, 'search agrees')
out('   integral points of E with X <= 10^6 (none with X <= -6: g\'(X) = 3X^2 - 9 > 0 for X <= -2, so g is increasing '
    'there, and g(-6) = -81 < 0): X in %s (= the 12 above)' % srch)

# ==========================================================================================================
out()
out('=' * 110)
out('(iv) TRANSFER: LMFDB, n(n+1)/2 = m(m+1)(2m+1)/6, and the model 3Y^2 = X^3 - X + 3')
out('=' * 110)
out('--- LMFDB model y^2 = x^3 - 9x + 81 = E (X = x, Y = y)')
with open(LMFDB_MD, encoding='utf-8') as fh:
    md = fh.read().replace(chr(0x2212), '-').replace(chr(0xB1), '+-')   # LMFDB text uses the Unicode minus sign
line = [ln for ln in md.splitlines() if ln.startswith('| 8604.a1 |') and '(y = ' in ln]
check(len(line) == 1, 'LMFDB line found')
xs_part, ys_part = line[0].split('|')[2].split('(y = ')
lx = [int(t) for t in re.findall(r'-?\d+', xs_part)]
ly = [int(t) for t in re.findall(r'\d+', ys_part)]
check(len(lx) == len(ly) == 12, '12 x-values with 12 |y|-values in lmfdb_data.md')
lm = sorted([(a, b) for a, b in zip(lx, ly)] + [(a, -b) for a, b in zip(lx, ly)])
check(lm == POINTS24, 'equal to the LMFDB list')
out('the 24 integral points (x, y) = (X, Y) = the 24 points listed by LMFDB (lmfdb_data.md: x in %s, both signs of y): '
    'True' % lx)
out()
out('--- n(n+1)/2 = m(m+1)(2m+1)/6 via x = 6m + 3, y = 18n + 9')
M_, N_, t_ = sp.symbols('M N t')
idt = sp.expand((18 * N_ + 9) ** 2 - ((6 * M_ + 3) ** 3 - 9 * (6 * M_ + 3) + 81)
                - 648 * (N_ * (N_ + 1) / 2 - M_ * (M_ + 1) * (2 * M_ + 1) / 6))
check(idt == 0, 'correspondence identity')
out('identity: (18n+9)^2 - ((6m+3)^3 - 9(6m+3) + 81) = 648 (T(n) - P(m)),  T(n) = n(n+1)/2, P(m) = m(m+1)(2m+1)/6  '
    '[exact]')
check(sp.expand((6 * t_ + 3) ** 3 - 9 * (6 * t_ + 3) + 81 - 81 * (8 * t_ * (t_ + 1) * (2 * t_ + 1) / 6 + 1)) == 0,
      'x = 6t+3 => x^3 - 9x + 81 = 81 (8 P(t) + 1)')
check(all((t * (t + 1) * (2 * t + 1)) % 6 == 0 for t in range(6)), 'P(t) is an integer')
out('so integer solutions (m,n) <-> integral points with x = 3 (mod 6) and y = 9 (mod 18); and x = 6t + 3 already '
    'forces y = 9 (mod 18):')
out('   x^3 - 9x + 81 = 81 (8 P(t) + 1) identically, with P(t) an integer (t(t+1)(2t+1) = 0 mod 6), so 9 | y and '
    '(y/9)^2 = 8 P(t) + 1 is odd.')
out('   Integral points with x = 3 (mod 6) -> (m, n) = ((x - 3)/6, (y - 9)/18):')
MN = []
for (x, y) in POINTS24:
    if x % 6 == 3:
        check(y % 18 == 9, 'y = 9 mod 18')
        mm, nn = (x - 3) // 6, (y - 9) // 18
        Pm, Tn = mm * (mm + 1) * (2 * mm + 1) // 6, nn * (nn + 1) // 2
        check(Pm == Tn, 'T(n) = P(m)')
        MN.append((mm, nn, Pm))
MN.sort()
check(len(MN) == 12, '12 pairs')
vals = sorted(set(c for _, _, c in MN))
check(vals == [0, 1, 55, 91, 208335], 'values')
check(sorted(set(a for a, _, _ in MN)) == [-1, 0, 1, 5, 6, 85], 'm values')
out('   x values used: %s (the other X: %s are not 3 mod 6)' % (sorted(set(x for x, _ in POINTS24 if x % 6 == 3)),
                                                              [X for X in ALLX if X % 6 != 3]))
for val in vals:
    out('   P(m) = T(n) = %-7d  (m, n) in %s' % (val, [(a, b) for a, b, c in MN if c == val]))
out('   the involution y -> -y of the curve is n -> -1 - n (T(-1-n) = T(n)); m <= -2 is impossible (P(m) < 0 <= T(n)); '
    'm = -1 gives P = 0.')
out('   values: {0, 1, 55, 91, 208335}; the numbers > 0 both triangular and square-pyramidal: 1, 55, 91, 208335.')
for val in vals[1:]:
    (mm, nn) = max((a, b) for a, b, c in MN if c == val)
    check(nn * (nn + 1) // 2 == val == mm * (mm + 1) * (2 * mm + 1) // 6 and mm > 0 and nn > 0, 'T_n = P_m')
    out('      %-7d = T_%d = n(n+1)/2 with n = %d,  = P_%d = m(m+1)(2m+1)/6 with m = %d' % (val, nn, nn, mm, mm))
srch_m = []
for mm in range(-1, 10 ** 6 + 1):
    D = 8 * (mm * (mm + 1) * (2 * mm + 1) // 6) + 1
    r = isqrt(D)
    if r * r == D:
        srch_m.append(mm)
check(srch_m == [-1, 0, 1, 5, 6, 85], 'search over m agrees')
out('   consistency (search, not proof): -1 <= m <= 10^6 with P(m) triangular (8 P(m) + 1 a square): m in %s' % srch_m)
out()
out('--- the model 3Y^2 = X^3 - X + 3 (Uchiyama, Tsukuba J. Math. 8 (1984)) via x = 3X, y = 9Y')
Xu_, Yu_ = sp.symbols('Xu Yu')
check(sp.expand((9 * Yu_) ** 2 - ((3 * Xu_) ** 3 - 9 * (3 * Xu_) + 81)
                - 27 * (3 * Yu_ ** 2 - (Xu_ ** 3 - Xu_ + 3))) == 0, 'E = 27 * model')
check(sp.expand(3 * (2 * N_ + 1) ** 2 - ((2 * M_ + 1) ** 3 - (2 * M_ + 1) + 3)
                - 24 * (N_ * (N_ + 1) / 2 - M_ * (M_ + 1) * (2 * M_ + 1) / 6)) == 0, 'model at (2m+1, 2n+1)')
out('identities: (9Y)^2 - ((3X)^3 - 9(3X) + 81) = 27 (3Y^2 - (X^3 - X + 3));  3(2n+1)^2 - ((2m+1)^3 - (2m+1) + 3) = '
    '24 (T(n) - P(m))  [exact]')
out('   so (X, Y) = (x/3, y/9); and (X, Y) = (2m+1, 2n+1) turns the model into T(n) = P(m) (X = x/3 = 2m+1, Y = y/9 = '
    '2n+1)')
check(sp.expand((3 * Xu_) ** 3 - 9 * (3 * Xu_) + 81 - 27 * (Xu_ ** 3 - Xu_ + 3)) == 0, 'x = 3k => rhs = 27 (k^3-k+3)')
check(all((y_ * y_) % 27 != 0 or y_ % 9 == 0 for y_ in range(27)), '27 | y^2 => 9 | y (residues mod 27)')
out('   integral points of the model <-> integral points of E with 3 | x: then y^2 = 27 (k^3 - k + 3) for x = 3k, so '
    '27 | y^2 and 9 | y (residues mod 27)')
UCH = []
for (x, y) in POINTS24:
    if x % 3 == 0:
        check(y % 9 == 0, '9 | y')
        Xm, Ym = x // 3, y // 9
        check(3 * Ym * Ym == Xm ** 3 - Xm + 3, 'on 3Y^2 = X^3 - X + 3')
        UCH.append((Xm, Ym))
UCH.sort()
check(len(UCH) == 18, '18 points on the model')
uX = sorted(set(a for a, _ in UCH))
check(uX == [-1, 0, 1, 3, 8, 11, 13, 171, 1704], 'X values on the model')
out('   x divisible by 3 among the 12: %s; the 18 integral points (X, Y) of 3Y^2 = X^3 - X + 3:' % [3 * a for a in uX])
out('   ' + ', '.join('(%d,%d)' % p for p in UCH))
oddX = sorted(set(a for a, _ in UCH if a % 2))
check(oddX == [2 * a + 1 for a in (-1, 0, 1, 5, 6, 85)], 'odd X = 2m+1')
check(all(b % 2 == 1 for a, b in UCH if a % 2), 'X odd => Y odd')
out('   odd X: %s = 2m + 1 for m in [-1, 0, 1, 5, 6, 85] (and then Y is odd, Y = 2n + 1): the same 12 pairs (m,n) as '
    'above' % oddX)
check(sorted(((a - 1) // 2, (b - 1) // 2) for a, b in UCH if a % 2) == sorted((a, b) for a, b, _ in MN),
      'same (m,n) as via x = 6m+3')

# ==========================================================================================================
out()
out('=' * 110)
out('CONCLUSION (conditional)')
out('=' * 110)
out('Assume (A) rank E(Q) = 2 for 8604.a1 [LMFDB/mwrank; only rank <= 2 is used, through (H)] and (B) Theorem A\': '
    'for each class')
out('(1,0), (0,1), (1,1) the coprime solutions of the field form G(m,n) in {1,4,9,36} are exactly +-core.KNOWN '
    '(01_coverage.py, 02_exclusion.py).')
out('Then, using also N1, N2 (index-discriminant formula, integral closure):')
out('   the coprime solutions of the class forms F_cl(m,n) in {1,4,9,36} are those of Theorem A (iii.a);')
out('   the integral points of y^2 = x^3 - 9x + 81 are exactly the 24 points x in %s, both signs of y;' % ALLX)
NONNEG = sorted((a, b) for a, b, _ in MN if a >= 0 and b >= 0)
check(NONNEG == [(0, 0), (1, 1), (5, 10), (6, 13), (85, 645)], 'solutions with m, n >= 0')
out('   n(n+1)/2 = m(m+1)(2m+1)/6 holds exactly for the 12 pairs (m,n) above, with values 0, 1, 55, 91, 208335;')
out('      with m, n >= 0 exactly (m,n) in %s;' % NONNEG)
out('   3Y^2 = X^3 - X + 3 has exactly the 18 integral points above (X in %s).' % uX)
out('Solving the trivial-class equation F0 = 3 (part ii) needs neither (A) nor (B) and is elementary.  Every class '
    'map, value and identity above is exact.')
out()
out('all %d checks passed  [done in %.1f s]' % (NCHECK[0], time.time() - T0))
