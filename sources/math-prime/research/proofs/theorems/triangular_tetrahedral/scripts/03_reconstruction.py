#!/usr/bin/env python3
"""03_reconstruction.py -- from the Thue solutions back to the integral points of 1611.a1, and on to
C(n,2) = C(m,3) (numbers both triangular and tetrahedral) and to Avanesov's 3y(y+1) = x(x+1)(x+2).

Curve 1611.a1 (LMFDB): y^2 + y = x^3 - 9x + 20.  Monic model E: Y^2 = X^3 - 144 X + 1296 (X = 4x, Y = 8y + 4).
K = Q(th), th^3 - 144 th + 1296 = 0.

(i)   THE REDUCTION.  Every item of the 2-descent that the reduction uses is read from
      research/proofs/targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt (section CURVE 1611.a1) and
      re-checked here by exact arithmetic: the maximal order O_K, the four classes and their distinctness, the
      lattices Xi_mu = {xi in K : mu xi^2 in O_K}, the forms A0, A1, A2, the conic, the completeness of its
      parametrisation, H = -A1(X(m,n)) = c F, the possible values of g, k = 1, X4 = A0(X(m,n)),
      Y6 = sqrt(N mu) N(xi), and Y6^2 = X4^3 - 144 X4 H^2 + 1296 H^3.  The chain of reasoning and the facts it
      uses are printed; the one non-elementary input is rank E(Q) = 2 (LMFDB, mwrank).
(ii)  THE TRIVIAL CLASS.  F0(m,n) = n (m^3 + 3 m^2 n - m n^2 + 3 n^3) = 6, solved completely and elementarily.
(iii) THE INTEGRAL POINTS, assuming that the coprime solutions of F(m,n) in {1,4,9,36} are +-core.KNOWN (proved
      separately): exactly 28 points.
(iv)  TRANSFER to the LMFDB model (22 points), to C(n,2) = C(m,3) and to 3y(y+1) = x(x+1)(x+2).

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
from core import KNOWN, RHS, Fval, descent_maps  # noqa: E402

T0 = time.time()
TARGET_DIR = os.path.normpath(os.path.join(HERE, '..', '..', '..', 'targets', 'rank2_signatures_2026-10'))
DESC = os.path.join(TARGET_DIR, 'weierstrass_2descent.out.txt')
LMFDB_MD = os.path.join(TARGET_DIR, 'lmfdb_data.md')
GA, GB = -144, 1296                        # g(t) = t^3 + GA t + GB;  E: Y^2 = g(X)

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
# arithmetic in K = Q[th]/(th^3 - 144 th + 1296); an element is [c0, c1, c2] = c0 + c1 th + c2 th^2
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
    """coordinates in the basis w0 = 1, w1 = th/6, w2 = th^2/36 of O_K"""
    return [Fr(al[0]), 6 * Fr(al[1]), 36 * Fr(al[2])]


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
    """(p, r): p prime, p not dividing 6, g(r) = 0 mod p, al(r) a non-residue mod p (al in Z[th])"""
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


# ==========================================================================================================
out('03_reconstruction.py -- from the Thue solutions back to the integral points of 1611.a1 (exact arithmetic)')
out('python %s, sympy %s' % (sys.version.split()[0], sp.__version__))
out('2-descent data read from research/proofs/targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt')
out('(section CURVE 1611.a1); every item used is re-checked below (a failing check stops the script, exit status 1).')
out()
with open(DESC, encoding='utf-8') as fh:
    TXT = fh.read()
i0, i1 = TXT.index('CURVE 1611.a1'), TXT.index('CURVE 8604.a1')
SEC = TXT[i0:i1]

# ==========================================================================================================
out('=' * 110)
out('(i) THE REDUCTION: every integral point of E gives coprime (m,n) with H(m,n) = g^2')
out('=' * 110)
out('--- (i.a) the curve, the field K and its maximal order')
mt = grab(r'monic model  E: Y\^2 = X\^3 \+ \((-?\d+)\) X \+ \((-?\d+)\),  map X = 4x, Y = 8y \+ 4', SEC,
          'the monic model')
check((int(mt.group(1)), int(mt.group(2))) == (GA, GB), 'monic model Y^2 = X^3 - 144 X + 1296')
xs, ys = sp.symbols('x y')
check(sp.expand((8 * ys + 4) ** 2 - ((4 * xs) ** 3 - 144 * (4 * xs) + 1296)
                - 64 * (ys ** 2 + ys - (xs ** 3 - 9 * xs + 20))) == 0, 'monic = 64 * LMFDB model')
out('E: Y^2 = X^3 - 144 X + 1296;  with X = 4x, Y = 8y + 4:  Y^2 - (X^3 - 144X + 1296) = '
    '64 (y^2 + y - x^3 + 9x - 20)  [exact]')
divs = sp.divisors(GB)
check(all(r ** 3 + GA * r + GB != 0 for r in divs + [-d for d in divs]), 'g has no rational root')
out('g(t) = t^3 - 144 t + 1296 has no rational root (rational-root test over the divisors of 1296): g is '
    'irreducible, K = Q(th) is a cubic field, E(Q)[2] = 0')
mt = grab(r'O_K = Z-span of  w0 = (.+),  w1 = (.+),  w2 = (.+)$', SEC, 'the basis of O_K')
WB = [pK(mt.group(k)) for k in (1, 2, 3)]
check(WB == [[1, 0, 0], [0, Fr(1, 6), 0], [0, 0, Fr(1, 36)]], 'O_K basis read from the file')
check(all(in_OK(kmul(WB[i], WB[j])) for i in range(3) for j in range(3)), 'O_K closed under multiplication')
e1, e2, e3 = charpoly(WB[1])
check((e1, e2, e3) == (0, -4, -6), 'lambda = th/6 is a root of T^3 - 4T + 6')
out('file: O_K = Z w0 + Z w1 + Z w2 with w0 = 1, w1 = th/6, w2 = th^2/36.  Re-checked: the structure constants '
    'w_i w_j are integral (an order);')
out('   lambda := w1 = th/6 has characteristic polynomial T^3 - 4T + 6, so w2 = lambda^2 and O_K = Z[lambda]: '
    'the ring of Avanesov\'s cubic lambda^3 - 4 lambda + 6 = 0 (his X = X_monic/6)')
tr = lambda al: charpoly(al)[0]                                      # noqa: E731
dO = det3([[tr(kmul(WB[i], WB[j])) for j in range(3)] for i in range(3)])
dZ = det3([[tr(kmul([int(i == 0), int(i == 1), int(i == 2)], [int(j == 0), int(j == 1), int(j == 2)]))
            for j in range(3)] for i in range(3)])
check(dZ == -4 * GA ** 3 - 27 * GB ** 2 == -33405696, 'disc Z[th]')
check(dO == -716 and Fr(dZ, dO) == 216 ** 2, 'disc O_K = -716, index 216')
check(int(grab(r'disc\(O_K\) = (-?\d+)', SEC, 'disc O_K').group(1)) == dO, 'disc O_K agrees with the file')
sqdiv = [k for k in range(1, isqrt(716) + 1) if 716 % (k * k) == 0]
check(sqdiv == [1, 2], 'square divisors of 716')
half = []
for c0 in (0, 1):
    for c1 in (0, 1):
        for c2 in (0, 1):
            if (c0, c1, c2) == (0, 0, 0):
                continue
            al = [Fr(c0 * WB[0][k] + c1 * WB[1][k] + c2 * WB[2][k]) / 2 for k in range(3)]
            half.append(all(Fr(c).denominator == 1 for c in charpoly(al)))
check(not any(half), 'no (c0 w0 + c1 w1 + c2 w2)/2 is integral')
out('disc(Z[th]) = -33405696, disc(O_K) = -716 = -2^2*179, [O_K : Z[th]] = 216.  Maximality (N1, N2): the index '
    'of O_K in the ring of integers')
out('   squares into 716, so it is 1 or 2; none of the 7 elements (c0 w0 + c1 w1 + c2 w2)/2, c in {0,1}^3 - 0, has '
    'an integral characteristic polynomial; hence O_K is the ring of integers of K.')
check(int(grab(r'\(D1 = (\d+)\)', SEC, 'D1').group(1)) == 6, 'D1 = 6 in the file')
out('the th-coordinate of a w0 + b w1 + c w2 is b/6: th-coordinates of elements of O_K lie in (1/6)Z  (D1 = 6)')

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
for cl in sorted(reps):
    R_ = reps[cl]
    check(onE(*R_['P']), 'representative point on E')
    check(knorm(R_['mu']) == R_['N'], 'N(mu) agrees with the file')
    p, r = R_['p'], R_['r']
    check((r ** 3 + GA * r + GB) % p == 0 and legendre(R_['mu'][0] - r, p) == -1, 'file witness re-checked')
    out('class %s: mu = %s from the point (%d,%d) on E (= %d*G1 + %d*G2 in the file); N(mu) = %d; file witness '
        '(p,r) = (%d,%d) re-checked' % (cl, kstr(R_['mu']), R_['P'][0], R_['P'][1], R_['ij'][0], R_['ij'][1],
                                         R_['N'], p, r))
MU = {(0, 0): [1, 0, 0], (1, 0): reps[(1, 0)]['mu'], (0, 1): reps[(0, 1)]['mu'], (1, 1): reps[(1, 1)]['mu']}
names = {(1, 0): 'mu10', (0, 1): 'mu01', (1, 1): 'mu11'}
tests = [((1, 0),), ((0, 1),), ((1, 1),), ((1, 0), (0, 1)), ((1, 0), (1, 1)), ((0, 1), (1, 1))]
for tcl in tests:
    al = MU[tcl[0]] if len(tcl) == 1 else kmul(MU[tcl[0]], MU[tcl[1]])
    wit = nonsquare_witness(al)
    check(wit is not None, 'non-square witness for %s' % (tcl,))
    out('   %-12s = %-28s is not a square in K: at the root r = %d of g mod p = %d its value is a non-residue'
        % ('*'.join(names[c] for c in tcl), kstr(al), wit[1], wit[0]))
out('   (if al in O_K were xi^2 with xi in K, xi would be integral (N2), so in O_K, which lies in Z[th][1/6]; the '
    'ring map Z[th][1/6] -> F_p, th -> r (p not dividing 6) would make al(r) a square mod p.)')
out('=> the classes [1], [mu10], [mu01], [mu11] are four DISTINCT elements of phi(E(Q)) (phi(O) = [1], and each mu '
    'is X0 - th of a point of E).')
# the line identity behind the homomorphism property, on a sample
sample = [(Fr(9), Fr(27)), (Fr(-15), Fr(9)), (Fr(-12), Fr(36)), (Fr(0), Fr(36)), (Fr(12), Fr(36)),
          (Fr(4), Fr(28)), (Fr(-8), Fr(44)), (Fr(48), Fr(324))]
nline = 0
for P1 in sample:
    for P2 in sample:
        S = ec_add(P1, P2)
        if S is None:
            continue
        lam = (3 * P1[0] ** 2 + GA) / (2 * P1[1]) if P1 == P2 else (P2[1] - P1[1]) / (P2[0] - P1[0])
        nu = P1[1] - lam * P1[0]
        lhs = kmul(kmul([P1[0], -1, 0], [P2[0], -1, 0]), [S[0], -1, 0])
        check(lhs == kmul([nu, lam, 0], [nu, lam, 0]), 'line identity')
        nline += 1
out('phi is a homomorphism (2-descent): for P1 + P2 + P3 = O on Y = lam X + nu, g(X) - (lam X + nu)^2 = '
    '(X - x1)(X - x2)(X - x3), and X = th gives')
out('   (x1 - th)(x2 - th)(x3 - th) = (lam th + nu)^2  [identity re-checked exactly on %d pairs of points]' % nline)
out('COMPLETENESS OF THE CLASS LIST (the non-elementary input): |phi(E(Q))| <= |E(Q)/2E(Q)| = 2^rank * |E(Q)[2]| '
    '= 2^rank (E(Q)[2] = 0 above;')
out('   E(Q) finitely generated by Mordell-Weil).  With rank E(Q) = 2 (LMFDB 1611.a1, computed by mwrank; only '
    'rank <= 2 is used) this is 4,')
out('   so phi(E(Q)) = {[1], [9 - th], [-15 - th], [-12 - th]}: every point P != O has X(P) - th = mu xi^2 with '
    'xi in K* for one of these four mu.')

# ----------------------------------------------------------------------------------------------------------
cls_pos = [(mt.start(), (int(mt.group(1)), int(mt.group(2))))
           for mt in re.finditer(r'^CLASS \((\d), (\d)\)', SEC, re.M)]
iextra = SEC.index('EXTRA for 1611.a1')
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
    # enumeration of (1/s) O_K / O_K
    S = []
    D0, D1_, D2 = 1296 * s * s, 216 * s * s, 36 * s * s
    mui = [int(c) for c in mu]
    for c0 in range(s):
        for c1 in range(s):
            for c2 in range(s):
                y = [36 * c0, 6 * c1, c2]                              # xi = y / (36 s)
                z = kmul(mui, kmul(y, y))                              # mu xi^2 = z / (1296 s^2)
                if z[0] % D0 == 0 and z[1] % D1_ == 0 and z[2] % D2 == 0:
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
    out('parametrisation (file = core.descent_maps): u = %s, v = %s, w = %s' % (fstr(Xu), fstr(Xv), fstr(Xw)))
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
    mt = grab(r'so k\^2 \| D1 = (\d+): k in \[(.*)\]', Bk, 'k')
    check(int(mt.group(1)) == 6 and plist('[' + mt.group(2) + ']') == [1], 'k = 1 in the file')
    check([k for k in range(1, 3) if 6 % (k * k) == 0] == [1], 'k^2 | 6 => k = 1')
    out('g := gcd(u,v,w) for coprime (m,n): T is unimodular, so gcd(m\',n\') = 1 and gcd(v,w) = |l| = %d|%s|; '
        'q = %s^2 mod %s gives gcd(q,%s) = 1;'
        % (abs(int(c_l)), xv, yv, xv, xv))
    out('   so g = gcd(q, %d) divides %d: g in %s (file, by CRT: %s)' % (abs(int(c_l)), abs(int(c_l)), gvals, gfile))
    out('k = 1: if (u,v,w) = k (u\',v\',w\'), then (X - th)/k^2 = mu (xi/k)^2 is in O_K, its th-coordinate -1/k^2 is in '
        '(1/6)Z, so k^2 | 6 and k = 1')
    out('   => (u,v,w) is primitive and (u,v,w) = +-X(m,n)/g for coprime (m,n)')
    # H = c F
    cH = sp.Rational(grab(r'with c = (\S+) and primitive', Bk, 'c').group(1))
    Fp = pexpr(grab(r'^\s+F\(m,n\) = (.+)$', Bk, 'F').group(1))
    if key != '00':
        check(sp.expand(Fval(key, m, n) - Fp) == 0, 'core F = file F')
    H = sp.expand(-A1.subs(subX, simultaneous=True))
    check(sp.expand(H - cH * Fp) == 0, 'H = -A1(X(m,n)) = c F')
    rhs_vals = sorted(set(sp.Rational(g * g) / cH for g in gvals))
    rhs_int = [int(r) for r in rhs_vals if r.is_integer]
    rfile = plist('[' + grab(r'so F\(m,n\) in \[(.*)\]', Bk, 'F values').group(1) + ']')
    check(rhs_int == rfile, 'right-hand sides agree with the file')
    out('A1(u,v,w) = -1 and (u,v,w) = +-X(m,n)/g give H(m,n) := -A1(X(m,n)) = g^2, with H = c F, c = %s,' % cH)
    out('   F(m,n) = %s;  so F(m,n) = g^2/c with g in %s: the integral values are F(m,n) in %s' % (fstr(Fp), gvals,
                                                                                               rhs_int))
    if key != '00':
        Dq = sp.Poly(sp.expand(Xu ** 2 - Fp), m, n)
        check(all(int(c) % 6 == 0 for c in Dq.coeffs()), 'u^2 = F mod 6 identically')
        for c in RHS:
            check(all(gcd(r, 6) == isqrt(c) for r in range(6) if (r * r - c) % 6 == 0), 'u^2 = c mod 6 => g = sqrt c')
        out('   moreover u^2 - F = 6 * (%s) identically, so u^2 = F (mod 6), and since g = gcd(u,6): every coprime '
            'solution of F(m,n) = c in {1,4,9,36} has g = sqrt(c)' % fstr(sp.expand((Xu ** 2 - Fp) / 6)))
    else:
        check(all(((mm ** 2 + 2 * mm * nn + 5 * nn ** 2) % 2 == 0) == ((mm - nn) % 2 == 0)
                  for mm in range(2) for nn in range(2)), 'u even iff m = n mod 2')
        out('   here H = (2/3) F0 = g^2 with g | 2: g = 1 would give F0 = 3/2, so F0(m,n) = 6 and g = 2, i.e. '
            'u = m^2 + 2mn + 5n^2 even, i.e. m = n (mod 2)')
    # X4, Y6
    X4 = sp.expand(A0.subs(subX, simultaneous=True))
    check(sp.expand(X4 - pexpr(grab(r'X4 = (.+)$', Bk, 'X4').group(1))) == 0, 'X4 = file')
    xiP = [sp.expand(e.subs(subX, simultaneous=True)) for e in xi_s]
    Y6 = sp.expand(s * det3(kmat(xiP)))
    check(sp.expand(Y6 - pexpr(grab(r'Y6 = (.+)$', Bk, 'Y6').group(1))) == 0, 'Y6 = file')
    if key != '00':
        check(sp.expand(X4 - X4_f(m, n)) == 0 and sp.expand(Y6 - Y6_f(m, n)) == 0, 'core X4, Y6 = recomputed')
    check(sp.expand(Y6 ** 2 - (X4 ** 3 - 144 * X4 * H ** 2 + 1296 * H ** 3)) == 0, 'Y6^2 = X4^3 - 144 X4 H^2 + ...')
    if key != '00':
        check(sp.expand(Y6 ** 2 - (X4 ** 3 - 144 * X4 * Fp ** 2 + 1296 * Fp ** 3)) == 0, 'identity with F')
    out('X = A0(u,v,w) = X4(m,n)/g^2 with X4 = A0(X(m,n)) = %s' % fstr(X4))
    out('Y^2 = N(X - th) = N(mu) N(xi)^2, Y = +-%d N(xi) = +-Y6(m,n)/g^3 with Y6 = %d N(X(m,n).b) = %s'
        % (s, s, fstr(Y6)))
    out('IDENTITY  Y6^2 = X4^3 - 144 X4 %s^2 + 1296 %s^3  holds in Z[m,n]  [exact]%s'
        % ('F' if key != '00' else 'H', 'F' if key != '00' else 'H',
           '' if key != '00' else '   (H = (2/3) F0; with F0 itself in place of H it does not hold)'))
    if key == '00':
        check(sp.expand(Y6 ** 2 - (X4 ** 3 - 144 * X4 * Fp ** 2 + 1296 * Fp ** 3)) != 0, 'not with F0')
    out('   so whenever H(m,n) = g^2:  (X,Y) = (X4/g^2, +-Y6/g^3) satisfies Y^2 = X^3 - 144X + 1296; X is an '
        'integer (A0 has integer coefficients, X(m,n)/g is integral)')
    out('   and Y is then an integer (Y^2 is).  Conversely X - th = mu xi^2 with xi = (X(m,n)/g).b, so the point '
        'lies in the class [mu].')
    # the file's table of known points of the class
    ftab = [(int(a), int(b), (int(c), int(d)), int(e), int(f)) for a, b, c, d, e, f in
            re.findall(r'^\s+(-?\d+)\s+(-?\d+)\s+\((-?\d+), (-?\d+)\)\s+(\d+)\s+(\d+)\s*$', Bk, re.M)]
    check(len(ftab) > 0, 'file table of the class')
    RED[key] = dict(cl=cl, mu=mu, s=s, bas=bas, Xp=Xp, X4=X4, Y6=Y6, H=H, F=Fp, cH=cH, gvals=gvals, ftab=ftab)

out()
out('--- (i.d) SUMMARY OF THE REDUCTION (chain of reasoning, facts used)')
for line in [
    '1. [class list; NON-ELEMENTARY INPUT] phi: E(Q) -> K*/K*^2, P -> [X(P) - th], is a homomorphism (line identity) '
    'with 2E(Q) in its kernel.',
    '   E(Q)[2] = 0 (g irreducible).  rank E(Q) = 2 (LMFDB, mwrank) and Mordell-Weil give |E(Q)/2E(Q)| = 4.  The four '
    'classes [1], [9-th],',
    '   [-15-th], [-12-th] are distinct (Legendre witnesses) and lie in phi(E(Q)); hence they are all of it.  Only '
    'rank <= 2 is used.',
    '2. [lattice; N1, N2] O_K = Z[th/6] is the ring of integers (index-discriminant formula + 7 half-elements).  For '
    'an integral point in class [mu],',
    '   X - th lies in O_K, so xi in Xi_mu = Z b1 + Z b2 + Z b3 (enumeration of (1/s)O_K/O_K, which uses integral '
    'closure), xi = u b1 + v b2 + w b3.',
    '3. [conic; elementary] A2(u,v,w) = 0, A1 = -1, A0 = X.  The conic A2 = 0 is smooth with the point P0 = (1,0,0); '
    'by the explicit inverse identity',
    '   every rational point is ~ X(m,n) for coprime integers (m,n) (unique up to sign).',
    '4. [primitivity; elementary] k^2 | 6 (th-coordinates of O_K in (1/6)Z) forces (u,v,w) primitive: (u,v,w) = '
    '+-X(m,n)/g, g = gcd of X(m,n).',
    '5. [equation; elementary] A1 = -1 gives H(m,n) = g^2.  Non-trivial classes: H = F and g | 6, so F(m,n) = g^2 in '
    '{1,4,9,36};',
    '   trivial class: H = (2/3) F0 and g | 2, so F0(m,n) = 6 with g = 2.  Then X = X4/g^2 and Y = +-Y6/g^3.',
    '6. [converse; elementary] each coprime (m,n) with H = g^2 gives an integral point of the class; in the three '
    'non-trivial classes u^2 = F (mod 6),',
    '   so EVERY coprime solution of F = c in {1,4,9,36} has g = sqrt(c) and gives an integral point.  The map '
    '+-(m,n) -> {(X,+-Y)} is a bijection',
    '   (xi is fixed up to sign by X, and P^1(Q) -> C(Q) is bijective by the inverse identity).',
    'FACTS USED:',
    '   external, non-elementary: rank E(Q) = 2 for 1611.a1 (LMFDB; only rank <= 2 is used), with the Mordell-Weil '
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
out('(ii) THE TRIVIAL CLASS (0,0): F0(m,n) = n (m^3 + 3 m^2 n - m n^2 + 3 n^3) = 6, solved completely')
out('=' * 110)
R0 = RED['00']
Cub = m ** 3 + 3 * m ** 2 * n - m * n ** 2 + 3 * n ** 3
check(sp.expand(R0['F'] - n * Cub) == 0, 'F0 = n * cubic (file form)')
out('F0 = %s = n * (%s)  (= the file\'s F of class (0, 0))' % (fstr(R0['F']), fstr(Cub)))
out('n divides F0(m,n) = 6, so n in {+-1, +-2, +-3, +-6}; for each n, m is an integer root of the monic cubic')
out('   m^3 + 3n m^2 - n^2 m + (3n^3 - 6/n) = 0, so (constant term != 0) m divides 3n^3 - 6/n (rational-root test):')
SOL0 = []
for nn in [1, -1, 2, -2, 3, -3, 6, -6]:
    c0 = 3 * nn ** 3 - 6 // nn
    check(6 % nn == 0 and c0 != 0, 'constant term non-zero')
    roots = sorted(d for d0 in sp.divisors(abs(c0)) for d in (d0, -d0)
                   if d ** 3 + 3 * nn * d ** 2 - nn ** 2 * d + c0 == 0)
    for r in roots:
        check(int((n * Cub).subs({m: r, n: nn})) == 6, 'F0 = 6')
        SOL0.append((r, nn))
    out('   n = %2d: m^3 + (%d) m^2 + (%d) m + (%d) = 0;  %d divisors tried;  integer roots m = %s'
        % (nn, 3 * nn, -nn ** 2, c0, 2 * len(sp.divisors(abs(c0))), roots))
SOL0.sort()
file0 = plist(grab(r'coprime solutions of F\(m,n\) in \[6\]: (\[.*\])$', SEC, 'trivial class solutions').group(1))
check(SOL0 == sorted(file0), 'solutions agree with the file')
check(all(gcd(a, b) == 1 for a, b in SOL0), 'all solutions coprime')
out('all integer solutions (automatically coprime: d^4 | 6 forces d = 1): %s  (= file)' % SOL0)
out('map with u = m^2+2mn+5n^2, v = -2mn-2n^2, w = -2n^2 (file), g = gcd(u,v,w); an integral point needs H = (2/3)*6 '
    '= 4 = g^2, i.e. g = 2:')
out('   %-10s %-16s %-3s %-22s %s' % ('(m,n)', '(u,v,w)', 'g', 'X = X4/4', 'Y = +-Y6/8'))
PTS = {}          # X -> dict(Y=|Y|, src=[...])
NOT_INTEGRAL = []
for (a, b) in SOL0:
    uu, vv, ww = [int(e.subs({m: a, n: b})) for e in R0['Xp']]
    g = gcd3(uu, vv, ww)
    X4v, Y6v = int(R0['X4'].subs({m: a, n: b})), int(R0['Y6'].subs({m: a, n: b}))
    Hv = int(R0['H'].subs({m: a, n: b}))
    check(Hv == 4, 'H = 4')
    Xq, Yq = Fr(X4v, 4), Fr(Y6v, 8)
    check(onE(Xq, Yq), 'point on E')
    if g == 2:
        check(Xq.denominator == 1 and Yq.denominator == 1, 'integral point')
        X, Y = int(Xq), abs(int(Yq))
        xi = [Fr(uu, g) * R0['bas'][0][k] + Fr(vv, g) * R0['bas'][1][k] + Fr(ww, g) * R0['bas'][2][k]
              for k in range(3)]
        check(kmul(xi, xi) == [X, -1, 0], 'X - th = xi^2 (trivial class)')
        PTS.setdefault(X, dict(Y=Y, src=[]))
        check(PTS[X]['Y'] == Y, 'same |Y|')
        PTS[X]['src'].append(('00', 6, (a, b)))
        hits = [r for r in R0['ftab'] if r[0] == X]
        check(len(hits) == 2 and all(r[2] in ((a, b), (-a, -b)) and r[3] == g and r[4] == 6 for r in hits),
              'file table of class (0, 0) agrees for X = %d' % X)
        note = 'integral point (%d, +-%d); X - th = xi^2 checked; = file table' % (X, Y)
    else:
        NOT_INTEGRAL.append(((a, b), g, Xq, Yq))
        note = 'g = %d != 2: H = 4 != g^2, NOT from an integral point; rational point (%s, +-%s) [non-integral]' % (
            g, Xq, abs(Yq))
    out('   %-10s %-16s %-3d %-22s %s' % ((a, b), (uu, vv, ww), g, Xq, note))
check(sorted(PTS) == [-8, 4, 40], 'trivial class gives X in {-8, 4, 40}')
check([s_[0] for s_ in NOT_INTEGRAL] == [(-7, 2), (7, -2)], 'only +-(7,-2) fails')
out('=> trivial class: X in {-8, 4, 40} (from +-(1,1), +-(1,-1), +-(3,-1)); the solution +-(7,-2) of F0 = 6 has g = 1 '
    'and gives no integral point.')
mt = grab(r'class \(0, 0\) \(reducible F = .*\n.*\n\s+integral points with phi\(P\) in this class: \(X,\|Y\|\) in '
          r'(\[.*?\])', SEC, 'trivial class points')
check(sorted((X, PTS[X]['Y']) for X in PTS) == plist(mt.group(1)), 'trivial class points agree with the file')

# ==========================================================================================================
out()
out('=' * 110)
out('(iii) ALL INTEGRAL POINTS, assuming the coprime solutions of F(m,n) in {1,4,9,36} are +-core.KNOWN')
out('=' * 110)
out('%-6s %-4s %-10s %-20s %-3s %-6s %-9s %s' % ('class', 'c', '+-(m,n)', '(u,v,w)', 'g', 'X', '|Y|',
                                                 'checks'))
for key in ('10', '01', '11'):
    Rk = RED[key]
    uvw_f, X4_f, Y6_f = descent_maps(key)
    for (a, b) in KNOWN[key]:
        rows = []
        for sg in (1, -1):
            mm, nn = sg * a, sg * b
            check(gcd(mm, nn) == 1, 'coprime')
            c = Fval(key, mm, nn)
            check(c in RHS, 'F in {1,4,9,36}')
            uu, vv, ww = uvw_f(mm, nn)
            g = gcd3(uu, vv, ww)
            check(c == g * g, 'F = g^2')
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
        out('%-6s %-4d %-10s %-20s %-3d %-6d %-9d F = g^2; X - th = mu xi^2; Y = %d N(xi); (-m,-n) the same'
            % (key, c, (a, b), uvw_, g, X, Y, Rk['s']))
        # agreement with the file's table of the class
        hits = [r for r in RED[key]['ftab'] if r[0] == X]
        check(len(hits) == 2 and all(r[2] in ((a, b), (-a, -b)) and r[3] == g and r[4] == c for r in hits),
              'file table of the class agrees for X = %d' % X)
ALLX = sorted(PTS)
EXPECT = [-15, -12, -8, 0, 4, 9, 12, 24, 40, 48, 108, 252, 420, 4953]
check(ALLX == EXPECT, 'the 14 X values')
POINTS28 = sorted([(X, PTS[X]['Y']) for X in ALLX] + [(X, -PTS[X]['Y']) for X in ALLX])
check(len(POINTS28) == 28 and all(PTS[X]['Y'] > 0 for X in ALLX), '28 points')
out()
out('trivial class (0,0) (part ii): X = -8, 4, 40 from +-(1,1), +-(1,-1), +-(3,-1) with c = F0 = 6, g = 2.')
# agreement with the file's table of integral points (and their classes from the Mordell-Weil decomposition)
seg = SEC[SEC.index('integral points (monic coordinates), decomposition'):SEC.index('grouped by class (X values')]
ftab = [(int(a), int(b), (int(c), int(d)), (int(e), int(f))) for a, b, c, d, e, f in
        re.findall(r'^\s+(-?\d+)\s+(-?\d+)\s+\((-?\d+), (-?\d+)\)\s+\((\d), (\d)\)\s*$', seg, re.M)]
check(sorted((a, b) for a, b, _, _ in ftab) == POINTS28, 'the 28 points = the file\'s list')
keycl = {'00': (0, 0), '10': (1, 0), '01': (0, 1), '11': (1, 1)}
check(all(keycl[PTS[a]['src'][0][0]] == cl for a, b, ij, cl in ftab), 'classes agree with the MW decomposition')
IJ = {(a, b): ij for a, b, ij, cl in ftab}
check(all(IJ[(X, -Y)] == (-IJ[(X, Y)][0], -IJ[(X, Y)][1]) for X, Y in POINTS28 if Y > 0), '(X,-Y) = -(X,Y)')
out('X values: %s; each with both signs of Y: %d integral points (X, Y):' % (ALLX, len(POINTS28)))
out('   ((m,n) and (-m,-n) give the same (u,v,w), X4, Y6 (forms of even degree); the two signs of Y come from +-xi)')
for X in ALLX:
    src = PTS[X]['src'][0]
    out('   X = %-5d Y = +-%-7d  from class %s, c = %-2d (m,n) = +-%-9s  file: (X,+|Y|) = %s G1 + %s G2'
        % (X, PTS[X]['Y'], src[0], src[1], src[2], IJ[(X, PTS[X]['Y'])][0], IJ[(X, PTS[X]['Y'])][1]))
out('agreement with the file: the same 28 points, and the class of each = (i mod 2, j mod 2) of its Mordell-Weil '
    'decomposition i G1 + j G2: True')
out()
out('Solutions of F = c that do NOT give integral points:')
out('   non-trivial classes: none (every listed solution has g = sqrt(c); by u^2 = F (mod 6) this holds for every '
    'coprime solution).')
for ((a, b), g, Xq, Yq) in NOT_INTEGRAL:
    out('   trivial class: (m,n) = %s, F0 = 6 but g = %d (H = 4 != 1): rational point (%s, +-%s), not integral'
        % ((a, b), g, Xq, abs(Yq)))
# consistency checks (not proofs): Thue box and integral-point search
out()
out('consistency checks (finite searches, not part of the argument):')
for key in ('10', '01', '11'):
    found = set()
    for mm in range(-300, 301):
        for nn in range(0, 301):
            if (nn == 0 and mm <= 0) or gcd(mm, nn) != 1:
                continue
            if Fval(key, mm, nn) in RHS:
                found.add((mm, nn))
    knownset = set()
    for (a, b) in KNOWN[key]:
        knownset.add((a, b) if (b > 0 or (b == 0 and a > 0)) else (-a, -b))
    check(found == knownset, 'box search = KNOWN for %s' % key)
    out('   class %s: coprime (m,n), |m|,|n| <= 300, F(m,n) in {1,4,9,36}: exactly +-%s (= core.KNOWN)'
        % (key, sorted(found)))
XMAX = 10 ** 6
srch = []
for X in range(-15, XMAX + 1):
    D = X ** 3 + GA * X + GB
    if D >= 0:
        r = isqrt(D)
        if r * r == D:
            srch.append(X)
check(-16 ** 3 + GA * (-16) + GB < 0 and srch == EXPECT, 'search agrees')
out('   integral points of E with X <= 10^6 (none with X <= -16: g is increasing there and g(-16) = -496 < 0): '
    'X in %s (= the 14 above)' % srch)

# ==========================================================================================================
out()
out('=' * 110)
out('(iv) TRANSFER: LMFDB model, C(n,2) = C(m,3), Avanesov\'s 3y(y+1) = x(x+1)(x+2)')
out('=' * 110)
out('--- LMFDB model y^2 + y = x^3 - 9x + 20 (X = 4x, Y = 8y + 4)')
kk = sp.symbols('k')
check(sp.expand(((4 * kk) ** 3 - 144 * 4 * kk + 1296) - 16 * (4 * kk ** 3 - 36 * kk + 81)) == 0, 'X = 4k')
out('the integral points (x,y) correspond to the integral points (X,Y) of E with 4 | X and Y = 4 (mod 8); for X = 4k, '
    'Y^2 = 16 (4k^3 - 36k + 81)')
out('   with 4k^3 - 36k + 81 odd, so Y = 4 * odd = 4 (mod 8): the condition 4 | X suffices and both signs of Y are kept')
LPTS = []
for (X, Y) in POINTS28:
    if X % 4 == 0:
        check(Y % 8 == 4, 'Y = 4 mod 8')
        x, y = X // 4, (Y - 4) // 8
        check(y * y + y == x ** 3 - 9 * x + 20, 'on the LMFDB model')
        LPTS.append((x, y))
LPTS.sort()
check(len(LPTS) == 22, '22 points')
out('X not divisible by 4: %s (no LMFDB point).  The 22 integral points (x, y):' % [X for X in ALLX if X % 4])
out('   ' + ', '.join('(%d,%d)' % p for p in LPTS))
with open(LMFDB_MD, encoding='utf-8') as fh:
    md = fh.read().replace(chr(0x2212), '-')            # LMFDB text uses the Unicode minus sign
line = [ln for ln in md.splitlines() if ln.startswith('| 1611.a1 |') and 'The page lists the points' in ln]
check(len(line) == 1, 'LMFDB line found')
lm = sorted((int(a), int(b)) for a, b in re.findall(r'\((-?\d+),(-?\d+)\)', line[0].split('The page lists the points')[1]))
check(lm == LPTS, 'equal to the LMFDB list')
out('   = the 22 points listed by LMFDB (lmfdb_data.md): True')
out()
out('--- C(n,2) = C(m,3) via x = 3(m-1), y = 9n - 5')
M_, N_ = sp.symbols('M N')
idt = sp.expand(((9 * N_ - 5) ** 2 + (9 * N_ - 5)) - ((3 * (M_ - 1)) ** 3 - 9 * 3 * (M_ - 1) + 20)
                - 162 * (N_ * (N_ - 1) / 2 - M_ * (M_ - 1) * (M_ - 2) / 6))
check(idt == 0, 'correspondence identity')
out('identity: [y^2 + y - (x^3 - 9x + 20)] at x = 3(m-1), y = 9n-5  =  162 (C(n,2) - C(m,3)),  C(n,2) = n(n-1)/2, '
    'C(m,3) = m(m-1)(m-2)/6  [exact]')
res27 = sorted(r for r in range(27) if (r * r + r - 20) % 27 == 0)
check(res27 == [4, 13, 22] and all(r % 9 == 4 for r in res27), 'y = 4 mod 9')
check(all((((3 * t) ** 3 - 9 * 3 * t + 20) - 20) % 27 == 0 for t in range(27)), 'x = 0 mod 3 => rhs = 20 mod 27')
out('so integer solutions (m,n) <-> integral points with 3 | x and 9 | y + 5; and 3 | x already forces y^2 + y = 20 '
    '(mod 27), i.e. y = 4 (mod 9)')
out('   (residues of y mod 27 with y^2 + y = 20: %s).  Integral points with 3 | x -> (m, n) = (x/3 + 1, (y+5)/9):' % res27)
MN = []
for (x, y) in LPTS:
    if x % 3 == 0:
        check((y + 5) % 9 == 0, 'y = 4 mod 9')
        mm, nn = x // 3 + 1, (y + 5) // 9
        Cm3, Cn2 = mm * (mm - 1) * (mm - 2) // 6, nn * (nn - 1) // 2
        check(Cm3 == Cn2, 'C(m,3) = C(n,2)')
        MN.append((mm, nn, Cm3))
MN.sort()
check(len(MN) == 16, '16 pairs')
vals = sorted(set(c for _, _, c in MN))
check(vals == [0, 1, 10, 120, 1540, 7140], 'values')
for val in vals:
    out('   C(m,3) = C(n,2) = %-5d  (m, n) in %s' % (val, [(a, b) for a, b, c in MN if c == val]))
out('   the involution y -> -1 - y of the curve is n -> 1 - n (C(1-n,2) = C(n,2)); m < 0 is impossible '
    '(C(m,3) < 0 <= C(n,2)).')
out('   values: {0, 1, 10, 120, 1540, 7140}; the numbers > 0 both triangular and tetrahedral: 1, 10, 120, 1540, 7140.')
for val in vals[1:]:
    (mm, nn) = max((a, b) for a, b, c in MN if c == val)
    j_, k_ = mm - 2, nn - 1
    check(k_ * (k_ + 1) // 2 == val == j_ * (j_ + 1) * (j_ + 2) // 6, 'T_k = Te_j')
    out('      %-5d = T_%d = k(k+1)/2 with k = n-1 = %d,  = Te_%d = j(j+1)(j+2)/6 with j = m-2 = %d'
        % (val, k_, k_, j_, j_))
srch_m = []
for mm in range(0, 10 ** 6 + 1):
    T = mm * (mm - 1) * (mm - 2) // 6
    r = isqrt(8 * T + 1)
    if r * r == 8 * T + 1:
        srch_m.append(mm)
check(srch_m == sorted(set(a for a, _, _ in MN)), 'search over m agrees')
out('   consistency (search, not proof): 0 <= m <= 10^6 with C(m,3) triangular: m in %s' % srch_m)
out()
out('--- Avanesov: 3y(y+1) = x(x+1)(x+2), model 6Y^2 = X^3 - 4X + 6 (X = 2x + 2, Y = 2y + 1)')
xa, ya, Xa, Ya = sp.symbols('xa ya Xa Ya')
check(sp.expand(6 * (2 * ya + 1) ** 2 - ((2 * xa + 2) ** 3 - 4 * (2 * xa + 2) + 6)
                - 8 * (3 * ya * (ya + 1) - xa * (xa + 1) * (xa + 2))) == 0, 'Avanesov model identity')
check(sp.expand((36 * Ya) ** 2 - ((6 * Xa) ** 3 - 144 * (6 * Xa) + 1296) - 216 * (6 * Ya ** 2 - (Xa ** 3 - 4 * Xa + 6)))
      == 0, 'monic = 216 * Avanesov')
check(sp.expand((ya + 1) * ya / 2 - (xa + 2) * (xa + 1) * xa / 6
                - (3 * ya * (ya + 1) - xa * (xa + 1) * (xa + 2)) / 6) == 0, 'C(y+1,2) - C(x+2,3)')
out('identities: 6(2y+1)^2 - ((2x+2)^3 - 4(2x+2) + 6) = 8 (3y(y+1) - x(x+1)(x+2));  (36 Y_A)^2 - ((6 X_A)^3 - '
    '144 (6 X_A) + 1296) = 216 (6 Y_A^2 - X_A^3 + 4 X_A - 6)')
out('   so X\' = 6 X_A = 12(x+1), Y\' = 36 Y_A = 36(2y+1) on E; and C(y+1,2) - C(x+2,3) = (3y(y+1) - x(x+1)(x+2))/6, '
    'i.e. (m,n) = (x+2, y+1)')
check(all((4 * t ** 3 - 4 * t + 3) % 3 == 0 and ((4 * t ** 3 - 4 * t + 3) // 3) % 2 == 1 for t in range(6)),
      'X\' = 12t => Y\'^2 = 1296 * odd')
check(sp.expand((12 * kk) ** 3 - 144 * 12 * kk + 1296 - 1296 * (4 * kk ** 3 - 4 * kk + 3) / 3) == 0, 'X = 12k')
out('   for X\' = 12k: Y\'^2 = 1296 (4k^3 - 4k + 3)/3 with (4k^3 - 4k + 3)/3 an odd integer, so 12 | X\' already gives '
    'Y\' = 36 * odd')
AV = []
for (X, Y) in POINTS28:
    if X % 12 == 0:
        check(Y % 36 == 0 and (Y // 36) % 2 == 1, 'Y = 36 * odd')
        x, y = X // 12 - 1, (Y // 36 - 1) // 2
        check(3 * y * (y + 1) == x * (x + 1) * (x + 2), 'Avanesov equation')
        AV.append((x, y))
AV.sort()
xsA = sorted(set(x for x, _ in AV))
check(xsA == [-2, -1, 0, 1, 3, 8, 20, 34], 'Avanesov x values')
check(sorted((x + 2, y + 1) for x, y in AV) == sorted((a, b) for a, b, _ in MN), 'same (m,n) as via LMFDB')
out('X\' divisible by 12 among the 14: %s; X_A = X\'/6 in %s; Avanesov solutions (x, y):'
    % ([X for X in ALLX if X % 12 == 0], [X // 6 for X in ALLX if X % 12 == 0]))
out('   ' + ', '.join('(%d,%d)' % p for p in AV))
out('   x in %s; (x+2, y+1) = the 16 pairs (m,n) of C(n,2) = C(m,3) above: True' % xsA)
out('   (so for x, y >= 0, Avanesov\'s x is the tetrahedral index j and y the triangular index k: Te_x = T_y)')

# ==========================================================================================================
out()
out('=' * 110)
out('CONCLUSION (conditional)')
out('=' * 110)
out('Assume (A) rank E(Q) = 2 for 1611.a1 [LMFDB/mwrank; only rank <= 2 is used] and (B) for each class (1,0), (0,1), '
    '(1,1)')
out('the coprime solutions of F(m,n) in {1,4,9,36} are exactly +-core.KNOWN.  Then, using also N1, N2 (index-'
    'discriminant formula, integral closure):')
out('   the integral points of Y^2 = X^3 - 144X + 1296 are exactly the 28 points X in %s, both signs of Y;' % ALLX)
out('   the integral points of y^2 + y = x^3 - 9x + 20 are the 22 points listed by LMFDB;')
out('   C(n,2) = C(m,3) holds exactly for the 16 pairs above, with values 0, 1, 10, 120, 1540, 7140;')
out('   3y(y+1) = x(x+1)(x+2) holds exactly for x in %s (16 pairs (x,y)).' % xsA)
out('Solving the trivial-class equation F0 = 6 (part ii) needs neither (A) nor (B) and is elementary.  Every class map, '
    'value and identity above is exact.')
out()
out('all %d checks passed  [done in %.1f s]' % (NCHECK[0], time.time() - T0))
