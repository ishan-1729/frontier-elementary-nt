"""verify_bruteforce.py -- independent sanity checks of the exclusion step 02_exclusion.py (brute force; NOT a proof).

Independent verifier (Claude subagent, Opus 5.5), 2026-10-15.  From core.py only DATA is used:
CLASSES[key]['F'], CLASSES[key]['BASIS'], CLASSES[key]['label'], UNITS, KNOWN (and core.Field.mt once, only to compare
it with the structure constants computed here).  All arithmetic is re-implemented here:
structure constants by sympy polynomial remainders, norms and inverses through the companion matrix of f (sympy,
exact), then plain Python integers.  The seeds are read from 01_coverage.out.txt and from the SEEDS block of
02_exclusion.py; the exponents, moduli and Proposition 17 data to compare with are read from 02_exclusion.out.txt.

Part 1  The orders: structure constants (compared with core.Field.mt), the forms of the task statement,
        N(n - m alpha) = F(m,n), N(eps1) = N(eps2) = 1, N(w) = c for the 16 seeds, the marker criterion
        (x2 = x3 = 0 <=> the element lies in Z + Z alpha), number of real roots (Sturm).
Part 2  EXACT brute force: eps1^a eps2^b w for ALL |a|, |b| <= 40 and all 16 orbits; every (a,b) whose two markers
        vanish is listed with its solution (m,n).  They must be exactly the exponents printed in
        02_exclusion.out.txt, and every pair +-(m,n) of core.KNOWN must occur exactly once over all 16 orbits.
Part 3  Modular scan of |a|, |b| <= 1000: both markers mod two primes near 2^30 (numpy), every hit confirmed exactly.
Part 4  Direct search of F(m,n) = c (c in {1,4,9,36}) over all integers n and |m| <= 10^5.  Complete in that strip:
        |F(m,n)| >= |n - r1 m| |n - r2 m| (Im alpha_3)^2 m^2, so a solution has n within 6/(|Im alpha_3| |m|) (+1 for
        rounding) of r1 m or r2 m (r1, r2 the real roots).  Compared with core.KNOWN.
Part 5  Independent re-check of the exclusion certificate with the moduli printed by 02_exclusion: for each q, the
        Hermite basis (h1,0), (-a1,b0) of Lambda_q (both vectors in Lambda_q, and h1*b0 = |G_q| with G_q enumerated),
        the candidate classes, Proposition 17 data (lambda, nu(J)) at each known zero; then L = cap Lambda_q (own CRT
        intersection) and EVERY class of Z^2/L is enumerated: it must be killed by some q (a marker != 0 mod q) or lie
        in the Lambda_q-class of a known zero at which Proposition 17 holds.  Also -1 not in <eps1,eps2> (some odd q
        with -1 not in G_q), and no other known zero lies in a Lambda_q-class where Proposition 17 holds.
Part 6  Adversarial test of Proposition 17 itself: planted coset systems with two zeros in one class (random seeds,
        markers = p-saturated integer rows orthogonal to both zeros); Proposition 17 must fail on every one.

A failed check prints CHECK FAILED and the script exits with status 1.
Usage: python verify_bruteforce.py [SCRIPTS_DIR]   (default: the folder of this file)
"""
import sys, os, re, ast, time, math
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE
sys.path.insert(0, SCRIPTS)
import core                      # DATA ONLY: CLASSES (F, BASIS, label), UNITS, KNOWN

T0 = time.time()
FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print('   CHECK FAILED: ' + msg)
    return bool(cond)


RHS = (1, 4, 9, 36)
# The forms as stated in the task (F(m,n) = a0 m^4 + a1 m^3 n + a2 m^2 n^2 + a3 m n^3 + a4 n^4).
TASK_FORMS = {'10': (36, 0, 0, -6, 1), '01': (-5, 14, 24, 2, 1), '11': (-8, 4, 24, 4, 1)}

# ======================================================================================================================
# Independent arithmetic in the order O
# ======================================================================================================================
tt = sp.symbols('t')


def _rat(c):
    if isinstance(c, Fr):
        return sp.Rational(c.numerator, c.denominator)
    return sp.Rational(c)


class Order:
    """O = Z omega_0 + ... + Z omega_3, omega_i = sum_j BASIS[i][j] alpha^j, f(alpha) = 0, f(t) = F(1,t)."""

    def __init__(self, key):
        self.key = key
        F = tuple(int(v) for v in core.CLASSES[key]['F'])
        assert F[4] == 1
        self.F = F
        a0, a1, a2, a3, _ = F
        self.f = sp.Poly(tt**4 + a3 * tt**3 + a2 * tt**2 + a1 * tt + a0, tt, domain='QQ')
        self.B = sp.Matrix([[_rat(c) for c in row] for row in core.CLASSES[key]['BASIS']])
        self.Binv = self.B.inv()
        Cm = sp.zeros(4, 4)                         # companion matrix: multiplication by alpha on 1, alpha, ...
        for j in range(3):
            Cm[j + 1, j] = 1
        for i, ai in enumerate((a0, a1, a2, a3)):
            Cm[i, 3] = -ai
        self.Comp = Cm
        om = [self._poly([self.B[i, j] for j in range(4)]) for i in range(4)]
        self.mt = [[None] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(4):
                r = (om[i] * om[j]).rem(self.f)
                co = self.from_power([r.coeff_monomial(tt**k) for k in range(4)])
                assert all(v.q == 1 for v in co), 'BASIS is not closed under multiplication'
                self.mt[i][j] = tuple(int(v) for v in co)

    @staticmethod
    def _poly(pc):
        return sp.Poly(sum(sp.Rational(pc[j]) * tt**j for j in range(4)), tt, domain='QQ')

    def to_power(self, x):
        return [sum(sp.Integer(x[i]) * self.B[i, j] for i in range(4)) for j in range(4)]

    def from_power(self, pc):
        return list(sp.Matrix([list(pc)]) * self.Binv)

    def mat_power(self, x):
        pc = self.to_power(x)
        M = sp.zeros(4, 4)
        P = sp.eye(4)
        for j in range(4):
            M += pc[j] * P
            P = P * self.Comp
        return M

    def norm(self, x):
        return self.mat_power(x).det()

    def inverse(self, x):
        pc = self.mat_power(x).inv() * sp.Matrix([1, 0, 0, 0])
        co = self.from_power(list(pc))
        assert all(v.q == 1 for v in co), 'inverse not in O'
        return tuple(int(v) for v in co)

    def mul(self, x, y):
        z0 = z1 = z2 = z3 = 0
        mt = self.mt
        for i in range(4):
            xi = x[i]
            if xi:
                row = mt[i]
                for j in range(4):
                    yj = y[j]
                    if yj:
                        c = xi * yj
                        r = row[j]
                        z0 += c * r[0]
                        z1 += c * r[1]
                        z2 += c * r[2]
                        z3 += c * r[3]
        return (z0, z1, z2, z3)

    def mulmod(self, x, y, M):
        z = self.mul(x, y)
        return (z[0] % M, z[1] % M, z[2] % M, z[3] % M)

    def pow(self, x, e, xinv):
        if e < 0:
            x, e = xinv, -e
        r, b = (1, 0, 0, 0), tuple(x)
        while e:
            if e & 1:
                r = self.mul(r, b)
            e >>= 1
            if e:
                b = self.mul(b, b)
        return r

    def powmod(self, x, e, M):
        assert e >= 0
        r, b = (1 % M, 0, 0, 0), tuple(v % M for v in x)
        while e:
            if e & 1:
                r = self.mulmod(r, b, M)
            e >>= 1
            if e:
                b = self.mulmod(b, b, M)
        return r

    def Fval(self, m, n):
        a0, a1, a2, a3, a4 = self.F
        return a0 * m**4 + a1 * m**3 * n + a2 * m**2 * n**2 + a3 * m * n**3 + a4 * n**4


def nu(x, p, cap):
    if x == 0:
        return cap
    v = 0
    while x % p == 0 and v < cap:
        x //= p
        v += 1
    return v


def pair_key(m, n):
    """Canonical representative of the pair +-(m,n)."""
    return max((m, n), (-m, -n))


# ======================================================================================================================
# Reading the files to compare with
# ======================================================================================================================
def read(name):
    with open(os.path.join(SCRIPTS, name), encoding='utf-8') as fh:
        return fh.read()


def seed_block(text, start):
    i = text.index(start)
    j = text.index('\n}', i)
    block = re.sub(r'#[^\n]*', '', text[i + len(start) - 1: j + 2])
    return ast.literal_eval(block)


LABEL2KEY = {core.CLASSES[k]['label']: k for k in core.CLASSES}


def parse_exclusion(text):
    orbits = []
    key = cur = None
    for line in text.splitlines():
        m = re.match(r'CLASS (\(\d,\d\))\s', line)
        if m:
            key = LABEL2KEY[m.group(1)]
            continue
        m = re.search(r'--- orbit c = (\d+), seed w = (\([^)]*\)); known solutions in it: (\[.*?\]) '
                      r'\(exponents (\[.*?\])\)', line)
        if m:
            cur = dict(key=key, c=int(m.group(1)), w=tuple(ast.literal_eval(m.group(2))),
                       known=[tuple(v) for v in ast.literal_eval(m.group(3))],
                       exps=[tuple(v) for v in ast.literal_eval(m.group(4))],
                       used=None, index=None, certs=[], Gq={}, done=None)
            orbits.append(cur)
            continue
        m = re.search(r'moduli used: (\[.*?\]);\s+final lattice index (\d+)', line)
        if m:
            cur['used'] = list(ast.literal_eval(m.group(1)))
            cur['index'] = int(m.group(2))
            continue
        m = re.search(r'class \((-?\d+), (-?\d+)\) contains the solution \((-?\d+), (-?\d+)\): Proposition 17 at '
                      r'q = (\d+) \(lambda = (\d+), nu\(J\) = (\d+), e = (\d+)\)', line)
        if m:
            g = [int(v) for v in m.groups()]
            cur['certs'].append(dict(cls=(g[0], g[1]), mn=(g[2], g[3]), q=g[4], lam=g[5], nuJ=g[6], e=g[7]))
            continue
        m = re.search(r'\bq = (\d+): \|G_q\| = (\d+), candidate classes mod Lambda_q: (\d+)', line)
        if m:
            cur['Gq'][int(m.group(1))] = (int(m.group(2)), int(m.group(3)))
            continue
        m = re.search(r'uncertified live classes left: (\d+)\s+->\s+(.*)$', line)
        if m:
            cur['done'] = (int(m.group(1)) == 0 and 'ORBIT DONE' in m.group(2))
    return orbits


print('verify_bruteforce.py -- independent sanity checks of 02_exclusion.py (brute force; NOT a proof)')
print(f'SCRIPTS = {SCRIPTS}')
print()

COV_OUT = read('01_coverage.out.txt')
EXC_PY = read('02_exclusion.py')
EXC_OUT = read('02_exclusion.out.txt')
SEEDS_COV = seed_block(COV_OUT, 'SEEDS = {')
SEEDS_EXC = seed_block(EXC_PY, 'SEEDS = {')
ORBITS = parse_exclusion(EXC_OUT)

# ======================================================================================================================
print('=' * 118)
print('PART 0. Files to compare with')
check(SEEDS_COV == SEEDS_EXC, 'SEEDS of 02_exclusion.py differ from those printed in 01_coverage.out.txt')
print(f'   SEEDS of 02_exclusion.py == SEEDS printed by 01_coverage.out.txt: {SEEDS_COV == SEEDS_EXC}')
nseeds = sum(len(v) for d in SEEDS_COV.values() for v in d.values())
print(f'   number of seeds: {nseeds};  orbits parsed from 02_exclusion.out.txt: {len(ORBITS)}')
check(nseeds == 16 and len(ORBITS) == 16, 'expected 16 seeds / 16 orbit sections')
order_cov = [(k, c, tuple(w)) for k in SEEDS_COV for c in SEEDS_COV[k] for w in SEEDS_COV[k][c]]
order_exc = [(o['key'], o['c'], o['w']) for o in ORBITS]
check(order_cov == order_exc, 'orbit list of 02_exclusion.out.txt differs from the seeds')
for o in ORBITS:
    check(o['used'] is not None and o['done'] is True, f'orbit {o["key"]} c={o["c"]} w={o["w"]} not reported DONE')
check('ALL ORBITS CERTIFIED' in EXC_OUT, '02_exclusion.out.txt does not end with ALL ORBITS CERTIFIED')
print('   all 16 orbit sections parsed; each reports ORBIT DONE; file says ALL ORBITS CERTIFIED:',
      'ALL ORBITS CERTIFIED' in EXC_OUT)

# ======================================================================================================================
print()
print('=' * 118)
print('PART 1. The orders (independent arithmetic)')
ORD, UNITS, UINV = {}, {}, {}
for key in ('10', '01', '11'):
    O = Order(key)
    ORD[key] = O
    check(O.F == TASK_FORMS[key], f'form of class {key} differs from the task statement')
    same_mt = (O.mt == [[tuple(v) for v in row] for row in core.Field(key).mt])
    check(same_mt, f'structure constants of {key} differ from core.Field')
    nreal = O.f.count_roots()          # real roots, exact (Sturm)
    B = O.B
    tri = (list(B.row(0)) == [1, 0, 0, 0] and list(B.row(1)) == [0, 1, 0, 0]
           and all(B[i, 2] == 0 for i in (0, 1, 3)) and all(B[i, 3] == 0 for i in (0, 1, 2))
           and B[2, 2] != 0 and B[3, 3] != 0)
    check(tri, f'basis of {key}: marker criterion does not apply')
    # N(n - m alpha) = F(m,n)
    nm_ok = all(O.norm((n, -m, 0, 0)) == O.Fval(m, n) for m in range(-3, 4) for n in range(-3, 4))
    check(nm_ok, f'N(n - m alpha) != F(m,n) for {key}')
    e1, e2 = (tuple(u) for u in core.UNITS[key])
    n1, n2 = O.norm(e1), O.norm(e2)
    check(n1 == 1 and n2 == 1, f'unit norms of {key}: {n1}, {n2}')
    e1i, e2i = O.inverse(e1), O.inverse(e2)
    check(O.mul(e1, e1i) == (1, 0, 0, 0) and O.mul(e2, e2i) == (1, 0, 0, 0), f'inverses of {key}')
    UNITS[key], UINV[key] = (e1, e2), (e1i, e2i)
    lab = core.CLASSES[key]['label']
    print(f'   class {lab}: F = {O.F} (= task statement: {O.F == TASK_FORMS[key]}); real roots of f (Sturm): {nreal};'
          f' structure constants == core.Field.mt: {same_mt}')
    print(f'      marker criterion (omega_0 = 1, omega_1 = alpha, alpha^2/alpha^3 only in omega_2/omega_3): {tri};'
          f'  N(n - m alpha) = F(m,n) on [-3,3]^2: {nm_ok}')
    print(f'      N(eps1) = {n1}, N(eps2) = {n2};  eps1^-1 = {e1i},  eps2^-1 = {e2i}')
    for c, ws in SEEDS_COV[key].items():
        for w in ws:
            nw = O.norm(w)
            check(nw == c, f'N(seed {w}) = {nw} != {c} in class {key}')
            print(f'      seed w = {str(tuple(w)):20s} N(w) = {nw}  (c = {c})')
    for (m, n) in core.KNOWN[key]:
        check(O.Fval(m, n) in RHS, f'KNOWN {(m, n)} of {key} has F = {O.Fval(m, n)}')

# ======================================================================================================================
print()
print('=' * 118)
R = 40
print(f'PART 2. EXACT brute force: eps1^a eps2^b w for all |a|, |b| <= {R} (big integers); zeros of the two markers')
FOUND = {}                 # (key, c, w) -> list of ((a,b), (m,n), element)
all_pairs = {}             # (class key, pair_key) -> list of orbits
t2 = time.time()
for key in ('10', '01', '11'):
    O = ORD[key]
    e1, e2 = UNITS[key]
    e1i, e2i = UINV[key]
    E1 = {0: (1, 0, 0, 0)}
    for a in range(1, R + 1):
        E1[a] = O.mul(E1[a - 1], e1)
        E1[-a] = O.mul(E1[-a + 1], e1i)
    for o in [o for o in ORBITS if o['key'] == key]:
        c, w = o['c'], o['w']
        W = {0: tuple(w)}
        for b in range(1, R + 1):
            W[b] = O.mul(W[b - 1], e2)
            W[-b] = O.mul(W[-b + 1], e2i)
        zeros = []
        for a in range(-R, R + 1):
            x = E1[a]
            for b in range(-R, R + 1):
                z = O.mul(x, W[b])
                if z[2] == 0 and z[3] == 0:
                    zeros.append(((a, b), z))
        rows = []
        for (ab, z) in zeros:
            pc = O.to_power(z)
            okz = (pc[2] == 0 and pc[3] == 0 and all(v.q == 1 for v in pc))
            n, m = int(pc[0]), -int(pc[1])
            Fv = O.Fval(m, n)
            check(okz and Fv == c and (n, -m) == (z[0], z[1]), f'zero {ab} of orbit {key} c={c} w={w}: F = {Fv}')
            rows.append((ab, (m, n), z))
            all_pairs.setdefault((key, pair_key(m, n)), []).append((c, w, ab))
        FOUND[(key, c, w)] = rows
        got = sorted(r[0] for r in rows)
        exp = sorted(o['exps'])
        same = (got == exp)
        check(same, f'orbit {key} c={c} w={w}: zeros {got} != exponents {exp} printed by 02_exclusion')
        known_here = sorted(pair_key(*mn) for mn in o['known'])
        sols_here = sorted(pair_key(*r[1]) for r in rows)
        check(known_here == sols_here, f'orbit {key} c={c} w={w}: solutions {sols_here} != known {known_here}')
        rep = lambda mn: mn if mn in core.KNOWN[key] else (-mn[0], -mn[1])
        sol_txt = ', '.join(f'(a,b) = {r[0]} -> element n - m alpha with (m,n) = {r[1]}, pair +-{rep(r[1])}, '
                            f'F = {O.Fval(*r[1])}' for r in rows) if rows else 'none'
        print(f'   class {core.CLASSES[key]["label"]} c = {c:2d} seed {str(w):20s}: zeros in the box: {sol_txt}')
        print(f'      02_exclusion.out.txt exponents {exp} -> identical: {same}')
keyed_found = sorted(kp for kp, lst in all_pairs.items() for _ in lst)
keyed_known = sorted((k, pair_key(m, n)) for k in core.KNOWN for (m, n) in core.KNOWN[k])
check(keyed_found == keyed_known, f'solution pairs found {keyed_found} != core.KNOWN {keyed_known}')
dup = {pk: lst for pk, lst in all_pairs.items() if len(lst) > 1}
check(not dup, f'a solution pair occurs more than once: {dup}')
print(f'   every pair of core.KNOWN occurs exactly once over the 16 orbits, and nothing else occurs: '
      f'{keyed_found == keyed_known and not dup}   [{time.time() - t2:.1f} s]')

# ======================================================================================================================
print()
print('=' * 118)
R2 = 1000
P1 = int(sp.prevprime(2**30))
P2 = int(sp.prevprime(P1))
print(f'PART 3. Modular scan |a|, |b| <= {R2}: both markers mod P1 = {P1} and P2 = {P2}; hits confirmed exactly')
t3 = time.time()


def modscan(O, e1, e1i, e2, e2i, w, R, P):
    n = 2 * R + 1

    def table(start, u, ui):
        T = np.zeros((n, 4), dtype=np.int64)
        x = tuple(v % P for v in start)
        T[R] = x
        for k in range(1, R + 1):
            x = O.mulmod(x, u, P)
            T[R + k] = x
        x = tuple(v % P for v in start)
        for k in range(1, R + 1):
            x = O.mulmod(x, ui, P)
            T[R - k] = x
        return T

    E = table((1, 0, 0, 0), e1, e1i)
    Wt = table(w, e2, e2i)
    mask = None
    for k in (2, 3):
        MT = np.array([[O.mt[i][j][k] % P for j in range(4)] for i in range(4)], dtype=np.int64)
        Rk = (E @ MT) % P                       # entries < 4 P^2 < 2^62
        Mk = (Rk @ Wt.T) % P
        mk = (Mk == 0)
        mask = mk if mask is None else (mask & mk)
    return mask


for o in ORBITS:
    key, c, w = o['key'], o['c'], o['w']
    O = ORD[key]
    e1, e2 = UNITS[key]
    e1i, e2i = UINV[key]
    mask = modscan(O, e1, e1i, e2, e2i, w, R2, P1) & modscan(O, e1, e1i, e2, e2i, w, R2, P2)
    hits = [(int(i) - R2, int(j) - R2) for (i, j) in np.argwhere(mask)]
    confirmed = []
    for (a, b) in hits:
        z = O.mul(O.mul(O.pow(e1, a, e1i), O.pow(e2, b, e2i)), w)
        if z[2] == 0 and z[3] == 0:
            confirmed.append((a, b))
    same = sorted(confirmed) == sorted(o['exps']) and sorted(hits) == sorted(confirmed)
    check(same, f'modular scan of orbit {key} c={c} w={w}: hits {hits}, confirmed {confirmed}, expected {o["exps"]}')
    print(f'   class {core.CLASSES[key]["label"]} c = {c:2d} seed {str(w):20s}: hits mod P1 and P2: {sorted(hits)};'
          f' exact zeros: {sorted(confirmed)};  = 02 exponents: {same}')
print(f'   [{time.time() - t3:.1f} s]')

# ======================================================================================================================
print()
print('=' * 118)
MB = 10**5
print(f'PART 4. Direct search: all integer solutions of F(m,n) = c, c in {RHS}, with |m| <= {MB} (every n)')
t4 = time.time()
for key in ('10', '01', '11'):
    O = ORD[key]
    a0, a1, a2, a3, _ = O.F
    rts = np.roots([1, a3, a2, a1, a0])
    real = sorted(float(r.real) for r in rts if abs(r.imag) < 1e-9)
    imag = max(abs(float(r.imag)) for r in rts)
    check(len(real) == 2, f'class {key}: expected 2 real roots')
    found = set()
    for n in range(-3, 4):                        # m = 0: n^4 = c
        if n**4 in RHS:
            found.add(pair_key(0, n))
    for m in range(-MB, MB + 1):
        if m == 0:
            continue
        rad = 6.0 / (imag * abs(m))
        m2 = m * m
        m3 = m2 * m
        m4 = m2 * m2
        for r in real:
            cen = r * m
            for n in range(math.floor(cen - rad) - 1, math.ceil(cen + rad) + 2):
                v = (((n + a3 * m) * n + a2 * m2) * n + a1 * m3) * n + a0 * m4
                if v in (1, 4, 9, 36):
                    found.add(pair_key(m, n))
    known = sorted(pair_key(m, n) for (m, n) in core.KNOWN[key])
    got = sorted(found)
    check(got == known, f'class {key}: direct search {got} != KNOWN {known}')
    print(f'   class {core.CLASSES[key]["label"]}: real roots {[round(r, 6) for r in real]}, |Im alpha_3| = {imag:.6f};'
          f' solutions (one of each +-pair): {[(m, n, O.Fval(m, n)) for (m, n) in got]}')
    print(f'      == core.KNOWN (up to sign): {got == known}')
print(f'   [{time.time() - t4:.1f} s]')

# ======================================================================================================================
print()
print('=' * 118)
print('PART 5. Independent re-check of the exclusion certificate (moduli as printed by 02_exclusion)')
t5 = time.time()
PREC = 30


def spf(q):
    return next(d for d in range(2, q + 1) if q % d == 0)


def kernel_data(O, e1, e2, e1i, q):
    one = (1 % q, 0, 0, 0)
    pw = [one]
    x = tuple(v % q for v in e1)
    while x != one:
        pw.append(x)
        x = O.mulmod(x, e1, q)
    h1 = len(pw)
    idx = {g: i for i, g in enumerate(pw)}
    y = tuple(v % q for v in e2)
    b0 = 1
    while y not in idx:
        y = O.mulmod(y, e2, q)
        b0 += 1
    a1 = idx[y]
    G = {one}
    frontier = [one]
    while frontier:
        nxt = []
        for g in frontier:
            for u in (e1, e2):
                h = O.mulmod(g, u, q)
                if h not in G:
                    G.add(h)
                    nxt.append(h)
        frontier = nxt
    in1 = O.powmod(e1, h1, q) == one
    in2 = O.mulmod(O.powmod(e1i, a1, q), O.powmod(e2, b0, q), q) == one
    minus = ((-1) % q, 0, 0, 0) in G
    return dict(h1=h1, b0=b0, a1=a1, pw=pw, G=len(G), basis_ok=(in1 and in2 and h1 * b0 == len(G)), minus=minus)


def lat_reduce(L, a, b):
    A, B, C = L
    k = b // C
    return ((a - k * B) % A, b - k * C)


def crt2(r1, m1, r2, m2):
    g = math.gcd(m1, m2)
    assert (r2 - r1) % g == 0
    l = m1 // g * m2
    m2g = m2 // g
    k = ((r2 - r1) // g * pow(m1 // g, -1, m2g)) % m2g if m2g > 1 else 0
    return (r1 + m1 * k) % l


def intersect(L1, L2):
    """Hermite form (A,B,C) of L1 cap L2, L_i = {(a,b): C_i | b, A_i | a - (b/C_i) B_i}."""
    A1, B1, C1 = L1
    A2, B2, C2 = L2
    C0 = C1 * C2 // math.gcd(C1, C2)
    g = math.gcd(A1, A2)
    D = (C0 // C1) * B1 - (C0 // C2) * B2
    t0 = g // math.gcd(g, D)
    C = t0 * C0
    A = A1 * A2 // g
    Bv = crt2((C // C1) * B1 % A1, A1, (C // C2) * B2 % A2, A2) % A
    return (A, Bv, C)


for o in ORBITS:
    key, c, w = o['key'], o['c'], o['w']
    O = ORD[key]
    e1, e2 = UNITS[key]
    e1i, e2i = UINV[key]
    zeros = FOUND[(key, c, w)]              # from Part 2: (exponent, (m,n), element)
    lab = core.CLASSES[key]['label']
    print(f'   class {lab} c = {c:2d} seed {w}: moduli used by 02 = {o["used"]}, known zeros {[(z[1], z[0]) for z in zeros]}')
    Lq, cand, p17 = {}, {}, {}
    for q in o['used']:
        p = spf(q)
        e = round(math.log(q, p))
        check(p**e == q and (p != 2 or e >= 2), f'q = {q} is not an admissible prime power')
        kd = kernel_data(O, e1, e2, e1i, q)
        check(kd['basis_ok'], f'Hermite basis of Lambda_{q} wrong for {key}')
        Gq02, nc02 = o['Gq'].get(q, (None, None))
        check(Gq02 == kd['G'], f'|G_{q}| = {kd["G"]} differs from 02 ({Gq02})')
        h1, b0, a1 = kd['h1'], kd['b0'], kd['a1']
        Lq[q] = (h1, (-a1) % h1, b0)
        tab = np.zeros((h1, b0), dtype=bool)
        y = tuple(v % q for v in w)
        for bb in range(b0):
            for aa in range(h1):
                r = O.mulmod(kd['pw'][aa], y, q)
                tab[aa, bb] = (r[2] == 0 and r[3] == 0)
            y = O.mulmod(y, e2, q)
        cand[q] = tab
        check(int(tab.sum()) == nc02, f'q = {q}: {int(tab.sum())} candidate classes, 02 printed {nc02}')
        # Proposition 17 at each known zero: x_k = l_k(eta_1 z)/q, y_k = l_k(eta_2 z)/q  (mod p^PREC)
        Mq = q * p**PREC
        eta1 = O.powmod(e1, h1, Mq)
        eta2 = O.mulmod(O.powmod(e1i, a1, Mq), O.powmod(e2, b0, Mq), Mq)
        iota = 1 if p == 2 else 0
        for (ab, mn, z) in zeros:
            u1, u2 = O.mulmod(eta1, z, Mq), O.mulmod(eta2, z, Mq)
            check(all(u[k] % q == 0 for u in (u1, u2) for k in (2, 3)), f'eta z not == z mod q at {q}')
            x1, x2, y1, y2 = u1[2] // q, u1[3] // q, u2[2] // q, u2[3] // q
            lam = min(nu(v, p, PREC) for v in (x1, x2, y1, y2))
            nJ = nu((x1 * y2 - x2 * y1) % p**PREC, p, PREC)
            ok = lam <= e - 1 - iota and nJ - lam <= e - 1 - iota and nJ < PREC
            p17[(q, pair_key(*mn))] = dict(lam=lam, nuJ=nJ, e=e, ok=ok, ab=ab)
        print(f'      q = {q:3d} (p = {p}, e = {e}): h1 = {h1}, b0 = {b0}, a1 = {a1}, |G_q| = {kd["G"]} = h1*b0: '
              f'{kd["basis_ok"]};  candidate classes {int(tab.sum())} (02: {nc02});  Prop. 17 at the known zeros: '
              + (', '.join(f'{mn}: lambda={d["lam"]}, nu(J)={d["nuJ"]}, ok={d["ok"]}'
                           for (qq, mn), d in p17.items() if qq == q) or 'no known zeros'))
    # consistency of Proposition 17 with the zeros: no OTHER known zero in a Lambda_q-class where it holds
    for (q, pk), d in p17.items():
        if d['ok']:
            zq = lat_reduce(Lq[q], *d['ab'])
            clash = [mn for (ab, mn, z) in zeros if ab != d['ab'] and lat_reduce(Lq[q], *ab) == zq]
            check(not clash, f'q = {q}: Prop. 17 holds at {pk} but the known zero(s) {clash} lie in the same class')
    # Proposition 17 data printed by 02
    for cert in o['certs']:
        d = p17.get((cert['q'], pair_key(*cert['mn'])))
        same = d is not None and (d['lam'], d['nuJ'], d['e']) == (cert['lam'], cert['nuJ'], cert['e']) and d['ok']
        check(same, f'Prop. 17 data of 02 for {cert} not reproduced ({d})')
    # L = cap Lambda_q and all its classes
    L = None
    for q in o['used']:
        L = Lq[q] if L is None else intersect(L, Lq[q])
    A, Bv, C = L
    check(A * C == o['index'], f'index of L = {A * C}, 02 printed {o["index"]}')
    av = np.tile(np.arange(A, dtype=np.int64), C)
    bv = np.repeat(np.arange(C, dtype=np.int64), A)
    killed = np.zeros(A * C, dtype=bool)
    covered = np.zeros(A * C, dtype=bool)
    for q in o['used']:
        Aq, Bq, Cq = Lq[q]
        aq = (av - (bv // Cq) * Bq) % Aq
        bq = bv % Cq
        killed |= ~cand[q][aq, bq]
        for (ab, mn, z) in zeros:
            if p17[(q, pair_key(*mn))]['ok']:
                za, zb = lat_reduce(Lq[q], *ab)
                covered |= (aq == za) & (bq == zb)
    unc = int((~killed & ~covered).sum())
    check(unc == 0, f'{unc} classes of Z^2/L neither killed nor covered by Proposition 17')
    zcls = []
    for (ab, mn, z) in zeros:
        za, zb = lat_reduce(L, *ab)
        i = zb * A + za
        zcls.append((mn, bool(killed[i]), bool(covered[i])))
        check(not killed[i] and covered[i], f'class of the known zero {mn} killed={killed[i]} covered={covered[i]}')
    print(f'      L = cap Lambda_q: Hermite (A,B,C) = {L}, index {A * C} (02: {o["index"]});  classes killed by some q: '
          f'{int(killed.sum())};  inside a Prop.-17 class of a known zero: {int(covered.sum())} '
          f'(of which also killed: {int((killed & covered).sum())});  neither: {unc}')
    for (mn, k_, c_) in zcls:
        print(f'      class of the known zero {mn}: killed = {k_}, inside a Prop.-17 class = {c_}')

# -1 not in <eps1, eps2>
print()
print('   -1 in <eps1,eps2>?  (an odd q with -1 not in G_q shows that -1 is not in <eps1,eps2>, so each pair +-z')
print('   meets <eps1,eps2> w in at most one element)')
for key in ('10', '01', '11'):
    O = ORD[key]
    e1, e2 = UNITS[key]
    e1i, _ = UINV[key]
    wit = None
    for q in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        kd = kernel_data(O, e1, e2, e1i, q)
        if not kd['minus']:
            wit = (q, kd['G'])
            break
    check(wit is not None, f'class {key}: no odd q found with -1 not in G_q')
    print(f'      class {core.CLASSES[key]["label"]}: q = {wit[0]} has -1 not in G_q (|G_q| = {wit[1]})' if wit
          else f'      class {key}: none found')
print(f'   [{time.time() - t5:.1f} s]')

# ======================================================================================================================
print()
print('=' * 118)
print('PART 6. Adversarial test of Proposition 17 itself in the same three orders: planted systems with two zeros')
print('   For a random seed v in O, a random (m0,n0) != 0 and u = eta_1^m0 eta_2^n0 v (eta_j as in 02, exact), take as')
print('   markers a p-saturated basis of {l in Z^4 : l.v = l.u = 0}.  The coset system (eta_1, eta_2; v) then has the')
print('   two zeros (0,0) and (m0,n0), so Proposition 17 must FAIL at v.  margin := max(lambda, nu(J) - lambda)')
print('   - (e - 1 - iota); margin <= 0 would refute Proposition 17, margin = 1 means the case sits at its threshold.')
t6 = time.time()
import random
RNG = random.Random(20261015)
QS6 = [3, 4, 5, 7, 8, 9, 11, 16, 25]
NTR = 30
BIG = 10**6


def p_saturate(rows, p):
    rows = [list(r) for r in rows]
    while True:
        for i in range(2):
            while any(rows[i]) and all(x % p == 0 for x in rows[i]):
                rows[i] = [x // p for x in rows[i]]
        dep = None
        for c1 in range(p):
            for c2 in range(p):
                if (c1, c2) != (0, 0) and all((c1 * a + c2 * b) % p == 0 for a, b in zip(rows[0], rows[1])):
                    dep = (c1, c2)
                    break
            if dep:
                break
        if dep is None:
            return rows
        c1, c2 = dep
        new = [(c1 * a + c2 * b) // p for a, b in zip(rows[0], rows[1])]
        if c1 % p:
            rows[0] = new
        else:
            rows[1] = new


ntests = nviol = 0
hist = {}
for key in ('10', '01', '11'):
    O = ORD[key]
    e1, e2 = UNITS[key]
    e1i, e2i = UINV[key]
    for q in QS6:
        p = spf(q)
        e = round(math.log(q, p))
        thr = e - 1 - (1 if p == 2 else 0)
        kd = kernel_data(O, e1, e2, e1i, q)
        h1, b0, a1 = kd['h1'], kd['b0'], kd['a1']
        eta = [O.pow(e1, h1, e1i), O.mul(O.pow(e1, -a1, e1i), O.pow(e2, b0, e2i))]
        etai = [O.pow(e1, -h1, e1i), O.mul(O.pow(e1, a1, e1i), O.pow(e2, -b0, e2i))]
        check(all(O.mul(eta[j], etai[j]) == (1, 0, 0, 0) for j in range(2)), 'eta inverses')
        check(all(all((eta[j][i] - (1 if i == 0 else 0)) % q == 0 for i in range(4)) for j in range(2)),
              f'eta not == 1 mod {q}')
        for _ in range(NTR):
            v = tuple(RNG.randint(-6, 6) for _ in range(4))
            m0, n0 = RNG.choice([(a, b) for a in range(-2, 3) for b in range(-2, 3) if (a, b) != (0, 0)])
            u = v
            for j, k in ((0, m0), (1, n0)):
                for _ in range(abs(k)):
                    u = O.mul(eta[j] if k > 0 else etai[j], u)
            Mx = sp.Matrix([list(v), list(u)])
            if Mx.rank() < 2:
                continue
            rows = []
            for vec in Mx.nullspace():
                den = sp.ilcm(*[sp.Rational(x).q for x in vec])
                rows.append([int(sp.Rational(x) * den) for x in vec])
            rows = p_saturate(rows, p)
            if not all(sum(l[i] * v[i] for i in range(4)) == 0 and sum(l[i] * u[i] for i in range(4)) == 0
                       for l in rows):
                check(False, 'planted markers not orthogonal')
                continue
            w1, w2 = O.mul(eta[0], v), O.mul(eta[1], v)
            X = [sum(l[i] * w1[i] for i in range(4)) for l in rows]
            Y = [sum(l[i] * w2[i] for i in range(4)) for l in rows]
            if not all(t % q == 0 for t in X + Y):
                check(False, 'eta v not == v mod q')
                continue
            X, Y = [t // q for t in X], [t // q for t in Y]
            lam = min(nu(t, p, BIG) for t in X + Y)
            nJ = nu(X[0] * Y[1] - X[1] * Y[0], p, BIG)
            ok = lam <= thr and nJ - lam <= thr
            ntests += 1
            if ok:
                nviol += 1
                check(False, f'Prop. 17 holds on a planted two-zero system: {key} q={q} v={v} (m0,n0)={(m0, n0)}')
            mg = max(lam - thr, nJ - lam - thr)
            mg = 'inf' if mg >= BIG // 2 else mg
            hist[mg] = hist.get(mg, 0) + 1
fin = sorted(k for k in hist if k != 'inf')
print(f'   {ntests} planted two-zero systems (3 orders x moduli {QS6}); Proposition 17 held on {nviol} of them')
print('   margin histogram: ' + ', '.join(f'{k}: {hist[k]}' for k in fin + (['inf'] if 'inf' in hist else [])))
check(ntests >= 500, 'too few planted tests')
print(f'   [{time.time() - t6:.1f} s]')

# ======================================================================================================================
print()
print('=' * 118)
if FAILS:
    print(f'CHECK FAILED: {len(FAILS)} failed check(s)')
    for f in FAILS:
        print('   ' + f)
else:
    print('ALL CHECKS PASSED: in |a|,|b| <= 40 (exact) and |a|,|b| <= 1000 (modular, hits exact) the zeros of the two')
    print('markers are exactly the exponents printed by 02_exclusion.out.txt; the direct search |m| <= 10^5 finds')
    print('exactly core.KNOWN; the exclusion certificate re-checks with independent code; Proposition 17 fails, as it')
    print('must, on every planted two-zero system.  (Brute force is not a proof.)')
print(f'[done in {time.time() - T0:.1f} s]')
if FAILS:
    sys.exit(1)
