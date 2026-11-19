"""02_exclusion.py -- every unit orbit of 01_coverage.py has exactly as many solutions as are known.

For an orbit seed w (norm c) the elements of the family <eps1, eps2> w are eps^(a,b) w, (a,b) in Z^2; their
marker coordinates x2, x3 vanish exactly at the solutions of the field-form equation G(m,n) = c (one element of
each pair +-(m,n): proof.md 2.2 and the start of 4; the facts used there, N(eps1) = N(eps2) = 1, the shape of the bases and -1 not in <eps1,eps2>, are
asserted below).  The script certifies, for each orbit, that the zero set in Z^2 is exactly the set of exponents
of the known solutions.

Certificate (exact, finite):
* For a prime power q, Lambda_q = {(a,b) : eps^(a,b) == 1 mod q} (lifting2.Kernel).  All members of a class
  (a0,b0) + Lambda_q are == eps^(a0,b0) w (mod q), so a class whose residue has a marker not == 0 (mod q) has no zero.
* SIEVE.  Moduli q_1, q_2, ... are applied in turn.  The live classes are classes modulo L = Lambda_{q_1} cap ...
  cap Lambda_{q_k}, and a class stays live only if its image modulo every Lambda_{q_i} is a candidate.
* CERTIFIED CLASS.  A live class containing the exponent of a known zero z is removed as soon as Proposition 17
  of lifting.binomial.multi holds at z for one of the moduli q_i used (the class lies inside the Lambda_{q_i}-class
  of z, and Proposition 17 says z is the only zero there).  Proposition 17 is applied to the coset system
  T_j = M(eta_j), eta_j = eps^(lambda_j) for a basis lambda_1, lambda_2 of Lambda_q, seed z (lifting2.CosetSystem).
* SIBLING CLASS.  A live class containing no known zero, but lying inside the Lambda_{q_i}-class of a known zero z
  at which Proposition 17 holds for q_i, is empty (that Lambda_{q_i}-class has the single zero z).
* The orbit is done when no live class remains: its zeros are then exactly the known ones, one per certified
  class.  Every Proposition 17 instance used (for certification or for a sibling class) is printed, and so are the
  data needed to check the certificate without re-running the greedy search: for each modulus used, the Hermite
  data (h1, a1, b0) of Lambda_q = Z(h1,0) + Z(-a1,b0); after each step, the Hermite triple (A,B,C) of the current
  lattice L; each certified class with the triple of the lattice it is taken modulo; each sibling class with its
  witness (q, zero).
* Every known solution of the class is located in exactly one orbit (asserted), so no solution is skipped.
Lattices of Z^2 are kept in Hermite form rows (A,0), (B,C) with 0 <= B < A:  (a,b) in L iff C | b and
A | a - (b/C) B.
"""
import sys, time
sys.path.insert(0, '.')
from core import CURVE, CLASSES, KNOWN, UNITS, Field, Embeddings, Fval
from lifting2 import Kernel, CosetSystem, unit_from_exponents, locate

T0 = time.time()
SEEDS = {   # from 01_coverage.out.txt
    '10': {1: [(-1, 0, 0, 0)], 4: [(1, 0, 2, -1)], 9: [(1, -2, 3, -1)], 36: [(-2, -2, 3, -1)]},
    '01': {1: [(-1, 0, 0, 0)], 4: [(1, 1, -2, 0)], 9: [(0, -1, 0, 0)], 36: [(3, -1, 0, 0)]},
    '11': {1: [(-1, 0, 0, 0)], 4: [(4, -1, 0, 0)], 9: [(0, 2, -3, 0)], 36: [(0, -1, 0, 0)]},
}
MODULI = [4, 8, 16, 32, 64, 3, 9, 27, 81, 5, 25, 7, 49, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
MAXG = 5000           # skip moduli whose |G_q| exceeds this (only to bound running time)


# ------------------------------------------------------------------------------------------------------
# two-dimensional lattices
# ------------------------------------------------------------------------------------------------------
class Lat:
    def __init__(self, gens):
        rows = [list(g) for g in gens if g[0] or g[1]]
        # Euclid on the second coordinate
        while True:
            nz = [r for r in rows if r[1] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda r: abs(r[1]))
            piv = nz[0]
            for r in nz[1:]:
                k = r[1] // piv[1]
                r[0] -= k * piv[0]
                r[1] -= k * piv[1]
            rows = [r for r in rows if r[0] or r[1]]
        rb = [r for r in rows if r[1] != 0]
        assert len(rb) == 1, 'not full rank'
        piv = rb[0]
        if piv[1] < 0:
            piv = [-piv[0], -piv[1]]
        A = 0
        for r in rows:
            if r[1] == 0:
                A = _gcd(A, r[0])
        assert A > 0, 'not full rank'
        self.A, self.C = A, piv[1]
        self.B = piv[0] % A

    def reduce(self, v):
        a, b = v
        k = b // self.C
        a, b = a - k * self.B, b - k * self.C
        return (a % self.A, b)

    def contains(self, v):
        return self.reduce(v) == (0, 0)

    def basis(self):
        return [(self.A, 0), (self.B, self.C)]

    def index(self):
        return self.A * self.C


def _gcd(a, b):
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def lattice_of_kernel(ker):
    return Lat([(ker.h1, 0), (-ker.a1, ker.b0)])


def image_and_intersection(L, M):
    """BFS of the image of L in Z^2/M.  Returns (list of preimage vectors, one per element of the image,
    i.e. coset representatives of L/(L cap M)), and the lattice L cap M (from Schreier generators)."""
    gens = L.basis()
    start = (0, 0)
    seen = {M.reduce(start): start}
    queue = [start]
    schreier = []
    while queue:
        v = queue.pop()
        for g in gens:
            for sgn in (1, -1):
                w = (v[0] + sgn * g[0], v[1] + sgn * g[1])
                key = M.reduce(w)
                if key in seen:
                    u = seen[key]
                    d = (w[0] - u[0], w[1] - u[1])
                    if d != (0, 0):
                        schreier.append(d)
                else:
                    seen[key] = w
                    queue.append(w)
    N = L.index() * M.index()                      # N Z^2 lies in both lattices
    inter = Lat(schreier + [(N, 0), (0, N)])
    # checks: the generators lie in both lattices, and [L : L cap M] equals the size of the image
    assert all(L.contains(d) and M.contains(d) for d in inter.basis())
    reps = list(seen.values())
    assert inter.index() == L.index() * len(reps), 'index mismatch in intersection'
    return reps, inter


# ------------------------------------------------------------------------------------------------------
print(f'02_exclusion.py ({CURVE}) -- zeros of the two markers on every unit orbit (exact)')
print('moduli tried, in order:', MODULI, f'(skipped when |G_q| > {MAXG})')
print()
GRAND = {}
for key in SEEDS:
    K = Field(key)
    Em = Embeddings(K)
    e1, e2 = UNITS[key]
    print('=' * 110)
    print(f'CLASS {CLASSES[key]["label"]}  field form G = {CLASSES[key]["F"]}  eps1 = {e1}  eps2 = {e2}')
    # facts behind "markers vanish exactly at the solutions" (proof.md 2.2): omega_0 = 1, omega_1 = beta, and the
    # (beta^2, beta^3)-block of omega_2, omega_3 is invertible; units of norm 1, so T_j lies in SL_4(Z)
    Bm = K.B
    assert Bm[0] == [1, 0, 0, 0] and Bm[1] == [0, 1, 0, 0]
    assert Bm[2][2] * Bm[3][3] - Bm[2][3] * Bm[3][2] != 0
    assert K.norm(e1) == 1 and K.norm(e2) == 1
    print('   checked: omega_0 = 1, omega_1 = beta, invertible (beta^2, beta^3)-block (markers x2 = x3 = 0 <=> element'
          ' in Z + Z beta); N(eps1) = N(eps2) = 1')
    kers = {}
    for q in MODULI:
        p = next(d for d in range(2, q + 1) if q % d == 0)
        e = 0
        t = q
        while t % p == 0:
            t //= p
            e += 1
        if p == 2 and e < 2:
            continue
        ker = Kernel(K, q)
        if ker.order > MAXG:
            continue
        kers[q] = (p, e, ker, lattice_of_kernel(ker))
    print('   moduli available (|G_q|):', ', '.join(f'{q}({kers[q][2].order})' for q in kers))
    # -1 is not in <eps1, eps2>: it is not even a residue of <eps1, eps2> modulo some q >= 3
    qneg = next((q for q in kers if kers[q][2].star is None), None)
    if qneg is None:
        for q in (19, 11, 47, 23, 29, 37, 43):
            if Kernel(K, q).star is None:
                qneg = q
                break
    assert qneg is not None
    print(f'   checked: -1 is not a residue of <eps1,eps2> mod {qneg}, so -1 is not in <eps1,eps2>; with the injectivity of'
          ' (a,b) -> eps^(a,b) w (cone certificate of 01, proof.md Remark 3.2(i)), each pair +-z meets the family'
          ' <eps1,eps2> w at most once')
    located = {mn: [] for mn in KNOWN[key]}     # known solution -> the orbits (c, w) it was located in
    for c, seeds in SEEDS[key].items():
        for w in seeds:
            known = [(m, n) for (m, n) in KNOWN[key] if Fval(key, m, n) == c]
            zeros = []          # (mn, exponent (a,b), element of the w-family)
            for (m, n) in known:
                z = K.from_mn(m, n)
                loc = locate(K, z, w, e1, e2, Em)
                if loc is None:
                    continue
                s, a, b = loc
                zf = z if s == 0 else tuple(-t for t in z)
                assert K.mul(K.mul(K.pow(e1, a), K.pow(e2, b)), w) == zf
                zeros.append(((m, n), (a, b), zf))
                located[(m, n)].append((c, w))
            print(f'   --- orbit c = {c}, seed w = {w}; known solutions in it: {[z[0] for z in zeros]} '
                  f'(exponents {[z[1] for z in zeros]})')
            # Proposition 17 data at each known zero for each modulus
            p17 = {}
            for q, (p, e, ker, Lq) in kers.items():
                etas = [unit_from_exponents(K, e1, e2, *ex) for ex in ker.eta_exponents('eps')]
                for (mn, ab, zf) in zeros:
                    cs = CosetSystem(K, p, e, zf, etas)
                    r = cs.prop17()
                    p17[(q, mn)] = r
            cands = {q: set(ker.candidates(w)) for q, (p, e, ker, Lq) in kers.items()}
            # sieve
            L = Lat([(1, 0), (0, 1)])
            live = [(0, 0)]
            used = []
            history = []
            p17_used = set()
            certified = {}
            sibling_list = []       # (step modulus, class representative, Hermite triple of its lattice, witness)
            cert_lat = {}           # certified class representative -> Hermite triple of its lattice
            for step in range(12):
                # greedy: choose the modulus that leaves the fewest uncertified live classes
                best = None
                for q, (p, e, ker, Lq) in kers.items():
                    if q in used:
                        continue
                    reps, inter = image_and_intersection(L, Lq)
                    newlive = []
                    for r in live:
                        for x in reps:
                            v = (r[0] + x[0], r[1] + x[1])
                            if ker.reduce(*v) in cands[q]:
                                newlive.append(inter.reduce(v))
                    # certification.  A live class v (mod inter) is
                    #  - certified if it contains a known zero z and Proposition 17 holds at z for a used modulus;
                    #  - empty (dropped) if, for some used modulus qq and known zero z with Proposition 17 at (qq, z),
                    #    v lies in the Lambda_qq-class of z but does not contain z: that class has the single zero z.
                    newcert = {}
                    rest = []
                    nsib = 0
                    sibs_this = []
                    mods = used + [q]
                    for v in newlive:
                        hit = [zz for zz in zeros if inter.contains((v[0] - zz[1][0], v[1] - zz[1][1]))]
                        if hit:
                            okq = [qq for qq in mods if p17[(qq, hit[0][0])]['ok']]
                            if okq:
                                newcert[v] = (hit[0][0], okq[0])
                            else:
                                rest.append(v)
                            continue
                        sib = [(qq, zz[0]) for qq in mods for zz in zeros
                               if p17[(qq, zz[0])]['ok'] and kers[qq][3].contains((v[0] - zz[1][0], v[1] - zz[1][1]))]
                        if sib:
                            nsib += 1
                            sibs_this.append((v, sib[0]))
                            continue
                        rest.append(v)
                    score = (len(rest), inter.index())
                    if best is None or score < best[0]:
                        best = (score, q, inter, rest, newcert, (len(newlive), len(newcert), nsib, len(rest)), sibs_this)
                    # (newcert and sibs_this record class representatives modulo inter)
                if best is None:
                    break
                (_, q, inter, rest, newcert, counts, sibs_best) = best
                p17_used.update(s_ for _, s_ in sibs_best)
                trip = (inter.A, inter.B, inter.C)
                for v, s_ in sibs_best:
                    sibling_list.append((q, v, trip, s_))
                for v in newcert:
                    cert_lat[v] = trip
                history.append((q,) + counts + (trip,))
                used.append(q)
                L = inter
                live = rest
                certified.update(newcert)
                if not live:
                    break
            done = not live
            # final consistency: each known zero certified exactly once
            cz = [v[0] for v in certified.values()]
            ok = done and sorted(cz) == sorted(z[0] for z in zeros)
            GRAND[(key, c, w)] = (ok, len(zeros))
            print(f'      moduli used: {used};  final lattice index {L.index()}')
            for (q, nl, nc, ns, nr, trip) in history:
                print(f'         after q = {q:3d}: {nl} candidate classes mod the new lattice L = {trip} (Hermite (A,B,C):'
                      f' rows (A,0), (B,C)); {nc} certified (contain a known zero, Prop. 17), {ns} empty (siblings of a'
                      f' certified zero), {nr} left')
            for v, (mn, qq) in certified.items():
                r = p17[(qq, mn)]
                print(f'      class {v} mod L = {cert_lat[v]} contains the solution {mn}: Proposition 17 at q = {qq} '
                      f'(lambda = {r["lam"]}, nu(J) = {r["nuJ"]}, e = {kers[qq][1]}) -> exactly one zero')
            for (qs, v, trip, (qq, mn)) in sibling_list:
                print(f'      sibling class {v} mod L = {trip} (step q = {qs}): it lies in the Lambda_{qq}-class of the'
                      f' solution {mn} but does not contain it -> empty by Proposition 17 at q = {qq}')
            for v, (mn, qq) in certified.items():
                p17_used.add((qq, mn))
            for (qq, mn) in sorted(p17_used):
                r = p17[(qq, mn)]
                pq, eq = kers[qq][0], kers[qq][1]
                bound = eq - 1 - (1 if pq == 2 else 0)
                print(f'      Prop. 17 instance used: q = {qq} (p = {pq}, e = {eq}), zero {mn}: lambda = {r["lam"]}, '
                      f'nu(J) = {r["nuJ"]}, delta = {r["delta"]} <= e-1-iota = {bound}: {r["ok"]}')
                assert r['ok']
            for q in used:
                kq = kers[q][2]
                print(f'      q = {q}: |G_q| = {kq.order}, Lambda_q = Z({kq.h1},0) + Z({-kq.a1},{kq.b0}), candidate classes'
                      f' mod Lambda_q: {len(cands[q])}')
            print(f'      uncertified live classes left: {len(live)}  ->  '
                  f'{"ORBIT DONE: zeros = the known solutions" if ok else "NOT DONE"}')
    bad = {mn: loc for mn, loc in located.items() if len(loc) != 1}
    print(f'   checked: each of the {len(located)} known solutions of the class lies in exactly one orbit: {not bad}')
    assert not bad, bad
print()
print('SUMMARY')
alldone = True
for (key, c, w), (ok, nz) in GRAND.items():
    print(f'   class {CLASSES[key]["label"]} c = {c:2d} seed {w}: {nz} known solution pair(s); certified: {ok}')
    alldone = alldone and ok
print('ALL ORBITS CERTIFIED' if alldone else 'CHECK FAILED: some orbit is not certified')
print(f'[done in {time.time() - T0:.1f} s]')
if not alldone:
    sys.exit(1)
