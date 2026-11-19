"""01_coverage.py -- elementary two-unit coverage: every element of norm c is +-eps1^a eps2^b w for a listed seed w.

Lemma (reduced representatives).  Let U be a finite subset of G = <eps1, eps2> (units of norm +1).  For an
element z of O with N(z) = c > 0, the set {g z : g in G, E(g z) <= E(z)} is finite (lattice points in an
ellipsoid), so it has an element w of least E.  For u in U, E(u w) >= E(w), and since s_v(u w) = f_v(u) s_v(w)
with f_v(u) = s_v(u),
        E(u w) - E(w) = sum_v s_v(w) (f_v(u) - 1) >= 0.                                        (*)
So s(w) lies in the cone C_U cut out by (*).  If f_v(u) <= f_v^+(u) (rational upper bounds), C_U lies in the
rational cone C_U^+.  Let P^+ = C_U^+ intersected with {r >= 0, r_1 + r_2 + r_3 = 1}.  If every vertex of P^+ has
all r_v > 0, then on P^+ the product r_1 r_2 r_3^2 is bounded below by its least value at a vertex, b_min
(a product of positive linear functions is log-concave on segments, so its minimum on a polygon is at a vertex).
Writing s(w) = E(w) r with r in P^+, c^2 = N(w)^2 = s_1 s_2 s_3^2 = E(w)^4 r_1 r_2 r_3^2 >= E(w)^4 b_min, hence
        E(w) <= T  for any rational T with T^4 >= c^2 / b_min.
The cone condition also forces eps1, eps2 to be multiplicatively independent: the log size vectors of U then span
the plane log s_1 + log s_2 + 2 log s_3 = 0 (proof.md, Remark 3.2(i)).

Enumeration.  E(x) = x^T G x with an interval Gram matrix G (core.Embeddings.gram).  With G_mid the midpoint
matrix and rho the largest radius, x^T G x >= x^T (G_mid - 4 rho I) x =: x^T Gt x for integer x.  The script lists
every x with x^T Gt x <= T exactly (Fincke-Pohst over Q), keeps those with N(x) = c, and merges them into orbits of
<-1, eps1, eps2> by exact unit quotients.  Every orbit of an element of norm c has a member in this list.
"""
import sys, time, math
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from core import CLASSES, RHS, KNOWN, UNITS, Field, Embeddings, Iv, Fval

T0 = time.time()
print('01_coverage.py -- reduced representatives for the unit orbits (exact; intervals with dyadic endpoints)')
print()

EXPONENTS = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)]


def polygon_vertices(cons):
    """Vertices of {r : r_1 + r_2 + r_3 = 1, g.r >= 0 for g in cons} (cons includes r_v >= 0), exactly."""
    verts = []
    n = len(cons)
    for i in range(n):
        for j in range(i + 1, n):
            A = [cons[i], cons[j], [Fr(1)] * 3]
            det = (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                   + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))
            if det == 0:
                continue
            # solve A r = (0, 0, 1) by Cramer
            rhs = [Fr(0), Fr(0), Fr(1)]
            r = []
            for k in range(3):
                Ak = [row[:] for row in A]
                for t in range(3):
                    Ak[t][k] = rhs[t]
                dk = (Ak[0][0] * (Ak[1][1] * Ak[2][2] - Ak[1][2] * Ak[2][1])
                      - Ak[0][1] * (Ak[1][0] * Ak[2][2] - Ak[1][2] * Ak[2][0])
                      + Ak[0][2] * (Ak[1][0] * Ak[2][1] - Ak[1][1] * Ak[2][0]))
                r.append(dk / det)
            if all(sum(g[t] * r[t] for t in range(3)) >= 0 for g in cons) and r not in verts:
                verts.append(r)
    return verts


def fincke_pohst(Q, T):
    """All integer x != 0 with x^T Q x <= T, for a rational positive definite 4x4 Q (exact)."""
    n = 4
    q = [[Fr(Q[i][j]) for j in range(n)] for i in range(n)]
    for i in range(n):
        assert q[i][i] > 0, 'Gram matrix not positive definite'
        for j in range(i + 1, n):
            q[j][i] = q[i][j]
            q[i][j] = q[i][j] / q[i][i]
        for k in range(i + 1, n):
            for l in range(k, n):
                q[k][l] -= q[k][i] * q[i][l]
    out = []
    x = [0] * n

    def rec(i, budget):
        center = -sum(q[i][j] * x[j] for j in range(i + 1, n))
        rad = math.sqrt(max(float(budget / q[i][i]), 0.0))
        lo = math.floor(float(center) - rad) - 1
        hi = math.ceil(float(center) + rad) + 1
        # exact certificate that the float-guided range misses nothing: t -> q_ii (t - center)^2 is convex, so
        # its integer solutions of "<= budget" form an interval around center; lo and hi lie outside it
        assert lo < center < hi
        assert q[i][i] * (lo - center) ** 2 > budget and q[i][i] * (hi - center) ** 2 > budget
        for xi in range(lo, hi + 1):
            t = xi - center
            used = q[i][i] * t * t
            if used <= budget:
                x[i] = xi
                if i == 0:
                    out.append(tuple(x))
                else:
                    rec(i - 1, budget - used)
        x[i] = 0

    rec(n - 1, Fr(T))
    return [v for v in out if any(v)]


SUMMARY = {}
for key in CLASSES:
    K = Field(key)
    Em = Embeddings(K)
    e1, e2 = UNITS[key]
    print('=' * 110)
    print(f'CLASS {CLASSES[key]["label"]}  F = {CLASSES[key]["F"]}  units eps1 = {e1}, eps2 = {e2}')
    assert K.norm(e1) == 1 and K.norm(e2) == 1
    print('   N(eps1) = N(eps2) = 1 (exact determinants)')
    # f-vectors of U
    cons = [[Fr(int(v == t)) for t in range(3)] for v in range(3)]          # r_v >= 0
    print('   U = {eps1^a eps2^b : a,b in {-1,0,1}, (a,b) != (0,0)};  f_v(u) = s_v(u) as intervals:')
    for (a, b) in EXPONENTS:
        u = K.mul(K.pow(e1, a), K.pow(e2, b))
        f = Em.svec(u)
        prod = f[0] * f[1] * f[2].sq()
        assert prod.lo <= 1 <= prod.hi
        cons.append([f[v].hi - 1 for v in range(3)])                        # relaxed: f_v^+ - 1
        print(f'      (a,b)=({a:2d},{b:2d})  f = ({float(f[0]):.6g}, {float(f[1]):.6g}, {float(f[2]):.6g})'
              f'   [interval widths < {float(max(x.width() for x in f)):.1e}]')
    verts = polygon_vertices(cons)
    assert verts, 'empty polygon'
    allpos = all(all(r > 0 for r in v) for v in verts)
    print(f'   relaxed polygon P+ has {len(verts)} vertices; all coordinates of all vertices > 0: {allpos}')
    assert allpos
    bmin = min(v[0] * v[1] * v[2] * v[2] for v in verts)
    for v in verts:
        print(f'      vertex ({float(v[0]):.6f}, {float(v[1]):.6f}, {float(v[2]):.6f})   r1 r2 r3^2 = {float(v[0]*v[1]*v[2]*v[2]):.6e}')
    print(f'   b_min = min r1 r2 r3^2 over vertices = {float(bmin):.6e};  E(w) <= (c^2/b_min)^(1/4) = {float(1/bmin)**0.25:.4f} sqrt(c)')
    # Gram matrix: rational lower form
    G = Em.gram()
    Gmid = [[G[i][j].mid() for j in range(4)] for i in range(4)]
    rho = max(G[i][j].hi - G[i][j].mid() for i in range(4) for j in range(4))
    Gt = [[Gmid[i][j] - (4 * rho if i == j else 0) for j in range(4)] for i in range(4)]
    print(f'   Gram matrix of E: max interval radius rho = {float(rho):.2e}; lower form Gt = Gmid - 4 rho I')
    SUMMARY[key] = {}
    for c in RHS:
        # least integer T with T^4 >= c^2 / bmin
        target = Fr(c * c) / bmin
        T = int(float(target) ** 0.25) + 1
        while Fr(T - 1) ** 4 >= target:
            T -= 1
        assert Fr(T) ** 4 >= target and (T == 1 or Fr(T - 1) ** 4 < target)
        pts = fincke_pohst(Gt, T)
        normc = [x for x in pts if K.norm(x) == c]
        # merge into orbits of <-1, eps1, eps2> by exact quotients
        orbits = []           # list of lists of elements
        import numpy as np
        fl = lambda x: (float(Em.sig1(x)), float(Em.sig2(x)))
        l1 = np.array([math.log(abs(v)) for v in fl(e1)])
        l2 = np.array([math.log(abs(v)) for v in fl(e2)])
        Mat = np.array([l1, l2]).T
        newunits = []
        for x in normc:
            placed = False
            for orb in orbits:
                w = orb[0]
                qt = K.divides_exactly(w, x)                # x / w
                if qt is None or K.norm(qt) != 1:
                    continue
                lq = np.array([math.log(abs(v)) for v in fl(qt)])
                ab = np.linalg.solve(Mat, lq)
                a, b = int(round(ab[0])), int(round(ab[1]))
                g = K.mul(K.pow(e1, a), K.pow(e2, b))
                if g == qt or g == tuple(-t for t in qt):
                    orb.append(x)
                    placed = True
                    break
                newunits.append(qt)
            if not placed:
                orbits.append([x])
        # the seed of each orbit is its listed member of least E (floats only choose among exact members)
        seeds = []
        for orb in orbits:
            orb.sort(key=lambda y: float(Em.E(y)))
            seeds.append((orb[0], len(orb)))
        print(f'   c = {c:2d}: T = {T} (T^4 >= c^2/b_min);  lattice points with Gt <= T: {len(pts)};  of norm c: {len(normc)};'
              f'  orbits of <-1,eps1,eps2>: {len(orbits)}')
        if newunits:
            print(f'      NOTE: {len(newunits)} unit quotients not in <-1,eps1,eps2> were met (orbits kept separate)')
        SUMMARY[key][c] = [s for s, _ in seeds]
        for (w, size) in seeds:
            print(f'      seed w = {w}   E(w) ~ {float(Em.E(w)):.4f}   ({size} listed members of this orbit)')
        if c == 1:
            # by-product: the unit group.  Lemma with c = 1 and c = -1 (same T, since only c^2 enters): every unit of
            # norm +1 is g*w with w in the c = 1 list, and every unit of norm -1 is g*w with N(w) = -1, E(w) <= T.
            one_orbit = len(orbits) == 1 and any(y in ((1, 0, 0, 0), (-1, 0, 0, 0)) for y in orbits[0])
            nminus = [y for y in pts if K.norm(y) == -1]
            print(f'      unit group: the units of norm +1 with E <= T form one orbit containing +-1: {one_orbit}; '
                  f'elements of norm -1 with E <= T: {len(nminus)}')
            assert one_orbit and not nminus
            print('      => O^x = {+-1} x <eps1> x <eps2>: eps1, eps2 are fundamental units of O (not needed below)')
print()
print('SEEDS = {')
for key in SUMMARY:
    print(f"    '{key}': {{" + ', '.join(f'{c}: {SUMMARY[key][c]}' for c in RHS) + '},')
print('}')
print(f'[done in {time.time() - T0:.1f} s]')
