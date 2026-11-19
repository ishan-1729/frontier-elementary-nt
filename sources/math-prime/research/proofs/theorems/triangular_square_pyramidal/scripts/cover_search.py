"""cover_search.py: how the box table of hypothesis_H.md, Section 2, was found (not part of the proof).

The search runs in floating point. The proof uses only the table, which 05_cover.py checks in exact rational
arithmetic; this script also re-checks the chosen (s, kappa) of every box exactly before writing the table.
Adapted from triangular_tetrahedral/scripts/cover_search_K.py, with the ring data of hring.py.

Setting. D = [-1/2,1/2] x [-1/2,1/2] x [0,1/2] in the coordinates of the ring's basis. A closed box with centre c and
half-widths h is accepted with a multiplier s in MULTS and an integer triple kappa when the Taylor bound
    B = |N(v)| + sum_i |g_i(v)| s h_i + sum_{i<=j} |q_ij(v)| s^2 h_i h_j + sum_m |n_m| (s h)^m,    v = s*c - kappa,
is < 1 (hlib.taylor_bound); then |N(s*xi - kappa)| <= B < 1 on the whole box.

Candidates. For s in MULTS and k in UNIT_K: kappa = U^(-k) * kappa' with kappa' in floor(U^k * s*c) + WINDOW^3,
where U is the unit of hring.py; since N(U) = 1, N(s*xi - kappa) = N(U^k*s*xi - kappa').
Split rules: 'long' (split the longest relative side) and 'reduce_sum' (split the axis whose halves have the least
sum of best bounds). The rule with fewer leaves is kept. On each leaf the table uses, among the accepted
candidates, the least s, then the least max |kappa_i|, then the least bound.
Usage: python cover_search.py   (writes TABLE, prints statistics)
"""
import itertools
import sys
import time
from fractions import Fraction as Fr

import numpy as np

import hring
from hlib import taylor_data, taylor_bound

sys.setrecursionlimit(20000)

# ---- configuration (the only card-specific part) ------------------------------------------------------------------
MULTS = (1,)                         # plain Euclidean division
UNIT = hring.EPS
UNIT_K = (0, 1, -1)
WINDOW = range(-3, 5)
TABLE = 'euclid_boxes.txt'
TITLE = "euclid_boxes.txt: the box table of Lemma 2.1 (hypothesis_H.md, Section 2)."
COORDS = "the coordinates (a,b,c) of a + b*l + c*l^2, l^3 = l - 3"
# --------------------------------------------------------------------------------------------------------------------

KMAX = 24
UNITF = 2.0**-KMAX
ROOT = ((-(1 << (KMAX - 1)), -(1 << (KMAX - 1)), 0), (1 << (KMAX - 1), 1 << (KMAX - 1), 1 << (KMAX - 1)))
ACCEPT = 1 - 1e-7
LETTERS = 'xyz'
CUBIC = hring.NORM_FORM
GRAD, QUAD = taylor_data(CUBIC)


def mats(u):
    return np.array(hring.matrix(u), dtype=np.int64)


MK = {0: np.eye(3, dtype=np.int64), 1: mats(UNIT), -1: mats(hring.inverse_unit(UNIT))}
WIN = np.array(list(itertools.product(WINDOW, repeat=3)), dtype=np.int64)


def fev(form, x, y, z):
    out = 0.0
    for m, c in form.items():
        out = out + float(c)*x**m[0]*y**m[1]*z**m[2]
    return out


def candidates(key):
    lo = np.array(key[0], dtype=float)*UNITF
    hi = np.array(key[1], dtype=float)*UNITF
    c = (lo + hi)/2
    h = (hi - lo)/2
    S, KAP, BB = [], [], []
    for s in MULTS:
        for k in UNIT_K:
            t = MK[k] @ (s*c)
            kp = np.floor(t).astype(np.int64) + WIN
            vp = t - kp
            kap = kp @ MK[-k].T
            v = s*c - kap
            x, y, z = v[:, 0], v[:, 1], v[:, 2]
            H = s*h
            B = np.abs(fev(CUBIC, vp[:, 0], vp[:, 1], vp[:, 2]))
            for i in range(3):
                B = B + np.abs(fev(GRAD[i], x, y, z))*H[i]
            for (i, j), q in QUAD.items():
                B = B + np.abs(fev(q, x, y, z))*H[i]*H[j]
            B = B + sum(abs(float(cf))*H[0]**m[0]*H[1]**m[1]*H[2]**m[2] for m, cf in CUBIC.items())
            S.append(np.full(len(B), s))
            KAP.append(kap)
            BB.append(B)
    return np.concatenate(S), np.vstack(KAP), np.concatenate(BB)


MEMO = {}


def min_bound(key):
    if key not in MEMO:
        S, KAP, B = candidates(key)
        MEMO[key] = float(B.min())
    return MEMO[key]


def halves(key, ax):
    lo, hi = list(key[0]), list(key[1])
    assert (lo[ax] + hi[ax]) % 2 == 0, "resolution exhausted"
    mid = (lo[ax] + hi[ax])//2
    hi1, lo2 = hi.copy(), lo.copy()
    hi1[ax] = mid
    lo2[ax] = mid
    return (tuple(lo), tuple(hi1)), (tuple(lo2), tuple(hi))


def rel_widths(key):
    w = [key[1][i] - key[0][i] for i in range(3)]
    w[2] *= 2
    return w


def greedy(rule):
    leaves = []

    def go(key, path):
        if min_bound(key) <= ACCEPT:
            leaves.append((path, key))
            return
        assert len(path) < 66, "too deep"
        w = rel_widths(key)
        if rule == 'long':
            ax = w.index(max(w))
        else:
            axes = [a for a in range(3) if 4*w[a] >= max(w)]
            ax = min(axes, key=lambda a: sum(min_bound(ch) for ch in halves(key, a)))
        lo_half, hi_half = halves(key, ax)
        go(lo_half, path + LETTERS[ax])
        go(hi_half, path + LETTERS[ax].upper())

    go(ROOT, '')
    return leaves


def choose(key):
    S, KAP, B = candidates(key)
    ok = np.nonzero(B <= ACCEPT)[0]
    order = sorted(ok, key=lambda i: (S[i], int(np.abs(KAP[i]).max()), B[i]))
    lo = [Fr(t, 1 << KMAX) for t in key[0]]
    hi = [Fr(t, 1 << KMAX) for t in key[1]]
    for i in order:
        kap = tuple(int(t) for t in KAP[i])
        Bx = taylor_bound(CUBIC, (GRAD, QUAD), lo, hi, int(S[i]), kap)
        if Bx < 1:
            return int(S[i]), kap, float(B[i]), Bx
    raise RuntimeError("no candidate passes the exact check")


def main():
    results = {}
    for rule in ('long', 'reduce_sum'):
        t0 = time.time()
        leaves = greedy(rule)
        results[rule] = leaves
        print(f"    {rule:10s} leaves {len(leaves):6d}   evaluations so far {len(MEMO):7d}   {time.time() - t0:7.1f} s"
              f"   depth {max(len(p) for p, _ in leaves)}", flush=True)
    rule = min(results, key=lambda r: (len(results[r]), max(len(p) for p, _ in results[r])))
    leaves = results[rule]
    print(f"kept: {rule} with {len(leaves)} boxes")
    rows = []
    for path, key in leaves:
        s, kap, Bf, Bx = choose(key)
        rows.append((path, key, s, kap, Bf, Bx))
    K = 0
    while any(any(t % (1 << (KMAX - K)) for t in key[0] + key[1]) for _, key, *_ in rows):
        K += 1
    sh = KMAX - K
    print(f"all box ends are multiples of 2^-{K}")
    print(f"largest number of splits on a path: {max(len(r[0]) for r in rows)}")
    for s in MULTS:
        print(f"    s = {s}: {sum(1 for r in rows if r[2] == s)} boxes")
    worst = max(rows, key=lambda r: r[5])
    print(f"largest exact bound B: {worst[5]} = {float(worst[5]):.9f} (box {worst[0]})")
    print(f"largest |kappa_i|: {max(max(abs(t) for t in r[3]) for r in rows)}")
    print(f"exact bound < 1 on every box: {all(r[5] < 1 for r in rows)}")
    with open(TABLE, 'w', newline='\n') as out:
        out.write(f"# {TITLE}\n")
        out.write("# Written by cover_search.py; read and checked by 05_cover.py. The proof uses only this table.\n")
        out.write(f"# D = [-1/2,1/2] x [-1/2,1/2] x [0,1/2] in {COORDS}.\n")
        out.write("# One line per box. The boxes are the leaves of a binary split tree of D, in depth-first order with\n")
        out.write("# the lower half first.\n")
        out.write("#   path   the splits from D to the box, one letter per split at the midpoint of an axis:\n")
        out.write("#          x, y, z = lower half along that axis, X, Y, Z = upper half\n")
        out.write(f"#   lo, hi the corners (lo_a lo_b lo_c hi_a hi_b hi_c), integers in units of 2^-{K}\n")
        out.write(f"#   s      the multiplier, in {set(MULTS)}\n")
        out.write("#   kappa  an integer triple (k1 k2 k3); on the box |N(s*xi - kappa)| <= B < 1 with the Taylor bound B\n")
        out.write(f"unit 2^-{K}\n")
        out.write("# path lo_a lo_b lo_c hi_a hi_b hi_c s k1 k2 k3\n")
        for path, key, s, kap, Bf, Bx in rows:
            ends = ' '.join(str(t >> sh) for t in key[0] + key[1])
            out.write(f"{path} {ends} {s} {kap[0]} {kap[1]} {kap[2]}\n")
    print(f"wrote {TABLE} ({len(rows)} boxes)")


if __name__ == '__main__':
    T0 = time.time()
    main()
    print(f"total time {time.time() - T0:.1f} s")
