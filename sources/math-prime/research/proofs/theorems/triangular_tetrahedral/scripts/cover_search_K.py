"""cover_search_K.py: how the box table dh_boxes.txt was found (not part of the proof).

The search runs in floating point. The proof uses only the table, which 05_dh_cover.py checks in exact rational
arithmetic. This script also re-checks the chosen (s, kappa) of every box exactly before writing the table.

Setting (hypothesis_H.md, Section 2). R = Z[l] with l^3 = 4l - 6; a triple (a,b,c) is a + b*l + c*l^2, and N is the
norm form of cubic.py. D = [-1/2,1/2] x [-1/2,1/2] x [0,1/2]. A closed box with centre c and half-widths h is
accepted with s in {1,2} and an integer triple kappa when the Taylor bound of 05_dh_cover.py (2.2),
    B = |N(v)| + sum_i |g_i(v)| s h_i + sum_{i<=j} |q_ij(v)| s^2 h_i h_j + sum_m |n_m| (s h)^m,    v = s*c - kappa,
is < 1; then |N(s*xi - kappa)| <= B < 1 on the whole box.

Candidates. For s in {1,2} and k in {-1,0,1}: kappa = EPS^(-k) * kappa' with kappa' in floor(EPS^k * s*c) + {-3..4}^3.
EPS = l^2 - 5l - 19 has norm 1, so N(s*xi - kappa) = N(EPS^k * s*xi - kappa'): the unit moves the search to where
the small values of N are. A candidate kappa can have coordinates of several hundred. In floating point |N(v)| is
evaluated as |N(EPS^k*s*c - kappa')|, which has small coordinates and is accurate, and a box is accepted when the
float bound is <= 1 - 1e-7.

Split rules tried. A box that is not accepted is halved at the midpoint of one axis. Side lengths are measured
relative to D (the third one doubled), and "thin" axes (relative side below 1/4 of the longest) are not split by
the reduce rules.
  long        split the longest side;
  reduce_sum  split the axis whose two halves have the least sum of best bounds;
  reduce_max  split the axis whose two halves have the least larger best bound;
  dp2, dp4    the fewest leaves over all midpoint split trees in which only axes with relative side >= 1/2
              (dp2) or >= 1/4 (dp4) of the longest are split, and among those the fewest leaves that need
              s = 2: an exact dynamic programme over all such boxes.
The rule with the fewest leaves (then the fewest s = 2 leaves) is kept. On each leaf the table uses, among the
accepted candidates, one with s = 1 if there is one, then the least max |kappa_i|, then the least bound.

Output: statistics on stdout, and the table dh_boxes.txt (format in its header).
"""
import itertools
import sys
import time
from fractions import Fraction as Fr

import numpy as np

from cubic import power, matrix, norm, EPS

sys.setrecursionlimit(10000)

KMAX = 24                       # box ends are integers in units of 2^-KMAX during the search
UNIT = 2.0**-KMAX
ROOT = ((-(1 << (KMAX - 1)), -(1 << (KMAX - 1)), 0), (1 << (KMAX - 1), 1 << (KMAX - 1), 1 << (KMAX - 1)))
ACCEPT = 1 - 1e-7
MK = {k: np.array(matrix(power(EPS, k)), dtype=np.int64) for k in (-1, 0, 1)}
WINDOW = np.array(list(itertools.product(range(-3, 5), repeat=3)), dtype=np.int64)
SK = [(s, k) for s in (1, 2) for k in (0, 1, -1)]
LETTERS = 'xyz'


def norm_f(x, y, z):
    return x**3 + 8*x*x*z - 4*x*y*y + 18*x*y*z + 16*x*z*z - 6*y**3 + 24*y*z*z + 36*z**3


def candidates(key):
    """All candidates of a box: arrays of s, k, kappa (n x 3) and float bounds B."""
    lo = np.array(key[0], dtype=float)*UNIT
    hi = np.array(key[1], dtype=float)*UNIT
    c = (lo + hi)/2
    h = (hi - lo)/2
    S, K, KAP, BB = [], [], [], []
    for s, k in SK:
        t = MK[k] @ (s*c)
        kp = np.floor(t).astype(np.int64) + WINDOW
        vp = t - kp                                  # EPS^k * s*c - kappa', small
        kap = kp @ MK[-k].T                          # kappa = EPS^-k * kappa'
        v = s*c - kap
        x, y, z = v[:, 0], v[:, 1], v[:, 2]
        h1, h2, h3 = s*h
        g1 = 3*x*x + 16*x*z - 4*y*y + 18*y*z + 16*z*z
        g2 = -8*x*y + 18*x*z - 18*y*y + 24*z*z
        g3 = 8*x*x + 18*x*y + 32*x*z + 48*y*z + 108*z*z
        q = (np.abs(3*x + 8*z)*h1*h1 + np.abs(-8*y + 18*z)*h1*h2 + np.abs(16*x + 18*y + 32*z)*h1*h3
             + np.abs(-4*x - 18*y)*h2*h2 + np.abs(18*x + 48*z)*h2*h3 + np.abs(16*x + 24*y + 108*z)*h3*h3)
        cub = (h1**3 + 8*h1*h1*h3 + 4*h1*h2*h2 + 18*h1*h2*h3 + 16*h1*h3*h3 + 6*h2**3 + 24*h2*h3*h3
               + 36*h3**3)
        B = (np.abs(norm_f(vp[:, 0], vp[:, 1], vp[:, 2])) + np.abs(g1)*h1 + np.abs(g2)*h2 + np.abs(g3)*h3
             + q + cub)
        S.append(np.full(len(B), s))
        K.append(np.full(len(B), k))
        KAP.append(kap)
        BB.append(B)
    return np.concatenate(S), np.concatenate(K), np.vstack(KAP), np.concatenate(BB)


MEMO = {}


def bounds(key):
    """(least float bound over the candidates with s = 1, least float bound over all candidates)."""
    if key not in MEMO:
        S, K, KAP, B = candidates(key)
        MEMO[key] = (float(B[S == 1].min()), float(B.min()))
    return MEMO[key]


def min_bound(key):
    return bounds(key)[1]


def halves(key, ax):
    lo, hi = list(key[0]), list(key[1])
    assert (lo[ax] + hi[ax]) % 2 == 0, "resolution 2^-KMAX exhausted"
    mid = (lo[ax] + hi[ax])//2
    hi1, lo2 = hi.copy(), lo.copy()
    hi1[ax] = mid
    lo2[ax] = mid
    return (tuple(lo), tuple(hi1)), (tuple(lo2), tuple(hi))


def rel_widths(key):
    w = [key[1][i] - key[0][i] for i in range(3)]
    w[2] *= 2
    return w


def splittable(key, ratio):
    w = rel_widths(key)
    return [ax for ax in range(3) if w[ax]*ratio >= max(w)]


def greedy(rule):
    """Leaves (path, key) of the tree built by a greedy rule, in depth-first order."""
    leaves = []

    def go(key, path):
        if min_bound(key) <= ACCEPT:
            leaves.append((path, key))
            return
        assert len(path) < 60
        if rule == 'long':
            w = rel_widths(key)
            ax = w.index(max(w))
        else:
            agg = sum if rule == 'reduce_sum' else max
            ax = min(splittable(key, 4), key=lambda a: agg(min_bound(ch) for ch in halves(key, a)))
        lo_half, hi_half = halves(key, ax)
        go(lo_half, path + LETTERS[ax])
        go(hi_half, path + LETTERS[ax].upper())

    go(ROOT, '')
    return leaves


def dp(ratio):
    """Over all trees that split only axes allowed by `ratio`: the fewest leaves, and among those the fewest
    leaves that need s = 2 (exact dynamic programme; costs are pairs compared lexicographically)."""
    best = {}

    def cost(key):
        if key in best:
            return best[key][0]
        b1, b = bounds(key)
        if b <= ACCEPT:
            best[key] = ((1, 0 if b1 <= ACCEPT else 1), None)
            return best[key][0]
        res = None
        for ax in splittable(key, ratio):
            lo_half, hi_half = halves(key, ax)
            c1, c2 = cost(lo_half), cost(hi_half)
            v = (c1[0] + c2[0], c1[1] + c2[1])
            if res is None or v < res[0]:
                res = (v, ax)
        best[key] = res
        return res[0]

    cost(ROOT)
    leaves = []

    def walk(key, path):
        n, ax = best[key]
        if ax is None:
            leaves.append((path, key))
            return
        lo_half, hi_half = halves(key, ax)
        walk(lo_half, path + LETTERS[ax])
        walk(hi_half, path + LETTERS[ax].upper())

    walk(ROOT, '')
    return leaves


def exact_bound(lo, hi, s, kap):
    """The bound B of 05_dh_cover.py (2.2), in exact rational arithmetic."""
    c = [(lo[i] + hi[i])/2 for i in range(3)]
    h1, h2, h3 = [s*(hi[i] - lo[i])/2 for i in range(3)]
    x, y, z = [s*c[i] - kap[i] for i in range(3)]
    g = (3*x*x + 16*x*z - 4*y*y + 18*y*z + 16*z*z, -8*x*y + 18*x*z - 18*y*y + 24*z*z,
         8*x*x + 18*x*y + 32*x*z + 48*y*z + 108*z*z)
    q = (abs(3*x + 8*z)*h1*h1 + abs(-8*y + 18*z)*h1*h2 + abs(16*x + 18*y + 32*z)*h1*h3
         + abs(-4*x - 18*y)*h2*h2 + abs(18*x + 48*z)*h2*h3 + abs(16*x + 24*y + 108*z)*h3*h3)
    cub = h1**3 + 8*h1*h1*h3 + 4*h1*h2*h2 + 18*h1*h2*h3 + 16*h1*h3*h3 + 6*h2**3 + 24*h2*h3*h3 + 36*h3**3
    return abs(norm((x, y, z))) + abs(g[0])*h1 + abs(g[1])*h2 + abs(g[2])*h3 + q + cub


def choose(key):
    """The preferred accepted candidate of a leaf, re-checked exactly."""
    S, K, KAP, B = candidates(key)
    ok = np.nonzero(B <= ACCEPT)[0]
    order = sorted(ok, key=lambda i: (S[i], int(np.abs(KAP[i]).max()), B[i]))
    lo = [Fr(t, 1 << KMAX) for t in key[0]]
    hi = [Fr(t, 1 << KMAX) for t in key[1]]
    for i in order:
        kap = tuple(int(t) for t in KAP[i])
        Bx = exact_bound(lo, hi, int(S[i]), kap)
        if Bx < 1:
            return int(S[i]), int(K[i]), kap, float(B[i]), Bx
    raise RuntimeError("no candidate passes the exact check")


def main():
    results = {}
    print("split rules (leaves, evaluations so far, seconds, largest number of splits on a path)")
    for rule in ('long', 'reduce_sum', 'reduce_max', 'dp2', 'dp4'):
        t0 = time.time()
        if rule.startswith('dp'):
            leaves = dp(int(rule[2:]))
        else:
            leaves = greedy(rule)
        results[rule] = leaves
        print(f"    {rule:10s} leaves {len(leaves):5d}   box evaluations so far {len(MEMO):6d}"
              f"   {time.time() - t0:6.1f} s   depth {max(len(p) for p, _ in leaves)}", flush=True)
    def n_s2(leaves):
        return sum(1 for _, key in leaves if bounds(key)[0] > ACCEPT)

    for rule in results:
        print(f"    {rule:10s} leaves that need s = 2: {n_s2(results[rule])}")
    rule = min(results, key=lambda r: (len(results[r]), n_s2(results[r]), max(len(p) for p, _ in results[r])))
    leaves = results[rule]
    print(f"kept: {rule} with {len(leaves)} boxes")

    rows = []
    for path, key in leaves:
        s, k, kap, Bf, Bx = choose(key)
        rows.append((path, key, s, k, kap, Bf, Bx))
    # the coarsest common unit 2^-K of all box ends
    K = 0
    while any(any(t % (1 << (KMAX - K)) for t in key[0] + key[1]) for _, key, *_ in rows):
        K += 1
    sh = KMAX - K
    print(f"all box ends are multiples of 2^-{K}")
    depth = max(len(r[0]) for r in rows)
    axis_depth = [max(r[0].lower().count(a) for r in rows) for a in LETTERS]
    print(f"largest number of splits on a path: {depth}; per axis (x, y, z): {axis_depth}")
    for s in (1, 2):
        for k in (0, 1, -1):
            n = sum(1 for r in rows if r[2] == s and r[3] == k)
            print(f"    s = {s}, kappa found with k = {k:2d}: {n:4d} boxes")
    print(f"s = 1: {sum(1 for r in rows if r[2] == 1)} boxes, s = 2: {sum(1 for r in rows if r[2] == 2)} boxes")
    worst = max(rows, key=lambda r: r[6])
    print(f"largest exact bound B: {worst[6]} = {float(worst[6]):.9f} (box {worst[0]})")
    print(f"largest |kappa_i|: {max(max(abs(t) for t in r[4]) for r in rows)}")
    print(f"exact bound < 1 on every box: {all(r[6] < 1 for r in rows)}")

    with open('dh_boxes.txt', 'w', newline='\n') as out:
        out.write("# dh_boxes.txt: the box table of Lemma 2.1(a) (hypothesis_H.md, Section 2).\n")
        out.write("# Written by cover_search_K.py; read and checked by 05_dh_cover.py. The proof uses only this table.\n")
        out.write("# D = [-1/2,1/2] x [-1/2,1/2] x [0,1/2] in the coordinates (a,b,c) of a + b*l + c*l^2, l^3 = 4l - 6.\n")
        out.write("# One line per box. The boxes are the leaves of a binary split tree of D, in depth-first order with\n")
        out.write("# the lower half first.\n")
        out.write("#   path   the splits from D to the box, one letter per split at the midpoint of an axis:\n")
        out.write("#          x, y, z = lower half along that axis, X, Y, Z = upper half\n")
        out.write(f"#   lo, hi the corners (lo_a lo_b lo_c hi_a hi_b hi_c), integers in units of 2^-{K}\n")
        out.write("#   s      1 or 2\n")
        out.write("#   kappa  an integer triple (k1 k2 k3); on the box |N(s*xi - kappa)| <= B < 1 with the Taylor bound B\n")
        out.write(f"unit 2^-{K}\n")
        out.write("# path lo_a lo_b lo_c hi_a hi_b hi_c s k1 k2 k3\n")
        for path, key, s, k, kap, Bf, Bx in rows:
            ends = ' '.join(str(t >> sh) for t in key[0] + key[1])
            out.write(f"{path} {ends} {s} {kap[0]} {kap[1]} {kap[2]}\n")
    print(f"wrote dh_boxes.txt ({len(rows)} boxes)")


if __name__ == '__main__':
    T0 = time.time()
    main()
    print(f"total time {time.time() - T0:.1f} s")
