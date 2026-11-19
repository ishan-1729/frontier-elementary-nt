"""05_cover.py: the division lemma of hypothesis_H.md, Section 2, checked exactly from the table TABLE.

Lemma 2.1.
  (a) For every real triple xi in D = [-1/2,1/2] x [-1/2,1/2] x [0,1/2] (coordinates in the basis of hring.py)
      there are s in MULTS and an integer triple kappa with |N(s*xi - kappa)| < 1.
  (b) (only when 2 is a multiplier) For each of the 7 classes xi in (1/2)R with xi not in R there is an integer
      triple kappa with 0 < |N(xi - kappa)| < 1.
Every rational triple is brought into D by subtracting nearest integers and, when its third coordinate is negative,
by negating it (N(-x) = -N(x)). That reduction and Corollary 2.2 are written in hypothesis_H.md; 2.7 below runs them.

Checks (exact rational arithmetic; sympy for the polynomial identities):
  2.1 the Taylor identity N(v+w) = N(v) + sum_i g_i(v) w_i + sum_{i<=j} q_ij(v) w_i w_j + N(w);
  2.2 the table: each box lies in D, has dyadic ends, an allowed s and an integer kappa, and, with v = s*centre - kappa
      and half-widths h, B = |N(v)| + sum_i |g_i(v)| s h_i + sum_{i<=j} |q_ij(v)| s^2 h_i h_j + sum_m |n_m| (s h)^m < 1;
  2.3 the boxes cover D: they are exactly the leaves of a binary split tree of D, replayed from the paths, and each
      cell of the grid of the table's unit lies in exactly one box; the volumes add up to vol(D) = 1/2;
  2.4 part (b), when 2 is a multiplier;
  2.5 remark, not used by (H): a point XS fixed modulo R by a unit U (U*XS - XS in R) whose coset XS + R has no element
      of norm below BOUND; the proof is the orbit-minimum descent on T2 = sigma1^2 + 2Q under U^(+-1), followed by an
      exact enumeration. Also the fixed points of +-U modulo R and their least norms (exploration);
  2.6 remark: Bernstein enclosures of N(s*xi - kappa) on each box (a sharper second certificate);
  2.7 confirmation, not used by the proof: random divisions 0 < |N(s*alpha - kappa*beta)| < |N(beta)|.
The table was found by cover_search.py; the search is not part of the proof.
Shared with the copy of this script in the other card that proves (H); only the configuration block differs.
"""
import itertools
import math
import random
import time
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

import hring
from hring import mul, norm, adj, sub, neg, ONE, NORM_FORM, F_COEFFS, BASIS_POW, f
from hlib import (check, taylor_data, taylor_bound, bernstein_coeffs, read_table, box_from_path, replay, form_eval,
                  bisect_root, real_forms, iv_eval, Iv, ldl3, coset_enumerate, cube_root_up)

# ---- configuration (the only card-specific part) ------------------------------------------------------------------
TABLE = 'euclid_boxes.txt'
MULTS = (1,)
ROOT_BRACKET = (-2, -1)              # f changes sign on it; its root is the real root of f
UNIT = hring.EPS                     # a unit of norm 1 with 0 < sigma1(UNIT) < 1
UNIT_NAME = 'EPS'
FIX_SIGN = 1                         # XS is fixed by FIX_SIGN * UNIT modulo R
XS = (Fr(-1, 3), Fr(-2, 9), Fr(4, 9))
BOUND = Fr(8, 9)                     # the least |N| on XS + R
C_TEXT = Fr(676515, 10**5)           # the rounded constant C of hypothesis_H.md (4.3') (compared in 2.5.3)
# --------------------------------------------------------------------------------------------------------------------

T_START = time.time()
HALF = Fr(1, 2)
D_LO = (-HALF, -HALF, Fr(0))
D_HI = (HALF, HALF, HALF)
CUBIC = NORM_FORM
TD = taylor_data(CUBIC)
GRAD, QUAD = TD

print("2.1 Taylor identity")
x, y, z, w1, w2, w3 = sp.symbols('x y z w1 w2 w3')
V, W = (x, y, z), (w1, w2, w3)
taylor = (form_eval(CUBIC, V) + sum(form_eval(GRAD[i], V)*W[i] for i in range(3))
          + sum(form_eval(q, V)*W[i]*W[j] for (i, j), q in QUAD.items()) + form_eval(CUBIC, W))
check("N(v+w) = N(v) + sum g_i(v) w_i + sum_(i<=j) q_ij(v) w_i w_j + N(w), with g and q as printed below",
      sp.expand(norm(tuple(V[i] + W[i] for i in range(3))) - taylor) == 0
      and sp.expand(form_eval(CUBIC, W) - norm(W)) == 0)
for i in range(3):
    print(f"    g_{i + 1} =", sp.expand(form_eval(GRAD[i], V)))
print("    q_ij:", ", ".join(f"q{i + 1}{j + 1} = {sp.expand(form_eval(q, V))}" for (i, j), q in QUAD.items()))
print("    sum of |n_m| =", sum(abs(co) for co in CUBIC.values()))

print(f"2.2 the table {TABLE}")
UNIT_EXP, ROWS = read_table(TABLE)
check(f"the table has a unit line (2^-{UNIT_EXP}) and {len(ROWS)} boxes", UNIT_EXP is not None and len(ROWS) > 0)
UNIT_W = Fr(1, 2**UNIT_EXP)
BOXES = []
for n, (path, ends, s, kap) in enumerate(ROWS, 1):
    lo = tuple(e*UNIT_W for e in ends[:3])
    hi = tuple(e*UNIT_W for e in ends[3:])
    tag = f"box {n} ({path})"
    if not (path and all(ch in 'xyzXYZ' for ch in path)):
        check(f"{tag}: the path is a word in x, y, z, X, Y, Z", False)
    if box_from_path(path, D_LO, D_HI) != (lo, hi):
        check(f"{tag}: the listed corners are those of the path", False)
    if not all(D_LO[i] <= lo[i] < hi[i] <= D_HI[i] for i in range(3)):
        check(f"{tag}: lies in D", False)
    if s not in MULTS:
        check(f"{tag}: s in {MULTS}", False)
    B = taylor_bound(CUBIC, TD, lo, hi, s, kap)
    if not B < 1:
        check(f"{tag}: B = {B} < 1", False)
    BOXES.append((path, lo, hi, s, kap, B))
check(f"all {len(BOXES)} boxes: the path is a word in x,y,z,X,Y,Z and gives the listed corners, the box lies in D, its"
      f" ends are multiples of 2^-{UNIT_EXP}, s in {set(MULTS)}, kappa is an integer triple, and B < 1", True)
for s in MULTS:
    print(f"    boxes with s = {s}: {sum(1 for b in BOXES if b[3] == s)}")
worst = max(BOXES, key=lambda b: b[5])
print(f"    largest B: {worst[5]} = {float(worst[5]):.9f}, on box {worst[0]} (s = {worst[3]}, kappa = {worst[4]})")
print(f"    largest |kappa_i|: {max(max(abs(t) for t in b[4]) for b in BOXES)}; longest path:"
      f" {max(len(b[0]) for b in BOXES)} splits")
edges = [Fr(1, 2), Fr(9, 10), Fr(99, 100), Fr(999, 1000), Fr(1)]
counts, prev = [], Fr(0)
for e in edges:
    counts.append(f"[{prev}, {e}): {sum(1 for b in BOXES if prev <= b[5] < e)}")
    prev = e
print("    number of boxes by B:", "; ".join(counts))

print("2.3 the boxes cover D")
PATHS = [b[0] for b in BOXES]
check("no path occurs twice", len(set(PATHS)) == len(PATHS))
order = replay('', PATHS)
check(f"the {len(PATHS)} paths are exactly the leaves of a binary split tree of D: every node is either a box of the"
      " table or is halved at the midpoint of one axis into two nodes, and the file lists the leaves in depth-first"
      " order", order == PATHS)
print("    By induction from the leaves, the boxes below a node cover that node's box (a box that is halved is the")
print("    union of its two halves). So the boxes cover D, and two of them meet at most in boundary points.")
nx, nz = 2**UNIT_EXP, 2**(UNIT_EXP - 1)
cells = np.zeros((nx, nx, nz), dtype=np.int64)
off = (nx//2, nx//2, 0)
for path, ends, s, kap in ROWS:
    cells[ends[0] + off[0]:ends[3] + off[0], ends[1] + off[1]:ends[4] + off[1], ends[2]:ends[5]] += 1
check(f"second argument: each of the {nx}*{nx}*{nz} = {cells.size} cells of side 1/{nx} that tile D lies in exactly"
      " one box", bool((cells == 1).all()))
vol = sum((b[2][0] - b[1][0])*(b[2][1] - b[1][1])*(b[2][2] - b[1][2]) for b in BOXES)
check(f"the volumes of the boxes add up to vol(D) = 1/2 (sum = {vol})", vol == HALF)

HALF_LIST = {}
if 2 in MULTS:
    print("2.4 part (b): the 7 classes xi in (1/2)R with xi not in R")
    for e in itertools.product((0, 1), repeat=3):
        if e == (0, 0, 0):
            continue
        xi = tuple(Fr(t, 2) for t in e)
        found = None
        for kap in itertools.product(range(-4, 5), repeat=3):
            n = norm(sub(xi, kap))
            if n != 0 and abs(n) < 1:
                key = (max(abs(t) for t in kap), sum(abs(t) for t in kap), abs(n))
                if found is None or key < found[0]:
                    found = (key, kap, n)
        check(f"xi = {tuple(str(t) for t in xi)}: kappa = {found[1]}, N(xi - kappa) = {found[2]}, and 0 < |N| < 1",
              found is not None and 0 < abs(found[2]) < 1)
        HALF_LIST[e] = found[1]
else:
    print("2.4 (no part (b): the only multiplier is 1)")

print(f"2.5 remark (not used by (H)): the coset of XS = {tuple(str(t) for t in XS)} needs |N| >= {BOUND}")
FU = tuple(FIX_SIGN*t for t in UNIT)
FU_INV = tuple(FIX_SIGN*t for t in hring.inverse_unit(UNIT))
d_plus, d_minus = sub(mul(FU, XS), XS), sub(mul(FU_INV, XS), XS)
check(f"{'' if FIX_SIGN == 1 else '-'}{UNIT_NAME}*XS - XS = {tuple(str(t) for t in d_plus)} and"
      f" {'' if FIX_SIGN == 1 else '-'}{UNIT_NAME}^-1*XS - XS = {tuple(str(t) for t in d_minus)} are integer triples,"
      " so these units map the coset C = XS + R into itself", all(t.denominator == 1 for t in d_plus + d_minus))
print("  2.5.1 the real root l1 of f, and the two real forms sigma1 and Q (identities checked in 04_ring.py)")
lo1, hi1 = bisect_root(f, ROOT_BRACKET[0], ROOT_BRACKET[1], 80)
check(f"l1 in [{float(lo1):.15f}, {float(hi1):.15f}], width 2^-80 (exact sign change)",
      f(lo1)*f(hi1) <= 0 and hi1 - lo1 <= Fr(1, 2**80))
L = sp.Symbol('L')
sigma1, Q, s_, p_ = real_forms(F_COEFFS, BASIS_POW, L)
fL = sp.Poly(f(L), L)
lv = Iv(lo1, hi1)
e_iv = iv_eval(sigma1(UNIT), L, lv)
check(f"e := sigma1({UNIT_NAME}) lies in [{float(e_iv.lo):.12g}, {float(e_iv.hi):.12g}] and 0 < e < 1"
      " (hypothesis_H.md writes E := Q(U) = 1/e)", 0 < e_iv.lo and e_iv.hi < 1)
print("  2.5.2 the scaled descent lemma under the unit")
print(f"    N({UNIT_NAME}) = 1 = e Q({UNIT_NAME}) gives Q({UNIT_NAME}) = 1/e; for y with r = sigma1(y)^2 and q = Q(y):")
print(f"    T2({UNIT_NAME} y) = e^2 r + 2q/e and T2({UNIT_NAME}^-1 y) = r/e^2 + 2eq; the sign does not change T2.")
es, r_, q_, tt, Es = sp.symbols('e r q t E', positive=True)
check("T2(U y) - T2(y) = (1 - e)(2q/e - (1 + e) r) and T2(U^-1 y) - T2(y) = (1 - e)((1 + e) r/e^2 - 2q)",
      sp.simplify(es**2*r_ + 2*q_/es - r_ - 2*q_ - (1 - es)*(2*q_/es - (1 + es)*r_)) == 0
      and sp.simplify(r_/es**2 + 2*es*q_ - r_ - 2*q_ - (1 - es)*((1 + es)*r_/es**2 - 2*q_)) == 0)
rho_hi_s, rho_lo_s = 2/(es*(1 + es)), 2*es**2/(1 + es)
check("with rho_hi = 2/(e(1+e)) and rho_lo = 2e^2/(1+e) (that is 2E^2/(E+1) and 2/(E(E+1)) for E = 1/e):"
      " r > rho_hi q gives T2(U y) < T2(y), and r < rho_lo q gives T2(U^-1 y) < T2(y)",
      sp.simplify((1 + es)*rho_hi_s - 2/es) == 0 and sp.simplify((1 + es)*rho_lo_s/es**2 - 2) == 0
      and sp.simplify(rho_hi_s.subs(es, 1/Es) - 2*Es**2/(Es + 1)) == 0
      and sp.simplify(rho_lo_s.subs(es, 1/Es) - 2/(Es*(Es + 1))) == 0)
psi = (tt + 2)**3/tt
check("psi(t) = (t+2)^3/t has psi'(t) = 2(t+2)^2(t-1)/t^2, and T2^3/N^2 = (r + 2q)^3/(r q^2) = psi(r/q)",
      sp.simplify(sp.diff(psi, tt) - 2*(tt + 2)**2*(tt - 1)/tt**2) == 0
      and sp.simplify((r_ + 2*q_)**3/(r_*q_**2) - psi.subs(tt, r_/q_)) == 0)
K3_e = 4*(1 + es + es**2)**3/(es**2*(1 + es)**2)
check("psi(rho_hi) = psi(rho_lo) = 4(1 + e + e^2)^3/(e^2(1 + e)^2) = 4(E^2 + E + 1)^3/(E^2(E + 1)^2) with E = 1/e",
      sp.simplify(psi.subs(tt, rho_hi_s) - K3_e) == 0 and sp.simplify(psi.subs(tt, rho_lo_s) - K3_e) == 0
      and sp.simplify(K3_e.subs(es, 1/Es) - 4*(Es**2 + Es + 1)**3/(Es**2*(Es + 1)**2)) == 0)
e_lo = e_iv.lo
RH = 2/(e_lo*(1 + e_lo))
RL = 2*e_lo**2/(1 + e_lo)


def psi_q(tv):
    return (tv + 2)**3/tv


K3 = max(psi_q(RL), psi_q(RH))
check(f"rho_hi <= RH := 2/(e-(1+e-)) = {float(RH):.6f} and rho_lo >= RL := 2e-^2/(1+e-) = {float(RL):.6e} (rho_hi"
      " decreases and rho_lo increases with e; e- is the lower end of the enclosure of e), and RL < 1 < RH",
      RL < 1 < RH)
check(f"K3 := max(psi(RL), psi(RH)) = psi(RL) = psi(RH) = 4(1 + e- + e-^2)^3/(e-^2(1 + e-)^2) = {float(K3):.6f}",
      psi_q(RL) == psi_q(RH) == 4*(1 + e_lo + e_lo**2)**3/(e_lo**2*(1 + e_lo)**2))
print("    psi decreases on (0, 1] and increases on [1, oo), so psi <= K3 on [RL, RH], which contains [rho_lo, rho_hi].")
print("    Lemma: for every y with N(y) != 0, T2(U y) < T2(y), or T2(U^-1 y) < T2(y), or T2(y)^3 <= K3 N(y)^2.")
print("    For y in C the set {U^k y : T2(U^k y) <= T2(y)} lies in C and is finite (T2 is positive definite and C is a")
print("    translate of the lattice R) and nonempty; its element y0 of least T2 has T2(U y0) >= T2(y0) and")
print("    T2(U^-1 y0) >= T2(y0), so T2(y0)^3 <= K3 N(y0)^2 = K3 N(y)^2.")
print("  2.5.3 a rational lower form for T2 and an exact enumeration of C")
T2 = sp.expand(sigma1((sp.Symbol('a'), sp.Symbol('b'), sp.Symbol('c')))**2
               + 2*Q((sp.Symbol('a'), sp.Symbol('b'), sp.Symbol('c'))))
syms = sp.symbols('a b c')
gram = [[sp.rem(sp.Poly(sp.expand(sp.Rational(1, 2)*sp.diff(T2, syms[k], syms[m])), L), fL).as_expr()
         for m in range(3)] for k in range(3)]
G_iv = [[iv_eval(gram[k][m], L, lv) for m in range(3)] for k in range(3)]
RHO = Fr(1, 2**20)
G_mid = [[Fr(round((G_iv[k][m].lo + G_iv[k][m].hi)/2*2**20), 2**20) for m in range(3)] for k in range(3)]
check("each entry of the Gram matrix G of T2 lies within rho = 2^-20 of the rational matrix G' (entries rounded)",
      all(G_mid[k][m] - RHO <= G_iv[k][m].lo and G_iv[k][m].hi <= G_mid[k][m] + RHO for k in range(3)
          for m in range(3)))
G_low = [[G_mid[k][m] - (3*RHO if k == m else 0) for m in range(3)] for k in range(3)]
dd = ldl3(G_low)
check(f"Glow := G' - 3 rho I has positive LDL pivots {[f'{float(t):.6f}' for t in dd[:3]]} (exact), and"
      " T2(y) >= Glow(y), since T2(y) - y^T G' y >= -rho (sum |y_i|)^2 >= -3 rho sum y_i^2", min(dd[:3]) > 0)
T_e = cube_root_up(K3*BOUND**2, 1000)
check(f"T_e = {T_e} = {float(T_e)} satisfies T_e^3 >= K3 BOUND^2", T_e**3 >= K3*BOUND**2)
fits = C_TEXT**3*BOUND**2 <= T_e**3
print(f"    with the rounded C = {float(C_TEXT)} of (4.3') in place of K3: C^3 BOUND^2 = {float(C_TEXT**3*BOUND**2):.4f}"
      f" {'<=' if fits else '>'} T_e^3 = {float(T_e**3):.4f}, so the rounded C {'would also' if fits else 'would not'}"
      " fit; the proof uses K3")
listC = coset_enumerate(G_low, XS, T_e)
check("every listed y is in C and has Glow(y) <= T_e",
      all(all((yv - xv).denominator == 1 for yv, xv in zip(yy, XS)) for yy in listC)
      and all(sum(G_low[k][m]*yy[k]*yy[m] for k in range(3) for m in range(3)) <= T_e for yy in listC))
normsC = sorted((abs(norm(yy)), yy) for yy in listC)
minC = normsC[0][0]
check(f"the {len(listC)} elements y of C with Glow(y) <= T_e have min |N(y)| = {minC} = BOUND, attained at"
      f" {[tuple(str(t) for t in yy) for n, yy in normsC if n == minC][:4]}", minC == BOUND)
print("    the next values of |N| on the list:", ", ".join(str(n) for n, _ in normsC[1:6]))
print(f"    Conclusion. If y in C had |N(y)| < {BOUND}, the element y0 of least T2 in its orbit (inside C) would have")
print(f"    Glow(y0)^3 <= T2(y0)^3 <= K3 N(y)^2 < K3 BOUND^2 <= T_e^3, so y0 would be on the list, where |N| >= BOUND:")
print(f"    impossible. So min over C of |N| = {BOUND}, and no kappa in R has |N(XS - kappa)| < {BOUND}.")
if BOUND >= 1:
    print("    Since BOUND > 1, XS has no translate of norm below 1: R is not norm-Euclidean, and the multiplier 2 is")
    print(f"    needed near XS. Indeed |N(2 XS - kappa)| < 1 for some kappa: min over |kappa_i| <= 6 is"
          f" {min(abs(norm(sub(tuple(2*t for t in XS), k))) for k in itertools.product(range(-6, 7), repeat=3))}.")
print("  2.5.4 the points fixed by +-U modulo R (exploration)")
T_1 = cube_root_up(K3, 1000)
for sg in (1, -1):
    Um = tuple(sg*t for t in UNIT)
    D1 = sub(Um, ONE)
    ND1 = norm(D1)
    gens = [tuple(Fr(t, ND1) for t in mul(adj(D1), ei)) for ei in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    gens = [tuple(t - math.floor(t) for t in gv) for gv in gens]
    classes = {(Fr(0), Fr(0), Fr(0))}
    frontier = list(classes)
    while frontier:
        new = []
        for v in frontier:
            for gv in gens:
                wv = tuple((p + q) - math.floor(p + q) for p, q in zip(v, gv))
                if wv not in classes:
                    classes.add(wv)
                    new.append(wv)
        frontier = new
    check(f"{'+' if sg == 1 else '-'}{UNIT_NAME}: {len(classes)} = |N({'' if sg == 1 else '-'}{UNIT_NAME} - 1)| classes"
          " of fixed points modulo Z^3", len(classes) == abs(ND1)
          and all(all(t.denominator == 1 for t in mul(D1, v)) for v in classes))
    mins = []
    for v in sorted(classes):
        if v == (0, 0, 0):
            continue
        lst = coset_enumerate(G_low, v, T_1)
        mins.append((min(abs(norm(yy)) for yy in lst), v))
    mins.sort(reverse=True)
    print(f"    the largest least |N| over the classes, from the elements with Glow <= T_1 = {T_1}. By the orbit argument")
    print("    with 1 in place of BOUND, a value below 1 is the exact least |N| of its class, and a class whose listed")
    print("    least is >= 1 has no element with |N| < 1 (its exact least is not computed here):")
    for m in mins[:4]:
        shown = f"{str(m[0]):>10s} = {float(m[0]):.6f}" if m[0] < 1 else f"{'>= 1':>21s}"
        print(f"        {shown}   class of {tuple(str(t) for t in m[1])}")
    below = [m[0] for m in mins if m[0] < 1]
    print(f"    nonzero classes with no element of norm below 1: {len(mins) - len(below)}; the other {len(below)} have"
          f" least |N| at most {max(below)}" + (f" (all values: {', '.join(str(v) for v in below)})"
                                                 if len(mins) <= 12 else ""))

print("2.6 remark: Bernstein enclosures (a second, sharper certificate; not needed for the proof)")
overall, corner_ok, bern_max = None, True, []
for path, lo, hi, s, kap, B in BOXES:
    bc = bernstein_coeffs(CUBIC, TD, lo, hi, s, kap)
    for e1, e2, e3 in itertools.product((0, 1), repeat=3):
        corner = (lo[0] if e1 == 0 else hi[0], lo[1] if e2 == 0 else hi[1], lo[2] if e3 == 0 else hi[2])
        if bc[3*e1][3*e2][3*e3] != norm(tuple(s*corner[i] - kap[i] for i in range(3))):
            corner_ok = False
    mx, where = max((abs(bc[I][J][K]), (I, J, K)) for I, J, K in itertools.product(range(4), repeat=3))
    bern_max.append(mx)
    if overall is None or mx > overall[0]:
        overall = (mx, where, path, lo, hi, s, kap)
check("on every box the corner Bernstein coefficients equal the values of N(s*xi - kappa) at the corners", corner_ok)
check(f"on every box max |b_IJK| < 1; the largest over all boxes is {overall[0]} = {float(overall[0]):.6f}",
      all(m < 1 for m in bern_max))
mx, (I, J, K), path, lo, hi, s, kap = overall
if all(t in (0, 3) for t in (I, J, K)):
    corner = tuple(lo[k] if e == 0 else hi[k] for k, e in enumerate((I//3, J//3, K//3)))
    print(f"    it is the corner coefficient at xi = {tuple(str(t) for t in corner)} of box {path} (s = {s},"
          f" kappa = {kap}), so the largest value of |N(s*xi - kappa)| over the boxes, each with its own (s, kappa) on")
    print(f"    its own closed box, is {mx} exactly.")
    others = [(p2, s2, k2, norm(tuple(s2*corner[i] - k2[i] for i in range(3))))
              for p2, lo2, hi2, s2, k2, B2 in BOXES
              if p2 != path and all(lo2[i] <= corner[i] <= hi2[i] for i in range(3))]
    print(f"    the corner lies in {len(others)} other boxes, with N(s*xi - kappa) ="
          f" {', '.join(str(v) for _, _, _, v in others)}")
else:
    print(f"    it is the coefficient {(I, J, K)} of box {path}: with the table's choices |N(s*xi - kappa)| <= {mx} on D.")

print("2.7 confirmation (not used by the proof): random divisions")
LEAF = {b[0]: (b[3], b[4]) for b in BOXES}
AXIS_OF = {}
for p in PATHS:
    for k in range(len(p)):
        AXIS_OF[p[:k]] = p[k].lower()


def locate(xi):
    node, lo, hi = '', list(D_LO), list(D_HI)
    while node not in LEAF:
        ch = AXIS_OF[node]
        ax = 'xyz'.index(ch)
        mid = (lo[ax] + hi[ax])/2
        if xi[ax] <= mid:
            hi[ax], node = mid, node + ch
        else:
            lo[ax], node = mid, node + ch.upper()
    return LEAF[node]


def div_step(alpha, beta):
    """s and kappa with 0 < |N(s alpha - kappa beta)| < |N(beta)|, as in Corollary 2.2; None if beta | alpha."""
    nb = norm(beta)
    xi = tuple(Fr(t, nb) for t in mul(alpha, adj(beta)))
    if all(t.denominator == 1 for t in xi):
        return None
    if HALF_LIST and all((2*t).denominator == 1 for t in xi):
        e = tuple(int(2*t) % 2 for t in xi)
        kb = HALF_LIST[e]
        return 1, tuple(int(xi[k] - Fr(e[k], 2)) + kb[k] for k in range(3)), 'b'
    r = tuple(math.floor(t + HALF) for t in xi)
    dlt = tuple(xi[k] - r[k] for k in range(3))
    sg = -1 if dlt[2] < 0 else 1
    s, kap = locate(tuple(sg*t for t in dlt))
    return s, tuple(s*r[k] + sg*kap[k] for k in range(3)), 'a'


rng = random.Random(20261015)
stats, ok, ntests = {}, True, 0
for _ in range(3000):
    beta = tuple(rng.randint(-20, 20) for _ in range(3))
    alpha = tuple(rng.randint(-60, 60) for _ in range(3))
    if beta == (0, 0, 0):
        continue
    ntests += 1
    res = div_step(alpha, beta)
    if res is None:
        stats['beta | alpha'] = stats.get('beta | alpha', 0) + 1
        continue
    s, kappa, how = res
    rem_n = abs(norm(sub(tuple(s*t for t in alpha), mul(kappa, beta))))
    if not 0 < rem_n < abs(norm(beta)):
        ok = False
    stats[f'{how}, s = {s}'] = stats.get(f'{how}, s = {s}', 0) + 1
check(f"{ntests} random pairs: 0 < |N(s alpha - kappa beta)| < |N(beta)| whenever beta does not divide alpha ({stats})",
      ok)
print(f"all checks of 05_cover.py passed ({time.time() - T_START:.1f} s)")
