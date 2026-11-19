"""hlib.py: helpers shared by the scripts of hypothesis_H.md (no checks of their own; a verbatim copy is in the other
card that uses it). Everything is exact: integers, Fractions and sympy. Floating point appears only in printing and in
the search scripts, never in a check.

Contents
- check(name, cond): print PASS/FAIL and stop at the first failure.
- Iv: closed intervals with rational ends (+, -, *, with enclosure of the exact result).
- bisect_root(f, lo, hi, n): a rational enclosure of a root of f from an exact sign change.
- Cubic forms as dicts {(i,j,k): coefficient}: evaluation, the Taylor data (gradient and second-order coefficients),
  the Taylor bound of a box, and the Bernstein coefficients of a box.
- Box tables: reading, replaying the split tree, and the cell count.
- ldl3 and coset_enumerate: all points of a coset z0 + Z^3 with y^T G y <= T, for a rational positive definite G.
"""
import itertools
import math
from fractions import Fraction as Fr


def check(name, cond):
    """Print PASS or FAIL and stop at the first failure."""
    print(("PASS  " if cond else "FAIL  ") + name)
    if not cond:
        raise SystemExit("check failed: " + name)


class Iv:
    """A closed interval with rational ends; +, - and * enclose the exact results."""

    def __init__(self, lo, hi=None):
        self.lo, self.hi = Fr(lo), Fr(lo if hi is None else hi)
        assert self.lo <= self.hi

    @staticmethod
    def of(o):
        return o if isinstance(o, Iv) else Iv(o)

    def __add__(self, o):
        o = Iv.of(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + (-Iv.of(o))

    def __rsub__(self, o):
        return Iv.of(o) - self

    def __mul__(self, o):
        o = Iv.of(o)
        p = (self.lo*o.lo, self.lo*o.hi, self.hi*o.lo, self.hi*o.hi)
        return Iv(min(p), max(p))

    __rmul__ = __mul__

    def __repr__(self):
        return f"[{float(self.lo):.12g}, {float(self.hi):.12g}]"


def bisect_root(f, lo, hi, n):
    """[lo, hi] with f(lo) and f(hi) of opposite signs, halved n times exactly; returns the final (lo, hi)."""
    lo, hi = Fr(lo), Fr(hi)
    flo = f(lo)
    assert flo*f(hi) < 0
    for _ in range(n):
        mid = (lo + hi)/2
        fm = f(mid)
        if fm == 0:
            return mid, mid
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return lo, hi


# ---- cubic forms --------------------------------------------------------------------------------------------------
MONO3 = [m for m in itertools.product(range(4), repeat=3) if sum(m) == 3]


def form_eval(form, v):
    """Value of the form {(i,j,k): c} at v."""
    return sum(c*v[0]**m[0]*v[1]**m[1]*v[2]**m[2] for m, c in form.items())


def taylor_data(cubic):
    """The gradient (three quadratic forms) and the second-order coefficients q_ij, i <= j (six linear forms) of a
    cubic form, as dicts, so that N(v + w) = N(v) + sum_i g_i(v) w_i + sum_{i<=j} q_ij(v) w_i w_j + N(w)."""
    grad = [dict() for _ in range(3)]
    for m, c in cubic.items():
        for i in range(3):
            if m[i]:
                mm = list(m)
                mm[i] -= 1
                grad[i][tuple(mm)] = grad[i].get(tuple(mm), 0) + c*m[i]
    quad = {}
    for i in range(3):
        for j in range(i, 3):
            q = {}
            for m, c in cubic.items():
                mm = list(m)
                if i == j:
                    if m[i] >= 2:
                        mm[i] -= 2
                        q[tuple(mm)] = q.get(tuple(mm), 0) + Fr(c*m[i]*(m[i] - 1), 2)
                else:
                    if m[i] and m[j]:
                        mm[i] -= 1
                        mm[j] -= 1
                        q[tuple(mm)] = q.get(tuple(mm), 0) + c*m[i]*m[j]
            quad[(i, j)] = q
    return grad, quad


def taylor_bound(cubic, tdata, lo, hi, s, kap):
    """B = |N(v)| + sum |g_i(v)| s h_i + sum |q_ij(v)| s^2 h_i h_j + sum |n_m| (s h)^m for the box [lo, hi], with
    v = s*centre - kappa and h the half-widths: then |N(s*xi - kappa)| <= B on the box (exact)."""
    grad, quad = tdata
    cen = [(lo[i] + hi[i])/2 for i in range(3)]
    H = [s*(hi[i] - lo[i])/2 for i in range(3)]
    v = tuple(s*cen[i] - kap[i] for i in range(3))
    return (abs(form_eval(cubic, v)) + sum(abs(form_eval(grad[i], v))*H[i] for i in range(3))
            + sum(abs(form_eval(q, v))*H[i]*H[j] for (i, j), q in quad.items())
            + sum(abs(c)*H[0]**m[0]*H[1]**m[1]*H[2]**m[2] for m, c in cubic.items()))


TB = [[Fr(math.comb(I, i), math.comb(3, i)) if i <= I else Fr(0) for i in range(4)] for I in range(4)]


def bernstein_coeffs(cubic, tdata, lo, hi, s, kap):
    """Bernstein coefficients b_IJK of P(u) = N(s*(lo + (hi - lo) u) - kappa) on [0,1]^3 (exact)."""
    grad, quad = tdata
    A = tuple(s*lo[i] - kap[i] for i in range(3))
    Dl = tuple(s*(hi[i] - lo[i]) for i in range(3))
    mono = {(0, 0, 0): form_eval(cubic, A)}
    for i in range(3):
        ei = [0, 0, 0]
        ei[i] = 1
        mono[tuple(ei)] = form_eval(grad[i], A)*Dl[i]
    for (i, j), q in quad.items():
        eij = [0, 0, 0]
        eij[i] += 1
        eij[j] += 1
        mono[tuple(eij)] = mono.get(tuple(eij), 0) + form_eval(q, A)*Dl[i]*Dl[j]
    for m, c in cubic.items():
        mono[m] = mono.get(m, 0) + c*Dl[0]**m[0]*Dl[1]**m[1]*Dl[2]**m[2]
    arr = [[[mono.get((i, j, k), 0) for k in range(4)] for j in range(4)] for i in range(4)]
    arr = [[[sum(TB[I][i]*arr[i][j][k] for i in range(I + 1)) for k in range(4)] for j in range(4)] for I in range(4)]
    arr = [[[sum(TB[J][j]*arr[I][j][k] for j in range(J + 1)) for k in range(4)] for J in range(4)] for I in range(4)]
    arr = [[[sum(TB[K][k]*arr[I][J][k] for k in range(K + 1)) for K in range(4)] for J in range(4)] for I in range(4)]
    return arr


# ---- box tables ---------------------------------------------------------------------------------------------------
AXES = 'xyz'


def box_from_path(path, d_lo, d_hi):
    lo, hi = list(d_lo), list(d_hi)
    for ch in path:
        ax = AXES.index(ch.lower())
        mid = (lo[ax] + hi[ax])/2
        if ch.islower():
            hi[ax] = mid
        else:
            lo[ax] = mid
    return tuple(lo), tuple(hi)


def read_table(fname):
    """(unit exponent, rows) with rows (path, ends (6 ints), s, kappa (3 ints))."""
    unit_exp, rows = None, []
    with open(fname, encoding='ascii') as fh:
        for line in fh:
            t = line.split()
            if not t or t[0].startswith('#'):
                continue
            if t[0] == 'unit':
                unit_exp = int(t[1].split('^-')[1])
                continue
            assert len(t) == 11, line
            rows.append((t[0], tuple(int(u) for u in t[1:7]), int(t[7]), tuple(int(u) for u in t[8:11])))
    return unit_exp, rows


def replay(prefix, paths):
    """The leaves below the node `prefix`, in depth-first order, after checking that they form a split tree."""
    if paths == [prefix]:
        return [prefix]
    if not paths:
        check(f"node '{prefix}' carries at least one box", False)
    if prefix in paths:
        check(f"the box '{prefix}' has no other box inside it", False)
    axes = {p[len(prefix)].lower() for p in paths}
    if len(axes) != 1:
        check(f"node '{prefix}' is split along one axis only (found {sorted(axes)})", False)
    ax = axes.pop()
    lower = [p for p in paths if p[len(prefix)] == ax]
    upper = [p for p in paths if p[len(prefix)] == ax.upper()]
    return replay(prefix + ax, lower) + replay(prefix + ax.upper(), upper)


# ---- lattice enumeration ------------------------------------------------------------------------------------------
def ldl3(G):
    """For a symmetric rational 3x3 matrix G: (d1, d2, d3, m12, m13, m23) with
    y^T G y = d1 (y1 + m12 y2 + m13 y3)^2 + d2 (y2 + m23 y3)^2 + d3 y3^2."""
    d1 = G[0][0]
    m12, m13 = G[0][1]/d1, G[0][2]/d1
    d2 = G[1][1] - d1*m12**2
    m23 = (G[1][2] - d1*m12*m13)/d2
    d3 = G[2][2] - d1*m13**2 - d2*m23**2
    return d1, d2, d3, m12, m13, m23


def coset_range(x0, alpha, beta, gamma):
    """All x in x0 + Z with alpha*(x + beta)^2 <= gamma (alpha > 0). The float estimate of the range is certified:
    the scan stops at points that fail the test, one on each side of the centre -beta."""
    if gamma < 0:
        return []
    centre, rad = -beta, math.sqrt(float(gamma/alpha))
    n_lo = math.floor(float(centre - x0) - rad) - 1
    n_hi = math.ceil(float(centre - x0) + rad) + 1
    while not (x0 + n_lo < centre and alpha*(x0 + n_lo + beta)**2 > gamma):
        n_lo -= 1
    while not (x0 + n_hi > centre and alpha*(x0 + n_hi + beta)**2 > gamma):
        n_hi += 1
    return [x0 + n for n in range(n_lo + 1, n_hi) if alpha*(x0 + n + beta)**2 <= gamma]


def coset_enumerate(G, z0, T):
    """All y in z0 + Z^3 with y^T G y <= T, for a rational positive definite G (exact)."""
    d1, d2, d3, m12, m13, m23 = ldl3(G)
    assert min(d1, d2, d3) > 0
    out = []
    for y3 in coset_range(z0[2], d3, 0, T):
        r3 = T - d3*y3**2
        for y2 in coset_range(z0[1], d2, m23*y3, r3):
            r2 = r3 - d2*(y2 + m23*y3)**2
            for y1 in coset_range(z0[0], d1, m12*y2 + m13*y3, r2):
                out.append((y1, y2, y3))
    return out


def cube_root_up(X, scale):
    """The least k/scale with (k/scale)^3 >= X (X > 0 rational)."""
    k = int(round(float(X)**(1/3)*scale)) - 2
    while Fr(k, scale)**3 < X:
        k += 1
    return Fr(k, scale)


# ---- the real place ---------------------------------------------------------------------------------------------
def real_forms(f_coeffs, basis_pow, L):
    """For f(t) = t^3 + f2 t^2 + f1 t + f0 with (f0, f1, f2) = f_coeffs and a symbol L standing for its real root,
    and a basis given by its coordinates in powers of t: functions sigma1(u) and Q(u) (sympy expressions in L), with
    the other two roots z, zbar the roots of t^2 - s t + p, s = -L - f2, p = L^2 + f2 L + f1 (so that
    f(t) = (t - L)(t^2 - s t + p) + f(L)). Q(u) = |u(z)|^2 = A^2 + s A B + p B^2, where u = A + B t mod t^2 - s t + p."""
    f0, f1, f2 = f_coeffs
    s = -L - f2
    p = L**2 + f2*L + f1

    def powc(u):
        return [sum(u[i]*basis_pow[i][k] for i in range(3)) for k in range(3)]

    def sigma1(u):
        c = powc(u)
        return c[0] + c[1]*L + c[2]*L**2

    def Q(u):
        c = powc(u)
        A, B = c[0] - c[2]*p, c[1] + c[2]*s
        return A**2 + s*A*B + p*B**2

    return sigma1, Q, s, p


def iv_eval(expr, L, lv):
    """Enclose a polynomial expression in the symbol L over the interval lv (Horner with intervals; exact ends)."""
    import sympy as sp
    poly = sp.Poly(sp.expand(expr), L)
    acc = Iv(0)
    for co in poly.all_coeffs():
        acc = acc*lv + Iv(Fr(int(sp.fraction(co)[0]), int(sp.fraction(co)[1])))
    return acc
