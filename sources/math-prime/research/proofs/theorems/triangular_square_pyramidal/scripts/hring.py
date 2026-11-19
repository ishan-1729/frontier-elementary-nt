"""hring.py: the ring R = Z[l], l^3 = l - 3, of hypothesis_H.md (definitions only; the checks are in 04_ring.py).

A triple (a, b, c) stands for a + b*l + c*l^2, where l = theta/3 satisfies l^3 = l - 3. So theta = 3l satisfies
theta^3 - 9 theta + 81 = 0, the cubic of E: y^2 = x^3 - 9x + 81 (8604.a1), and X - theta is the triple (X, -3, 0).
Every function works on integers, on Fractions and on sympy expressions.
The independent verifiers (verify_H_*.py) must not import this module.
"""
from hlib import check  # noqa: F401  (re-exported for the scripts)

NAME = "R = Z[l], l^3 = l - 3"
F_COEFFS = (3, -1, 0)                # f(t) = t^3 + 0 t^2 - t + 3: (f0, f1, f2) with f = t^3 + f2 t^2 + f1 t + f0
BASIS_POW = ((1, 0, 0), (0, 1, 0), (0, 0, 1))   # the basis 1, l, l^2 in powers of l


def f(t):
    return t**3 - t + 3


def mul(u, v):
    """Product of two triples, using l^3 = l - 3 and l^4 = l^2 - 3l."""
    a, b, c = u
    d, e, g = v
    return (a*d - 3*b*g - 3*c*e,
            a*e + b*d + b*g + c*e - 3*c*g,
            a*g + b*e + c*d + c*g)


def add(u, v):
    return tuple(x + y for x, y in zip(u, v))


def sub(u, v):
    return tuple(x - y for x, y in zip(u, v))


def neg(u):
    return tuple(-x for x in u)


def scal(r, u):
    return tuple(r*x for x in u)


def matrix(u):
    """Columns are u*1, u*l, u*l^2 (rows: coordinates)."""
    a, b, c = u
    return ((a, -3*c, -3*b),
            (b, a + c, b - 3*c),
            (c, b, a + c))


NORM_FORM = {(3, 0, 0): 1, (2, 0, 1): 2, (1, 2, 0): -1, (1, 1, 1): 9, (1, 0, 2): 1, (0, 3, 0): -3, (0, 1, 2): 3,
             (0, 0, 3): 9}


def norm(u):
    """N(a,b,c) = det matrix(u) = a^3 + 2a^2c - ab^2 + 9abc + ac^2 - 3b^3 + 3bc^2 + 9c^3 (checked in 04_ring.py)."""
    a, b, c = u
    return a**3 + 2*a**2*c - a*b**2 + 9*a*b*c + a*c**2 - 3*b**3 + 3*b*c**2 + 9*c**3


def trace(u):
    a, b, c = u
    return 3*a + 2*c


def adj(u):
    """Adjugate triple: u * adj(u) = (N(u), 0, 0)."""
    a, b, c = u
    return (a*a + 2*a*c - b*b + 3*b*c + c*c,
            -a*b - 3*c*c,
            -a*c + b*b - c*c)


def inverse_unit(u):
    n = norm(u)
    assert n in (1, -1)
    return scal(n, adj(u))


def power(u, k):
    """u^k for an integer k; negative k only for units."""
    if k < 0:
        u = inverse_unit(u)
        k = -k
    r = (1, 0, 0)
    for _ in range(k):
        r = mul(r, u)
    return r


def divides(d, u):
    """u/d as an integer triple if d divides u in R, else None (d nonzero; integer triples)."""
    n = norm(d)
    q = mul(u, adj(d))
    if all(x % n == 0 for x in q):
        return tuple(x // n for x in q)
    return None


ONE = (1, 0, 0)
LAM = (0, 1, 0)             # l = theta/3
THETA = (0, 3, 0)
EPS = (-1, 1, 1)            # l^2 + l - 1, norm 1 (Uchiyama's unit)
EPS_INV = (2, -2, 1)        # l^2 - 2l + 2
P3A = (0, 1, 0)             # l, norm -3
P3B = (-1, 1, 0)            # l - 1, norm -3
P3C = (1, 1, 0)             # l + 1, norm -3; P3A*P3B*P3C = l^3 - l = -3
P239 = (-1, 0, 3)           # 3l^2 - 1, norm 239; f'(theta) = 3 theta^2 - 9 = 9*P239


def x_minus_theta(X):
    """The triple of X - theta = X - 3l."""
    return (X, -3, 0)


def cofactor(X):
    """The triple of X^2 + theta X + theta^2 - 9, so that (X - theta) * cofactor = X^3 - 9X + 81."""
    return (X*X - 9, 3*X, 9)
