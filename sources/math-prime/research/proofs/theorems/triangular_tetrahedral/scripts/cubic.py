"""Definitions shared by the scripts of hypothesis_H.md (no checks here).

A triple (a, b, c) stands for a + b*l + c*l^2, where l = theta/6 satisfies l^3 = 4l - 6.
So theta^3 - 144 theta + 1296 = 0 for theta = 6l, and X - theta is the triple (X, -6, 0).
Every function works on integers, on fractions and on sympy expressions.
The independent verifiers (verify_*.py) must not import this module.
"""


def mul(u, v):
    """Product of two triples, using l^3 = 4l - 6 and l^4 = 4l^2 - 6l."""
    a, b, c = u
    d, e, f = v
    return (a*d - 6*b*f - 6*c*e,
            a*e + b*d + 4*b*f + 4*c*e - 6*c*f,
            a*f + b*e + c*d + 4*c*f)


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
    return ((a, -6*c, -6*b),
            (b, a + 4*c, 4*b - 6*c),
            (c, b, a + 4*c))


def norm(u):
    """N(a,b,c) = det matrix(u) (checked in 04_ring.py)."""
    a, b, c = u
    return (a**3 + 8*a**2*c - 4*a*b**2 + 18*a*b*c + 16*a*c**2
            - 6*b**3 + 24*b*c**2 + 36*c**3)


def adj(u):
    """Adjugate triple: u * adj(u) = (N(u), 0, 0)."""
    a, b, c = u
    return (a*a + 8*a*c - 4*b*b + 6*b*c + 16*c*c,
            -a*b - 6*c*c,
            -a*c + b*b - 4*c*c)


def inverse_unit(u):
    """Inverse of a triple of norm +1 or -1."""
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
    """Return u/d as an integer triple if d divides u in Z[l], else None (d nonzero, integer triples)."""
    n = norm(d)
    q = mul(u, adj(d))
    if all(x % n == 0 for x in q):
        return tuple(x // n for x in q)
    return None


ONE = (1, 0, 0)
LAM = (0, 1, 0)            # l = theta/6
EPS = (-19, -5, 1)         # l^2 - 5l - 19, norm 1 (Avanesov's epsilon)
EPS_INV = (95, -101, 40)   # 40l^2 - 101l + 95 (Avanesov's epsilon_0)
PI2 = (-4, 1, 1)           # l^2 + l - 4, norm -2; PI2^3 = -2*EPS
P3A = (-1, 1, 0)           # l - 1, norm -3, vanishes at l = 1 mod 3
P3B = (5, 2, 0)            # 2l + 5, norm -3, vanishes at l = -1 mod 3
P3C = (3, -3, 1)           # l^2 - 3l + 3, norm 3, vanishes at l = 0 mod 3
P179 = (107, -114, 45)     # 45l^2 - 114l + 107, norm 179; PI2^2 * P179 = 3l^2 - 4


def x_minus_theta(X):
    """The triple of X - theta = X - 6l."""
    return (X, -6, 0)


def check(name, cond):
    """Print PASS or FAIL and stop at the first failure."""
    print(("PASS  " if cond else "FAIL  ") + name)
    if not cond:
        raise SystemExit("check failed: " + name)
