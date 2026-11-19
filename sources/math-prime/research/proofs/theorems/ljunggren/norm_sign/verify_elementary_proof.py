"""Exact arithmetic checks supporting elementary_proof.md.

Run with Python 3 and SymPy. The document proves its polynomial lemmas for
all N; the finite checks here are supplementary, not a substitute for that
proof, and do not invoke any Diophantine-equation solver.
"""

from fractions import Fraction
from math import comb, factorial, prod

import sympy as sp


def quartic(u, v):
    return u**4 - 956*u**3*v - 6*u*u*v*v + 956*u*v**3 + v**4


def generalized_binomial(a, k):
    return prod((a-j for j in range(k)), start=Fraction(1)) / factorial(k)


def gaussian_mul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def gaussian_power(a, k):
    result = (1, 0)
    for _ in range(k):
        result = gaussian_mul(result, a)
    return result


def approximant(n):
    real = imag = 0
    for k in range(n+1):
        coefficient = (
            (-1)**k * 4**n * 2**k
            * generalized_binomial(Fraction(n)-Fraction(1, 4), k)
            * comb(2*n-k, n)
        )
        assert coefficient.denominator == 1
        z = gaussian_mul(gaussian_power((0, 1), k),
                         gaussian_power((239, 1), n-k))
        real += coefficient.numerator*z[0]
        imag += coefficient.numerator*z[1]
    return real, imag


def main():
    m, n, u, v, w, x = sp.symbols("m n u v w x")
    q = m**4+4*m**3*n-6*m*m*n*n-4*m*n**3+n**4
    assert sp.expand(quartic(3*m+2*n, -2*m+3*n)-13**4*q) == 0
    assert sp.expand((3*m+2*n)**2+(-2*m+3*n)**2
                     -13*(m*m+n*n)) == 0
    assert sp.expand(quartic(-v, u)-quartic(u, v)) == 0
    assert sp.expand(quartic(u+v, v-u)+4*quartic(u, v)) == 0
    assert sp.expand(quartic(u-v, u+v)+4*quartic(u, v)) == 0
    for m0 in range(16):
        for n0 in range(16):
            if (m0-n0) % 2:
                assert int(q.subs({m: m0, n: n0})) % 16 == 1
    print("PASS: quartic transforms, norms, and opposite-parity congruence")

    # Direct continuation of the screenshot's r,s,t,u branch.
    d0, e0, r0, h0 = sp.symbols("d0 e0 r0 h0")
    inherited_square = (5*r0**4+28*r0**3*h0+60*r0*r0*h0*h0
                        +56*r0*h0**3+20*h0**4)
    centered_leg = 169*e0*e0+478*d0*e0+338*d0*d0
    assert sp.expand(169*inherited_square.subs({r0: e0+2*d0, h0: e0+d0})
                     -centered_leg**2-4*d0*d0*e0*e0) == 0
    assert sp.expand((u*u-v*v-478*u*v)**2
                     -8*169**2*u*u*v*v-quartic(u, v)) == 0
    assert sp.expand((e0*e0+2*d0*d0)**2-2*(2*d0*e0)**2
                     -(e0*e0-2*d0*d0)**2) == 0
    print("PASS: centered Vierzahlensatz continuation and coefficient-956 form")

    r = Fraction(5, 12)
    slope = -956+12*r+2868*r*r+4*r**3
    assert slope == Fraction(-195607, 432) and slope < -450
    # sqrt(2)-1 < 5/12, checked by squaring positive quantities.
    assert Fraction(17, 12)**2 > 2
    assert Fraction(114244, 450) < 254
    print("PASS: rational divided-difference and quartic-error constants")

    # Sharper positive narrow-cone closure, independent of the root bracket.
    assert -956+Fraction(12, 10)+Fraction(2868, 100)+Fraction(4, 1000) == -Fraction(231529, 250)
    assert -Fraction(231529, 250) < -900
    assert Fraction(28561, 900) < 32
    assert 4*64**10*4000**20 == 200**21*51200**9
    assert 956*53 < 51200 < 956*54
    assert 6*53 < 956 and 956**2 > 28561
    assert sp.expand(quartic(u, v)
                     -(u*u*(u*(u-956*v)-6*v*v)+956*u*v**3+v**4)) == 0
    print("PASS: sharper positive-quartic height and integer-gap identity")

    for degree in range(11):
        polynomial_b = sum(
            sp.binomial(degree+sp.Rational(1, 4), degree-j)
            * sp.binomial(degree-sp.Rational(1, 4), j) * x**j
            for j in range(degree+1)
        )
        polynomial_a = sp.expand(x**degree*polynomial_b.subs(x, 1/x))
        remainder = sp.expand(w*polynomial_b.subs(x, w**4)
                              -polynomial_a.subs(x, w**4))
        quotient, residual = sp.div(remainder, (w-1)**(2*degree+1), w)
        assert residual == 0
        assert all(c > 0 for c in sp.Poly(quotient, w).all_coeffs())
        k_product = prod(16*j*j-1 for j in range(1, degree+1))
        assert k_product * 4**degree * factorial(2*degree) == factorial(4*degree+1)
        assert quotient.subs(w, 1) == sp.Rational(k_product,
                                                  factorial(2*degree+1))
        assert quotient.subs(w, 1) == sp.Rational(comb(4*degree+1, 2*degree), 4**degree)
        assert quotient.subs(w, 1) <= 4**degree
    print("PASS: finite polynomial identities and positive coefficients, N=0..10")

    h_values = [approximant(degree) for degree in range(13)]
    assert h_values[0] == (1, 0) and h_values[1] == (1912, 2)
    for degree in range(12):
        s0, t0 = h_values[degree]
        s1, t1 = h_values[degree+1]
        determinant = s0*t1-s1*t0
        expected = Fraction(
            (-1)**degree*2**(2*degree+1)
            * prod(16*j*j-1 for j in range(1, degree+1)),
            factorial(degree)*factorial(degree+1),
        )
        assert determinant == expected and determinant != 0
        assert abs(s0) <= 4000**degree
    print("PASS: integer approximants and explicit nonzero determinants, N=0..11")

    assert 3840*59 < 478**2
    assert 4000**10 < 59**21
    assert 200**21 > 64000**10
    assert 127000**10 < 500000**9
    assert 508**10 * 4000**20 < 200**21 * 500000**9
    print("PASS: approximation and height-bound integer comparisons")

    a, b, c, d = 191, 182597, 765, 731344
    assert quartic(b, a) == 6396269108732
    assert quartic(d, c) == -401982685574639
    assert b*c-a*d == 1
    assert b+d == 913941 > 500000
    assert 254*d == 185761376 < 600**3
    print("PASS: rational bracket, determinant-one exclusion, U < 600")

    assert Fraction(1)-Fraction(599, 956) > Fraction(1, 3)
    assert 150*10**3 > 114244
    assert 450*7**3*Fraction(947, 956) > 114244
    expected = {
        (3, 1): (-22916, 22972),
        (4, 1): (-57199, 57521),
        (5, 1): (-114244, 115196),
        (5, 2): (-200719, 200801),
        (6, 1): (-199679, 201841),
        (6, 2): (-366656, 367552),
    }
    sector = {(u0, v0) for u0 in range(1, 7) for v0 in range(1, u0+1)
              if v0*v0+2*u0*v0-u0*u0 <= 0}
    assert sector == set(expected)
    for (u0, v0), values in expected.items():
        assert (quartic(u0, v0), quartic(u0, -v0)) == values
    assert quartic(13, 0) == 28561
    assert quartic(5, 1) == -114244
    # A weaker sector |v| <= u would invalidate the terminal classification.
    assert quartic(3, -2) == 28561
    assert (3+abs(-2))**2 > 2*3**2
    assert 239**2+1 == 2*13**4
    print("PASS: exhaustive six-row terminal table and reconstructed solutions")
    print("All exact checks passed. Universal arguments are separately proved in Lean.")


if __name__ == "__main__":
    main()
