# A quintic descent for squares in the thirteenth power sum

Let \(T=x(x+1)/2\) for a positive integer \(x\), and let
\[
R(T)=960T^5-2800T^4+4592T^3-4720T^2+2764T-691.
\]
Define \(g(W)=60^5 R(W/60)/960\), let \(w\) be a root, and put \(K=\mathbb Q(w)\). The polynomial \(g\) is monic and irreducible.

**Theorem.** \(S_{13}(x)\) is a square if and only if there is a \(z\in\mathcal O_K\) satisfying
\[
             (60-w)z^2=60T-w.
\]
Thus every positive integral solution lifts to the same quadratic cover as \(T=1,z=\pm1\). The [relative-unit proof](../16_relative_units/README.md) determines all its triangular integral parameters with positive endpoint. The class-group computation below is unconditional.

## The norm and its ideal

The identity \(105S_{13}(x)=T^2R(T)\) follows by taking the finite difference in \(x\) and checking the value at 0. Set \(\gamma=(60T-w)/(60-w)\). Its norm is
\[
N_{K/\mathbb Q}(\gamma)=g(60T)/g(60)=R(T)/105=S_{13}(x)/T^2.
\]
A square in \(K\) therefore implies an integral square on the original side, since a rational square which is an integer is an integral square. For the forward direction suppose this norm is a rational square.

The order index \([\mathcal O_K:\mathbb Z[w]]\) is \(2^7 3^6 5^6\), and \(g(60)=2^4 3^5 5^5 7\). For a prime \(p\notin\{2,3,5,7\}\) a prime ideal dividing \(60T-w\) corresponds, by Dedekind factorization, to the linear factor \(W-60T\) of \(g\bmod p\). It is unique over that residue and has residue degree one, even when the factor is repeated. Its valuation equals the rational \(p\)-valuation of the norm of the numerator. The denominator has no prime factor above \(p\), so the valuation is nonnegative and even.

For the four other primes enumerate the triangular residues \(a\bmod p^4\) for which \(105R(a)\) is square modulo \(p^4\), including zero as a square. Computing \(h(h+1)/2\) for \(h\) in a complete interval of length \(2p^4\) enumerates every triangular residue. For each surviving \(a\) and each \(\mathfrak p\mid p\), the verifier checks
\[
v_{\mathfrak p}(60a-w)=v_{\mathfrak p}(60-w)
  < e_{\mathfrak p}\bigl(v_p(60)+4\bigr).
\]
The strict inequality makes this valuation constant for every lift \(T\equiv a\bmod p^4\).

| \(p\) | Admissible triangular residues | Valuations at primes over \(p\), in verifier order |
|---:|---:|---|
| 2 | 8 | 3, 0, 1 |
| 3 | 4 | 2, 0, 1, 1 |
| 5 | 125 | 1, 1, 3 |
| 7 | 343 | 0, 0, 1 |

It follows that \(\gamma\) is an algebraic integer and \((\gamma)=\mathfrak a^2\) for an integral ideal. PARI computes class number one for \(K\). The verifier uses the isomorphic reduced quintic and requires the full \(\mathtt{bnfcertify}\) check to succeed. Consequently \(\gamma=u v^2\) for some unit \(u\) and some \(v\in\mathcal O_K\). [PARI/GP number-field reference](https://pari.math.u-bordeaux.fr/dochtml/ref/General_number_fields.html)

## The unit is a square

The signature is \((1,2)\). Dirichlet's unit theorem gives eight unit square classes, since the only roots of unity are \(\pm1\). Use two units and \(-1\) from PARI in the original \(w\) presentation, keeping the units in factored form. For each prime \(\mathfrak p\mid3\), evaluate the quadratic Hilbert symbol \((\,\cdot\,,3)_{\mathfrak p}\). Across all four admissible residues \(a\bmod81\), compare the four symbols for \((60a-w)(60-w)\) against the four symbols for each of the eight unit products. Only the identity product matches. The product of numerator and denominator represents the same square class as their quotient. The strict valuation inequality above implies that replacing \(a\) by an arbitrary lift changes the numerator by a factor in \(1+\mathfrak p\), a local square at these odd primes.

In particular the seven nonempty products all differ from the identity at some prime over 3, since \(T=1\) is included in the test. The three computed units give independent classes, hence all eight classes by the signature. This last check verifies the unit basis modulo squares independently of its computation in the original presentation. We conclude \(u\in\mathcal O_K^{*2}\), proving the theorem.

## The remaining curve

Write \(\theta=w/60\) and \(f(t)=105R(t)\). The theorem also puts each square solution on the single quartic over \(K\)
\[
V^2=(1-\theta)\frac{f(t)}{t-\theta},\qquad t=T,\qquad V=\frac{105 y}{Tz}
\]
where \(S_{13}(x)=y^2\). It has the point \((1,105)\). Indeed \(f(T)=(105y/T)^2\) and \(z^2=(T-\theta)/(1-\theta)\). The four remaining roots of the separable quintic are distinct, so the quartic is a genus-one curve. Stage 16 determines its square-sum endpoints via the integral norm equation.

A separate check on the genus-two curve \(C:Y^2=f(t)\) gives \(\#J(\mathbb F_{11})=214\) and \(\#J(\mathbb F_{13})=157\) for \(J=\operatorname{Jac}(C)\). These are good primes. Prime-to-\(p\) torsion injects under reduction, and the two counts are coprime, with \(11\nmid157\) and \(13\nmid214\). Thus \(J(\mathbb Q)\) has trivial torsion. The divisor \([(1,105)-\infty]\) is nonzero, since a degree-one function would imply genus zero. The rank is at least one. The upper bound below is independent of the integral square-class calculation.

Run the [PARI/GP verifier](verify.gp) with \(\mathtt{gp\ -q\ -s\ 256M\ verify.gp}\). The final class-group certification can take several minutes. The [odd-prime argument](../13_cyclotomic_cofactor/README.md) deals separately with every target exponent divisible by an odd prime.

## An upper bound on the Jacobian rank

The [second PARI/GP verifier](rank.gp) gives \(1\leq\operatorname{rank}J(\mathbb Q)\leq3\). For the upper bound change coordinates to \(u=60t\) to obtain \(v^2=105g(u)\), then to the monic curve with polynomial \(F(X)=105^5g(X/105)\). Its root field is \(K\). The descent map from \(J(\mathbb Q)/2J(\mathbb Q)\) into \(K^*/K^{*2}\) has square norm and sends an affine point minus infinity to \(105(u-w)\). The only primes whose squares divide the discriminant of \(F\) are \(2,3,5,7\). Away from these primes all local images have even valuations by Stoll's Proposition 4.6. [Stoll, *Implementing 2-descent for Jacobians of hyperelliptic curves*, §4](https://www.impan.pl/shop/publication/transaction/download/product/83397)

There are thirteen prime ideals of \(K\) above \(2,3,5,7\). Together with the certified trivial class group and the three unit square classes forming a basis, this gives an ambient square-class group of dimension 16. The verifier factors the constituents of each candidate generator to check even valuation outside these thirteen primes. It also checks independence of the sixteen generators by their joint images under the real sign and local characters. This establishes that the tested group contains every possible global image.

For \(p=2,3,5\), form the vector of valuations modulo two and quadratic Hilbert symbols in each factor of \(K\otimes\mathbb Q_p\). The following \(u\)-coordinates support local points, as checked by the ordinary square criterion in \(\mathbb Q_p\).

| \(p\) | \(u\) giving independent local classes | Dimension |
|---:|---|---:|
| 2 | \(-1740,-1620,-57,-46\) | 4 |
| 3 | \(-1680,-1560,-56\) | 3 |
| 5 | \(-1740,-1620\) | 2 |

These dimensions equal \(\dim J(\mathbb Q_p)/2J(\mathbb Q_p)\) by Stoll's Lemma 4.4. The measured vectors therefore span the image under the selected characters of every local divisor class. A global class must have positive real norm and lie in all three measured local spans. Of the \(65{,}536\) products of the sixteen generators, the real sign leaves \(32{,}768\). The successive local tests leave \(256,16,8\). Consequently the image of \(J(\mathbb Q)/2J(\mathbb Q)\) has dimension at most three. The known infinite-order divisor gives the lower bound. Determining whether the rank is one would require a further argument.

The [paired descent](paired.md) combines this curve with a second genus-two quotient and extends the single square-class result to rational triangular points.
