# An effective endpoint bound for squares at odd indices

Write \(S_k(x)=\sum_{j=1}^x j^k\). For every odd \(k\ge7\) and every positive integral solution of \(S_k(x)=y^2\),
\[
                 \log\log(3x)\ \ll\ k^2\log(2k).
\]
The implied constant is absolute and effective. The same bound applies to every even target exponent. At \(k=19\) the argument below gives **\(\log(2x+1)<10^{350}\)** on every global square class. The known finite exclusions leave a large interval before this upper bound. [Numerical specialization](nineteen.md) · [Finite checks](verify.py)

For an odd prime target \(p\), Notebook II’s binomial inversion supplies a second construction in a degree-\(p\) field. Together the two arguments give \(\log\log(3x)\ll k^2\log^2(2k)\) for **all** targets at every \(k\ge6\). [Odd prime and all-target proofs](odd_prime.md)

## The centered source polynomial

Put \(N=2x+1,\ 2m=k+1\). Reflection and translation give the monic polynomial
\[
 F(s)=2^{2m}\bigl(B_{2m}((1+\sqrt s)/2)-B_{2m}\bigr),\qquad
 F(N^2)=2^{2m}(2m)y^2 .
\]
The coefficients before the constant subtraction are \(\binom{2m}{2j}2^{2j}B_{2j}(1/2)\) at \(s^{m-j}\). For \(j>0\), \(2^{2j}B_{2j}(1/2)=(2-2^{2j})B_{2j}\). These are the editorial Bernoulli identities in Andrews and Berndt's Lost IV §8.4. Our [simple-root theorem](../17_all_indices/roots.md) uses these coefficients, together with a real-root count and a multiple-root bound, to prove that \(F\) has a simple nonreal root \(z\) at every index in question. This supplies exactly the conjugate pair needed below. [Source image, physical p. 210](../../../notebooks/lost-4/pages/0210.webp)

Let \(a=(k+2)!\), \(P=aF\), \(b=a2^{2m}(2m)\), and \(H=k\log(2k)\). By von Staudt–Clausen every Bernoulli denominator through \(2m\) divides \(a\). Euler’s even-zeta evaluation, Entry 25(i), gives the factorial bound \(|B_j|\le 4 j!\) and hence \(\log\|P\|_\infty+\log|ab|\ll H\). [Notebook I source and coefficient bounds](../notebook-methods/deep/regular-1/chapters/ch05-entry25.md)

Write \(f\) for the primitive irreducible polynomial of \(z\), with positive leading coefficient. Gauss's lemma gives \(\operatorname{lc}(f)\mid a\). The Mahler measure bound
\(\|f\|_\infty\le 2^m M(f)\le 2^m M(P)\le 2^m\sqrt{m+1}\|P\|_\infty\)
shows that every conjugate of \(z\) has modulus at most \(\exp(O(H))\). In particular \(az\) is integral and \(z\ne0\).

## Uniform ramification in a single root field

Set \(K=\mathbb Q(z)\), \(h=[K:\mathbb Q]\le m\), \(L=K(\sqrt{N^2-z})\), and \(d=[L:\mathbb Q]\le 2m\). Lemma 3.4 of Bérczes et al. for a single root gives \(\log|D_K|\ll m^2+m\log\|f\|_\infty\ll kH\). We bound the relative discriminant uniformly in \(N\). [Bérczes et al., Lemmas 3.1 and 3.4](https://arxiv.org/pdf/2310.09704)

For a prime ideal \(\mathfrak p\) of \(K\) over \(p\nmid2ab\), divide \(P(X)=(X-z)G(X)\) over the local integers. If \(v_{\mathfrak p}(N^2-z)>0\) and \(P'(z)\) is a local unit, then \(G(N^2)\equiv G(z)=P'(z)\pmod{\mathfrak p}\). Thus \(v_{\mathfrak p}(N^2-z)=v_{\mathfrak p}(b y^2)\) is even. At odd residue characteristic, even valuation yields a split or unramified quadratic extension. At the remaining primes away from \(2ab\) the relative discriminant exponent is at most one. Their product of norms divides \(|\operatorname{Res}(f,P')|\), since \(\operatorname{lc}(f)\) is a \(p\)-unit and \(v_p(\operatorname{Res}(f,P'))=\sum_{\mathfrak p\mid p}f_{\mathfrak p}v_{\mathfrak p}(P'(z))\). This counts only the affected prime ideals. Simplicity makes this resultant a nonzero integer. Hadamard's bound on the Sylvester matrix gives \(\log|\operatorname{Res}(f,P')|\ll kH\).

Locally multiply the radicand by a square to have valuation zero or one. The polynomial discriminant of this integral quadratic presentation bounds the relative discriminant exponent by \(2v_{\mathfrak p}(2)+1\). Thus the primes over \(2ab\) cost at most \(4^h\operatorname{rad}(2ab)^h\) in the discriminant norm. The tower formula and Lemma 3.1 yield
\[
  \log |D_L|\ll kH,\qquad R_L\le |D_L|^{1/2}(\log^*|D_L|)^{d-1},\qquad
  \log^* R_L\ll kH.
\]
Here \(\log^* v=\max(1,\log v)\). This also holds when the quadratic extension is trivial.

## Two close conjugates of one algebraic integer

Fix the embedding \(z\in\mathbb C\) above and choose \(s=\sqrt{N^2-z}\) with positive real part. Set \(A=a(N+s)\). Both \(A\) and \(A'=a(N-s)\) are algebraic integers and \(AA'=a^2 z\). If \(L/K\) is quadratic, this is the relative norm. If \(L=K\), the integral absolute norm of \(A\) divides that of \(a^2z\). Consequently \(0\le\log|\mathrm N_{L/\mathbb Q}A|\ll kH\) in either case. At all embeddings \(|\sigma A|\le a(N+\sqrt{N^2+\max_{\tau}|\tau z|})\), so its absolute logarithmic height satisfies \(h(A)\ll\log N+H\).

In the chosen embedding, \(|N+s|\ge N\) and \(|s-N|\le |z|/N\). The imaginary part of \(s\) is nonzero. Ordinary complex conjugation therefore gives
\[
             0<|A/\bar A-1|\le 2|z|/N^2.                     \tag{1}
\]
This works even when \(s\in K\). Both numerator and denominator are conjugates of the same balanced element and the same list of units in the next step.

## Balancing by units

Write \(r=r_1+r_2-1\) for the unit rank of \(L\). For each archimedean place take the coordinate \(w_v(u)=d_v\log|\sigma_v u|\), with \(d_v=1,2\) as usual. Delete one coordinate of a unit logarithm. The resulting lattice in \(\mathbb R^r\) has covolume \(R_L\). If \(r=0\), \(L\) has a nonreal embedding and is imaginary quadratic. There \(|A|^2=|\mathrm N A|\) and \(|A|\ge aN\), already giving the conclusion.

For \(r>0\), take independent log-unit vectors at the successive minima of the unit cube and denote their sup norms by \(L_i\), their units by \(u_i\), and the column matrix by \(U\). Minkowski's second theorem gives \(\prod L_i\le R_L\). Also \(|\det U|\ge R_L\) by the sublattice index. For a unit, \(dh(u_i)\le r L_i\) because the sum of the positive full coordinates equals the larger of the positive and negative sums among the retained coordinates. Voutier's lower bound consequently gives \(L_i\ge 2/(r(\log(3d))^3)\). Hence with \(E=\max(1,r(\log(3d))^3/2)^r\),
\[
                      \prod_i\max(1,L_i)\le R_L E,\qquad
                      \log E\ll k\log(2k).
\]
Only independence is needed here. [Minkowski's minima and unit heights](https://arxiv.org/abs/1508.01969) · [Voutier bound, Lemma 3.8](https://arxiv.org/pdf/2310.09704)

Subtract \(\log|\mathrm N A|/(r+1)\) from every coordinate \(w_v(A)\) and delete the same last coordinate, obtaining \(w\). Integrality and the norm give \(\|w\|_\infty\ll d(\log N+H)+kH\), since the total of the negative logarithms in absolute value is at most the positive total \(dh(A)\). Solve \(Uc=w\) over the reals, round each \(c_i\) to an integer \(b_i\), and set \(\gamma=A\prod u_i^{-b_i}\). From Cramer's rule,
\[
 |b_i|\le 1+r!\|w\|_\infty E .
\]
Every centered coordinate of \(\gamma\) is bounded by \(r\max_i L_i/2\) in the retained rows and \(r^2\max_i L_i/2\) in the omitted row. It is integral, and its absolute norm equals that of \(A\). These estimates prove
\[
 \log\max(1,h(\gamma))\ll kH,\qquad
 \prod_i\max(1,h(u_i))\le R_L E,\qquad
 \log B\ll k\log(2k)+\log\log(3N),
\]
where \(B=\max(3,|b_1|,\ldots,|b_r|)\).

## A complex linear form

Work in the compositum \(L\bar L\), of degree \(D\le d^2\). Apply the archimedean form of Proposition 3.12 of Bérczes et al. to \(\gamma/\bar\gamma\) and \(u_i/\bar u_i\), with exponents \(1,b_i\). Their product is \(A/\bar A\). Keep all \(r+1\ge2\) displayed factors, including any with zero exponent or height zero. Proposition 3.12 permits such factors, and (1) gives its required nonvanishing. The product of the height parameters \(\max(h(\alpha),V(D))\) is at most \(\exp(O(kH))\) by \(h(v/\bar v)\le2h(v)\) and the preceding product bound. The proposition's factor \(12(16eD)^{3(r+1)+2}(\log^*D)^2\) contributes \(\exp(O(k\log(2k)))\). Passing from its normalized archimedean valuation to the usual complex modulus preserves this lower bound up to an absolute factor. [Proposition 3.12](https://arxiv.org/pdf/2310.09704)

Together with (1), this gives \(t=\log(3N)\le\exp(O(kH))(\log t+kH)\). For \(t\ge 1\), \(\log t\le\sqrt t\). Dividing by \(\sqrt t\) and absorbing \(kH\) proves the asserted effective bound on \(\log t\). The regulator used throughout belongs to the field of degree at most \(k+1\). The larger compositum enters only in the last linear form. The [even-index estimate](../47_even_height/README.md) supplies the complementary square bound for even indices.
