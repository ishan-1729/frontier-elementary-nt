# A shorter proof of the remaining positive quartic case

**Verification status: PASS.** The complete classification and all adopted steps
have been checked in Lean. The full build and all six audits exited 0: 70
mathematical modules and 1,720 local mathematical declarations. The independent
new-proof audit covers all302 added declarations and excludes the old error
construction. See [the verification receipt](../../../../../lean/FINAL_AUDIT.md).

The key discovery is that the same determinant controls two parts of the
proof: the sign of a polynomial norm, and the separation of distinct integer
fractions. This replaces the earlier error-factorization and coefficient
majorization argument. It retains the rational approximants and the elementary
construction that establishes their integrality.

## 1. The remaining arithmetic target

The Vierzahlensatz reduction has already reduced the nonterminal positive
branch to integers satisfying

\[
V>0,\qquad U>956V,\qquad 0<\Phi(U,V)\le28561,
\]
where
\[
\Phi(U,V)=U^4-956U^3V-6U^2V^2+956UV^3+V^4.
\]

The separate arithmetic branch argument classifies the right-hand side1
and excludes169. The argument below excludes the entire displayed interval,
and therefore the remaining right-hand side28561.

Let \(\xi\) be the small positive root of
\[
f(x)=1-956x-6x^2+956x^3+x^4.
\]
The established elementary root bounds give
\[
0<\xi<1/956.
\]
On \([0,1/8]\), the difference quotient of \(f\) is less than−900.
Consequently the hypothetical quartic pair satisfies
\[
0<\xi-V/U<\frac{28561}{900U^4}<\frac{32}{U^4}.
\tag{1}
\]
This comparison is a finite polynomial identity and inequalities; it does not
use the approximation construction.

## 2. An integer recurrence

Temporarily let \(a\) be a real parameter. Define polynomial pairs by
\[
S_0=1,\quad S_1=8a,\qquad T_0=0,\quad T_1=2,
\]
and, for either \(Y=S\) or \(Y=T\), with \(n\ge1\),
\[
n(n+1)Y_{n+1}=8an(2n+1)Y_n+4(16n^2-1)Y_{n-1}.
\tag{2}
\]

At \(a=239\), these are the already constructed integer sequences. Their
integrality is a separate proved input: division by \(n(n+1)\) in(2) alone
would not establish it. The finite construction recursively solves
\[
h(1-h^2)(1-4ah-6h^2+4ah^3+h^4)=z(1+h^2)^4
\]
for integer coefficients, with \(h_0=0\); each new coefficient occurs with
coefficient1. The accompanying integral quotients give precisely \(S_n,T_n\).
The construction and recurrence identification are proved in
[IntegralCoefficientConstruction.lean](../../../../../lean/Ljunggren/IntegralCoefficientConstruction.lean)
and [NormSignSequenceBridge.lean](../../../../../lean/Ljunggren/NormSignSequenceBridge.lean).
Only this integrality construction is retained from the coefficient method.

## 3. The replacement for the long error calculation

For fixed \(n\ge1\), write \(s=S_n,t=T_n\), and set
\[
\begin{aligned}
W&=s^2+t^2,\\
A&=s^4-6s^2t^2+t^4,\qquad B=4st(s^2-t^2),\\
F&=A-aB,\qquad\Psi=aA+B,\\
c&=(16n^2-1)/(2n),\qquad
\Delta=S_{n-1}T_n-S_nT_{n-1}.
\end{aligned}
\]

The recurrence makes \(\Delta\) a nonzero constant in \(a\), of sign
\((-1)^{n-1}\). A direct induction also gives
\[
(1+a^2)s'=nas+t/4+cS_{n-1},\qquad
(1+a^2)t'=nat-s/4+cT_{n-1}.
\]
Primes here denote formal polynomial derivatives.

The decisive finite identity is
\[
\boxed{F'\Psi-F\Psi'=4c\Delta W^3.}
\tag{3}
\]
It follows directly from the two displayed derivative equations. It involves
neither a binomial remainder nor an infinite estimate.

For \(a\ge1\), the positive recurrence gives \(s>0,t>0,4at\le s\),
so \(\Psi>0\). Moreover,
\[
\deg W=2n,\qquad \deg\Psi=4n+1,\qquad \deg F\le4n,
\]
and \(W,\Psi\) have positive leading coefficients.
The nonzero right-hand side of(3) implies \(F\ne0\).
If \(m=\deg F\), comparison of leading terms in(3) gives
\[
m=2n,\qquad
\operatorname{sign}(\operatorname{lc}F)=-\operatorname{sign}\Delta.
\]
Indeed, the leading coefficient on the left is
\((m-4n-1)\operatorname{lc}(F)\operatorname{lc}(\Psi)\), and its prefactor
is strictly negative because \(m\le4n\).

At any real zero \(r\ge1\) of \(F\), identity(3) instead gives
\(\operatorname{sign}F'(r)=\operatorname{sign}\Delta\).
Thus every possible zero would cross in the direction opposite to the
eventual leading sign. Taking the largest zero shows that none exists.
Therefore
\[
\boxed{(-1)^nF(a)>0\qquad(a\ge1).}
\tag{4}
\]
The case \(n=0\) is immediate from \(F=1\).

The full finite derivation is in
[the norm-sign lemma](elementary_norm_sign_lemma.md).
The Lean proof also checks the equivalent identity
\[
(1+a^2)WF'-\{2(1+a^2)W'+aW\}F=4c\Delta\Psi,
\]
and uses its leading coefficients and the same largest-zero argument.
The root-sign argument uses ordinary continuity of real polynomials.

Specializing(4) at239 puts consecutive fractions \(T_n/S_n\) on opposite
sides of \(\xi\), because \(f\) is strictly decreasing on the interval
containing those fractions. In particular, the errors
\(E_n=\xi S_n-T_n\) satisfy
\[
E_nE_{n+1}<0.
\tag{5}
\]

## 4. The error bound now follows from the determinant

From this point onward \(a=239\). Write(2) as
\(Y_{n+1}=\alpha_nY_n+\beta_nY_{n-1}\), with \(n\ge1\).
Its coefficients satisfy
\[
2868\le\alpha_n<3824,\qquad0<\beta_n<64.
\]
Thus direct induction gives
\[
0<S_n\le4000^n,\qquad S_{n+1}\ge1912\,2868^n.
\]
For \(D_n=S_nT_{n+1}-S_{n+1}T_n\), subtraction gives
\[
D_0=2,\qquad D_n=-\beta_nD_{n-1}.
\]
Hence \(D_n\ne0\) and \(|D_n|\le2\,64^n\).

Opposite error signs in(5) imply
\[
|E_n|S_{n+1}\le|D_n|.
\]
For \(n\ge1\), since \(40\cdot64<2868\),
\[
|E_n|\le\frac1{956}\left(\frac{64}{2868}\right)^n
<\frac1{956\,40^n}.
\tag{6}
\]
For \(n=0\), the same strict bound is \(\xi<1/956\).

## 5. An integer cannot lie strictly between0 and1

Choose the least \(N\) such that \(U\le400\,40^N\).
Since \(U>956\), we have \(N\ge1\), and
\(U>400\,40^{N-1}\).

Because \(D_N\ne0\), one of
\(UT_N-VS_N\) and \(UT_{N+1}-VS_{N+1}\) is a nonzero integer.
For its index \(j\), equations(1) and(6) give
\[
1\le|UT_j-VS_j|
<\frac{400}{956}+\frac{32\,4000^{N+1}}{U^3}
<\frac12+\frac{32\,4000^{N+1}}{U^3}.
\]
Therefore
\[
U^3<64\,4000^{N+1}.
\tag{7}
\]
If \(N\ge2\), however,
\[
(400\,40^{N-1})^3\ge64\,4000^{N+1}.
\]
Equality holds at \(N=2\), and each subsequent step multiplies the ratio
by16. This contradicts(7). Thus \(N=1\) and \(U\le16000\).

Now use the explicit first integer difference
\[
UT_1-VS_1=2(U-956V)>0.
\]
The same argument, with \(S_1=1912\), gives
\[
U^3<64\cdot1912<956^3,
\]
contradicting \(U>956\). The remaining28561 case is excluded.

## 6. Relation to the proposed arithmetic descent

The first actual integer differences are
\[
\begin{aligned}
R_0&=-V,\\
R_1&=2k,\qquad k=U-956V,\\
R_2&=-6w,\qquad w=5V-956k,\\
R_3&=4(21k-4780w).
\end{aligned}
\]
Thus the first Euclidean remainders previously extracted from the
Vierzahlensatz equations are exactly the first terms used here.
While three consecutive \(R\)'s alternate in sign, their positive magnitudes
satisfy \(59|R_n|<|R_{n-1}|\). This ensures a finite sign exit for every
integer input; it alone does not exclude a quartic solution.
The new polynomial sign argument supplies the uniform root control that
was missing from that attempted descent.

## 7. Scope

The argument is elementary in its stated ingredients: finite polynomial
identities, a proved integer-coefficient construction, polynomial root signs,
positive recurrences, and separation of integers. It is still an approximation
proof. No historical novelty claim is made.

For the complete equation, both signed quartic sectors are handled using the
same root bracket and base40 estimate. The full-sector comparison forces
the selected index to be at most2, then gives \(U<32000\) and \(U<79\).
The previously checked finite sector arithmetic completes the classification
of positive solutions of the original equation:
\((x,y)=(1,1),(239,13)\). The screenshot branch gives
\((r,s,t,u)=(8,4,5,8)\).

The assembled statements are in
[CompleteNormSign.lean](../../../../../lean/Ljunggren/CompleteNormSign.lean).
The proof-dependency audit is
[NormSignProofAudit.lean](../../../../../lean/NormSignProofAudit.lean).
