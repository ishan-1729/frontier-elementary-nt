# An elementary norm-sign proof of Ljunggren’s equation

The only positive integer solutions of \(x^2+1=2y^4\) are \((1,1)\) and \((239,13)\). The proof below retains both signs produced by the arithmetic reduction and recovers the original variables at the end.

The new ingredient is a finite polynomial sign argument: the determinant of two neighboring recurrence pairs determines the direction in which their quartic norm could cross zero. Its leading coefficient has the opposite sign, so no such crossing occurs. This gives exact root brackets; their determinant then gives the entire error estimate. The earlier error-factorization and coefficient-majorization argument is no longer needed. Integrality still comes from the proved elementary integer-coefficient construction. The method remains rational approximation and integer separation, rather than a solution-preserving descent. No historical novelty is claimed.

**Verification status:** the new construction, actual root brackets, base-40 termination, and complete original-equation and screenshot classification statements have compiled in Lean. The combined build and all six proof-dependency audits passed, covering70 mathematical modules and1,720 local mathematical declarations. The main statements are in [CompleteNormSign.lean](../../../../../lean/Ljunggren/CompleteNormSign.lean).

The self-contained replacement lemma is [the finite norm-sign proof](elementary_norm_sign_lemma.md). Its use in the remaining positive quartic case, including the sharper first-index contradiction, is explained in [the positive-case proof](elementary_norm_sign_proof.md). This document supplies the full original reduction, both signed sectors, and recovery of the variables.

There is also an independent [arithmetic proof of the entire \(13\mid y\) branch](elementary_arithmetic_branch.md), using a seventh-power identity, coprime factor splitting, small congruences, and a consecutive-square gap. It handles the centered values 1 and 169. That branch result is a partial arithmetic replacement; the uniform norm-sign argument below closes the remaining case and the full equation.

## 1. Two coprime factorizations

The equation implies that \(x\) and \(y\) are odd. The case \(y=1\) gives \(x=1\), so assume \(y>1\).

Put
\[
A=\frac{x+1}{2},\qquad B=\frac{x-1}{2}.
\]
Then
\[
A-B=1,\qquad \gcd(A,B)=1,\qquad A^{2}+B^{2}=y^{4}.
\]

We use the primitive Pythagorean parametrization twice. Recall that this follows from the same coprime-factor principle as the Vierzahlensatz: for positive integers \(r,s,t\), if
\[
r^{2}+s^{2}=t^{2},\qquad \gcd(r,s)=1,
\]
with \(r\) odd and \(s\) even, then
\[
\frac{t+r}{2}\frac{t-r}{2}=\left(\frac{s}{2}\right)^{2}.
\]
The two factors on the left are coprime, so each is a square. Thus
\[
r=c^{2}-d^{2},\qquad s=2cd,\qquad t=c^{2}+d^{2}.
\]

Applied to \(A^{2}+B^{2}=y^{4}\), this gives
\[
\{A,B\}=\{c^{2}-d^{2},\,2cd\},
\qquad c^{2}+d^{2}=y^{2},
\]
where \(c,d\) are coprime and of opposite parity. Since \(A-B=1\),
\[
c^{2}-d^{2}-2cd=\pm1.
\]

Parametrizing \(c^{2}+d^{2}=y^{2}\) gives
\[
\{c,d\}=\{m_{0}^{2}-n_{0}^{2},\,2m_{0}n_{0}\},
\qquad y=m_{0}^{2}+n_{0}^{2},
\]
where \(m_{0},n_{0}\) are again coprime and of opposite parity.

Define
\[
Q(m,n)=m^{4}+4m^{3}n-6m^{2}n^{2}-4mn^{3}+n^{4}.
\]
Both assignments of the second pair of Pythagorean legs must be retained. Explicitly,
\[
\begin{aligned}
c=m_0^2-n_0^2,\ d=2m_0n_0
&\Longrightarrow c^2-d^2-2cd=Q(m_0,-n_0),\\
c=2m_0n_0,\ d=m_0^2-n_0^2
&\Longrightarrow c^2-d^2-2cd=-Q(m_0,n_0).
\end{aligned}
\]
Thus, for a suitable choice \(n=\pm n_{0}\) and \(m=m_{0}\),
\[
Q(m,n)=\pm1,\qquad y=m^{2}+n^{2}.
\]

The square completion
\[
Q(m,n)=(m^2+2mn-n^2)^2-8m^2n^2
\]
shows that \(Q(m,n)\) is a square modulo \(4\). Since a square is congruent
to \(0\) or \(1\) modulo \(4\), the value \(-1\) is impossible. This
argument works for all integers \(m,n\), without a parity assumption, and
is verified by `quartic_ne_neg_one` in Lean. We have reduced the original equation to
\[
\boxed{Q(m,n)=1,\qquad y=m^{2}+n^{2}.}
\tag{1}
\]

## 2. A change of variables centered at the known solution

Set
\[
U_{0}=3m+2n,\qquad V_{0}=-2m+3n,
\]
and define
\[
\Phi(U,V)
=
U^{4}-956U^{3}V-6U^{2}V^{2}+956UV^{3}+V^{4}.
\]

Direct expansion gives
\[
\Phi(3m+2n,-2m+3n)=13^{4}Q(m,n),
\]
and
\[
(3m+2n)^{2}+(-2m+3n)^{2}=13(m^{2}+n^{2}).
\]
Therefore
\[
\Phi(U_{0},V_{0})=28561,\qquad U_{0}^{2}+V_{0}^{2}=13y.
\tag{2}
\]

We also have
\[
\Phi(-V,U)=\Phi(U,V),
\tag{3}
\]
and
\[
\Phi(U+V,V-U)=\Phi(U-V,U+V)=-4\Phi(U,V).
\tag{4}
\]

Using rotations \((U,V)\mapsto (-V,U)\), arrange
\[
U>0,\qquad |V|\le U.
\]
If
\[
|V|>(\sqrt2-1)U,
\]
use \((U,V)\mapsto (U+V,V-U)\) when \(V>0\), and
\((U,V)\mapsto (U-V,U+V)\) when \(V<0\).
The new first coordinate is positive, and the new ratio of absolute coordinates is
\[
\frac{U-|V|}{U+|V|}<\sqrt2-1.
\]
The squared norm doubles under this transformation.

Thus every original solution produces integers \(U,V\) satisfying
\[
\boxed{
U>0,\qquad |V|\le(\sqrt2-1)U,
}
\tag{5}
\]
together with one of the two systems
\[
\boxed{
\begin{array}{c|c}
\Phi(U,V)&U^{2}+V^{2}\\ \hline
28561&13y\\
-114244&26y.
\end{array}}
\tag{6}
\]

We will prove that the only possibilities in this sector are
\[
\boxed{(U,V)=(13,0)\quad\text{or}\quad(5,1).}
\tag{7}
\]

## 3. Solutions must approximate a particular quartic irrational

Put
\[
f(t)=\Phi(1,t)
=1-956t-6t^{2}+956t^{3}+t^{4}.
\]

We first define its small positive root using radicals, so no existence theorem from analysis is needed. Let
\[
\tau=\sqrt{239^{2}+1},\qquad
\delta=\tau-239,\qquad
\xi=\frac{\delta}{1+\sqrt{1+\delta^{2}}}.
\]
Then
\[
\delta^{2}+478\delta=1,
\qquad
\delta=\frac{2\xi}{1-\xi^{2}}.
\]
Substituting the second identity into the first gives
\[
f(\xi)=0.
\]
Moreover,
\[
0<\delta=\frac1{\tau+239}<\frac1{478},
\]
so
\[
\boxed{0<\xi<\frac1{956}.}
\tag{8}
\]

For distinct real \(s,t\),
\[
\begin{aligned}
\frac{f(t)-f(s)}{t-s}
={}&-956-6(t+s)+956(t^{2}+ts+s^{2})\\
&+t^{3}+t^{2}s+ts^{2}+s^{3}.
\end{aligned}
\]
When \(|s|,|t|\le R:=5/12\), this is at most
\[
-956+12R+2868R^{2}+4R^{3}
=-\frac{195607}{432}<-450.
\]
Thus \(f\) is strictly decreasing on this interval, and
\[
|f(t)|>450|t-\xi|
\tag{9}
\]
whenever \(t\ne\xi\) lies in that interval.

Since \(\sqrt2-1<5/12\), equations (5)–(6) imply
\[
\boxed{
\left|\xi-\frac VU\right|
<
\frac{|\Phi(U,V)|}{450U^{4}}
<
\frac{254}{U^{4}}.
}
\tag{10}
\]

Our task is therefore to show that a rational number cannot approximate \(\xi\) this closely when \(U\) is large.

## 4. The integer recurrence and its norm sign

For a parameter \(a\), define polynomial pairs by
\[
S_0=1,\quad S_1=8a,\qquad T_0=0,\quad T_1=2,
\]
and, for \(n\ge1\) and \(Y=S,T\),
\[
n(n+1)Y_{n+1}
=8an(2n+1)Y_n+4(16n^2-1)Y_{n-1}.
\]
These are polynomials over the rationals. Their values at an integer parameter are integers, but that assertion requires a proof: the division in this recurrence alone does not supply it.

The retained construction recursively solves
\[
h(1-h^2)(1-4ah-6h^2+4ah^3+h^4)=z(1+h^2)^4,\qquad h(0)=0,
\]
with integer coefficients. At each degree the new coefficient of \(h\) occurs with coefficient 1; all other contributions involve earlier coefficients. The inverse of \(1+h^2\), whose constant coefficient is 1, also has integer coefficients. Set
\[
\mathcal S(z)=\left(\frac{h}{1+h^2}\right)',\qquad
\mathcal T(z)=\left(\frac{h^2}{1+h^2}\right)',
\]
where these primes mean formal differentiation in \(z\). Their integer coefficients have precisely the displayed initial values and recurrence. This identification is proved by finite coefficient algebra in [IntegralCoefficientConstruction.lean](../../../../../lean/Ljunggren/IntegralCoefficientConstruction.lean); [NormSignSequenceBridge.eval_integer](../../../../../lean/Ljunggren/NormSignSequenceBridge.lean) identifies them with the parameter-polynomial sequences for every index. The construction and recurrence identification are detailed in Sections 3–4 of [the earlier real recurrence manuscript](../../../../../Archive/ljunggren_working_2026_09/elementary_real_recurrence_proof.md). Its error construction is not used here.

Now fix \(n\ge1\), write \(s=S_n,t=T_n\), and put
\[
\begin{aligned}
W&=s^2+t^2,\\
A&=s^4-6s^2t^2+t^4,\qquad B=4st(s^2-t^2),\\
F&=A-aB,\qquad \Psi=aA+B,\\
c&=\frac{16n^2-1}{2n},\qquad
\Delta=S_{n-1}T_n-S_nT_{n-1}.
\end{aligned}
\]
The recurrence makes \(\Delta\) a nonzero constant polynomial of sign \((-1)^{n-1}\). A simultaneous finite induction gives
\[
(1+a^2)s'=nas+t/4+cS_{n-1},\qquad
(1+a^2)t'=nat-s/4+cT_{n-1},
\]
where primes now mean formal differentiation in the parameter \(a\). Expanding these identities gives
\[
\boxed{F'\Psi-F\Psi'=4c\Delta W^3.}
\]
The [norm-sign lemma](elementary_norm_sign_lemma.md) supplies the full induction and expansion; the actual identity is checked as `NormSignActual.actual_wronskian`.

For \(a\ge1\), the positive recurrence gives \(s>0,t>0,4at\le s\). Hence \(W>0\) and \(\Psi>0\). The polynomial degrees satisfy
\[
\deg W=2n,\qquad \deg\Psi=4n+1,\qquad \deg F\le4n,
\]
with positive leading coefficients for \(W\) and \(\Psi\). Since the right side of the boxed identity is nonzero, \(F\ne0\). If \(m=\deg F\), the leading term on its left has degree \(m+4n\) and coefficient
\[
(m-4n-1)\operatorname{lc}(F)\operatorname{lc}(\Psi)\ne0.
\]
The right has degree \(6n\). Thus \(\deg F=2n\), and its leading sign is opposite to the sign of \(\Delta\).

At every hypothetical real zero \(r\ge1\), the same identity gives
\[
F'(r)\Psi(r)=4c\Delta W(r)^3.
\]
The zero is simple, and its derivative has the sign of \(\Delta\). Taking the largest such zero gives a contradiction: immediately to its right \(F\) has that sign, while sufficiently far to the right it has the opposite leading sign. Continuity would force a larger zero. Therefore
\[
(-1)^nF(a)>0\qquad(a\ge1,n\ge1).
\]
At \(n=0\), \(F=1\), so the same conclusion holds separately.

Specialize now to \(a=239\), and suppress this parameter in \(S_n,T_n\). Every ratio \(T_n/S_n\) is in \([0,1/956]\). Because \(f\) is strictly decreasing there, the norm sign gives
\[
\boxed{(-1)^n(\xi S_n-T_n)>0\qquad(n\ge0).}
\]
This is the strict consecutive root bracket, proved without any prior error estimate.

At this parameter, write the recurrence as
\[
Y_{n+1}=\alpha_nY_n+\beta_nY_{n-1},\qquad
\alpha_n=\frac{1912(2n+1)}{n+1},\quad
\beta_n=\frac{4(16n^2-1)}{n(n+1)},\qquad n\ge1.
\]
Its coefficients satisfy \(2868\le\alpha_n<3824\) and \(0<\beta_n<64\). Thus
\[
0<S_n\le4000^n,\qquad S_{n+1}\ge1912\,2868^n.
\]
The upper bound follows, for example, from
\(3824\cdot4000+64<4000^2\) and the two initial values. For
\(D_n=S_nT_{n+1}-S_{n+1}T_n\), subtraction gives
\[
D_0=2,\qquad D_n=-\beta_nD_{n-1}\quad(n\ge1),\qquad
0<|D_n|\le2\,64^n.
\]
Since consecutive fractions strictly bracket \(\xi\),
\[
|\xi S_n-T_n|<\frac{|D_n|}{S_{n+1}}
\le\frac1{956}\left(\frac{64}{2868}\right)^n.
\]
For \(n\ge1\), use \(40\cdot64<2868\); at \(n=0\), use (8). We have established, for every \(n\ge0\),
\[
\boxed{
S_n,T_n\in\mathbb Z,\quad
0<S_n\le4000^n,\quad
|\xi S_n-T_n|<\frac1{956\,40^n},\quad
D_n\ne0.
}
\tag{11}
\]
The actual certificate is `NormSignCertificate.strong_error`, together with `S_growth` and `determinant_ne_zero`. Integrality has been proved at every index before any index is selected.

## 5. Two indices and the first pair suffice

Choose the least \(n\ge0\) with \(U\le400\cdot40^n\). If \(n\ge1\), minimality gives \(U>400\cdot40^{n-1}\). The determinant condition in (11) makes one of UT_n−VS_n and UT_{n+1}−VS_{n+1} nonzero. For its index j, equation (10) implies
\[
1\le|UT_j-VS_j|
<\frac{400}{956}+\frac{254\,4000^{n+1}}{U^3}
<\frac12+\frac{254\,4000^{n+1}}{U^3}.
\]
Consequently
\[
U^3<508\,4000^{n+1}.
\tag{12}
\]
If n≥3, then
\[
\frac{400^3\,40^{3(n-1)}}{508\,4000^{n+1}}
=\frac{16^{n-1}}{127}
\ge\frac{256}{127}>1,
\]
contradicting the preceding two bounds. Therefore n≤2. If n=0, U≤400; otherwise (12) gives
\[
U^3<508\cdot4000^3<32000^3.
\]
In all cases U<32000.

The first pair can be halved to (956,1), and (11) at n=1 gives
\[
|1-956\xi|<1/76480.
\]
The exceptional possibility U=956V is excluded directly:
\[
\Phi(956V,V)=-4569679V^4.
\]
Since U>0 would force V≥1, its absolute value exceeds 114244. Thus U−956V is a nonzero integer. Equations (10)–(11) now give
\[
\begin{aligned}
1&\le|U-956V|\\
&\le U|1-956\xi|+956|U\xi-V|\\
&<\frac{32000}{76480}+\frac{956\cdot254}{U^3}
<\frac12+\frac{956\cdot254}{U^3}.
\end{aligned}
\]
Hence
\[
U^3<508\cdot956=485648<79^3,
\qquad \boxed{U<79.}
\tag{13}
\]
There is no large-height calculation or rational-endpoint argument.

## 6. The remaining small values

For direct comparison with the checked small-sector lemma, (13) implies U≤599.

If \(V=0\), equation (6) gives
\[
U^{4}=28561=13^{4},
\]
so
\[
(U,V)=(13,0).
\]

Suppose \(V\ne0\). From (8) and (13),
\[
|V-\xi U|
\ge |V|-\xi U
>
1-\frac{599}{956}
>
\frac13.
\]
Using (9),
\[
|\Phi(U,V)|
>
450U^{3}|V-\xi U|
>
150U^{3}.
\]
Since \(|\Phi(U,V)|\le114244\), this gives \(U\le9\).

Now the stronger bound
\[
|V-\xi U|>1-\frac9{956}=\frac{947}{956}
\]
shows that \(U\ge7\) is impossible:
\[
450\cdot7^{3}\cdot\frac{947}{956}>114244.
\]
Therefore
\[
U\le6.
\]

The sector restriction \(|V|\le(\sqrt2-1)U\) leaves exactly the following nonzero possibilities:
\[
\begin{array}{c|c|r|r}
U&|V|&\Phi(U,|V|)&\Phi(U,-|V|)\\ \hline
3&1&-22916&22972\\
4&1&-57199&57521\\
5&1&\mathbf{-114244}&115196\\
5&2&-200719&200801\\
6&1&-199679&201841\\
6&2&-366656&367552
\end{array}
\]
Thus the only nonzero possibility is
\[
(U,V)=(5,1).
\]

We have proved (7).

## 7. Recovering \(x\) and \(y\)

Return to the norm relations in (6).

For \((U,V)=(13,0)\),
\[
13y=U^{2}+V^{2}=169,
\]
so \(y=13\), and then
\[
x^{2}=2\cdot13^{4}-1=239^{2}.
\]
Thus \(x=239\).

For \((U,V)=(5,1)\),
\[
26y=U^{2}+V^{2}=26,
\]
so \(y=1\), and then \(x=1\).

Both pairs satisfy the original equation. Therefore
\[
\boxed{
x^{2}+1=2y^{4},\quad x,y>0
\quad\Longleftrightarrow\quad
(x,y)=(1,1)\ \text{or}\ (239,13).
}
\qquad\square
\]


## 8. The positive endpoint and recovery of \(r,s,t,u\)

The same certificate gives the stronger statement used in the Vierzahlensatz branch:
\[
V>0,\quad U>478V,\quad 0<\Phi(U,V)\le28561
\quad\Longrightarrow\quad\text{a contradiction}.
\]
Indeed, the secant bound first gives \(U>956V\) and
\(0<\xi-V/U<32/U^4\). Since \(U>956\), the least \(N\) such that \(U\le400\,40^N\) is at least 1. Integer separation now gives \(U^3<64\,4000^{N+1}\). For \(N\ge2\),
\[
(400\,40^{N-1})^3=64\,4000^{N+1}\,16^{N-2},
\]
contradicting minimality and the strict upper bound, including the equality of scales at \(N=2\). Hence \(N=1\) and \(U\le16000\). The particular nonzero integer
\(UT_1-VS_1=2(U-956V)\) then gives
\[
U^3<64\cdot1912<956^3,
\]
contrary to \(U>956\). The full details and hypotheses are in [the positive-case proof](elementary_norm_sign_proof.md). This is `NormSignProof.positive_quartic_exclusion`; its specialization `remaining_28561_exclusion` needs only \(V>0\), \(956V<U\), and \(\Phi(U,V)=28561\).

For clarity, the screenshot theorem retains its exact source hypotheses. The variables \(p,q,a,r,s,t,u,y\) are positive integers, \(\gcd(s,t)=1\), and
\[
\begin{aligned}
rs&=7p-a-6,& rt&=4(a-1),\\
tu&=7p-a+2,& su&=2(a+5),\\
p^2-2q^2&=-1,& y^2&=p^2+2pq+2q^2.
\end{aligned}
\]
The [Vierzahlensatz normalization](../../../../../Archive/ljunggren_working_2026_09/derivation_r_s_t_u.md) gives positive integers \(R,H\) with
\[
r=8R,\quad u=8H,\quad s=R+3H,\quad t=3R+2H,\quad R^2-2H^2=-1,
\]
\[
p=R^2+4RH+2H^2,\quad q=R^2+2RH+2H^2,\quad
a=5R^2+4RH+2H^2,
\]
and
\[
y^2=5R^4+28R^3H+60R^2H^2+56RH^3+20H^4.
\]
The established arithmetic reduction sends every nonterminal normalized pair to the positive quartic conditions just excluded. The remaining pair is \(R=H=1\), and positivity selects \(y=13\). Substitution recovers
\[
\boxed{(r,s,t,u)=(8,4,5,8),\qquad (p,q,a,y)=(7,5,11,13).}
\]
The normalization is `Rstu.ScreenshotData.normalized_coordinates`; the reduction and terminal recovery are assembled in `NormSignProof.normalized_rstu_classification` and `NormSignProof.screenshot_classification`.

## Verification and scope

[CompleteNormSign.lean](../../../../../lean/Ljunggren/CompleteNormSign.lean) has compiled cleanly. Its `ljunggren_classification`, `ljunggren_iff`, and `ljunggren_nat_iff` classify the original equation from positivity and the equation alone. The root bracket, coefficient integrality, and transformed classification are proved internally, not left as hypotheses of these final statements. `quartic_classification` covers both values \(28561,-114244\), either sign of \(V\), and the integer sector condition \((U+|V|)^2\le2U^2\). `positive_quartic_exclusion` covers the full positive interval stated in Section 8, and `screenshot_classification` retains precisely the source hypotheses displayed there.

The certificate (11) is proved in [NormSignCertificate.lean](../../../../../lean/Ljunggren/NormSignCertificate.lean). The bounds \(U<32000\) and \(U<79\) are `FullNormSignTermination.height_32000_from_strong_certificate` and `quartic_height_79_from_strong_certificate` in [FullNormSignTermination.lean](../../../../../lean/Ljunggren/FullNormSignTermination.lean). The remaining finite sector checks are in [Bounds.lean](../../../../../lean/Ljunggren/Bounds.lean), and the exact norm rows and original-variable recovery are in [SectorReduction.lean](../../../../../lean/Ljunggren/SectorReduction.lean).

The human sign proof uses the Wronskian identity checked by `NormSignActual.actual_wronskian`. The generic Lean sign mechanism uses the equivalent polynomial identity
\[
(1+a^2)WF'-\{2(1+a^2)W'+aW\}F=4c\Delta\Psi,
\]
with the same degree comparison and largest-zero argument. Both identities and the actual all-index brackets have been checked. The real ingredients are ordinary polynomial continuity, a finite root set, and leading-term signs; no generalized binomial estimate, complex-number error construction, Jacobi symbol, quadratic reciprocity, or prior classification is used as an input. The existing integer-coefficient construction remains part of the proof; the collaborator's optional alternative integrality proposal is not adopted here.

The authorized [ChatGPT collaborator](https://chatgpt.com/c/6aafeb87-357c-83ee-a6dd-6af6b789f062) independently reviewed the sign, bracketing, and base-40 termination argument. Its [preserved review](../../../../../Archive/ljunggren_working_2026_09/collaborator_norm_sign_review.md) distinguishes that endorsement from optional proposals. The combined build and all six recursive proof-dependency audits passed; [NormSignProofAudit.lean](../../../../../lean/NormSignProofAudit.lean) is the audit source. The [verification receipt](../../../../../lean/FINAL_AUDIT.md) records the project's published audit status.

The immediately preceding full real-coefficient manuscript is preserved byte-for-byte in [its archive](../../../../../Archive/ljunggren_working_2026_09/elementary_proof_real_coefficient_archive.md). The earlier binomial proof remains in [its separate archive](../../../../../Archive/ljunggren_working_2026_09/elementary_proof_binomial_archive.md).
