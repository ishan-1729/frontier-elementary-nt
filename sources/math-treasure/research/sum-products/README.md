# Sum-products: recovered statements and proof draft

T026/R09, 2026-10-25. **Independently reviewed written proof (R13, 2026-10-26).**
See [review scope](review.md); this is not formal verification. Source images
and all six local text notes were personally read; see `sources.json` for
current paths, catalog identities, byte counts and SHA-256 hashes. These are
derived proofs of the owner's recorded questions, not claims of new priority.
The Capture image credits the binary coefficient observation to Gatsby.
The two Mysterious README variants differ only by numbered link markers;
both were retained and compared using Math/VERSIONS.html (V-41f56c7f74).

## Indexed subset generating functions

Let a_1,...,a_n be positive integers, with repeated values allowed, and
u_i,v_i in a commutative ring R. Subsets below are subsets of **indices**.
Distributivity, choosing one term from each factor, proves

    ∏(u_i+v_i y x^{a_i})
      = Σ_{I⊆{1,...,n}} (∏_{i∉I}u_i)(∏_{i∈I}v_i) y^{|I|}x^{Σ_{i∈I}a_i}.

Thus the coefficient of x^k y^m is the weighted count of m-element
representations of k. This repairs the source's ambiguous set X_n notation:
if equal a_i are collapsed into a literal set, the indexed formula changes.
For u_i=v_i=1 write P_{m,n}(x)=[y^m]∏(1+yx^{a_i}). Separating subsets
according as n belongs gives

    P_{m,n}=P_{m,n-1}+x^{a_n}P_{m-1,n-1}.

Initial conditions: P_{0,0}=1 and P_{m,0}=0 for m≠0; also P_{m,n}=0
for m<0 or m>n. Consequently f_k(m,n)=f_k(m,n-1)+f_{k-a_n}(m-1,n-1),
with negative k counts zero. Summing m gives the owner's univariate identity.

In the original polynomial ∏_{j=0}^{10}(x^{2^j}+2^j), an exponent
selects its unique binary subset. Since 2047−2012=35=1+2+32, the omitted
positions are 0,1,5 and the coefficient is 2^0·2^1·2^5=64. All other
positions contribute coefficient 1. This proves the recovered screenshot
problem, without purporting to recover the missing Kainui problem.

## Completeness, uniqueness, and infinite products

For a nondecreasing infinite sequence of positive integers, every nonnegative
integer is a sum of finitely many distinct indexed terms iff a_1=1 and

    a_j ≤ 1+S_{j-1},  S_j=Σ_{i=1}^j a_i   (j≥2).

Necessity: if a_1>1, 1 is missed. If a_j>S_{j-1}+1, earlier terms cannot
reach S_{j-1}+1 and all later terms exceed it. Sufficiency: inductively,
subset sums of the first j terms cover [0,S_j], since the two integer
intervals [0,S_{j-1}] and [a_j,a_j+S_{j-1}] touch or overlap. S_j≥j tends
to infinity. This is the classical completeness criterion, not a new result.

Coverage **and** uniqueness of all finite indexed subset representations
force a_j=2^{j-1}. Indeed the covered initial interval makes a_j≤S_{j-1}
a duplicate representation (singleton j versus earlier indices); completeness
then forces a_j=S_{j-1}+1. Starting at 1 proves the assertion inductively.
Conversely the disjoint adjacent intervals prove uniqueness for powers of 2.
Completeness alone allows 1,1,1,... and hence does not imply uniqueness or
locally finite representation counts.

For the unweighted formal product ∏(1+yx^{a_i}) in Z[y][[x]], require
local finiteness: for each M only finitely many i have a_i≤M. Then every
coefficient stabilizes and equals the finite representation count. For
positive integers this condition is also necessary: infinitely many terms
below M imply infinitely many equal terms at some exponent ≤M, so its
coefficient of y cannot stabilize. No infinite formal product of the
unnormalized constants u_i is being asserted.

An exact sufficient analytic hypothesis for the normalized weighted product
∏(1+c_i y z^{a_i}) is Σ_i |c_i|ρ^{a_i}<∞ for a fixed ρ>0. On |z|≤ρ,
|y|≤Y<∞ the absolute subset expansion is bounded by
∏(1+Y|c_i|ρ^{a_i})≤exp(YΣ|c_i|ρ^{a_i}). Tail estimates prove uniform
convergence and absolute summability of the expansion, permitting regrouping
by degree and cardinality. On open disks where these bounds hold locally,
the product is analytic. Local finiteness alone does not imply this analytic
condition: multiplicities can grow faster than exponentially. For binary
weights and c_i=1, |z|<1 gives convergence and the formal/analytic identity
∏_{j≥0}(1+z^{2^j})=1/(1−z), from finite geometric sums and their limit.
No boundary-value convergence at |z|=1 is claimed.

## Alpha: finite mixed radix identity

For integers r_i≥2, 1≤i≤N, let Q_i=r_1...r_i and Q=Q_N. In Z[t],

    ∏_{i=1}^N Σ_{d=0}^{r_i−1} t^{d Q/Q_i} = Σ_{k=0}^{Q−1}t^k.

Successive Euclidean division starting with the last radix gives each integer
0≤k<Q a unique expansion Σd_i Q/Q_i: the last digit is k mod r_N;
dividing the remainder by r_N reduces to N−1 digits. Expanding the product
therefore counts each exponent exactly once. With r_i=r and t a chosen
r^N-th root of x, this is the source Alpha formula, PROVIDED x^{d/r^i}
means t^{d r^{N-i}}. For positive real x the usual positive powers work;
for nonzero complex x a single chosen logarithm and x^q=exp(q Log x)
works. Independently chosen roots need not work. At x=0 use the polynomial
interpretation t=0 and t^0=1. No quotient division by t−1 is needed.

A concrete trigonometric product consequence, without a missing Brilliant
statement, follows by putting t=e^{2iθ} and summing geometric progressions:

    ∏_{i=1}^N sin(r_i Qθ/Q_i)/sin(Qθ/Q_i) = sin(Qθ)/sin θ,

where all displayed denominators are nonzero. Phases cancel since
Σ(r_i−1)Q/Q_i=Q−1; removable cases follow by limits. This is one proven
connection to finite sine products, not an attribution to Jack Lam's solution.

## Omega: Dirichlet inverse proof

Let n≥2 and a,b be integers with 0≤b≤a≤n. These integer restrictions
are reconstructed hypotheses: the Omega image does not explicitly state
integrality of a,b. Its auxiliary cosine formula also omits the necessary
finite integer k range supplied below. Extending the displayed Omega claim
to real a,b fails: n=2, a=b=1/2 gives left side 1/2 and right side −3/4.
Set θ_m=mπ/n. Then

    Σ_{m=1}^{n−1} sin(aθ_m)sin(bθ_m)/(1+cos θ_m)
       = (−1)^{a+b} b(n−a).

For interior indices define L with diagonal 2 and adjacent entries −1.
Its inverse is H_{ij}=min(i,j)(n−max(i,j))/n: extend H_{0j}=H_{nj}=0;
it is linear on each side of j, and its downward slope jump at j is 1.
Hence 2H_{ij}−H_{i−1,j}−H_{i+1,j}=δ_{ij}, proving LH=I directly.
With D_{ii}=(−1)^i, A=DLD has diagonal 2 and adjacent entries +1,
and inverse K=DHD.

The vectors v_m(i)=sin(iθ_m) are eigenvectors of A with eigenvalues
λ_m=2+2cos θ_m>0, by the sine addition formula and zero boundary values.
They have pairwise distinct eigenvalues, so symmetry of A gives orthogonality.
Their squared norms are n/2, since
Σ_{i=1}^{n−1}cos(2imπ/n)=−1 (sum the nontrivial n-th roots of unity).
There are n−1 such nonzero orthogonal vectors, hence they form a basis and

    K_{ab}=(2/n)Σ_m sin(aθ_m)sin(bθ_m)/[2+2cos θ_m].

Multiplication by n gives the claim. Boundary a=n or b=0 yields zero on
both sides and was not inferred from an interior inverse entry.

## Omega: corrected cosine polynomial and Chebyshev recurrence

For integer 0≤k≤2n the actual formula is

    T(k)=Σ_{m=1}^{n−1} cos(kθ_m)/(1+cos θ_m)
        =(−1)^k [k²/2−nk+(4n²−1)/12]−1/4.

Extend it by period 2n, not by the same unbounded quadratic. The proposed
reconstruction with constants (2n²+1)/6 and −1/2 is false: n=2,k=1
would give −1/2 instead of 0. The source image specifies only unnamed
constants, so this correction preserves the source statement.

Here is a full recurrence proof also interpreting the Chebyshev request.
Define C_k=Σ_{m=1}^{n−1}cos(kmπ/n). Pair m with n−m for odd k, and use
the geometric sum for even k. Thus C_0=n−1, C_k=0 for odd k, and
C_k=−1 for even 0<k<2n. Chebyshev's recurrence
T_{k+1}(z)+T_{k−1}(z)=2zT_k(z), with T_k(cos θ)=cos(kθ), implies

    T(k+1)+2T(k)+T(k−1)=2C_k.

On the cyclic indices modulo 2n the proposed F(k) satisfies this recurrence:
for 1≤k≤2n−1 its quadratic second difference gives −(−1)^k−1=2C_k;
at 0, using F(2n−1)=F(1), it gives 2n−2. Also F(2n)=F(0).
The kernel of this cyclic recurrence is spanned by (−1)^k: the general
homogeneous solution is (A+Bk)(−1)^k and periodicity forces B=0.
To fix that final ambiguity, the cosine sum has zero alternating mean:
Σ_{k=0}^{2n−1}(−1)^k cos(kmπ/n)=0 for 1≤m<n, by geometric summation.
The same is true of F because

    Σ_{k=0}^{2n−1}(k²/2−nk)=−2n(4n²−1)/12,
    Σ_{k=0}^{2n−1}(−1)^k=0.

Thus F=T. This determines all four constants, including the required
finite integer range. Finally 2sin u sin v=cos(u−v)−cos(u+v) gives
[T(a−b)−T(a+b)]/2=(−1)^{a+b}b(n−a), an independent recurrence route
to Omega with the necessary integer domain made explicit.

## Beatty assertion and unresolved bridges

The complete-sequences note's Beatty assertion is valid for positive
irrational α,β with 1/α+1/β=1 (in particular both exceed 1). The number
of terms floor(jα)≤N is floor((N+1)/α), since (N+1)/α is nonintegral.
The analogous two counts sum to N, since their arguments are nonintegers
whose sum is N+1. Therefore the total multiplicity at each positive integer
N, obtained by subtracting counts at N−1, is exactly one. This proves both
disjointness and coverage. Positivity cannot be dropped silently.

OPEN: Kainui's actual linked identity; Monster Limit 2; the specific Brilliant
family-of-products and Jack Lam problem/solution; the requested broader Euler
reflection connections; prime/quadratic-form, Goldbach and Landau bridges;
and a specified general complete-sequence multiplication with associativity
or group/ring/module properties. Completeness and a subset generating function
alone prove none of those claims. The [2] Sum To Product note is an organization
request and video pointer, not a stated Rogers–Ramanujan theorem to discharge.
No rule on arbitrary representations was invented.

## Attribution, search scope, and reproducibility

The local notes motivate these reconstructions. The subset expansion, binary
and mixed-radix mechanisms, completeness criterion, Beatty partition and
finite difference inverse are classical mechanisms; novelty is not claimed.
A limited web search on 2026-10-25 used “complete sequences Brown criterion
1961 note complete sequences” and “Knuth Fibonacci multiplication 1988
Fibonacci number system multiplication”. Brown's own later paper states the
criterion and cites his 1961 paper: [Brown, 1978, introduction](https://www.fq.math.ca/Scanned/16-1/brown1.pdf).
Only that introductory attribution is imported; the proof above is independent.
The [Knuth PDF](https://www.cs.umb.edu/~rvetro/vetroBioComp/compression/FIBO/KnuthFibb.pdf)
was identified but fetching timed out in this worker; no theorem from it is
used. Peeranswer and Brilliant family-of-products were inaccessible through
the web tool; other linked Brilliant/video contents were not read here.
This is not an exhaustive equivalence or priority search for Omega or Alpha.

Run `python3 experiments/sum-products/verify.py` and
`git diff --check -- research/sum-products experiments/sum-products`.
The verifier uses exact integers/Fractions, enumerates specified finite controls,
and compares with its saved record; it does not replace the general proofs.
