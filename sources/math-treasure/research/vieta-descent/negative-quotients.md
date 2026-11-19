# All quotients in the negative-q family

**Status: checked written proof (25 October 2026); see [review](negative-quotients-review.md).**
This completes the varying-quotient obligation for the negative-denominator-offset family in Ishan Singh's unfinished negative-q notes. Integral reflections and the owner programme are recovered source work. Finiteness and effective classification have a classical antecedent in Mills's more general rational quadratic problem; this is a self-contained, explicit bound and implementation, not a claim of novelty. See [source audit](source-audit.md). The previously checked [signed-level proof](signed-levels.md) supplies positivity, descent, uniqueness of the reduced seed and the matrix orbit normal form used below.

## 1. A sharp product bound

Let a,b,c,N be positive integers, a|c and b|c, and suppose

F(x,y)=ax²+by²−cxy=−N,  x,y positive integers.

Write Δ=c²−4ab. A solution forces Δ>0 by completing the square. A reduced solution satisfies 2ax≤cy and 2by≤cx. Multiplication and the equation give the exact identity

(cx−2by)(cy−2ax)
= (c²+4ab)xy−2c(ax²+by²)
= 2cN−Δxy.

Both factors on the left are nonnegative at a reduced point, so

xy ≤ 2cN/Δ = (2N/c)·d/(d−4) ≤ 10N/c,

where d=c²/(ab)=(c/a)(c/b) is an integer at least 5. The last inequality is equivalent to d≥5. The sharper rational inequality retains the coefficient dependence; the uniform constant 10 is sharp. Indeed (a,b,c,N,x,y)=(1,5,5,1,2,1) satisfies F=−1, is reduced (4≤5 and 10≤10), and has xy=2=10N/c. The first inequality is equality precisely when at least one reduced inequality is equality. Equality in the uniform bound additionally requires d=5. No claim of a uniform 10 holds here without the integral-reflection hypotheses.

## 2. Eliminate the unknown quotient from the search bound

Fix arbitrary positive integers a,b,p,q. We classify all positive integer triples (x,y,k) with positive denominator and

(ax²+by²)/(abp·xy−q)=k.                                      (1)

There is no assumption q≤abp and no coprimality assumption. Setting c=kabp and N=kq turns (1) into F=−N, with a|c and b|c. The signed-level theorem therefore gives a unique reduced seed in each positive reflection orbit, keeping k fixed during descent. At that seed the preceding bound becomes

xy ≤ 10kq/(kabp)=10q/(abp).

Thus B=floor(10q/(abp)) bounds the product of seed coordinates **independently of k**. This cancellation closes the fixed-k limitation for precisely family (1).

The following terminating algorithm returns all quotients and all their ordered reduced seeds:

```
B = (10*q) // (a*b*p)
for x = 1,...,B:
    for y = 1,...,B//x:
        D = a*b*p*x*y - q
        if D <= 0: continue
        numerator = a*x*x + b*y*y
        if numerator % D != 0: continue
        k = numerator // D
        c = k*a*b*p
        if 2*a*x <= c*y and 2*b*y <= c*x:
            retain (k,x,y)
```

Every retained k is positive because numerator and D are positive. The enumeration is finite, with exactly Σ_{x=1}^B floor(B/x) tested ordered pairs, interpreted as zero when B=0. All divisions and comparisons are integral and exact. The algorithm neither scans nor guesses an upper bound for k.

Under the owner’s original restriction q≤abp, B≤10. Consequently at most
Σ_{x=1}^{10} floor(10/x)=27 ordered candidate pairs suffice for every fixed
parameter tuple in the original source regime. The endpoint product 10 is
attained there: (a,b,p,q,k,x,y)=(1,5,1,5,1,5,2) has numerator and denominator
45, is reduced with 2ax=cy, and has xy=10.

**Completeness proof.** Any original solution has Δ>0. At every positive solution the two integral reflections X(x,y)=(kbp·y−x,y) and Y(x,y)=(x,kap·x−y) preserve the equation and positivity, because their root products are (by²+kq)/a and (ax²+kq)/b. If the point is not reduced, its unique decreasing reflection strictly decreases x+y. This reaches a reduced seed, and the uniform bound places it inside the finite loop. The denominator there is positive: the preserved equation gives D=(ax²+by²)/k>0 at every orbit point. Thus every original solution has a quotient appearing in the output and lies in the orbit of a retained seed. Conversely, each retained tuple satisfies (1), and every word in X,Y preserves positivity, the numerator identity and the positive denominator. Consequently the output gives exactly the possible quotients, and the complete solution set is the union of their positive reflection orbits. Uniqueness of the reduced seed is the maximal-height argument of the signed-level proof; different retained seeds for the same k represent distinct orbits. ∎

Explicitly put m=kbp, ℓ=kap and

A=XY=[[mℓ−1,−m],[ℓ,−1]],  X=[[-1,m],[0,1]].

For each retained seed z, all its solutions are A^j z and XA^j z, j∈Z. These matrices have the signed-level recurrence formulas already proved. Negative powers of A are integral because det A=1. No positivity filter or further denominator filter is needed on these genuine negative-level orbits. Different words may represent the same point when a reflection fixes it; only the reduced seed, not its word parameter, is unique.

## 3. Complete small examples and exclusions

For a=b=p=1:

| q | All possible k | Reduced ordered seeds by k |
|---|---|---|
| 1 | 5 | 5: (1,2), (2,1) |
| 2 | 4, 10 | 4: (2,2); 10: (1,3), (3,1) |
| 3 | 3, 4, 8, 13, 17 | Listed in the exact verification record |

These are whole-quotient classifications because the finite loop is proved exhaustive, not merely because a coordinate box was searched. The q=2 seed (2,2) must be retained: a gcd-one filter would incorrectly erase quotient 4. More generally scaling coordinates by t and q by t² preserves the quotient and reducedness, illustrating why primitive normalization must track q.

The uniform constant remains sharp within (1): (a,b,p,q,k,x,y)=(1,5,1,1,1,2,1) has denominator 9 and numerator 9. If B=0 the classification is empty by the proved bound; for instance (a,b,p,q)=(4,4,1,1). No promise that B>0 implies existence is made.

The family with a positive denominator offset, or a general varying cross coefficient unrelated to the level, is outside this cancellation argument. Positive-level duplicate orbit representatives and nonintegral-reflection lattices retain their previously stated limitations. No novelty conclusion follows from these computations or from the limited source comparison.

## 4. Executable classifier and finite diagnostics

Run, for example:

```
python3 experiments/vieta-descent/negative_quotients.py 1 1 1 2
python3 experiments/vieta-descent/verify.py
```

The first command emits all quotients and reduced seeds for the supplied positive integers a,b,p,q. The second recomputes and compares the stored deterministic record. Intentional regeneration uses `--write-record`. The record's `negative_quotients` section gives exact ranges, counts, complete example seeds and limitations.

Diagnostics enumerate a,b=1,…,6, p=1,…,4 and q=1,…,24 (3,456 tuples). For each tuple they execute the entire proved product-bound algorithm; independently compare every found quotient's seeds against the fixed-level rectangle; descend every integral-quotient positive-denominator solution with x,y=1,…,80 and regenerate it by reversing the descent; and verify A^j and XA^j for j=−8,…,8 at every retained seed, including denominator positivity. Scaling checks use a=b=p=1,q=1,…,8,t=1,…,4. Existing positive, negative and zero-level diagnostics are retained, with the product identity additionally checked at every existing negative reduced seed. These finite diagnostics are not substitutes for the universal proof.
