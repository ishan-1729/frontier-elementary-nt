# All positive quotients for an unrestricted positive offset

T026/R197, 2026-10-26. **Independently approved in R198; see [review record](review.md).** This resolves the mathematical scope of ITEM-VIETA-POSITIVE-OFFSET as an explicit reconstruction beyond the owner's q≤abp manuscripts; it does not claim that those manuscripts already prove the unrestricted result.

Fix positive integers a,b,p. Seek positive integers x,y,k such that

    (a x²+b y²)/(abp xy+q)=k.                               (1)

Initially q may be any positive real number. Any solution gives q=(a x²+b y²)/k−abp xy∈Q. Thus irrational q has no solutions. From now on write q=r/s>0 in lowest terms, r,s positive integers. Equation (1) gives C=kq=a x²+b y²−kabp xy∈Z, so s|k. In particular all integer q>0, including q>abp, and k=1 are covered.

## Fixed-k boundary descent, with every sign retained

For fixed positive integer k put c=kabp, m=c/a=kbp, l=c/b=kap and F(x,y)=a x²+b y²−cxy. The [reviewed positive-level boundary theorem](../vieta-descent/README.md) applies whenever a,b,c,C are positive integers, a|c and b|c. Every solution of (1) supplies precisely these hypotheses; if s∤k there are none. We use its arbitrary-positive-C theorem, not the smaller axis barrier. The complete argument needed here follows.

Set X(x,y)=(my−x,y), Y(x,y)=(x,lx−y). Direct substitution shows that both are integral involutions preserving F. They preserve gcd, since adding a multiple of one coordinate to the other preserves their common divisors. For x,y>0 with F=C, put x'=my−x,y'=lx−y. Then

    a xx'=b y²−C,           b yy'=a x²−C.

If x'≥x and y'≥y, adding the resulting inequalities would give −2C≥0, impossible. At least one flip strictly decreases its changed coordinate. Choose one such flip. If it would produce a negative coordinate, **stop before applying it**. Otherwise perform it, strictly decreasing the integer height x+y; if a coordinate becomes zero, stop on that axis. At any remaining positive point another decreasing choice exists. Thus this procedure terminates after at most its initial height many executed flips. No assertion that all flips stay positive is needed or true in this generality.

An axis endpoint is (u,0) or (0,v), with u or v positive. If an X flip would instead produce (−u,v), u,v≥1, the current point is (mv+u,v); if a Y flip would produce (u,−v), its current point is (u,lu+v). In every case the nonnegative boundary representation is

    a u²+cuv+b v²=C,          (u,v)≠(0,0).                 (2)

Conversely an axis representation gives positive witnesses Y(u,0)=(u,lu) or X(0,v)=(mv,v). A mixed representation u,v≥1 gives **both** positive witnesses

    z_X=(mv+u,v)=X(−u,v),      z_Y=(u,lu+v)=Y(u,−v).       (3)

Expansion proves their F-value is C, and their gcd is gcd(u,v). Both orientations are retained; identifying them would require an orbit argument not provided or needed here. Reverse the finite list of executed involutions to recover any original positive solution from its endpoint seed. This proves boundary coverage and its converse, with no uniqueness assertion and no primitive-seed restriction. It reproduces the applicable reviewed proof rather than importing a negative-level reduction.

## A finite mixed list independent of k

Substitute c=kabp,C=kq into (2):

    k(q−abp uv)=a u²+b v².                                 (4)

For mixed u,v≥1 the right side is positive, so uv<q/(abp). In reduced rational notation this is s abp uv<r, or exactly

    uv≤B,       B=floor((r−1)/(s abp)).                     (5)

The integer r−1 is essential at a strict endpoint. Enumerate all ordered pairs u=1,…,B and v=1,…,floor(B/u), computing

    k=s(a u²+b v²)/(r−s abp uv).                            (6)

Its denominator is a positive integer. Keep precisely integral values; they are automatically positive. Each retained value has both positive witnesses (3), so every retained quotient is genuine. Conversely every mixed endpoint of every solution appears in this list by (4)–(5). The number of ordered candidates is exactly Σ_(u=1)^B floor(B/u), zero when B=0. For B≥1, u,v≤B and the positive denominator in (6) is at least 1, hence

    k≤s(a+b)B².                                            (7)

This is a bound on mixed candidates, not on the axis families or all possible quotients. Neither coprimality of a,b nor primitivity of u,v is used.

## Two explicit infinite axis families

For either weight w=a (x-axis) or w=b (y-axis), the axis condition is w z²=kq, or k=s w z²/r. Since gcd(r,s)=1, its integrality is equivalent to r|w z². Define

    h_w=product_(prime P|r) P^ceil(max(v_P(r)−v_P(w),0)/2),
    K_w=s w h_w²/r.                                        (8)

Prime valuations prove r|w z² exactly when h_w|z. They also prove K_w is a positive integer. Consequently all axis quotients and coordinates are exactly

    z=h_w t,       k=K_w t²,        t=1,2,… .                (9)

For r=1 the product is 1. This includes k=1 when K_w=t=1. Each family gives positive original solutions by the corresponding axis witness above, for every t. Thus these are infinite families of distinct quotient values, regardless of whether a particular fixed-k orbit is finite. Formula (6) and formula (9) also imply s|k: in (6), reducing k(r−s abp uv)=s(a u²+b v²) modulo s gives kr≡0, and gcd(r,s)=1.

**Classification theorem.** All possible positive integer quotients in (1) are exactly

    {K_a t²:t≥1} ∪ {K_b t²:t≥1} ∪ {integral values in (6)}. (10)

Necessity follows by the terminating boundary descent and (4)–(9). Sufficiency follows from the explicit positive witnesses. This proves the unrestricted all-k result, not just a finite experiment or a fixed-k reduction. Quotient values may be deduplicated as numbers. Keep both labelled axis families and every mixed pair/orientation with its gcd label when constructing solutions; overlaps do not justify deleting seed information.

If q≤abp then B=0 and only axis families remain, recovering the sufficient regime of the owner manuscripts. Beyond it mixed families may occur. Even when B>0, a given mixed quotient can coincide with an axis quotient or another mixed quotient; the theorem does not assert disjointness or canonical orbit representatives.

## Gcd labels and normalization

Every orbit preserves g=gcd(x,y). Axis witnesses have g=z; mixed witnesses have g=gcd(u,v), which may exceed 1. Dividing coordinates x=gx_0,y=gy_0 in (1) yields

    (a x_0²+b y_0²)/(abp x_0y_0+q/g²)=k.

The offset changes to q/g². Retaining q unchanged would lose the original equation; we do not normalize away gcd labels in the classifier.

## Fully evaluated examples

For a=b=p=1,q=3, r=3,s=1,B=2. Both weights have h=3,K=3, so both labelled axis families give k=3t². The only ordered mixed candidates are (u,v)=(1,1),(1,2),(2,1). Formula (6) gives respectively k=2/(3−1)=1, 5/(3−2)=5 and 5/(3−2)=5. Therefore the complete quotient set is

    {3t²:t≥1} ∪ {1,5}.

The k=1 mixed witnesses are (2,1),(1,2); its finite orbit must not be described as an infinite family of distinct positive pairs. At k=5 the two ordered pairs provide witnesses (11,2),(1,7) and (7,1),(2,11). Axis examples start with k=3 and points (3,9),(9,3). All are checked directly by (1).

For a=b=p=1,q=8,k=2, the mixed pair (2,2) satisfies k=(4+4)/(8−4)=2 and yields (6,2),(2,6), both of gcd 2. Removing nonprimitive mixed pairs loses these positive solutions. The other mixed pairs at this k are exactly (1,3),(3,1): (2) becomes (u+v)²=16, so positive integral u,v must sum to 4. The axis equation z²=16 also supplies both axes with z=4; their positive witnesses have gcd 4, distinct from the mixed gcd-2 witnesses. This last distinction shows why mere quotient membership cannot replace the labelled seed list. Dividing (6,2) by 2 changes q to 2: (3,1) has quotient 2 for q=2, but quotient 10/11 for q=8.
