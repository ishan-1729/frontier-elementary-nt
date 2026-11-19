# Review record: mordell_k_minus17

All reviewers below are Claude subagents of the same model family as the author. So this is in-session review, not
an independent audit. After round 2 the card was WRITTEN (Theorem A) and PARTIAL (Theorem B, under hypothesis (H)).
Round 3 reviewed the proof of (H) in [hypothesis_H.md](hypothesis_H.md); since then both theorems are WRITTEN
(in-session review).

## Round 1: scripts and independent verification (2026-10-15)
Two workflows ran, with ten agents (five per curve; this record covers the 10404.f2 agents, and the literature
agent). Each worked alone, wrote its own script, and reported in a structured form. The author wrote
[core.py](scripts/core.py), [01_coverage.py](scripts/01_coverage.py) and [02_exclusion.py](scripts/02_exclusion.py)
from the template's, with new data; [lifting2.py](scripts/lifting2.py) is a verbatim copy.

| Agent | Task | Verdict | What it found |
|---|---|---|---|
| fields | wrote [00_fields.py](scripts/00_fields.py) | 158 checks pass | Fields, maximal orders (also by sympy's round_two), markers, the substitution lemma for classes (0,1) and (1,1), units, known solutions. The three fields are pairwise non-isomorphic, and the resolvent cubic of each is a shift and scaling of y³ − 17, so its field is Q(17^{1/3}). In class (0,1), O = Z[β] and β = −ε₁ is a unit. |
| reconstruction | wrote [03_reconstruction.py](scripts/03_reconstruction.py) | 6435 checks pass | The 2-descent re-derived. Classes (1,0), (0,1): u² ≡ F_cl (mod 4) and g = √c automatically. Class (1,1): g = √c unless 3 \| c and 3 \| n, which no solution of Theorem A meets. The trivial class gives x = 8 only; its solution ±(1,1) of F₀ = 9 gives the rational point (−8/9, ±109/27). 16 points. The template's integrality argument for X and its bound on g do not apply here, and are replaced. |
| coverage | wrote [verify_coverage.py](scripts/verify_coverage.py), no project imports | 210 checks pass | Own Sturm and Krawczyk enclosures, own polygon by clipping. The linear bound of Remark 3.2(ii) is used where the enumeration is feasible; elsewhere a subdivided linear bound, which gives a log-free weak form of step 4 (π ≥ 0.9·b_min; wording corrected in round 2). One orbit per seed of 01 (17); the unit pairs are fundamental. |
| exclusion | wrote [verify_exclusion.py](scripts/verify_exclusion.py), no project imports | all checks pass | Every number printed by 02 reproduced, on all 17 orbits; an order-free check of all 10,554,462 classes of the largest lattice leaves none unaccounted. The λ = 1 instance (8, (1,2)) checked against Proposition 17's statement and proof. Negative controls bite. |
| brute force and audit | wrote [verify_bruteforce.py](scripts/verify_bruteforce.py) | all checks pass | Exact box \|a\|,\|b\| ≤ 200 and modular box ≤ 1000: only the certified zeros. Direct search \|m\| ≤ 10⁶ with rigorous windows: only the known solutions of G and of F_cl. 1200 planted two-zero systems refused by Proposition 17. |
| literature | wrote [literature_check.md](literature_check.md) and the computations of its Appendix A, [lit_checks_17.py](scripts/lit_checks_17.py) | the script's checks pass | (H) is a published theorem, Cassels (1950, §22); the published proofs (Nagell 1930 onwards) take the cubic route or elliptic logarithms; the quartic two-unit route was not found in print; Bennett–Ghadermarzi cite Ljunggren's 1963 paper with a wrong title, year and pages; two target screens omitted Nagell and Hemer. The README and proof.md §6.1 then cited Cassels, and the screens were corrected. |

**Fixes made after round 1.**
1. *Printing.* [01_coverage.py](scripts/01_coverage.py) printed three positive vertex coordinates of class (1,1) as
   0.000000 (the least is 3.37·10⁻⁹). It now prints them in exponent form, with the least coordinate. Positivity was
   always checked exactly.
2. *A checkable certificate.* [02_exclusion.py](scripts/02_exclusion.py) now prints, for each modulus used, the basis
   of Λ_q; after each step, the Hermite form of the current lattice; each certified class with its lattice; and each
   sibling class with its witness. It also asserts that every known solution lies in exactly one orbit.
3. *Wording.* 02's line "each pair ±z meets the family at most once" now cites the injectivity from the cone
   certificate as well. 01's docstring states the lemma for c ≠ 0.
4. *proof.md.*
   - The template's "λ = 0 at every zero used" is false here: one instance, (8, (1,2)) in class (1,0), has λ = 1 and
     ν(J) = 2. It is valid (both hypotheses hold with equality, and the strict inequalities of the proof read 4 < 5),
     and the basis remark is restated: a positive verdict of Proposition 17 never depends on the basis.
   - Remark 3.2(ii) now says where verify_coverage.py uses the linear bound and where it uses the subdivided one.
   - Table 4.5 lists the moduli of the class (1,1) orbits in seed order.
   - The margins are stated exactly: all eleven instances have δ = e − 1 − ι.
   - README: the bound 7.11√c was rounded to 7.2√c; now 7.12√c.
5. *Re-runs.* After the changes, 01 and 02 were re-run (only the printed format changed, and the new lines), and the
   three verifiers were re-run against the new outputs. Their parsers of 02's output were adapted to the new line
   format.

**Facts recorded by the verifiers.**
- *Phantoms.* In orbit (0,1), c = 1, −1 ≡ ε₁⁸¹ (mod 27): the three classes left after q = 27 are the phantoms of the
  three known zeros. They die as siblings at q = 7 and q = 19, where −1 is not a residue.
- *The same-modulus sibling rule is needed.* In class (1,1), ε₁² = 1 − 3β, so the two zeros of the orbit c = 1 share a
  Λ₃-class; a cross-modulus sibling rule would accept an incomplete certificate.
- *Moduli.* In orbit (0,1), c = 1, [27, 7, 19] suffices; in orbit (1,1), c = 1, [9, 5, 4] suffices.
- *No element of norm 4 in class (0,1)*, also for a reason independent of the enumeration: O = Z[β] and
  g ≡ t⁴ + t + 1 (mod 2) is irreducible, so 2 is inert.

**Process notes.** Several agents ran one read-only `git log`, which their briefs forbade; nothing changed. The author
converted some files from CRLF to LF line endings while the agents worked; the agents checked that the content was
unchanged.

All scripts whose inputs changed were re-run after the fixes, and all pass (outputs in [scripts/](scripts/)).

## Round 2: referees on the text (2026-10-15)
Three referees worked alone on this card, each with one lens: coverage (§§2–3 with core.py, 00 and 01), exclusion
(§4 with lifting2.py and 02, re-run) and global (§§0, 1, 5–7, the README, the literature check, this record and the
shared files changed for the card). triangular_square_pyramidal had three referees with the same briefs; the six ran in
two workflows. Each re-ran the scripts it needed and wrote its own scratch checks, among them:
- an 80-digit recomputation of table 3.3, an own enumeration of the lattice points, and 680 greedy descents;
- a re-implementation of §4 from the tables of proof.md alone, Proposition 17 in five bases, and a direct check of
  every class of the final lattices of all 17 orbits;
- an own re-derivation of §1 and §5, and an enumeration of the binary cubic forms of discriminant −1836;
- page images of Cassels (1950) and Hemer (1956), and the reviews and records of the other sources.

**Result.** 29 findings: 0 critical, 1 major, 28 minor. All three verdicts were "accept with minor revisions". No
finding changes Theorem A′, A or B or any certificate, and all were fixed.

**The major finding: what Cassels's proof of (H) uses.**
- The card said that Cassels's proof of the class list uses his Lemma 6, which rests on class field theory. The
  global referee read pp. 249–251 and 261–262: at D = −17, §22 takes the list "by the argument of § 12", and §12 uses
  only ideal factorization (his Lemma 7), the class group and the units. Lemma 6 serves other D in the paper, for
  example D = +17 (p. 262).
- The author checked this on the page images and the text layer. proof.md §6.1 and §6.5, the README, the literature
  check §3 and GENERALIZATIONS §4 item 7 are corrected. What an elementary proof of (H) has to replace is ideal
  factorization in Q(∛17), with h = 1 and the unit η = 18 + 7θ.

**Minor findings, by group (all fixed).**
- *Notation and wording in §§2–4*, as on the sister card: Γ for the unit group, a₀..a₃ named, T₁ defined, Remark 3.2(i)
  completed, the logarithm-free route written out, §0.3's list corrected, r and (s, s') in §4, and q₁, …, q_s.
- *Remark 3.2(ii) overstated what verify_coverage.py proves.* Its subdivision gives only a weak form of step 4,
  π ≥ β with β ≥ 0.9·b_min, and its enumeration covers the resulting bounds. Reworded here and in the round-1 table.
- *Rounding.* 7.11·√c was rounded down from 7.1128·√c; the bounds are now rounded up everywhere, and the README
  agrees.
- *Field data.* p-maximality is now stated per class (2, 3 and 17; 3 and 17; 2 and 17, with 3 ∥ disc O). The mod-4
  congruence of 2.3 notes −3 ≡ 1. ε₂ of class (1,1) has its closed form 10 − 32β − β²/2 + 3β³/2.
- *The basis remark and the margins.* J changes by a sign modulo p^{e+min(λ,e)} (the referee gave a counterexample to
  the general claim, in a synthetic algebra); the conclusion stands. The λ-inequality has margin 2 at three instances.
- *Printed evidence.* 01 prints the Sturm counts, root widths, pivots and certified ranges. The bullets of §4.5 on
  phantoms, unnecessary moduli and ε₁² = 1 − 3β cite verify_exclusion.out.txt (parts 0–2, control (a6)). The
  phantoms of orbit (1,1), c = 1 are now mentioned.
- *Lemma 4.4 and the steps of 02*, as on the sister card.
- *The README's sibling rule* now says "for the same q".
- *(H) and the rank.* The CATALOGUE row, GENERALIZATIONS §3.2 item 6, the method card's Limits, the README status
  line and §0.3 now say that (H) is Cassels's theorem, which is not elementary, and that rank 2 is also proved by
  Cassels and, unconditionally, by Gebel–Pethő–Zimmer (1998).
- *Attribution.*
  - The cubic route is now credited to the algebraic proofs "whose method is known": Hemer by its review, Nagell by
    inference.
  - The four cubic equations are marked as computed here.
  - "Two-index argument" is now "two-exponent Skolem-type (unit) argument".
  - The "Not new" list gains the quartic route itself (Bennett–Ghadermarzi; KANT) and Skolem and Ljunggren for the
    two-index mechanism.
- *The project's own mordell_negative_k*, finiteness by the one-unit cubic route, is now cited. The open questions
  include a complete cubic-route proof, which needs no (H).
- *Details.*
  - §1.3: A₀ has denominators 3 only in class (1,1) and the trivial class.
  - §5: the extra solution is the projective point (1:0) of F_cl.
  - §0.1 and the Conventions of the literature check explain the slug.
  - §6.2: phantoms and the same-modulus rule occurred for 1611.a1 already.
  - §7 lists lit_checks_17.py and the descent script.
  - The parametrization is 𝐗(m,n), and N1's subring is B.
- *Literature check.*
  - §5 marks its corrections of the target files as addressed, and its line references are corrected.
  - §3 says where the Cassels scan was read: in the session's scratch folder, not stored in the repository.
- *This record.*
  - The round-1 table gained the literature agent's row.
  - "A shift of y³ − 17" is now "a shift and scaling".
  - verify_bruteforce.py's docstring now says |m| ≤ 10⁶, as its code and output do.

**Process notes.** As on the sister card, the referee workflows were resumed once after the author's context ended, and
the author applied fixes while referees were still reading. Referees fetched page scans and records into their scratch
folders; none of it is stored in the repository.

After the round-2 fixes every script of the card was re-run, and all pass.

## Round 3: the proof of (H) (2026-10-15 to 2026-10-26)
The same workflow as on the sister card ([triangular_square_pyramidal/review.md](../triangular_square_pyramidal/review.md),
round 3): for this card three verifiers, which wrote their own scripts and imported nothing from the card, and a
referee for hypothesis_H.md, the proof.md changes (§§0–1 with Lemma 1.3′, §§5–7), the README and the shared files; one
literature agent for both cards.

| Agent | Task | Verdict | What it found |
|---|---|---|---|
| units | wrote [verify_H_units.py](scripts/verify_H_units.py) (Sections 1 and 4) | 86 checks pass | The product formula twice, the norm as determinant and resultant, disc(1, θ, ω) = −867; θ₁ to width 2⁻²¹⁰; the descent identities symbolically; exactly 275 triples with T₂ ≤ C, each decided exactly in Q(θ₁), none within 1/20 of C; (4.1) a second time from a fundamental domain, without the descent lemma; a search \|a\|, \|b\|, \|c\| ≤ 40 finds exactly ±1, ±η; planted errors caught. |
| cover | wrote [verify_H_cover.py](scripts/verify_H_cover.py) (Section 2) | 72 checks pass | The 795 boxes by its own Bernstein expansion (no subdivision needed) and by a second bound, \|N\| = \|σ₁\|·Q with both factors convex; disjoint interiors (315615 pairs) and volume ½; the seven half-lattice values; the reduction to D on 11789 points; 3248 divisions; Remark 2.5 with its own diagonal lower form (303 elements, least 1115/1028 attained once); planted errors caught. |
| classes | wrote [verify_H_classes.py](scripts/verify_H_classes.py) (Section 5) | 70 checks pass | The identities in X; explicit Bézout identities for the three pairs of {ρ, σ, θ}; its own enumeration of the 32 products, of which exactly four have square norm; the identities of (5.7); its own characters at 23, 89 and 73; all 8 X with \|X\| ≤ 10⁶ in their classes, with explicit S; (5.4)–(5.6) replayed on each point; 13 planted errors caught. Also (−2 − θ)(−1 − θ)(4 − θ) = (5 + θ)². |
| referee | the whole card | accept with minor revisions | 12 findings: 1 major, 11 minor. Re-ran 04–07 (identical outputs) and wrote its own checks of the ring, the table (cell cover and a convexity bound), Remark 2.5 and the units by a fundamental domain, and [Ξ_μ : R] = 1, 1, 3, 4. Read Cassels 1950 (§11, §22, Table 2 and its notes), Cavallar–Lemmermeyer (the row "−867 N 1115/1028", p. 17) and LMFDB 3.1.867.1. |
| literature | both cards | — | [literature_check.md](literature_check.md) §7: Cioffari (1979) excluded m = 17 from the norm-Euclidean pure cubic fields; Lemmermeyer's 1989 point for m = 17 is −ξ*, with M(K, −ξ*) = 1115/1028; Cavallar–Lemmermeyer record M = 1115/1028; no published division algorithm for this ring, and no ideal-free proof of (H). Identification scripts [lit_H_lemmermeyer.py](scripts/lit_H_lemmermeyer.py), [lit_H_lmfdb.py](scripts/lit_H_lmfdb.py). |

**Result.** No finding changes Theorem H, Theorem B or any certificate. The major finding is the one of the sister
card: proof.md §6.1 cited "review.md, round 3" before it existed, and this file's header still called Theorem B PARTIAL.
This section is that record. All findings were fixed.

**Fixes made after round 3.**
1. *The record.* This section and the header; the verifiers listed in hypothesis_H.md (Checks) and in proof.md §7.
2. *K₃ in Remark 2.5.* The text attributed K₃ ≤ 3782604 to (4.3′) with its rounded C = 155.8144, which gives only
   C³ = 3782881.8; with that, T_e = 164.482 would not suffice (C³·(1115/1028)² = 4450269 > T_e³ = 4449950). The text now
   uses the exact threshold K₃ = 4(E² + E + 1)³/(E²(E + 1)²) = 3782603.77…, as 05 always did, and says that the rounded
   C would not fit; 05 prints that comparison. This is the pitfall of the 1611.a1 card's round 3.
3. *The coset C.* Only −η and −η⁻¹ map C onto itself; η maps C onto −C, where |N| takes the same values.
4. *The constants of (4.3).* The text deferred to the sister card, whose numbers are for ε. It now states h, ρ_hi,
   ρ_lo and the printed bounds E < 971.9445, ρ_hi < 1941.8911, ρ_lo > 2.114951·10⁻⁶. With these rounded bounds the old
   B = 77905209/10⁶ failed B³·ρ_lo ≥ 1 by about 10⁻⁷ (06 had used the exact enclosure of E); 06 now computes A and B
   from the printed bounds, and B is 77905221/10⁶. C = h(A) is unchanged.
5. *Section 1.* (1.6) gives A and B explicitly (04 checks the reduction), and says that Q is only semidefinite on
   real triples and positive on nonzero rational ones because σ₁Q = N ≠ 0. The discriminant is −867, and 04 now
   checks it; 867 is its absolute value. (1.1) uses (d, e, g).
6. *Section 2.* The Taylor coefficients are printed; the cell count is the cover argument; (2.1)(b) says how a general
   half-lattice point reduces to the seven; (2.2) and (2.3) are written out; (2.4) names the box of the sharp maximum
   and its five neighbours at that corner (05 prints them).
7. *Section 4.* (4.4) defines T_low and the inclusion it gives; (4.5) is written out; (4.6) says that η is
   fundamental by (4.1).
8. *Cassels.* (4.6): Table 2 (p. 270) lists, for D = 17, h = 1 and 18 − 7δ with δ = ∛17 = −θ, which is η; the entry
   has no star, which by note 2 (p. 271) marks a fundamental unit; only §11 (p. 249) speaks of an odd power. The author
   checked these pages in the referee's copy. (5.8): his λ and λη are ρσ and ηρσ in some order, since his definition
   fixes λ only up to η.
9. *Section 5.* (5.2) N(α) ≠ 0 by (1.4); (5.4) is written out with this card's symbols, where "word for word" was not
   literal; (5.5) Y ≠ 0 and N(R₀) ≠ 0; (5.6) S is written out, with the reason it lies in R, and the points that
   realize the four candidates are named; (5.7) defines the characters for this basis, including ω, and restricts the
   constancy argument to nonzero values.
10. *Section 6.* The values N(μ) are named for each μ, since the order differed from that of Theorem H.
11. *Scripts.* 05 printed 1163/1028 as the least |N| on the class of ξ*, the least on a truncated list; it now prints
    "≥ 1" for such classes (the exact least, 1115/1028, is proved in 2.5.3). 07 cites (5.7), not 5.6, and says that S is
    unique up to sign.
12. *Wording.* "Unique factorization" is replaced by the divisibility theory of Sections 2–3 (header, proof.md §6.5,
    README). proof.md's header says that every check is exact and floating point only seeds searches; the README's
    step 1 cites Lemma 1.3′.
13. *Provenance.* hypothesis_H.md §7 records the row "−867 N 1115/1028" of Cavallar–Lemmermeyer, Lemmermeyer's 1989
    point, Cioffari's exclusion and Selmer's tables. literature_check.md has a new §7. The method card
    norm.triples.multiplier names −867 next to −716; norm.triples covers this ring's basis (1, θ, ω).

**Process notes.** As on the sister card: the workflow was resumed on 2026-10-26 and all nine agents ran again. The
15-October classes verifier (75 checks) and units verifier (76 checks) of this card, both passing, were overwritten by
the 26-October ones; their files are kept in
[Archive/h_review_round3_2026-10-15/](../../../../Archive/h_review_round3_2026-10-15/). Referees downloaded page scans
into their scratch folders; none of it is stored in the repository.

After the round-3 fixes every (H) script of the card was re-run: 04–07, cover_search (the table regenerates byte for
byte) and the three verifiers. All pass.
