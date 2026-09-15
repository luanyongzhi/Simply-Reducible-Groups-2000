# Independently reproducible small-dimensional exclusions

This is a working audit note, not a claim that all order-1024 groups have been classified. The C++17 source `mf_small.cpp` and its output `mf_small_results.txt` supply an exact finite check. They do not use the supplied group census or GAP.

## Conventions and criteria

Let V=F2^m. Encode a scalar alternating form by the upper-triangular entries in lexicographic order (12,13,...,1m,23,...). An r-dimensional subspace S of Alt(V) is the dual image of the commutator map beta. Its elements are the scalar forms beta_lambda, lambda in W*.

The multiplicity-free criterion used in this audit is

    rank A + rank B + rank(A+B)
      = 2 [m - dim(rad A intersect rad B)]

for all distinct nonzero A,B in S. The proof of equivalence is to be supplied in the manuscript's character-theoretic reduction. In the enumeration, including pairs involving zero would give a tautology, so excluding them has no effect.

Choose a basis A_1,...,A_r of S. Write scalar quadratic refinements as

    q_i(v) = sum_{j<k} (A_i)_{jk} v_j v_k + sum_j t_{ij} v_j.

For each lambda in F2^r and v in rad(sum_i lambda_i A_i), ambivalence is precisely

    sum_{i,j} lambda_i v_j t_{ij}
      = sum_{j<k} (sum_i lambda_i A_i)_{jk} v_j v_k.

These are affine linear equations in mr unknown bits. Gaussian elimination returns -1 for inconsistency or the dimension of the solution affine space. No enumeration of all 2^(mr) refinements is needed.

## Normalization and exhaustive enumeration

Every nonzero alternating form on F2^5 has rank 2 or 4. Any multiplicity-free subspace of dimension at least 2 contains a rank-2 form: if two independent forms both have rank 4, the criterion gives rank(A+B)<=2; their nonzero sum must have rank 2. GL(5,2) is transitive on rank-2 forms. Thus it suffices to inspect subspaces containing the specific form e1* wedge e2*.

Subspaces containing that vector are in bijection with subspaces of the quotient 9-dimensional coordinate space. The program enumerates unique reduced row echelon representatives in that quotient, by pivot set and allowed free entries. It prunes a row extension as soon as a triple violates the rank criterion. The number of unvisited full RREF matrices below a rejected branch is accumulated exactly. Every run checks that the number rejected this way plus the number of surviving matrices equals the Gaussian binomial coefficient.

The rank criterion is symmetric in A,B,A+B. When a subspace S is extended by x, every new triple consists of one member of S and two members of x+S. Thus checking pairs (a,x+s), with a in S\{0} and s in S, is exhaustive.

## Exact results

| m | r | Enumeration domain | Total subspaces | MF survivors | Admitting an ambivalent refinement |
|---|---|---|---:|---:|---:|
| 3 | 3 | All | 1 | 1 | 0 |
| 4 | 4 | All | 651 | 0 | 0 |
| 5 | 4 | Containing e12 | 788035 | 1691 | 3 |
| 5 | 5 | Containing e12 | 3309747 | 84 | 0 |

The three (5,4) subspaces admitting refinements are the full stars a wedge V*, with a=e1*, e2*, or e1*+e2*. Each has refinement affine dimension 5, hence 32 refinements before stabilizer orbits. These are three different subspaces in the normalized search, NOT three GL-orbits; the stars form one GL-orbit. The familiar three group types arise only after the stabilizer acts on the 32 quadratic refinements.

All 84 (5,5) MF survivors have 15 rank-2 members. These 15 members, together with zero, form the unique full-star hyperplane a wedge V*. The remaining 16 forms have rank 4. The code independently checks the common wedge factor of the rank-2 forms and the inconsistency of every refinement system.

Consequences:

* There is no class-2 SR group with dim(G/G')=4 and dim G'=4. In particular the draft's positive-looking claim of 546840 different four-dimensional subspaces in Alt(F2^4) cannot be correct: there are only [6 choose 4]_2=651.
* The order-1024 stem split (m,r)=(4,6) is excluded directly by the rank criterion: the full Alt(F2^4) contains A=e12+e34 and B=e13+e24, both rank 4, with rank(A+B)=2 and common radical zero. The criterion would require 10=8.
* The order-1024 stem split (m,r)=(5,5) has no SR group. This follows from the complete normalized enumeration, independently of the order-512 census.

## A short hand proof for the final ambivalence obstruction

Here is a general useful lemma, which explains the 84 discarded cases.

**Full-star obstruction.** Let V=F2 e0 direct-sum U with dim U>=2. If S contains the full star e0* wedge U*, no ambivalent quadratic system can have S strictly larger than that star.

Proof. The star components can be combined into a U-valued quadratic map. Their polar map is (t,u),(s,v) |-> t v+s u. A quadratic refinement has the form

    q_star(t,u)=t u+L(u)+t w.

Ambivalence at (0,u) says a(L(u))=0 whenever a(u)=0, for every functional a in U*. Hence L(u) is in span(u) for every u; linearity and dim U>=2 force L=delta I, delta in F2. Therefore

    q_star(t,u)=t u+delta u+t w.

An additional scalar form can be modified by a star form so that its polar is an alternating form omega on U alone. Since it is outside the star, omega is nonzero. Its refinement is

    q_omega(t,u)=Q_omega(u)+c t+ell(u).

For each u in U, put a=omega(u,-). Then (1,u) belongs to the radical of omega+e0* wedge a. Ambivalence for that scalar form gives

    Q_omega(u)+c+ell(u)+omega(u,w)=0

for all u, because omega(u,u)=0. Taking u=0 gives c=0. The displayed identity then makes Q_omega linear, contradicting its nonzero polar omega. QED.

## Reproduction

    g++ -O3 -std=c++17 mf_small.cpp -o mf_small
    ./mf_small > mf_small_results.txt

The local execution completed in approximately 0.02 seconds (machine dependent). This is possible because the criterion rejects whole RREF branches early; the source does not literally visit millions of completed matrices. Timings are not part of the mathematical certificate.

The surviving MF-space counts are counts of distinct coordinate subspaces in the stated normalization. They must not be reported as numbers of group isomorphism types. Nothing in this check resolves the (6,4), (7,3), or (8,2) order-1024 strata.
