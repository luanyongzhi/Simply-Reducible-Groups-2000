# Classification manuscript: independent review for the September revision

Reviewed on 2026-09-08. This review concerns the current `output/SRgroup2000.tex`, the prior independent audits in `audit_math/`, and the uploaded Chapter XXVI of Berkovich–Kazarin–Zhmud'. It does not rerun the expensive classification. The new parent-generation logs are being audited separately.

## Overall assessment

The current class-two classification has a sound mathematical chain: central splitting; the quadratic-map dictionary; exact rank/radical and ambivalence criteria; exhaustive quotient-polar coverage; faithful graph equivalence and full stabilizer actions; and the two small exclusions plus the four-pencil refinement calculation. No new mathematical defect was found in these arguments on rereading the final manuscript and the prior independent code audits.

The numerical conclusion is a computer-assisted classification. The proof of coverage is mathematical; the finite orbit totals are supplied by exact enumeration and its execution certificates. It should not be described as an unaudited invariant-matching argument, but neither should numerical counts be presented as though deduced entirely by hand.

The major revision needed is expository: state the classification outcome and the remaining certification assumptions clearly, then arrange every retained theorem and algorithm around the census. The old manuscript's qualifications about the 581 higher-class groups should be updated only to the extent supported by the newly supplied `scriptG` source and complete logs.

## Exact book citations verified from the uploaded PDF

The supplied chapter is Chapter XXVI, *Simply reducible groups*, beginning at printed p.419.

1. **Direct products:** Exercise 1.8 and its solution, printed pp.426–427 (PDF pages 8–9). Recommended citation:
   `\cite[Chapter~XXVI, Exercise~1.8]{BKZ19}`.
   This is explicitly an exercise with a proof in the book, rather than a numbered theorem. The manuscript may retain a short proof for the combined quotient/direct/central-product proposition, while crediting this known direct-product result.

2. **Odd-order impossibility:** Exercise 1.1 and its solution, printed p.421 (PDF page 3). The solution states the stronger fact that the identity is the only real element of a finite odd-order group. Recommended citation:
   `\cite[Chapter~XXVI, Exercise~1.1]{BKZ19}`.
   This supports an early introductory sentence: a nontrivial SR group has even order because every element is conjugate to its inverse, whereas an odd-order group has no nonidentity element with that property. The trivial group is SR and must be explicitly exempted.
   Do not call Exercise 1.1 an “odd-order theorem” with an invented theorem number.

3. **Elementary center and abelianization:** Exercise 1.7 and its solution, printed p.426 (PDF page 8). It proves that both `Z(G)` and `G/G'` have exponent at most two and also explicitly uses quotient closure. Recommended citation:
   `\cite[Chapter~XXVI, Exercise~1.7]{BKZ19}`.

The existing bibliography entry `BKZ19` already identifies the second edition of *Characters of finite groups*, Vol.2, De Gruyter, 2019, series volume64. No additional bibliographic entry is needed.

## Mathematical checks and assumptions that must remain explicit

### Central splitting and parameter range

The special factor in `G = S x C_2^j` is unique up to isomorphism by Krull–Remak–Schmidt cancellation, and “special” is not the same as “directly indecomposable.” In particular `D_8 x D_8` is special. Do not remove this distinction when simplifying the prose.

For the **special** target of order1024, define `V=S/S'`, `W=S'=Z(S)`, `m=dim V`, `r=dim W` before listing parameter pairs. Since `Phi(S)=S'`, `m` is the minimal number of generators. The conditions `m+r=10`, `r>=1`, and surjectivity `r<=binom(m,2)` imply `4<=m<=9`. The six preliminary pairs are `(4,6),(5,5),(6,4),(7,3),(8,2),(9,1)`. The scalar pair `(9,1)` is impossible because a nondegenerate alternating scalar form has even rank. This explanation is needed before the five computational branches, and only concerns special factors; the 95 non-special groups are recovered by padding the 95 class-two groups of order512.

### Exact SR criterion

The character-support argument and the tensor coefficient

`2^((rho_lambda+rho_mu+rho_(lambda+mu))/2 + dim(R_lambda intersection R_mu) - m)`

are correct even before ambivalence is imposed. The scalar phases on character supports may have order four; the argument uses their linear-character sum and does not assume they are signs. Multiplicity-freeness depends only on the polar tensor. Reality/ambivalence is the separate affine condition `lambda(q(v))=0` for `v in R_lambda`.

Checking one pair per two-dimensional subspace of `W*` is sufficient. The cases `lambda=0`, `mu=0`, and `lambda=mu` are automatic. The final draft should use the unambiguous quantifier “distinct nonzero lambda, mu” instead of relying on a chained inequality.

Retain the simple formula

`# involutions = 2^r * |q^(-1)(0)| - 1`.

It follows because every fiber over `v in V` contains `2^r` elements with square `q(v)`. This formula is needed by the final pencil orbit-separation certificate. It can stand directly after the quadratic dictionary; the Fourier/Gauss-sum/Clifford discussion is unnecessary for the classification proof.

### Small exclusions

For `(4,6)`, the two displayed rank-four forms have a rank-two sum, giving the contradiction `4+4+2 != 8`; this is a direct proof.

For `(5,5)`, the normalization to contain `e12` is exhaustive because an MF space of dimension at least two cannot have every nonzero form of rank four. The raw quotient-Grassmannian count is `[9 choose4]_2=3309747`. The exact enumeration yields84 MF coordinate spaces, all with a full-star hyperplane and inconsistent affine reality equations. The number84 is neither a count of groups nor a count of GL-orbits. The earlier independent Python construction already verifies the structure and reality obstruction for all84; no expensive rerun is needed.

The full-star parameter calculation and the “no proper ambivalent extension” lemma are therefore part of the classification proof and must survive the scope reduction. Their discussion as an independent infinite family can be shortened.

### Pencil branch:18

For a stem pencil, the three radical-equation row spaces have zero total dependence, so every stem pencil has an ambivalent refinement. In the MF case the affine solution set has size `2^m`. This is an existence theorem, not an empirical assertion.

The canonical alternating-pencil blocks have nonnegative rank defect. Defect zero permits only two-dimensional scalar-colored symplectic blocks and three-dimensional singular blocks. At `m=8` this gives exactly four polar types, including two types with the same rank multiset. Their refinement counts are4,3,6,5, for18.

For the last type, the15 verified tensor automorphisms yield subgroup orbits of sizes48,48,48,96,16. They give at most5 full orbits. The pairs `(q|_X=0, #involutions)` are distinct and intrinsic, giving at least5. Thus there is no need to assert that the15 matrices generate the full stabilizer. Keep this upper/lower-bound explanation in the main text; put the finite generator/orbit procedure in its referenced pseudocode appendix.

The ordinary source `Scharlau76` is a useful direct citation for the alternating-pencil canonical decomposition. Falcone–Vaccaro explicitly attributes the complete alternating case to Scharlau on p.56 and discusses the indecomposable modules and alternating doubling on pp.57–59. The already-cited DD18 reference supports the singular decomposition. An updated wording can cite all three without attributing the full regular classification to DD18 alone.

Primary source checked online:
- Falcone–Vaccaro (2004), article record: https://dml.cz/handle/10338.dmlcz/132949
- Full paper: https://dml.cz/bitstream/handle/10338.dmlcz/132949/ActaOlom_43-2004-1_5.pdf

### Quotient coverage and exact graph equivalence

The quotient argument correctly uses **all** relevant class-two SR parents of order512, including non-special parents. The quotient is `P=S/L`, with `L` a line in `S'`; it preserves `P/P'=S/S'`. Replacing that quotient with `P/Z(P)` or restricting parents to stem groups would lose valid targets.

Enumerating polar extensions independently of the particular parent's quadratic refinement is safe because all admissible quadratic refinements of each target polar space are subsequently enumerated. A parent stabilizer subgroup may be proper: it only leaves extra representatives. Its generators must be verified invertible and to preserve the parent polar space, and global target equivalence is still essential.

The five-color graph is exact. Line incidence over `F_2` forces linearity of V- and W-point permutations, while the tensor edge records beta; q-value edges record the full quadratic pair. Exact canonical graph bytes are used for equivalence. A SHA hash is only an identifier; it is not a substitute for exact comparison. Distinct PC codes likewise do not prove nonisomorphism.

For refinement orbits, the full **target** stabilizer is needed unless an independent invariant separation argument proves a smaller subgroup sufficient. The graph automorphism group supplies that full stabilizer for `(6,4)` and `(7,3)`. The affine action `q -> D^(-1) q A` is correct and retains its quadratic correction term.

The prior independent audits inspected the actual quotient-vector code, affine solver, graph encoding, automorphism validation, and orbit traversal. They checked that all successful-parent orbit masses match exhaustive totals, including both m7 parent-cover paths. The recorded counts are5 polar types with28 q-orbits and8 polar types with80 q-orbits.

## Scope changes recommended for the classification paper

Remove derived/Fitting-length discussion, diagonal Clifford families, collision/epicenter/norm-circle constructions, and the modulo-two character-ring section. They are not needed for the census proof as currently structured. Preserve them in the earlier source or a separate future-paper notes file rather than implying those results were false.

Retain only the following dependencies from the class-two structural material:
- central splitting and exact order512-to1024 non-special correspondence;
- special-group quadratic dictionary and involution formula;
- exact MF and ambivalence test;
- full-star formula and extension obstruction for the `(5,5)` certificate;
- pencil canonical types and18 refinements;
- quotient coverage and exact graph equivalence;
- finite construction of the221 class-two representatives and their reconstruction certificate.

The abstract should state the census objective and result concisely. Technical pairs and separate orbit counts belong at the start of the order1024/class-two section. Until the new higher-class evidence establishes its full certificate, an unqualified global classification assertion would still overstate the verified conclusion. The final qualification should reflect the new log audit, not the old screenshots.

GAP version numbers can be removed from narrative prose as requested; retain actual executable/package versions in bibliography or archived run metadata. Replace bare mentions of important script filenames in proof paragraphs with references to their appendix algorithms, while the reproducibility appendix maps algorithms to exact files.

The long count table should remain, with `f(n)` defined as the number of **SR-group isomorphism classes** of order n. It is different from `NumberSmallGroups(n)`, which counts all group isomorphism classes. A finite table by itself does not prove a universal formula; separate any proven recurrence from observed patterns or conjectures.
