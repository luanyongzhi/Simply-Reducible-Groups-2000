The `(m,r)=(8,2)` quadratic refinements: 18 representatives

The enumeration below is conditional only on the separate alternating-pencil
normal-form argument giving the four polar types listed here. It completely
settles the quadratic-refinement and isomorphism step for each of those four
types. No generic group-isomorphism search is used in the proof.

Let `W=F_2^2`, with nonzero vectors `w1,w2,w3=w1+w2`. A regular block of
colour `wa` is a two-dimensional symplectic plane whose commutator takes
values in the line `<wa>`. The three regular polar types have respectively
`(3,1,0)`, `(2,2,0)`, `(2,1,1)` planes of the three colours. Their rank profiles
are `[2,6,8]`, `[4,4,8]`, `[4,6,6]`.

For a regular polar type, the colour subspaces are exactly the radicals of the
three nonzero scalar forms. They are therefore preserved, up to the
permutations induced by `GL(W)`, by every polar equivalence. The reality
condition forces the quadratic map on a colour subspace to take values in its
colour line. It is a nonsingular scalar quadratic form there, determined up to
isometry by its Arf invariant. Conversely the independent isometries on colour
subspaces, and all permutations of subspaces of equal dimension, are actual
polar equivalences. Thus the complete counts of refinements up to equivalence
are:

| Polar type | Complete Arf data | Count |
|---|---|---:|
| R310 | Ordered pair of signs on dimensions 6 and 2 | 4 |
| R220 | Unordered pair of signs on two spaces of dimension 4 | 3 |
| R211 | One sign on dimension 4 and an unordered pair on dimensions 2 and 2 | 6 |

The final polar type, denoted S1S1R1, has basis
`(a1,b1,c1,a2,b2,c2,p,q)` and nonzero commutators
`beta(ai,bi)=w1`, `beta(ai,ci)=w2` for `i=1,2`, and `beta(p,q)=w1`.
Its three scalar forms are `B1,B2,B12=B1+B2`; their ranks are `6,4,6`.
Their radicals are

```
R2  = <b1,b2,p,q>,
R1  = <c1,c2>,
R12 = <b1+c1,b2+c2>.
```

Because `B2` is the unique rank-four form, the subspace
`X=R2 intersect (R1+R12)=<b1,b2>` is intrinsic under the full pencil-equivalence
group, including interchange of `B1` and `B12`. Consequently whether
`q|X=0` is a group-isomorphism invariant. The number of involutions is another
such invariant.

`classify_last_pencil.g` first solves all reality constraints and obtains
exactly 256 quadratic refinements. It defines 15 explicit invertible linear
changes `(T,C)`, checks the full tensor identities

`sum_b C[b,a] T B_b T^transpose = B_a` for `a=1,2`,

and then performs an exhaustive breadth-first orbit traversal of the 256
refinements under `Q'(ei)=q(T ei) C` (the script uses row vectors). This yields:

| Orbit | Size | `q|X=0` | Number of involutions |
|---:|---:|---|---:|
| 1 | 48 | true | 383 |
| 2 | 48 | false | 383 |
| 3 | 48 | false | 255 |
| 4 | 96 | false | 191 |
| 5 | 16 | true | 255 |

The orbit sizes sum to 256. Each orbit is an orbit of a verified subgroup of
the full tensor stabilizer, so all its refinements define isomorphic groups.
The five invariant pairs are different, so no two subgroup orbits can merge
under the full stabilizer. This proves exactly five full orbits without needing
to prove that the 15 generators generate the entire stabilizer.

The last polar type and R211 have the same rank profile. They cannot be
equivalent: the intrinsic subspace `R2 intersect (R1+R12)` has dimension 2 for
S1S1R1 and dimension 0 for R211. The other two regular types have different
rank profiles. Therefore these four lists do not overlap and contain
`4+3+6+5=18` pairwise nonisomorphic groups.

Validation performed with GAP 4.15.1 and SmallGrp 1.5.4:

* `test_sr_class2.g`: 79 comparisons with the definition using exact ordinary
  characters, 79 extraction/reconstruction `IdGroup` comparisons for orders
  8,16,32,64, and all four refinements of the dimension-two symplectic form.
* `audit_uploaded_groups.g`: every class-two row of the uploaded 2-power
  lists through order 512 passes the new exact SR criterion. Class-two counts
  are `2,2,7,10,20,42,95`; stem counts are `2,0,5,3,10,22,53`.
* `build_m8r2_representatives.g`: all 18 constructed groups have order 1024,
  nilpotency class 2, centre and derived subgroup of order 4. All satisfy the
  form criterion and both Wigner moment identities, computed independently
  from actual group elements and conjugacy classes.
* The generated representative file was read back, and all 18 saved pc codes
  were reconstructed using `PcGroupCode(code,1024)` and passed the SR and centre
  checks again. The CSV contains 18 rows and 18 distinct pc codes. Code
  distinctness by itself is not being used as an isomorphism proof.

Reproduction, from this directory:

```
gap -q -b test_sr_class2.g
gap -q -b classify_last_pencil.g
gap -q -b build_m8r2_representatives.g
```

The reusable core is `sr_class2.g`. The last two commands generate
`last_pencil_orbits.csv`, `last_pencil_representatives.g`,
`sr1024_class2_m8r2.csv`, and `sr1024_class2_m8r2_representatives.g`.
The four-polar-type classification should be proved in the paper before stating
that 18 exhausts the entire `(8,2)` stratum.
