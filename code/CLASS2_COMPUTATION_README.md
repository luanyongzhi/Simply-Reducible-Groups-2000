Order 1024, nilpotency class 2: reproducible classification

The new stem counts are 18 for `(dim G/G', dim G')=(8,2)`, 80 for `(7,3)`,
and 28 for `(6,4)`. The `(5,5)` and `(4,6)` strata are empty by the separate
finite form checks. The 126 new stem groups, together with the 95 groups
obtained by padding smaller stem groups with central elementary abelian
factors, give 221 class-two groups. The completeness argument depends on the
previous order-512 SR census, whose positive class-two entries were checked
again in this work.

The group construction and all final positive-record checks use GAP 4.15.1.
The installed SmallGrp version used for source data and regression tests is
1.5.4. Finite graph canonization uses `pynauty==2.8.8.1` (nauty 2.8.8);
small form and parent-extension traversal use C++17. These are implementation
choices: the classification is a finite linear-algebra and group-orbit
calculation, not a call to the unavailable `SmallGroup(1024,i)` catalogue.

Reconstruct a result in GAP, from the directory containing the files:

```gap
Read("sr1024_special_representatives.g");
Length(SR2_special_representatives);  # 126
R := SR2_special_representatives[1];
G := PcGroupCode(R.code,1024);
Size(G);                           # 1024
NilpotencyClassOfGroup(G);          # 2
```

`sr1024_special_representatives.csv` contains the same 126 records, their
dimensions, family labels, centre/derived orders, involution and conjugacy
class counts, two Wigner moments, and reconstructible pc codes. Their
pairwise nonisomorphism follows from the complete form-orbit classification;
distinct pc codes alone are not used as an isomorphism certificate.

The core GAP library is `sr_class2.g`. Given a list `B` of alternating matrices,
`SR2_PolarData(B)` computes scalar ranks and radicals,
`SR2_MFPolar(P)` tests the exact multiplicity-free rank/intersection criterion,
and `SR2_RealityAffine(P)` solves all reality conditions as one affine linear
system. `SR2_ForEachRefinement(B,callback)` streams all SR refinements of that
fixed polar tensor. `SR2_Group(B,Q)` constructs the pc group and
`SR2_Extract(G)` extracts its data. The core does not independently enumerate
polar tensors or remove isomorphic refinements.

First run the GAP checks:

```text
gap -q -b test_sr_class2.g
gap -q -b audit_uploaded_groups.g
```

The first test exhausts groups of orders 8, 16, 32 and 64, checks the 79
eligible class-two form encodings against exact ordinary-character tensor
products, and verifies all 79 extraction/reconstruction identities using
`IdGroup`. The second audits every uploaded class-two row at orders through
512; it gives class-two counts `2,2,7,10,20,42,95` and stem counts
`2,0,5,3,10,22,53`.

For the `(8,2)` pencil stratum:

```text
gap -q -b classify_last_pencil.g
gap -q -b build_m8r2_representatives.g
```

The first command verifies 15 explicit tensor automorphisms, partitions all
256 refinements of the singular pencil into five orbits, and separates those
orbits by intrinsic invariants. The second adds 13 representatives from the
three regular pencil types and checks all 18 with actual group-element
Wigner moments. `pencil_m8r2_certificate.md` explains the complete refinement
classification and its relation to the pencil normal-form theorem.

For the two remaining stem dimensions, generate the covering polar lists
using `small_forms_notes.md` and `extensions_notes.md` in the accompanying
small-form computation directory. The relevant inputs to the graph stage are
`candidate_m6_nets.txt` (356 coordinate nets) and
`subgroup_m7_net_representatives.txt` (44 coordinate nets covering the parent
extension orbits). With those input paths adjusted to their location:

```text
python -m pip install pynauty==2.8.8.1
python canonical_class2.py candidate_m6_nets.txt class2_m6r4_orbits.json
python canonical_class2.py subgroup_m7_net_representatives.txt class2_m7r3_orbits.json
```

The results are five polar orbits with `6,7,6,6,3` quadratic orbits and eight
polar orbits with `7,6,10,6,18,12,12,9` quadratic orbits. Canonicalization uses
the entire coloured graph, not scalar rank/Gauss fingerprints. Exact
canonical bytes are compared internally; SHA-256 digests are recorded only
to identify the archived certificates. Each graph automorphism is checked
to induce linear maps on all V and W points and to preserve the full polar
tensor. Its action on quadratic refinements is computed as
`q'(v)=S^{-1}q(Tv)`. Full quadratic-pair graph canonization also separates
the reported stabilizer orbits.

`m7_fullparent_crosscheck.json` records an independent second cover: the 24
coordinate nets produced using complete parent stabilizers give exactly the
same eight polar classes as the 44-net subgroup cover. The extension programs
also check that orbit sizes sum to the exhaustive accepted extension counts.
`twin_tripwire_certificate.json` checks three known, distinct order-512
groups which share polar data: their quadratic-pair canonical certificates
are distinct. This specifically tests that the graph stage retains q.

To regenerate and verify the combined 126 stem representatives:

```text
python export_special_to_gap.py class2_m6r4_orbits.json class2_m7r3_orbits.json --output verify_special_representatives.g
gap -q -b verify_special_representatives.g
```

The log must end with `COMPLETE 126`. Every record is checked for its order,
nilpotency class, centre and derived sizes, the exact form criterion, both
independent Wigner moment identities computed from elements and conjugacy
classes, and a saved-pc-code reconstruction. The program aborts on any failed
check. The representative file was additionally read back as GAP source.

The graph reduction is exact for the following elementary reason. Over F2,
each projective line has points `{x,y,x+y}`, and `beta(x,y)` is the same for
each pair of different points on that line. The coloured graph records all
such lines in V, their beta values in W, and all lines in W; W zero has its
own colour. An isomorphism of these graphs therefore gives additive, hence
linear, point bijections on V and W preserving beta, and conversely every
tensor equivalence gives such an isomorphism. Adding the edge from each
nonzero V point to its q-value encodes the entire quadratic refinement.

The class-at-least-three completion is now documented by the original
539-job checkpoint and the census audit included in this archive. All expected
jobs are present and all 581 positive codes match the CSV. Together with the
221 class-two groups and the elementary abelian group, this gives 803.
See the classification manuscript and ../validation/revision_census_audit.md.

Primary software references:

* [GAP/SmallGrp manual](https://gap-packages.github.io/smallgrp/doc/chap1.html)
* [Pynauty source and API](https://github.com/pdobsan/pynauty)
* [Nauty documentation](https://users.cecs.anu.edu.au/~bdm/nauty/)
