# Portable computation: scalar-polar candidate covers

These programs supply complete candidate covers. Use the separate canonical-graph and quadratic-orbit classifier to obtain group-isomorphism counts.

Requirements: GAP 4.15.1 with SmallGrp, and a GCC/Clang C++17 compiler supporting unsigned 128-bit integers. No optional GAP packages are needed beyond SmallGrp for extraction. The full-parent orbit path also uses stabilizer generators produced by the accompanying graph-classification scripts. The elementary-subgroup path generates its own certified subgroup generators and does not require graph software.

## Files and order of use

1. `mf_small.cpp`: independent small-dimensional exclusion and validation. See `small_forms_notes.md` for proofs and exact counts.
2. `extract_parents.g`: obtains the 95 class-2 parents from the supplied `SR_groups_results512.csv` IDs, using `sr_class2.g`. The IDs are embedded in the script after independent class-2 filtering. It produces `parents512_forms.txt`; a copy of this exact input is supplied.
3. `extend_parents.cpp`: direct enumeration of polar extensions, without orbit pruning.
4. `parent_subgroup_orbits.cpp`: orbit pruning, either with self-generated verified elementary subgroups or supplied full parent-polar stabilizers.
5. Pass the resulting covers to the separate canonical-graph classifier. Its standard inputs are `candidate_m6_nets.txt` and `subgroup_m7_net_representatives.txt`; `fullparent_m7_net_representatives.txt` provides the independently reduced second m=7 cover.

All commands below run from the directory containing these files. Put `sr_class2.g` in that directory for portable GAP extraction.

    g++ -O3 -std=c++17 mf_small.cpp -o mf_small
    ./mf_small > mf_small_results.txt

    gap -q -b -A extract_parents.g > extract_parents.log

    g++ -O3 -std=c++17 extend_parents.cpp -o extend_parents
    ./extend_parents 6 < parents512_forms.txt > extend_m6_results.txt
    ./extend_parents 7 < parents512_forms.txt > extend_m7_results.txt

    g++ -O3 -std=c++17 parent_subgroup_orbits.cpp -o parent_subgroup_orbits
    ./parent_subgroup_orbits 6 < parents512_forms.txt > subgroup_m6_results.txt
    ./parent_subgroup_orbits 7 < parents512_forms.txt > subgroup_m7_results.txt

When the full stabilizer files have been generated:

    ./parent_subgroup_orbits 6 parent_polar_auts_m6.txt < parents512_forms.txt > fullparent_m6_results.txt
    ./parent_subgroup_orbits 7 parent_polar_auts_m7.txt < parents512_forms.txt > fullparent_m7_results.txt

## Expected checks

* Extraction: GAP version 4.15.1; 95 parents.
* Small exclusion: all 651 (4,4) spaces fail MF; among 3309747 normalized (5,5) spaces, 84 pass MF and none admits an ambivalent refinement.
* m=6 direct enumeration: 36855 quotient vectors; 358 surviving parent-extension hits; 356 distinct coordinate nets.
* m=7 direct enumeration: 4194296 quotient vectors; 290490 surviving parent-extension hits; 290488 distinct coordinate nets.
* Elementary-subgroup covers: m=6 gives 11 parentwise orbits, 9 coordinate representatives, orbit-size sum 358; m=7 gives 45 parentwise orbits, 44 coordinate representatives, sum 290490.
* Full-parent covers: m=6 gives the same 11/9/sum358; m=7 gives 25 parentwise orbits, 24 coordinate representatives, sum290490.

The 356 and 290488 coordinate-net counts are not group-isomorphism counts. Neither are the parentwise orbit counts. Exact global equivalence and quadratic stabilizer orbits are indispensable.

For mathematical completeness, input formats, actions, and interpretation, see `extensions_notes.md`. The large direct candidate files can be regenerated rapidly and need not be permanently stored if disk space is limited. Keep the source, the 95-parent input, small representative covers, stabilizer generators, and computation logs.
