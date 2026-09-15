# The order-1024 class-two catalogue

`sr1024_class2.csv` contains **221 isomorphism-class representatives**: 95 groups obtained as `K x C2` from the complete class-two order-512 census, and 126 special representatives constructed from the classified quadratic maps. The special counts for `(m,r)=(6,4),(7,3),(8,2)` are respectively 28, 80, and 18.

The complete class-two split is:

| `(m,r)` | Special | Non-special | Total |
|---|---:|---:|---:|
| `(6,4)` | 28 | 3 | 31 |
| `(7,3)` | 80 | 51 | 131 |
| `(8,2)` | 18 | 33 | 51 |
| `(9,1)` | 0 | 8 | 8 |
| Total | 126 | 95 | 221 |

The exclusion of `(9,1)` in the special classification does not exclude its eight non-special groups.

The export is a new GAP reconstruction and validation of those representatives. It does not merely filter an existing CSV, and it does not repeat the expensive orbit classification. Its completeness and absence of duplicate isomorphism types use the classification proof and certificates in the paper; presentation codes and numerical invariants are not used to establish isomorphism classes.

## Minimal files

Place these four input/code files together in one directory:

1. `export_sr1024_class2.g` — the entry point;
2. `sr_class2.g` — construction, extraction, and exact class-two SR criterion;
3. `sr1024_special_representatives.g` — all 126 classified `(B,Q)` representatives with saved presentation codes;
4. `SR_groups_results512.csv` — the complete 317-row order-512 SR census.

Only GAP and its SmallGrp package are needed for this export. ANUPQ, Python, C++, and graph-canonicalization packages are not required. The original export was run with GAP 4.15.1 and SmallGrp 1.5.4; its archived execution log is supplied. The 14 September revision only changes the canonical input filename and file lookup paths; GAP was unavailable for a new execution of that revision. The script also recognizes the canonical census in `../data/`, so it runs directly from this archive's `code/` directory. Historical fallback locations remain supported.

Run from that directory:

```sh
gap -q -b export_sr1024_class2.g < /dev/null > export_sr1024_class2.log 2>&1
```

The script writes `sr1024_class2.csv` only after every representative passes. The log must end with the `COMPLETE 221` message. A shell exit status alone is not sufficient evidence that an interactive GAP run reached its endpoint.

Alternatively, in a GAP session in the same directory:

```gap
Read("export_sr1024_class2.g");
```

To choose a different location for the order-512 census, set the optional input variable before reading the script:

```gap
SR1024Class2_Input512 := "SR_groups_results512.csv";;
Read("export_sr1024_class2.g");
```

## Schema and reconstruction

The first eleven columns exactly match the order and names in `sr1024_class3plus.csv`:

`No,ParentOrder,ParentId,Step,Class,Exponent,CenterSize,DerivedLength,NrConjClasses,HasC2DirectFactor,Code`.

Here `ParentOrder`, `ParentId`, and `Step` refer to the lower exponent-2 central-series parent. Since every group in this file has class two and `Phi(G)=G'`, this parent is `G/G' = C2^m`; its order is `2^m`, its SmallGrp identifier is `ParentId`, and `Step=r`. These fields describe the canonical parent even though the export constructs groups from forms or direct products, rather than running ANUPQ.

The appended fields are:

| Field | Meaning |
|---|---|
| `Order` | Always 1024 |
| `DerivedSize` | Order of the derived subgroup |
| `Involutions` | Number of nonidentity elements of order two |
| `m`, `r` | `m=log2(|G/G'|)`, `r=log2(|G'|)`; thus `m+r=10` |
| `IsSpecial` | Whether `Z(G)=G'=Phi(G)` |
| `Source` | `class2_padding` or `class2_special` |
| `SourceID` | For padding: the original order-512 SmallGrp identifier; for special groups: the record number in the 126-representative file |
| `Family` | Construction label identifying the polar/quadratic representative or the direct-product construction |
| `WignerS2`, `WignerS3` | Exact square-root moments recomputed from the group |

`No` is a row number, not a SmallGrp identifier at order 1024. Each row represents the group reconstructed in GAP by:

```gap
G := PcGroupCode(code,1024);;
```

where `code` is the entire exact integer in the `Code` column. Conversely, `CodePcGroup(G)` encodes a pc presentation. **Do not convert `Code` to a floating-point number**, and do not interpret distinct codes as a proof of nonisomorphism.

For a special representative the full defining form data are also available:

```gap
Read("sr_class2.g");;
Read("sr1024_special_representatives.g");;
R := SR2_special_representatives[1];;
G := SR2_Group(R.B,R.Q);;
```

For a non-special representative with `SourceID=i`, its alternative construction is:

```gap
G := DirectProduct(SmallGroup(512,i),CyclicGroup(2));;
```

## Checks performed by the export

For each group the exporter verifies order 1024, nilpotency class two, the exact form SR criterion, the special/direct-factor flag, the independent Wigner equalities, and reconstruction from the saved pc presentation. The Wigner moments are computed from actual group elements, independently of the form calculation. The exporter requires 221 rows, 95 padding groups, and the special split 28/80/18 before publishing the CSV.

A separate comparison against the previous complete 803-row catalogue is recorded in `sr1024_class2_comparison.json`: the 221 entries agree on presentation code and every shared group invariant, after matching `(Source,SourceID)`. This is a consistency check between catalogue exports, not a substitute for the orbit-classification proof.

The archived original clean run reached `COMPLETE 221` in 10.665 seconds with process exit status zero. Its log contains no syntax errors or warnings.
