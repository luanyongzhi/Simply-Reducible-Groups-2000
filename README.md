# Reproducibility files for the classification of finite simply reducible groups

This repository accompanies the manuscript revision of 14 September 2026,
*Classification of simply reducible groups*. It contains the completed census (7889 nontrivial isomorphism types, including
803 at order 1024), the generating programs, exact representative data, and
computation certificates. The trivial group has f(1)=1; other odd orders have
zero entries. This distribution was assembled on 15 September 2026 from
the checked manuscript reproduction package; it does not represent a new
run of the large catalogue scans or descendant generation.

## Start here

- Browse the published catalogues in [`data/`](data/). For all order-1024
  representatives, use [`SR_groups_results1024.csv`](data/SR_groups_results1024.csv).
- [`PAPER_FILE_INDEX.md`](PAPER_FILE_INDEX.md) maps every explicitly named
  code, CSV, and execution-evidence file in the paper to its repository path.
- Use the quick Python checks below to inspect the saved evidence and arithmetic
  without starting group generation.
- For the full class-two workflow, see
  [`CLASS2_COMPUTATION_README.md`](code/CLASS2_COMPUTATION_README.md) and
  [`COMPUTE_README.md`](code/COMPUTE_README.md).
- Keep the `code/`, `data/`, and `validation/` directories together. Run each
  command from the working directory specified below. The `Code` CSV column
  must retain arbitrary-precision integers.

Requirements depend on the task: Python 3's standard library for the quick
checks; GAP and SmallGrp for reconstruction and exports; ANUPQ for descendant
generation; a C++17 compiler with unsigned 128-bit integer support for the
extension programs; and `pynauty==2.8.8.1` for the graph classification.
The workflow documents record the versions used for the archived computation.
`SHA256SUMS.txt` lists checksums of the delivered files.

## Repository contents

| Files | Role |
|---|---|
| `data/SR_groups_results2_2000.csv`, `data/SR_groups_results512.csv`, `data/SR_groups_results1536.csv` | the three disjoint SmallGroups catalogues, with 6210, 317 and 559 rows. |
| `data/sr1024_class3plus.csv`, `data/sr1024_class2.csv`, `data/SR_groups_results1024.csv`, `data/f_of_n.csv` | the high-class, class-two, combined order-1024 catalogues and full counting table. The first two lists are parts of the third, not additional contributions. |
| `code/sr_search*.g` | the original catalogue-search algorithms cited in the paper. The upper-range alternative overlaps the other scans; do not add its output twice. |
| `code/scriptG.g`, `code/scriptG_reproduce.g` | the descendant-generation source and reproduction variant. Their mathematical job domains agree. Use the reproduction variant for a new run; see below. |
| `code/scriptG_results.FINAL.bak` | the completed machine-readable checkpoint. Its DONE records establish which jobs finished; its SRHIT records carry the accepted presentation codes. |
| `code/makeCSV.g` | reconstructs the accepted groups from the checkpoint, repeats the exact SR test, and exports the high-class CSV. |
| `code/scriptG_console.FINAL.bak` | Supporting execution evidence. The human-readable console transcript is not an input to `makeCSV.g`, but records the run cited in the paper and is hashed by the audit. |
| `code/audit_census_checkpoint.py` and `validation/census_checkpoint_audit.json` | verifies the expected parent jobs, completion records, accepted records and CSV correspondence without generating descendants. |
| `code/read_sr1024_class3plus.g` | a short exact CSV reader and group reconstruction helper. |
| `code/README_sr1024_class2.md`, the form-classification sources and their input/certificate files | construction and verification of all 221 class-two representatives. The ordered full computation is documented in `code/CLASS2_COMPUTATION_README.md` and `code/COMPUTE_README.md`. |
| `code/counting_structures.py`, `code/abelian_hall_counts.py`, `code/verify_progressions.py` | exact arithmetic checks and derived counting tables for the arithmetic section. |
| Other `.log`, `.json`, `.txt`, `.md` files in `code/` and `validation/` | Retain as computation certificates or historical execution evidence. Consult the computation READMEs before removing any representative input file. Typesetting-only logs are optional for public distribution. |

`data/f_of_n.csv` has no header: its two columns are `n,f(n)` and it
contains 1000 rows for the even orders 2 through 2000. The representative
catalogues have header rows. `data/f_table_body.tex` supplies the same count
table in the LaTeX format used by the paper.

`SRHIT` is **not a separate file**. It is the first field of an accepted-result
line in the checkpoint:

```text
SRHIT k parent_id step presentation_code
```

The parent order is `2^k`; the target order is 1024 and `step=10-k`.
A completion record is `DONE k parent_id step accepted_count descendant_count`.
The delivered checkpoint has 539 completed jobs, 581 accepted records and no
missing expected jobs. The 581 accepted records agree with the high-class CSV
in row order, parent metadata and exact code.

## Canonical order-512 input and provenance

There is one canonical order-512 CSV, `data/SR_groups_results512.csv`.
All executable references use that filename. The descendant scripts first
look in their working directory, then in `../data/`; the class-two exporter
also searches `../data/`. Thus the delivered directory layout needs no extra
copy of the input file.

The two historical alias CSV copies were removed. The source modifications
are confined to the filename, its lookup paths, and the associated checkpoint
audit/documentation. Existing archived execution logs are preserved unchanged;
they record the historical run, not a fresh GAP run of this revision. The
reproduction variant already changed the two ANUPQ `Workspace` option names
to the documented `PqWorkspace` in the earlier package; this revision preserves
that choice. `validation/repository_revision_audit.json` records source changes,
removed duplicate files and unchanged checkpoint/CSV hashes.

## Reconstruct a group from its CSV code

The `Code` column consists of arbitrary-precision presentation integers.
Read the entire value as an integer or string, never as a floating-point
number. `No` is only a catalogue row number; it is not a SmallGroups ID.
In a GAP session, the first high-class row can be reconstructed directly:

```gap
a := 23429787501007350394706814058502769492323970612862;;
G := PcGroupCode(a,1024);;
Size(G);                         # 1024
NilpotencyClassOfGroup(G);        # 3
```

For any row, run from this package's root:

```gap
Read("code/read_sr1024_class3plus.g");;
rows := SR1024ReadClass3Plus("data/sr1024_class3plus.csv");;
Length(rows);                    # 581
G := SR1024GroupFromClass3Row(rows[1]);;
SR1024CheckSR(G);                 # true
```

Replace `rows[1]` by `rows[i]` to select the row labelled `No=i`. The helper
checks the schema and exact integer fields; reconstruction checks the order
and nilpotency class. `SR1024CheckSR` computes the two exact Wigner moment
identities from actual group elements. It is an optional check of the selected
group, not a completeness or isomorphism test.

`PcGroupCode(a,1024)` reconstructs the presentation encoded by `a`; conversely,
`CodePcGroup(G)` encodes a pc presentation. Distinct codes do not prove that
groups are nonisomorphic. These functions are documented in the official
[GAP Reference Manual, section 46.9](https://docs.gap-system.org/doc/ref/chap46.html).
The absence of order 1024 from the standard catalogue is documented in the
[SmallGrp manual](https://gap-packages.github.io/smallgrp/doc/chap1.html).

The new helper has been inspected statically and its row parser checked
against all 581 CSV rows. GAP is unavailable in the revision environment,
so the helper has not been executed here. The archived `makeCSV.g` run and
class-two export already reconstructed the delivered codes using the same
GAP `PcGroupCode` operation.

## Quick checks without new group generation

Run these commands from the package root. They require only Python's standard
library and finish quickly. Do not use Python's `-O` option, because the
checkpoint audit uses assertions.

```sh
cd validation
python3 ../code/audit_census_checkpoint.py
cd ..
python3 code/counting_structures.py data/f_of_n.csv --class2-audit code/audited_class2.csv --class2-1024 data/sr1024_class2.csv --output-dir validation
python3 code/abelian_hall_counts.py data/f_of_n.csv --output validation/abelian_hall_counts.json
python3 code/verify_progressions.py data/f_of_n.csv --output validation/progression_checks.json
```

The checkpoint audit checks 539 jobs and 581 accepted groups. The arithmetic
checks cover all 500 twice-odd orders; all 250 four-times-odd orders; all 107
available orders `4p^e`; all 204 four-times-squarefree-odd orders; all 193
applicable nilpotent-number orders; and all 51 orders `8p` with prime `p>=5`.
They also perform independent formal-Euler inversion checks. Consequences
include 382 directly indecomposable groups at order 1024, of which 91 have
class two. These arithmetic checks use the completed census as input; they
do not independently certify that census by rerunning its searches.

## Re-export the high-class CSV from the saved checkpoint

Use a separate working directory to preserve the delivered CSV and checkpoint:

```sh
mkdir -p work/class3_export
cp code/makeCSV.g work/class3_export/
cp code/scriptG_results.FINAL.bak work/class3_export/scriptG_results.log
cd work/class3_export
gap -q -b makeCSV.g < /dev/null > makeCSV.log 2>&1
```

Expected output is `sr1024_class3plus.csv` with 581 rows, per-step counts
`[474,105,2,0,0,0,0]`, and 221 groups with a C2 direct factor. This reconstructs
and retests supplied presentations; it does not rerun ANUPQ. Require the
completion summary in the log; a shell exit status alone does not establish
that an interactive GAP script reached the end.

## Re-export class two and the combined catalogue

From `code/`:

```sh
gap -q -b export_sr1024_class2.g < /dev/null > export_sr1024_class2.log 2>&1
gap -q -b build_complete_1024.g < /dev/null > build_complete_1024.log 2>&1
```

The class-two exporter must reach `COMPLETE 221`: 126 special groups, split
28/80/18, and 95 products with C2. The combined export must give 803 groups.
These commands regenerate output copies in `code/`; the published CSV files
remain in `data/`. GAP with SmallGrp is sufficient for these exports. ANUPQ
and graph software are needed only for the corresponding full generation or
orbit-classification stages.

## Rerun descendant generation only when intended

A new generation run can be long and resource intensive. From the package root:

```sh
mkdir -p work/descendants
cp code/scriptG_reproduce.g work/descendants/
cp data/SR_groups_results512.csv work/descendants/
cd work/descendants
gap -q -b scriptG_reproduce.g < /dev/null > scriptG_console.log 2>&1
```

Install GAP, SmallGrp, ANUPQ and its dependencies first. The output checkpoint
is `scriptG_results.log`. The original checkpoint format writes accepted
records before a job's DONE line. After interruption within a job, remove
that job's provisional accepted records before restarting it, or introduce
transactional per-job output. The supplied completed checkpoint has no such
repeated records. Full orbit generation is documented separately in the
computation READMEs; the present revision did not repeat either large search.
