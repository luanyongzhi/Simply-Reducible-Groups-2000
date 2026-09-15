#!/usr/bin/env python3
"""Audit the completed SR census checkpoint without generating descendants.

Run from the archive's code directory:
    python3 audit_census_checkpoint.py

Or specify an unpacked archive root explicitly:
    python3 /path/to/audit_census_checkpoint.py --root /path/to/archive

Inputs are found in ROOT/data, ROOT/code, or ROOT itself. The successful
result is written as census_checkpoint_audit.json in the current directory.
Only the Python standard library is required. No GAP process is started.
"""
from pathlib import Path
import argparse
import csv, collections, hashlib, json

if not __debug__:
    raise RuntimeError("Run this audit without Python's -O option; assertions are required.")
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1],
                    help='unpacked archive root (default: parent of the script directory)')
ROOT = parser.parse_args().root.resolve()

def locate(name):
    for directory in (ROOT / 'data', ROOT / 'code', ROOT):
        path = directory / name
        if path.is_file():
            return path
    raise FileNotFoundError(f"Missing {name}: searched data/, code/, and root under {ROOT}")
def read_csv(path):
    with path.open(newline='') as stream:
        return [{k.strip(): v.strip() for k, v in row.items()} for row in csv.DictReader(stream)]
def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
main = read_csv(locate('SR_groups_results2_2000.csv'))
p512 = read_csv(locate('SR_groups_results512.csv'))
p1536 = read_csv(locate('SR_groups_results1536.csv'))
target = read_csv(locate('sr1024_class3plus.csv'))
assert len(p512) == 317
keys = [(int(r['Order']), int(r['SmallGroup_ID'])) for r in main+p512+p1536]
assert len(keys) == len(set(keys)) == 7086
expected = set()
for row in main + p512:
    order = int(row['Order'])
    if order not in (8,16,32,64,128,256,512): continue
    k = order.bit_length()-1
    if int(row['CenterSize']) == order: continue  # abelian SR = elementary abelian
    expected.add((k, int(row['SmallGroup_ID']), 10-k))
lines = [line.split() for line in locate('scriptG_results.FINAL.bak').read_text().splitlines()]
assert all(row[0] in ('DONE','SRHIT') for row in lines)
done = [tuple(map(int,row[1:])) for row in lines if row[0]=='DONE']
hits = [tuple(map(int,row[1:])) for row in lines if row[0]=='SRHIT']
assert all(len(row)==5 for row in done)
assert all(len(row)==4 for row in hits)
assert len(done) == len({row[:3] for row in done}) == 539
assert {row[:3] for row in done} == expected
assert len(hits) == len(set(hits)) == 581
perjob = collections.Counter(row[:3] for row in hits)
assert all(perjob[row[:3]] == row[3] for row in done)
assert set(perjob) <= expected
assert [(2**k, i, s, str(a)) for k,i,s,a in hits] == [
    (int(row['ParentOrder']),int(row['ParentId']),int(row['Step']),row['Code'])
    for row in target]
result = {
    'catalogue_records': {'other_orders':len(main),'512':len(p512),'1536':len(p1536),'total':len(keys)},
    'canonical_order512_records':len(p512),
    'job_roster_exact':True,
    'jobs':len(done),
    'positive_jobs':len(perjob),
    'zero_output_jobs':len(done)-len(perjob),
    'hits':len(hits),
    'descendants':sum(row[4] for row in done),
    'checkpoint_csv_recordwise_match':True,
    'per_step':[
        {'step':s,'parent_order':2**(10-s),'jobs':sum(row[2]==s for row in done),
         'descendants':sum(row[4] for row in done if row[2]==s),
         'sr':sum(row[2]==s for row in hits)} for s in range(1,8)],
    'sha256': {name:sha256(locate(name)) for name in ('scriptG.g','makeCSV.g','scriptG_console.FINAL.bak','scriptG_results.FINAL.bak','SR_groups_results512.csv','sr1024_class3plus.csv')},
}
(Path.cwd()/'census_checkpoint_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
