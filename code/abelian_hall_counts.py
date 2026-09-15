"""Exact multiset counts for SR groups with abelian odd Hall and Sylow 2.
No external modules. This verifies arithmetic against the supplied census;
it does not independently identify the Sylow/Hall subgroups of census groups.
"""
from collections import defaultdict
from math import comb
import argparse, csv, hashlib, json
from pathlib import Path
from counting_structures import dihedral_multisets
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("census", type=Path)
parser.add_argument("--output", type=Path, default=Path("abelian_hall_counts.json"))
args = parser.parse_args()
source = args.census
census = dict((int(n), int(f)) for n, f in csv.reader(source.open()))
LIMIT = max(census)
MAX_E = LIMIT.bit_length()
part = [0] * (MAX_E + 1)
part[0] = 1
for d in range(1, MAX_E + 1):
    for j in range(d, MAX_E + 1):
        part[j] += part[j-d]
def abelian_count(n):
    out, p = 1, 2
    while p*p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        out *= part[e]
        p += 1
    if n > 1:
        out *= part[1]
    return out
# Nontrivial odd abelian factors; each costs one C2 in Dih(A).
counts = {(0,1): 1}
for d in range(3, LIMIT//2+1, 2):
    colors = abelian_count(d)
    updated = defaultdict(int)
    for (k,m), value in counts.items():
        r, power = 0, 1
        while 2**(k+r)*m*power <= LIMIT:
            updated[k+r,m*power] += value*comb(colors+r-1,r)
            r += 1
            power *= d
    counts = dict(updated)
# Arbitrarily many trivial factors Dih(1)=C2.
baseline = defaultdict(int)
for (k,m), value in counts.items():
    while 2**k*m <= LIMIT:
        baseline[2**k*m] += value
        k += 1
assert all(baseline[n] <= f for n,f in census.items())
for n in census:
    k, m = 0, n
    while m % 2 == 0:
        k, m = k + 1, m // 2
    assert baseline[n] == dihedral_multisets(k, m), n
assert all(baseline[2*m] == abelian_count(m) for m in range(1,LIMIT//2+1,2))
def h(m):
    conv=sum(abelian_count(d)*abelian_count(m//d) for d in range(1,m+1) if m%d==0)
    root=int(m**.5)
    return (conv+(abelian_count(root) if root*root==m else 0))//2
assert all(baseline[4*m] == h(m) for m in range(1,LIMIT//4+1,2))
# Exact tensor multiplicity for the minimal even circuit over C3^4.
from itertools import product
signs=[s for s in product((1,-1),repeat=4) if s[0]*s[1]*s[2]*s[3]==1]
def orbit(v):
    return {tuple(v[i]*s[i]%3 for i in range(4)) for s in signs}
lam=(1,1,1,0); mu=(2,0,0,1); nu=(0,1,1,1)
mult=sum(tuple((lam[i]+w[i])%3 for i in range(4)) in orbit(nu) for w in orbit(mu))
assert mult==2
result={
 'bound':LIMIT,
 'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
 'subclass_total_nontrivial':sum(baseline[n] for n in census),
 'census_total':sum(census.values()),
 'all_baselines_at_most_census':True,
 'multiset_product_equals_cycle_index_at_every_order':True,
 'twice_odd_formulas_checked':len(range(1,LIMIT//2+1,2)),
 'four_times_odd_formulas_checked':len(range(1,LIMIT//4+1,2)),
 'minimal_even_circuit':{'group_order':648,'lambda_orbit':len(orbit(lam)),'mu_orbit':len(orbit(mu)),'nu_orbit':len(orbit(nu)),'tensor_multiplicity':mult},
 'examples':{str(n):{'B':baseline[n],'f':census[n]} for n in [6,18,36,60,72,108,120,180,216,300,360,648,840,1024,1500,1620,2000]},
 'by_v2':{str(k):{'subclass':sum(baseline[n] for n in census if n%(2**k)==0 and n%(2**(k+1))!=0),'census':sum(f for n,f in census.items() if n%(2**k)==0 and n%(2**(k+1))!=0)} for k in range(1,11)}
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
