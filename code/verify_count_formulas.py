#!/usr/bin/env python3
"""Exact integer checks of the arithmetic formulas against the census table."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def partition_number(e):
    values = [1] + [0] * e
    for part in range(1, e + 1):
        for target in range(part, e + 1):
            values[target] += values[target - part]
    return values[e]


def abelian_count(m):
    total = 1
    p = 2
    while p * p <= m:
        exponent = 0
        while m % p == 0:
            m //= p
            exponent += 1
        if exponent:
            total *= partition_number(exponent)
        p += 1
    return total  # any remaining prime has exponent 1


def dihedral_pair_count(m):
    ordered = sum(abelian_count(d) * abelian_count(m // d)
                  for d in range(1, m + 1) if m % d == 0)
    square_root = math.isqrt(m)
    diagonal = abelian_count(square_root) if square_root ** 2 == m else 0
    assert (ordered + diagonal) % 2 == 0
    return (ordered + diagonal) // 2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('counts_csv', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    counts = {int(n): int(f) for n, f in csv.reader(args.counts_csv.open())}
    twice_odd = sorted(n for n in counts if n % 4 == 2)
    four_odd = sorted(n for n in counts if n % 8 == 4)
    assert all(counts[n] == abelian_count(n // 2) for n in twice_odd)
    differences = [dict(order=n, f=counts[n], h=dihedral_pair_count(n // 4),
                        excess=counts[n] - dihedral_pair_count(n // 4))
                   for n in four_odd if counts[n] != dihedral_pair_count(n // 4)]
    assert differences == [dict(order=300, f=6, h=5, excess=1),
                           dict(order=1500, f=11, h=10, excess=1),
                           dict(order=1620, f=21, h=20, excess=1)]
    assert all(counts[2*n] >= f for n, f in counts.items() if 2*n in counts)
    result = dict(source_sha256=hashlib.sha256(args.counts_csv.read_bytes()).hexdigest(),
                  twice_odd_orders_checked=len(twice_odd),
                  four_odd_orders_checked=len(four_odd),
                  four_odd_records=sum(counts[n] for n in four_odd),
                  abelian_hall_records=sum(dihedral_pair_count(n//4) for n in four_odd),
                  exceptions=differences, all_available_doubling_inequalities_pass=True)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
