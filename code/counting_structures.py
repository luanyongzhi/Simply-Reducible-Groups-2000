#!/usr/bin/env python3
"""Exact arithmetic consequences of an SR census; Python standard library only.

No group or character computations occur here. The indecomposable extraction
uses direct-product closure and Remak--Krull--Schmidt. The lower bound counts
products of an arbitrary SR 2-group and generalized dihedral groups over
nontrivial abelian odd groups. Input CSV consists of order,count rows, no header.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path


def indecomposable_counts(f):
    """Triangular inversion of the formal Dirichlet Euler product."""
    limit = len(f) - 1
    reconstructed = [0] * (limit + 1)
    reconstructed[1] = 1
    indecomposable = [0] * (limit + 1)
    for order in range(2, limit + 1):
        number = indecomposable[order] = f[order] - reconstructed[order]
        if number < 0:
            raise ValueError(f"negative indecomposable count at {order}: {number}")
        if not number:
            continue
        old = reconstructed[:]
        power, multiplicity = order, 1
        while power <= limit:
            multiset_count = math.comb(number + multiplicity - 1, multiplicity)
            for cofactor in range(1, limit // power + 1):
                reconstructed[cofactor * power] += old[cofactor] * multiset_count
            power *= order
            multiplicity += 1
    assert reconstructed == f
    return indecomposable


def dirichlet_convolution(left, right):
    limit = len(left) - 1
    result = [Fraction(0)] * (limit + 1)
    for a in range(1, limit + 1):
        if left[a]:
            for b in range(1, limit // a + 1):
                if right[b]:
                    result[a * b] += left[a] * right[b]
    return result


def verify_with_formal_log(f, indecomposable):
    """Independently check log F = sum_{r,n>=2} d(n)/(r n^(rs))."""
    limit = len(f) - 1
    base = list(map(Fraction, f))
    base[1] = 0
    power = base
    computed = [Fraction(0)] * (limit + 1)
    for exponent in range(1, limit.bit_length()):
        coefficient = Fraction((-1) ** (exponent + 1), exponent)
        for n in range(2, limit + 1):
            computed[n] += coefficient * power[n]
        power = dirichlet_convolution(power, base)
    expected = [Fraction(0)] * (limit + 1)
    for n in range(2, limit + 1):
        power, exponent = n, 1
        while power <= limit:
            expected[power] += Fraction(indecomposable[n], exponent)
            power *= n
            exponent += 1
    assert computed == expected


@lru_cache(None)
def partition_number(n):
    values = [1] + [0] * n
    for part in range(1, n + 1):
        for total in range(part, n + 1):
            values[total] += values[total - part]
    return values[n]


@lru_cache(None)
def abelian_count(n):
    answer, p = 1, 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        answer *= partition_number(exponent)
        p += 1
    return answer  # residual prime contributes p(1)=1


@lru_cache(None)
def divisors(n):
    return tuple(d for d in range(1, n + 1) if n % d == 0)


@lru_cache(None)
def dihedral_multisets(k, m, allow_trivial=True):
    """Cycle-index recurrence: k D_k(m)=sum_i sum_{d^i|m} a(d)D_{k-i}(m/d^i)."""
    if k == 0:
        return int(m == 1)
    total = 0
    for d in divisors(m):
        if d == 1 and not allow_trivial:
            continue
        power = d
        for cycle_length in range(1, k + 1):
            if m % power == 0:
                total += abelian_count(d) * dihedral_multisets(
                    k - cycle_length, m // power, allow_trivial
                )
            power *= d
            if power > m:
                break
    assert total % k == 0
    return total // k


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("census", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--class2-audit", type=Path,
                        help="Optional audited_class2.csv covering orders through 512")
    parser.add_argument("--class2-1024", type=Path,
                        help="Optional sr1024_class2.csv covering all 221 class-two types")
    args = parser.parse_args()
    records = [(int(n), int(v)) for n, v in csv.reader(args.census.open())]
    limit = max(n for n, _ in records)
    f = [0] * (limit + 1)
    f[1] = 1
    for n, value in records:
        f[n] = value
    assert set(n for n, _ in records) == set(range(2, limit + 1, 2))
    indecomposable = indecomposable_counts(f)
    verify_with_formal_log(f, indecomposable)
    rows = []
    for n in range(2, limit + 1, 2):
        k, m = 0, n
        while m % 2 == 0:
            k, m = k + 1, m // 2
        baseline = dihedral_multisets(k, m)
        strengthened = sum(
            f[2 ** (k - j)] * dihedral_multisets(j, m, False)
            for j in range(k + 1)
        )
        assert 0 <= baseline <= strengthened <= f[n]
        generalized_dihedral = abelian_count(n // 2)
        assert generalized_dihedral <= f[n]
        assert baseline == sum(dihedral_multisets(j, m, False) for j in range(k + 1))
        rows.append((n, f[n], indecomposable[n], f[n] - indecomposable[n],
                     baseline, strengthened, f[n] - strengthened, generalized_dihedral))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "sr_counting_structures.csv").open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow(("order", "sr_count", "directly_indecomposable_count",
                         "directly_decomposable_count", "dihedral_product_count",
                         "two_group_dihedral_product_count", "outside_product_subclass",
                         "all_abelian_generalized_dihedral_count"))
        writer.writerows(rows)
    certificate = {
        "input_sha256": hashlib.sha256(args.census.read_bytes()).hexdigest(),
        "limit": limit,
        "formal_euler_reconstruction": True,
        "independent_formal_log_check": True,
        "all_lower_bounds_valid": True,
        "sr_nontrivial_total": sum(f) - 1,
        "indecomposable_nontrivial_total": sum(indecomposable),
        "decomposable_nontrivial_total": sum(f) - 1 - sum(indecomposable),
        "pure_two_power_counts": [
            {"order": 2 ** k, "sr": f[2 ** k], "indecomposable": indecomposable[2 ** k]}
            for k in range(1, limit.bit_length())
        ],
        "dihedral_product_subclass_total": sum(row[4] for row in rows),
        "two_group_dihedral_product_subclass_total": sum(row[5] for row in rows),
        "dihedral_baseline_exact_orders": sum(row[1] == row[4] for row in rows),
        "strengthened_bound_exact_orders": sum(row[1] == row[5] for row in rows),
        "generalized_dihedral_bound_exact_orders": sum(row[1] == row[7] for row in rows),
    }
    if args.class2_audit or args.class2_1024:
        if not (args.class2_audit and args.class2_1024):
            parser.error("supply both --class2-audit and --class2-1024")
        lower = list(csv.DictReader(args.class2_audit.open()))
        upper = list(csv.DictReader(args.class2_1024.open()))
        counts = Counter(int(row["Order"]) for row in lower)
        counts.update(int(row["Order"]) for row in upper)
        class2_limit = max(counts)
        at_most_two = [0] * (class2_limit + 1)
        for k in range(class2_limit.bit_length()):
            at_most_two[2 ** k] = counts[2 ** k] + 1
        indec_two = indecomposable_counts(at_most_two)
        verify_with_formal_log(at_most_two, indec_two)
        certificate["class2_inputs_sha256"] = {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [args.class2_audit, args.class2_1024]
        }
        certificate["class_at_most_two_independent_formal_log_check"] = True
        certificate["class2_directly_indecomposable"] = [
            {"order": 2 ** k, "class2_count": counts[2 ** k],
             "class2_indecomposable": indec_two[2 ** k] if k >= 2 else 0}
            for k in range(1, class2_limit.bit_length())
        ]
        assert indec_two[1024] == 91
        special_1024 = sum(row["IsSpecial"].lower() == "true" for row in upper)
        assert special_1024 == 126
        assert special_1024 - indec_two[1024] == 2 * 10 + math.comb(6, 2)
        certificate["special_1024"] = {
            "total": special_1024, "indecomposable": indec_two[1024],
            "decomposable": special_1024 - indec_two[1024],
            "product_8_times_128": 20, "product_32_times_32": 15,
        }
    (args.output_dir / "counting_structures_certificate.json").write_text(
        json.dumps(certificate, indent=2) + "\n"
    )
    print(json.dumps(certificate, indent=2))


if __name__ == "__main__":
    main()
