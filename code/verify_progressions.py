#!/usr/bin/env python3
"""Check the proved arithmetic progressions against an existing SR census.

Uses only the Python standard library. These are checks of numerical
consequences of the proofs, not independent group enumerations.
"""
import argparse
import csv
import hashlib
import json
from math import gcd, isqrt, prod
from pathlib import Path

from counting_structures import abelian_count, divisors, partition_number


def factor(n):
    answer = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            answer[p] = answer.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        answer[n] = 1
    return answer


def h(m):
    root = isqrt(m)
    numerator = sum(abelian_count(d) * abelian_count(m // d) for d in divisors(m))
    numerator += abelian_count(root) if root * root == m else 0
    assert numerator % 2 == 0
    return numerator // 2


def four_prime_power(e):
    numerator = sum(partition_number(i) * partition_number(e - i) for i in range(e + 1))
    numerator += partition_number(e // 2) if e % 2 == 0 else 0
    assert numerator % 2 == 0
    return numerator // 2


def is_nilpotent_number(m):
    psi = prod(p ** i - 1 for p, a in factor(m).items() for i in range(1, a + 1))
    return gcd(m, psi) == 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("census", type=Path)
    parser.add_argument("--output", type=Path, default=Path("progression_checks.json"))
    args = parser.parse_args()
    with args.census.open() as handle:
        records = [(int(n), int(f)) for n, f in csv.reader(handle)]
    f = dict(records)
    limit = max(f)
    assert len(records) == len(f)
    assert set(f) == set(range(2, limit + 1, 2))
    checked = {key: [] for key in ["twice_odd", "four_prime_power", "four_squarefree_odd",
                                  "four_nilpotent_number", "eight_prime_at_least_five"]}
    exceptions = []
    for m in range(1, limit // 2 + 1, 2):
        assert f[2 * m] == abelian_count(m), ("twice_odd", m)
        checked["twice_odd"].append(2 * m)
    for m in range(1, limit // 4 + 1, 2):
        baseline = h(m)
        delta = f[4 * m] - baseline
        assert delta >= 0
        if delta:
            exceptions.append({"order": 4 * m, "odd_part": m, "h": baseline, "e": delta})
        fac = factor(m)
        if len(fac) <= 1:
            e = next(iter(fac.values()), 0)
            assert f[4 * m] == four_prime_power(e), ("four_prime_power", m)
            checked["four_prime_power"].append(4 * m)
        if all(e == 1 for e in fac.values()):
            expected = 2 ** (len(fac) - 1) if m > 1 else 1
            assert f[4 * m] == expected, ("four_squarefree_odd", m)
            checked["four_squarefree_odd"].append(4 * m)
        if is_nilpotent_number(m):
            assert f[4 * m] == baseline, ("four_nilpotent_number", m)
            checked["four_nilpotent_number"].append(4 * m)
    for p in range(5, limit // 8 + 1, 2):
        if factor(p) == {p: 1}:
            assert f[8 * p] == 3, ("eight_prime", p)
            checked["eight_prime_at_least_five"].append(8 * p)
    if limit == 2000:
        assert [(x["order"], x["e"]) for x in exceptions] == [(300, 1), (1500, 1), (1620, 1)]
        assert f[24] == 4
        assert (f[6], f[18], f[12], f[20], f[60], f[300], f[588]) == (1, 2, 1, 1, 2, 6, 5)
    sequence = [four_prime_power(e) for e in range(11)]
    assert sequence == [1, 1, 3, 5, 11, 18, 34, 55, 95, 150, 244]
    result = {
        "input_sha256": hashlib.sha256(args.census.read_bytes()).hexdigest(),
        "limit": limit,
        "all_applicable_checks_passed": True,
        "checked_order_counts": {key: len(value) for key, value in checked.items()},
        "checked_orders": checked,
        "four_times_odd_exceptions_in_census": exceptions,
        "four_prime_power_sequence_e_0_through_10": sequence,
        "constructive_bounds_beyond_census_not_exact_counts": {
            str(2 * q * q * (q + 1)): h(q * q * (q + 1) // 2) + 1 for q in (13, 17)
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "checked_orders"}, indent=2))


if __name__ == "__main__":
    main()
