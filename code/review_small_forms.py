"""Independent audit checkpoint: construct all 84 star-plus-symplectic spaces.

Uses high-pivot elimination, explicit sets of radical vectors, and a separate
construction rather than the C++ RREF traversal. No GAP dependency.
"""
from itertools import combinations

M = 5
PAIRS = tuple(combinations(range(M), 2))


def wedge(x, y):
    return sum(((((x >> i) & 1) * ((y >> j) & 1)) ^
                (((x >> j) & 1) * ((y >> i) & 1))) << k
               for k, (i, j) in enumerate(PAIRS))


def form_rows(f):
    a = [0] * M
    for k, (i, j) in enumerate(PAIRS):
        if f >> k & 1:
            a[i] ^= 1 << j
            a[j] ^= 1 << i
    return a


def rank(a):
    basis = {}
    for x in a:
        while x:
            p = x.bit_length() - 1
            if p in basis:
                x ^= basis[p]
            else:
                basis[p] = x
                break
    return len(basis)


def radical(f):
    rows = form_rows(f)
    return {v for v in range(1 << M)
            if all((v & a).bit_count() % 2 == 0 for a in rows)}


def span(basis):
    s = [0]
    for x in basis:
        s += [y ^ x for y in s]
    return s


def q_solution_dimension(basis):
    forms = span(basis)
    nvars = M * len(basis)
    eq = {}
    for lam, f in enumerate(forms):
        for v in radical(f):
            lhs = sum(v << (M * i) for i in range(len(basis)) if lam >> i & 1)
            rhs = sum(1 for k, (i, j) in enumerate(PAIRS)
                      if f >> k & 1 and v >> i & 1 and v >> j & 1) % 2
            while lhs:
                p = lhs.bit_length() - 1
                if p in eq:
                    old, value = eq[p]
                    lhs ^= old
                    rhs ^= value
                else:
                    eq[p] = lhs, rhs
                    break
            else:
                if rhs:
                    return None
    return nvars - len(eq)


def main():
    all_spaces = set()
    per_line = []
    for line in (1, 2, 3):  # three covector lines contained in <e1,e2>
        pivot = (line & -line).bit_length() - 1
        complement = [1 << i for i in range(M) if i != pivot]
        star = [wedge(line, x) for x in complement]
        found = 0
        for coeff in range(64):
            w = 0
            for k, (i, j) in enumerate(combinations(range(4), 2)):
                if coeff >> k & 1:
                    w ^= wedge(complement[i], complement[j])
            if rank(form_rows(w)) != 4:
                continue
            basis = star + [w]
            s = span(basis)
            assert len(set(s)) == 32 and 1 in s
            rad = {f: radical(f) for f in s}
            ranks = {f: M - (len(rad[f]).bit_length() - 1) for f in s}
            for a, b in combinations(s[1:], 2):
                t = len(rad[a] & rad[b]).bit_length() - 1
                assert ranks[a] + ranks[b] + ranks[a ^ b] == 2 * (M - t)
            assert q_solution_dimension(basis) is None
            all_spaces.add(tuple(sorted(s)))
            found += 1
        per_line.append(found)
    assert len(all_spaces) == 84
    print("Independent star-family checkpoint: 3 lines x", per_line,
          "=", len(all_spaces), "distinct spaces; all MF; all q systems inconsistent")
    # Known surviving full star gives a check against falsely rejecting every q.
    star = [wedge(1, 1 << i) for i in range(1, 5)]
    assert q_solution_dimension(star) == 5
    print("Positive affine-system control: full 4-dimensional star has q dimension 5")


if __name__ == "__main__":
    main()
