#!/usr/bin/env python3
"""Exact polar/quadric orbit classification via coloured incidence graphs.

Requires pynauty==2.8.8.1. Input: m r affine_dimension B1 ... Br per line,
where B_a is the bit mask of upper-triangular coefficients in lexicographic
(i,j), i<j order. This classifies the supplied candidate universe only;
its exhaustive coverage is the responsibility of the producer's theorem.

Vertices have five ordered colours: nonzero V points; zero W point;
nonzero W points; V projective lines; W projective lines. A V-line is joined
to its three points and beta(x,y); a W-line is joined to its three points.
Since beta is constant on a two-dimensional F2 subspace, this is exact.
Addition-line preservation forces point permutations to be linear, so graph
isomorphism equals GL(V)xGL(W) tensor equivalence. For a quadratic pair, add
one edge from each nonzero V point to its q-value in W.
"""
from __future__ import annotations
import argparse, csv, hashlib, itertools, json, time
from collections import Counter
from functools import lru_cache
from pathlib import Path
import pynauty


@lru_cache(None)
def pairs(m):
    return tuple(itertools.combinations(range(m), 2))


@lru_cache(None)
def lines(m):
    return tuple((x, y, x ^ y) for x in range(1, 1 << m)
                 for y in range(x + 1, 1 << m) if y < (x ^ y))


def wedge(x, y, m):
    return sum(((((x >> i) & 1) & ((y >> j) & 1)) ^
                 (((x >> j) & 1) & ((y >> i) & 1))) << k
               for k, (i, j) in enumerate(pairs(m)))


@lru_cache(None)
def line_wedges(m):
    return tuple(wedge(x, y, m) for x, y, z in lines(m))


def apply_beta(B, z):
    return sum(((b & z).bit_count() & 1) << a for a, b in enumerate(B))


def graph(m, r, B, q_values=None):
    nv, nw = (1 << m) - 1, 1 << r
    vl, wl = lines(m), lines(r)
    base, wbase = nv + nw, nv + nw + len(vl)
    n = wbase + len(wl)
    adj = {i: [] for i in range(n)}
    def edge(x, y):
        adj[x].append(y); adj[y].append(x)
    for k, ((x, y, z), bz) in enumerate(zip(vl, line_wedges(m))):
        node = base + k
        for v in (x - 1, y - 1, z - 1, nv + apply_beta(B, bz)):
            edge(node, v)
    for k, triple in enumerate(wl):
        for w in triple:
            edge(wbase + k, nv + w)
    if q_values is not None:
        for x in range(1, 1 << m):
            edge(x - 1, nv + q_values[x])
    colors = [set(range(nv)), {nv}, set(range(nv + 1, nv + nw)),
              set(range(base, wbase)), set(range(wbase, n))]
    return pynauty.Graph(n, directed=False, adjacency_dict=adj,
                        vertex_coloring=[c for c in colors if c])


def solve(rows, n):
    """Return one solution plus a kernel basis; rows include RHS at bit n."""
    pivots = {}
    for row in rows:
        while row & ((1 << n) - 1):
            col = (row & ((1 << n) - 1) & -row).bit_length() - 1
            if col in pivots:
                row ^= pivots[col]
            else:
                pivots[col] = row
                break
        else:
            if row >> n:
                return None
    def back(x, homogeneous):
        for col in sorted(pivots, reverse=True):
            row = pivots[col]
            bit = ((row & x).bit_count() & 1)
            if not homogeneous:
                bit ^= (row >> n) & 1
            x ^= bit << col
        return x
    return back(0, False), [back(1 << col, True) for col in range(n) if col not in pivots]


def scalar_rows(mask, m):
    rows = [0] * m
    for k, (i, j) in enumerate(pairs(m)):
        if mask >> k & 1:
            rows[i] ^= 1 << j; rows[j] ^= 1 << i
    return rows


@lru_cache(None)
def square_wedges(m):
    return tuple(sum((((v >> i) & 1) & ((v >> j) & 1)) << k
                     for k, (i, j) in enumerate(pairs(m))) for v in range(1 << m))


def reality_affine(m, r, B):
    eqs = []
    for lam in range(1, 1 << r):
        mask = 0
        for a, b in enumerate(B):
            if lam >> a & 1: mask ^= b
        rad = solve(scalar_rows(mask, m), m)[1]
        for v in rad:
            row = sum(lam << (i * r) for i in range(m) if v >> i & 1)
            rhs = (mask & square_wedges(m)[v]).bit_count() & 1
            eqs.append(row | rhs << (m * r))
    return solve(eqs, m * r)


def affine_elements(solution):
    first, basis = solution
    out = [first]
    for v in basis:
        out += [x ^ v for x in out]
    return out


def quadratic_values(m, r, B, Q):
    qbasis = [(Q >> (i * r)) & ((1 << r) - 1) for i in range(m)]
    out = []
    for v, sq in enumerate(square_wedges(m)):
        value = apply_beta(B, sq)
        for i in range(m):
            if v >> i & 1: value ^= qbasis[i]
        out.append(value)
    return out


def linear_image(x, basis):
    value = 0
    for i, y in enumerate(basis):
        if x >> i & 1: value ^= y
    return value


def stabilizer_actions(m, r, B, g):
    gens, order1, order2, _, _ = pynauty.autgrp(g)
    nv = (1 << m) - 1
    actions = []
    for perm in gens:
        T = [perm[(1 << i) - 1] + 1 for i in range(m)]
        S = [perm[nv + (1 << a)] - nv for a in range(r)]
        for x in range(1, 1 << m):
            assert linear_image(x, T) == perm[x - 1] + 1
        for w in range(1 << r):
            assert linear_image(w, S) == perm[nv + w] - nv
        invS = [0] * (1 << r)
        for w in range(1 << r): invS[linear_image(w, S)] = w
        for i, j in pairs(m):
            assert apply_beta(B, wedge(T[i], T[j], m)) == linear_image(
                apply_beta(B, wedge(1 << i, 1 << j, m)), S)
        actions.append((T, invS))
    return actions, [order1, order2]


def classify_quadratics(m, r, B, g):
    solution = reality_affine(m, r, B)
    assert solution is not None
    candidates = set(affine_elements(solution))
    actions, group_order = stabilizer_actions(m, r, B, g)
    todo = set(candidates); orbits = []; pair_certificates = []
    while todo:
        Q = min(todo); todo.remove(Q); orbit = [Q]
        for current in orbit:
            values = quadratic_values(m, r, B, current)
            for T, invS in actions:
                other = sum(invS[values[T[i]]] << (i * r) for i in range(m))
                assert other in candidates
                if other in todo:
                    todo.remove(other); orbit.append(other)
        values = quadratic_values(m, r, B, Q)
        cert = pynauty.certificate(graph(m, r, B, values))
        pair_certificates.append(cert)
        orbits.append(dict(Q=Q, size=len(orbit), involutions=(1 << r)*values.count(0)-1,
                           pair_certificate_sha256=hashlib.sha256(cert).hexdigest()))
    assert sum(o['size'] for o in orbits) == len(candidates)
    # Full pair graph canonicalization independently separates stabilizer orbits.
    assert len(set(pair_certificates)) == len(orbits)
    return dict(affine_dimension=len(solution[1]), candidate_count=len(candidates),
                stabilizer_generators=len(actions), graph_aut_order_scientific=group_order,
                quadratic_orbits=orbits)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    start = time.monotonic(); unique = {}; membership = []; rep_data = []
    for rownum, line in enumerate(args.input.read_text().splitlines(), 1):
        if not line.strip() or line.startswith('#'): continue
        m, r, dim, *B = map(int, line.split())
        assert len(B) == r
        g = graph(m, r, B)
        cert = pynauty.certificate(g)
        key = (m, r, cert)  # Exact canonical bytes, never hash-only deduplication.
        if key not in unique:
            idx = len(rep_data) + 1; unique[key] = idx
            rec = dict(polar_orbit=idx, m=m, r=r, B=B, supplied_affine_dimension=dim,
                       graph_vertices=g.number_of_vertices,
                       polar_certificate_sha256=hashlib.sha256(cert).hexdigest(),
                       first_input_line=rownum)
            rec.update(classify_quadratics(m, r, B, g))
            assert rec['affine_dimension'] == dim
            rep_data.append(rec)
            print('NEW', idx, 'm', m, 'r', r, 'Qorbits', len(rec['quadratic_orbits']),
                  'seconds', round(time.monotonic()-start, 3), flush=True)
        membership.append(dict(input_line=rownum, polar_orbit=unique[key]))
        if rownum % 100 == 0:
            print('PROGRESS', rownum, 'polar_orbits', len(unique), flush=True)
    result = dict(complete_for_supplied_input=True, input=str(args.input),
                  input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  pynauty_version=pynauty.__version__, elapsed_seconds=time.monotonic()-start,
                  polar_orbits=rep_data, membership=membership,
                  total_group_orbits=sum(len(p['quadratic_orbits']) for p in rep_data))
    args.output.write_text(json.dumps(result, indent=2))
    print('COMPLETE polar_orbits', len(rep_data), 'group_orbits', result['total_group_orbits'],
          'seconds', round(result['elapsed_seconds'], 3), flush=True)


if __name__ == '__main__':
    main()
