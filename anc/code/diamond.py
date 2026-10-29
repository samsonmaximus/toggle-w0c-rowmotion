"""Check, for crystallographic root posets, the hypothesis of Lemma 5.3
(two distinct upper covers of a root always have a common upper cover) and
that the squares span the cycle space of the Hasse diagram (Corollary 5.4,
Lemma 5.5).  Exact arithmetic: roots are integer vectors, and the rank of
the square cycles is computed modulo the prime 2^61 - 1 (full rank modulo a
prime implies full rank over Q).

Usage: python3 diamond.py
"""
import itertools

P = 2**61 - 1


def cartan(t, n):
    """A[i][j] = <alpha_j, alpha_i^vee> (Bourbaki numbering)."""
    A = [[2 if i == j else 0 for j in range(n)] for i in range(n)]

    def link(i, j, aij=-1, aji=-1):
        A[i][j], A[j][i] = aij, aji

    if t in "ABC":
        for k in range(n - 2):
            link(k, k + 1)
        if n >= 2:
            if t == "A":
                link(n - 2, n - 1)
            elif t == "B":   # alpha_n short
                link(n - 2, n - 1, -1, -2)
            else:            # C: alpha_n long
                link(n - 2, n - 1, -2, -1)
    elif t == "D":
        for k in range(n - 2):
            link(k, k + 1)
        link(n - 3, n - 1)
    elif t == "E":
        for a, b in [(1, 3), (3, 4), (4, 5), (2, 4)] + [(k, k + 1) for k in range(5, n)]:
            link(a - 1, b - 1)
    elif t == "F":           # alpha_1, alpha_2 long; alpha_3, alpha_4 short
        link(0, 1)
        link(1, 2, -1, -2)
        link(2, 3)
    elif t == "G":           # alpha_1 short
        link(0, 1, -3, -1)
    return A


def positive_roots(A):
    n = len(A)
    simple = [tuple(int(i == j) for j in range(n)) for i in range(n)]
    roots, frontier = set(simple), list(simple)
    while frontier:
        new = []
        for r in frontier:
            for i in range(n):
                c = sum(r[j] * A[i][j] for j in range(n))
                s = tuple(r[j] - (c if j == i else 0) for j in range(n))
                if all(x >= 0 for x in s) and any(s) and s not in roots:
                    roots.add(s)
                    new.append(s)
        frontier = new
    return roots


def plus(r, *idx):
    return tuple(r[j] + sum(1 for i in idx if i == j) for j in range(len(r)))


def rank_mod_p(rows):
    pivots, rank = {}, 0
    for row in rows:
        r = {k: v % P for k, v in row.items() if v % P}
        while r:
            c = min(r)
            if c not in pivots:
                pivots[c] = r
                rank += 1
                break
            pr = pivots[c]
            f = r[c] * pow(pr[c], P - 2, P) % P
            for k, v in pr.items():
                r[k] = (r.get(k, 0) - f * v) % P
                if r[k] == 0:
                    del r[k]
    return rank


def check(t, n, expected_N):
    A = cartan(t, n)
    R = positive_roots(A)
    assert len(R) == expected_N, (t, n, len(R))
    idx = {r: k for k, r in enumerate(sorted(R))}
    failures = 0
    for g in R:
        ups = [i for i in range(n) if plus(g, i) in idx]
        for s, u in itertools.combinations(ups, 2):
            if plus(g, s, u) not in idx:
                failures += 1
    edges = {}
    for g in R:
        for i in range(n):
            if plus(g, i) in idx:
                edges[(idx[g], idx[plus(g, i)])] = len(edges)
    squares = []
    for g in R:
        for s, u in itertools.combinations(range(n), 2):
            a, b, c = plus(g, s), plus(g, u), plus(g, s, u)
            if a in idx and b in idx and c in idx:
                squares.append({edges[(idx[g], idx[a])]: 1, edges[(idx[a], idx[c])]: 1,
                                edges[(idx[b], idx[c])]: -1, edges[(idx[g], idx[b])]: -1})
    cyc = len(edges) - len(R) + 1
    rk = rank_mod_p(squares)
    ok = failures == 0 and rk == cyc
    print(f"{t}{n}: roots {len(R)}, edges {len(edges)}, cycle space dim {cyc}, "
          f"squares {len(squares)}, rank {rk}, diamond failures {failures} -> {'OK' if ok else 'FAIL'}")
    return ok


cases = ([("A", n, n * (n + 1) // 2) for n in range(1, 10)]
         + [("B", n, n * n) for n in range(2, 9)]
         + [("C", n, n * n) for n in range(3, 9)]
         + [("D", n, n * (n - 1)) for n in range(4, 9)]
         + [("E", 6, 36), ("E", 7, 63), ("E", 8, 120), ("F", 4, 24), ("G", 2, 6)])
results = [check(*c) for c in cases]
print("ALL OK" if all(results) else "SOME FAILED")
