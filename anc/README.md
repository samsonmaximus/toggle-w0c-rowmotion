# Ancillary files: "Toggling in w0(c) root order is conjugate to rowmotion: resolution of AIM Conjecture 3.2 and Thomas's conjecture"

Author: Samson Fraser (samsonmfraser@gmail.com).
Requirements: Python 3.11 with sympy, networkx and python-flint; a C++17 compiler (g++).

## `code/` — main checks
Compile the C++ checkers with `g++ -O2 -o fastcheck fastcheck.cpp` (and the same for `squares` and `omega`).
Type letters are A–G; E needs n = 6, 7 or 8, F needs n = 4, G needs n = 2. An unknown type is rejected.

| Claim in the paper | Command | Expected output |
|---|---|---|
| Theorem 6.1, types E6, E7, E8, F4, G2: every w0(c) orientation has a potential (exact integer check) | `./fastcheck E 6`, `E 7`, `E 8`, `F 4`, `G 2` | "without potential 0" for 720, 5040, 40320, 24 and 2 Coxeter words. Logs: `fastcheck_big.log`, `fastcheck_FG_E78.log` |
| Remark 6.2: the same check for A_n (n ≤ 9) and B_n, C_n, D_n (n ≤ 8) | `./fastcheck A 9`, `./fastcheck B 8`, ... | "without potential 0" (`fastcheck_big.log`) |
| Squares span the cycle space and have parallel opposite edges (exact rank mod 2^31−1) | `./squares F 4` (and similarly A–G) | "SQUARES GENERATE", "non-parallel squares total 0" (`squares_exact.log`) |
| Cover rule (Corollary 4.2, Lemma 4.4) and the Ω_c sign rule (Lemma 4.3) | `./omega B 4`, `./omega E 8`, `./omega F 4`, ... | "cover-rule violations 0", "pairs beta<beta' with Omega<0: 0" |
| H3 root poset = lex-greater representatives of Cuntz–Stump Fig. 2 (32 order ideals) | `python3 h3cs.py` | "OK: lex-greater representatives reproduce h3_poset.json" |
| H3: K_c is conjugate to Row for all words, and has order 10 at an exact rational point (so Row² ≠ id) | `python3 h3run.py` | "orientation ROW K_c order: {10}" for every word (about 1 minute) |
| H3: randomized identity test modulo four large primes | `python3 modcheck.py 10` | "all 25 random points fixed: 25" for each of the four primes |
| I2(m): order lcm(2, m) | `python3 i2m.py` | rowmotion order equals lcm(2, m) for m = 3, ..., 14 |
| Floating-point sanity check of rowmotion orders (Section 8): A_n (n ≤ 5), B_n (n ≤ 6) | `python3 roworder.py` | order 2(n+1) in type A_n, 2n in type B_n |

Data: `h3_poset.json` holds the 15 D6 roots of the H3 root poset and its covers. `h3_map.json` holds the labelling by H3 roots (Table 1 of the paper; it indexes `roots.positive_roots(cartan('H3',3))`) and the covers.

## `h3_certificate/` — deterministic H3 certificate
`python3 h3_sym_half.py` computes Row^5(x) and Row^{-5}(x) as reduced fractions in Z[x_1..x_15] (FLINT, gcd-normalised, boundary labels 1). It checks that they are identical, so Row^10 = id as rational maps, and that Row^5 ≠ id. This takes about 12 minutes; the run log is `h3_sym_half.log`. The poset is read from `h3_iso.json` (key `('lexgreater', 'N^6 (D6 order)')`), the same poset as `code/h3_poset.json`.
