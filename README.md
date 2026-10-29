# Toggling in w₀(c) root order is conjugate to rowmotion

**Resolution of AIM Conjecture 3.2 and Thomas's conjecture**

Samson Fraser · samsonmfraser@gmail.com

📄 **Paper:** [`toggle-w0c-conjugacy.pdf`](toggle-w0c-conjugacy.pdf) (version 1.2, 29 October 2026, 15 pages)

🔗 **Zenodo:** [doi:10.5281/zenodo.23602944](https://doi.org/10.5281/zenodo.23602944) (always points to the newest version)

---

## The problem in brief

Take a partially ordered set and a subset closed downward. *Toggling* an element adds it to the subset or removes it, whenever the result is still closed downward. Toggling every element once, from the top down, gives **rowmotion**, a map studied since the 1970s.

Toggles also act on numerical labels, where each toggle becomes a rational map. In that *birational* setting it is hard to tell whether a product of toggles ever returns to where it started. Even for rowmotion this is known only for special posets.

Two open problems ask about one particular product of toggles. Fix a Coxeter element *c*, and toggle the positive roots of a root system in the order given by the *c*-sorting word of the longest element w₀. Call the result **K_c** (Thomas's *Krew_c* in the combinatorial setting).

- **Thomas's conjecture** (2013; Conjecture 4.4.5 in N. Williams's thesis *Cataland*): for every finite Weyl group, Krew_c is conjugate to rowmotion in the toggle group.
- **AIM Conjecture 3.2** (problem list of the 2015 American Institute of Mathematics workshop on Dynamical Algebraic Combinatorics): birational K_c has order 2h in types A, B, H₃ and I₂(m), where h is the Coxeter number.

## What the paper proves

1. **Thomas's conjecture holds in every finite crystallographic type**, for combinatorial, piecewise-linear and birational toggles. The same conjugacy holds in types I₂(m) and H₃ (for an explicit labelling of the H₃ root poset).
2. **AIM Conjecture 3.2 is resolved.** K_c returns to the identity after 2h steps in all four types. The literal claim "order 2h" is false in some of them:

   | Type | Exact order of K_c |
   |---|---|
   | Aₙ (n ≥ 2) | 2h = 2(n+1) |
   | Bₙ | h = 2n |
   | H₃ | h = 10 |
   | I₂(m) | lcm(2, m): 2h for m odd, h for m even |

3. **A consequence for noncrossing partitions** (Corollary 1.4): in every finite Weyl group, Krew_c has the same orbit structure as the Kreweras complement on noncrossing partitions, and an explicit bijection carries one to the other. In type A this is close to the result of Dequêne, Frieden, Iraci, Schreier-Aigner, Thomas and Williams, whose bijection is more direct.

## How the proof works

1. **A toggle order is a set of arrows.** Toggles of elements that are not neighbours commute. So a toggle product depends only on which end of each edge of the Hasse diagram is toggled first. Draw that as an arrow on each edge.
2. **Rotating the order is conjugation.** Moving the first toggle to the end conjugates the product by that toggle. In the arrow picture, it turns a source into a sink.
3. **Pretzel's theorem** (1986) says when two arrow patterns are related by such moves: when they have the same circulation (arrows forward minus arrows backward) around every cycle. Rowmotion's pattern has circulation zero, so it suffices to show the same for K_c.
4. **Which way each arrow points.** Dyer's reflection orders and a sign rule of Reading and Speyer decide the direction of every edge from a skew-symmetric form Ω_c.
5. **Small squares are enough.** Opposite edges of each small square in the root poset point the same way, so every square has circulation zero. The squares generate all cycles: in simply laced types because any two upper covers of a root have a common upper cover, in types B and C through a planar grid model, and in F₄ and G₂ by short direct arguments.

The exact orders then follow from known periodicity results for birational rowmotion (Grinberg–Roby in type A, Johnson–Liu in type B) and from direct computation in types I₂(m) and H₃.

## What uses a computer

**For crystallographic types the proof uses no computer.** Only type H₃ is computational. Its root poset is not built from simple-root steps, so the conjugacy is checked on a cycle basis, and Row¹⁰ = id is certified by an exact symbolic computation.

Everything else in the code folder is an independent check of results that are proved by hand.

## Reproducing the computations

Requirements: Python 3.11 with `sympy`, `networkx` and `python-flint`, and a C++17 compiler.

```bash
cd anc/code
python3 diamond.py                                # squares span the cycles, all types: ends "ALL OK"
g++ -O2 -o fastcheck fastcheck.cpp && ./fastcheck E 8   # a potential for all 40320 orderings, about 15 s
python3 h3run.py                                  # H3 conjugacy, about 1 minute
cd ../h3_certificate && python3 h3_sym_half.py    # H3 certificate Row^10 = id, about 12 minutes
```

[`anc/README.md`](anc/README.md) lists every computational claim in the paper with its command and expected output. Run logs are included next to the scripts.

To build the paper: `pdflatex main.tex` (run it twice). It needs no `.bib` file.

## Repository layout

| Path | Contents |
|---|---|
| `toggle-w0c-conjugacy.pdf` | The paper |
| `main.tex` | LaTeX source |
| `anc/README.md` | Map from each computational claim to the command that checks it |
| `anc/code/` | Python and C++ checks, with run logs |
| `anc/h3_certificate/` | Exact symbolic certificate for type H₃ |
| `CITATION.cff` | Citation metadata (GitHub's "Cite this repository" button) |

## Version history

- **1.2** (29 October 2026)
  - The proof for crystallographic types no longer uses a computer. Types A, D and E use a new argument that the squares span all cycles (Lemmas 5.2–5.3), and F₄ is a short hand check. The earlier computer searches are kept as independent checks (Remark 6.2).
  - New Corollary 1.4 on noncrossing partitions.
  - New Table 1 showing where each type is proved.
  - Notation fixes and acknowledgements.
  - These changes follow comments from Hugh Thomas.
- **1.1** (23 October 2026): code and certificates added to the Zenodo record.
- **1.0** (23 October 2026): first public version.

## How to cite

```bibtex
@misc{fraser2026toggling,
  author    = {Fraser, Samson},
  title     = {Toggling in $w_0(c)$ root order is conjugate to rowmotion:
               resolution of {AIM} {C}onjecture 3.2 and {T}homas's conjecture},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23602944},
  note      = {Preprint}
}
```

## License

[CC BY 4.0](LICENSE), for the paper and the code.
