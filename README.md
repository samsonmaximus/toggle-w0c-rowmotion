# Toggling in w₀(c) root order is conjugate to rowmotion: resolution of AIM Conjecture 3.2 and Thomas's conjecture

**Samson Fraser** (samsonmfraser@gmail.com)

Preprint : [`toggle-w0c-conjugacy.pdf`](toggle-w0c-conjugacy.pdf)
Zenodo DOI: *[10.5281/zenodo.23603681](https://doi.org/10.5281/zenodo.23603681)*

## Abstract

Toggling is a simple operation on a partially ordered set: it adds an element to a downward-closed subset, or removes one, whenever the result is still downward closed. Toggling every element in turn, from the top down, gives rowmotion, a map studied since the 1970s. Toggles also act on numerical labels, where they become rational maps, and there it is hard to tell whether a product of toggles ever returns to where it started. Even for rowmotion, this is known only for special posets.

In 2015, the problem list of the American Institute of Mathematics workshop on Dynamical Algebraic Combinatorics conjectured that one such product has order 2h on the root posets of types A, B, H₃ and I₂(m), where h is the Coxeter number. This product toggles the positive roots in an order determined by a Coxeter element. We resolve the conjecture. We prove that the product is conjugate to rowmotion in every finite crystallographic type and in types H₃ (for an explicit labelling of its roots) and I₂(m); this also proves a 2013 conjecture of H. Thomas. In the four types of the conjecture, it follows that the product returns to the identity after 2h steps. Its exact order is 2h in types Aₙ (n ≥ 2) and I₂(m) with m odd, but only h in types Bₙ, H₃ and I₂(m) with m even. The exceptional types and H₃ rely on exact computer calculations.

## Contents

| Path                       | What it is                                                                                          |
| -------------------------- | --------------------------------------------------------------------------------------------------- |
| `toggle-w0c-conjugacy.pdf` | The paper                                                                                           |
| `main.tex`                 | LaTeX source (self-contained, no .bib file)                                                         |
| `anc/README.md`            | Every computational claim in the paper, with the command that reproduces it and the expected output |
| `anc/code/`                | Python and C++ checks, with run logs                                                                |
| `anc/h3_certificate/`      | Exact symbolic certificate that Row¹⁰ = id in type H₃                                               |

## Building and checking

Paper: `pdflatex main.tex` (run it twice).

Computations need Python 3.11 with `sympy`, `networkx` and `python-flint`, and a C++17 compiler. See [`anc/README.md`](anc/README.md) for the commands. Most checks take seconds; the H₃ certificate takes about 12 minutes.

## License

[CC BY 4.0](LICENSE) for the paper and the code.
