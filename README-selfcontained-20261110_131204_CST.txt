Subrings and cotypes of truncated polynomial rings
Sarthak Chimni and Ramin Takloo-Bighash
Self-contained revision — 10 November 2026

MAIN FILES
chimni-TB2-selfcontained-20261110_131204_CST.tex
chimni-TB2-selfcontained-20261110_131204_CST.pdf
chimni-TB2-selfcontained-review-20261110_131204_CST.txt

The PDF has 37 pages. The original three draft references are omitted;
the remaining bibliography contains third-party contextual sources.
The manuscript includes all proof-critical finite computations.

BUILD
Use a standard LaTeX installation with the packages listed in the preamble:
  latexmk -pdf -interaction=nonstopmode -halt-on-error chimni-TB2-selfcontained-20261110_131204_CST.tex
No external bibliography, figure, or source file is required.

REPRODUCE CHECKS
The following scripts use Python 3 and its standard library only:
  python3 printed-enumeration-20261110_131204_CST.py
  python3 verify_chimni_20261110_131204_CST.py
  python3 verify_symbolic_20261110_131204_CST.py

The first script reproduces the two consecutive listings printed in
Appendix A. The second independently checks finite coefficients by
truncated-polynomial products, triangular division, and linear lifting.
The third checks exact quartic rational polynomial identities and finite
comparisons. The recorded outputs give the precise tested ranges.
Passing finite tests does not replace the mathematical proofs.

SOURCE SHA256
3b9abc80be08a8da6d89d70731e68de2f695bb405284b3dc724e484b244fae91  chimni-TB2-selfcontained-20261110_131204_CST.tex
