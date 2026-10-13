#!/usr/bin/env python3
"""
num_spurious.py -- position of the spurious (round-off / pseudospectral)
branch of the double-precision Chebyshev collocation matrices of
num_spectra.py, as a function of N.  For each case the right edge of the set
of computed eigenvalues that are NOT within 1e-3 of a closed-form eigenvalue
is printed, together with K_N = (d-1)/2 - s - edge  (s -> -1 for (u_t)^2).
Runtime: a few seconds.  python3 code/num_spurious.py
"""
import os
import sys
from fractions import Fraction as Fr

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from num_spectra import opmat, formula  # noqa: E402

CASES = [("MEMS q=2, d=1", 1, 0, Fr(-2, 3), Fr(4, 9), False),
         ("MEMS q=2, d=3", 3, 0, Fr(-2, 3), Fr(4, 9), False),
         ("u^3, d=3", 3, 0, Fr(1), Fr(6), False),
         ("e^u, d=1", 1, 0, Fr(0), Fr(2), False),
         ("(u_t)^2, d=1", 1, 0, Fr(0), Fr(0), True),
         ("(u_t)^2, d=3", 3, 0, Fr(0), Fr(0), True)]
NS = (20, 40, 80, 120)

if __name__ == "__main__":
    print("right edge of the spurious set (l = 0), double precision; K_N = (d-1)/2 - s - edge")
    for (name, d, l, s, V0, der) in CASES:
        fl = [float(x) for x in formula(s, V0, der, l, remin=-200)]
        edges = []
        for N in NS:
            A, _ = opmat(N, d, l, float(s), float(V0), der)
            ev = np.linalg.eigvals(A)
            spur = [e for e in ev if min(abs(e - f) for f in fl) > 1e-3]
            edges.append(max(e.real for e in spur) if spur else np.nan)
        thr = (d - 1) / 2 - (-1.0 if der else float(s))
        print("  %-14s (d-1)/2-s = %5.2f | edge N=%s: %s | K_N: %s"
              % (name, thr, "/".join(map(str, NS)),
                 " ".join("%6.2f" % e for e in edges),
                 " ".join("%5.2f" % (thr - e) for e in edges)))
