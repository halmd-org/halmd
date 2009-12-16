#!/usr/bin/python
from pylab import *
import sys
if len(sys.argv) < 2:
    sys.exit("Usage: halmd.py [OUTPUT] ...")
fig = figure(figsize=(1, 1))
ax = axes((0, 0, 1, 1), frameon=False)
lj = lambda r: 4 * (pow(r, -12) - pow(r, -6))
x = linspace(0.969, 2.95, 1000)
ax.plot(x, lj(x), color="#00545c", lw=2, alpha=0.8)
ax.plot(x, lj(x), color="#00545c", lw=4, alpha=0.6)
ax.plot(x, lj(x), color="#00545c", lw=6, alpha=0.4)
ax.plot(x, lj(x), color="#00545c", lw=8, alpha=0.2)
ax.plot(x, lj(x), color="#00545c", lw=10, alpha=0.1)
H = array(((0, 0, 1), (2, 1, 1)))
A = array(((0,), (2,)))
L = array(((0, 0, 0), (0, 1, 2)))
M = array(((0, 0, 1), (0, 2, 2)))
D = array(((0, 1, 1), (2, 1, 2)))
ax.scatter(H[0] * 0.3 + 1.35, H[1] * 0.3 + 0.35, s=120, color="red", alpha=0.2)
ax.scatter(H[0] * 0.3 + 1.35, H[1] * 0.3 + 0.35, s=40, color="red", alpha=0.8)
ax.scatter(A[0] * 0.3 + 1.95, A[1] * 0.3 + 0.35, s=120, color="green", alpha=0.2)
ax.scatter(A[0] * 0.3 + 1.95, A[1] * 0.3 + 0.35, s=40, color="green", alpha=0.8)
ax.scatter(L[0] * 0.3 + 2.55, L[1] * 0.3 + 0.35, s=120, color="darkblue", alpha=0.2)
ax.scatter(L[0] * 0.3 + 2.55, L[1] * 0.3 + 0.35, s=40, color="darkblue", alpha=0.8)
ax.scatter(M[0] * 0.3 + 1.95, M[1] * 0.3 - 0.95, s=120, color="#00545c", alpha=0.2)
ax.scatter(M[0] * 0.3 + 1.95, M[1] * 0.3 - 0.95, s=40, color="#00545c", alpha=0.8)
ax.scatter(D[0] * 0.3 + 2.55, D[1] * 0.3 - 0.95, s=120, color="grey", alpha=0.2)
ax.scatter(D[0] * 0.3 + 2.55, D[1] * 0.3 - 0.95, s=40, color="grey", alpha=0.8)
ax.axis((0.8, 3.1, -1.15, 1.15))
for tick in ax.xaxis.get_major_ticks():
    tick.set_visible(False)
for tick in ax.yaxis.get_major_ticks():
    tick.set_visible(False)
for fn in sys.argv[1:]:
    savefig(fn, transparent=True)
