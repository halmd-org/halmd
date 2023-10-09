#!/usr/bin/env python3

# use arbitrary precision floating-point maths with 18 digits
import mpmath as mp
mp.dps=18

# set potential parameters
eps    = mp.mpf(.5)
sigma  = mp.mpf(.75)
rmin   = mp.mpf(2)
B      = mp.mpf(1.5)
cutoff = mp.mpf(8) * sigma

# define truncated and shifted Morse potential U(r)
potential = lambda r : eps / (2 * B * B - 1) * (mp.exp(-2 * (r - rmin) / sigma * B) - 2 * B * B * mp.exp(-(r - rmin) / sigma / B))
potential_trunc = lambda r : potential(r) - potential(cutoff)

# fval = -U'(r) / r
fval = lambda r : eps / (B - .5 / B) * (mp.exp(-2 * (r - rmin) / sigma * B) - mp.exp(-(r - rmin) / sigma / B)) / sigma / r

# output is suitable to be pasted into the .cpp file of the potential unit test
for r in (.25, .5, .75, 1., 2., 3., 5., 6.):
    print("      , {{{{{0}, {1}, {2}}}}}".format(r, mp.nstr(fval(r), 16), mp.nstr(potential_trunc(r), 16)))

    # verify force calculation by doing the derivative numerically
    dr = 1e-6
    tolerance = 1e-4
    assert(mp.fabs(fval(r) - (-(potential(r + dr) - potential(r)) / dr / r)) < tolerance * max(mp.fabs(fval(r)), eps / sigma**2))

# generate plot of the potential U(r)
from numpy import linspace
from matplotlib.pyplot import axhline, plot, show, ylim

x = linspace(0, float(cutoff), num=100)
y = [potential_trunc(r) for r in x]
plot(x, y)
axhline(y=0, color='k', lw=0.5)
ylim(-1.5 * float(eps), 5 * float(eps))
show(block=True)
