#!/usr/bin/env python

# use arbitrary precision floating-point maths with 18 digits
import mpmath as mp
mp.dps=18

# set potential parameters
eps    = mp.mpf(.25)
sigma  = mp.mpf(2)
rmin   = mp.mpf(1.5)
B      = mp.mpf(1.41)
cutoff = mp.mpf(5) * sigma

# define truncated and shifted Morse potential U(r)
potential = lambda r : eps / (2 * B * B - 1) * (mp.exp(-2 * (r - rmin) / sigma / B) - 2 * B * B * mp.exp(-(r - rmin) / sigma / B))
potential_trunc = lambda r : potential(r) - potential(cutoff)

# fval = -U'(r) / r
fval = lambda r : 2 * eps / (2 * B * B - 1) * (mp.exp(-2 * (r - rmin) / sigma / B) / B - B * mp.exp(-(r - rmin) / sigma / B)) / sigma / r

# output is suitable to be pasted into the .cpp file of the potential unit test
for r in (.25, .5, .75, 1., 2., 5., 10.):
    print("      , {{{{{0}, {1}, {2}}}}}".format(r, mp.nstr(fval(r), 16), mp.nstr(potential_trunc(r), 16)))

