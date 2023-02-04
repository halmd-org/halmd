#!/usr/bin/env python

# use arbitrary precision floating-point maths with 18 digits
import mpmath as mp
mp.dps=18

# set potential parameters
eps    = mp.mpf(.5)
sigma  = mp.mpf(2)
rmin   = mp.mpf(.75)
cutoff = mp.mpf(5) * sigma

# define truncated and shifted Morse potential U(r)
potential = lambda r : eps * ((1 - mp.exp(-r/sigma + rmin))**2 - 1)
potential_trunc = lambda r : potential(r) - potential(cutoff)

# fval = -U'(r) / r
fval = lambda r : -2 * eps * (1 - mp.exp(-r/sigma + rmin)) * (mp.exp(-r/sigma + rmin) / sigma) / r

# output is suitable to be pasted into the .cpp file of the potential unit test
for r in (.25, .5, .75, 1., 2., 5., 10.):
    print("      , {{{{{0}, {1}, {2}}}}}".format(r, mp.nstr(fval(r), 16), mp.nstr(potential_trunc(r), 16)))

