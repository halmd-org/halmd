#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2014 Felix Höfling
#
# This file is part of HALMD.
#
# HALMD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# define and parse command line arguments
import argparse
parser = argparse.ArgumentParser(description="Generate and plot reference data for modulation functions")
parser.add_argument("--kappa", type=float, required=True, help="inverse penetration depth")
parser.add_argument("--width", type=float, help="width of liquid slab")
parser.add_argument("--offset", type=float, default=0, help="interface position")
parser.add_argument("--samples", type=float, nargs='+', default=(-5, -1, 0, 0.1, 1, 9),
                    help="sample points along z-axis")
parser.add_argument("modulation", choices=("exponential", "exponential_slab", "catenary"),
                    help="name of modulation function")
args = parser.parse_args()

import sys
from pylab import *

def unit_step(x):
    return (sign(x) + 1) / 2

def cond(x, y):
    return (x - 1) * unit_step(y) + 1

def exponential(kappa, offset):
    return lambda x: cond(exp(-kappa * (x - offset)), kappa * (x - offset))

def exponential_slab(kappa, width, offset):
    return lambda x: cond(exp(-kappa * (x - offset)) * unit_step(kappa * (sign(kappa) * width - x + offset)), kappa * (x - offset))

def catenary(kappa, width, offset):
    return lambda x: cond(cosh(kappa * (x - offset)) / cosh(kappa * width / 2), width / 2 - abs(x - offset));

if args.modulation == "exponential":
    f = exponential(args.kappa, args.offset)
elif args.modulation == "exponential_slab":
    f = exponential_slab(args.kappa, args.width, args.offset) if args.width else None
elif args.modulation == "catenary":
    f = catenary(args.kappa, args.width, args.offset) if args.width else None

if not f:
    sys.exit("missing argument: width")

samples = array(args.samples)
for z in samples:
    print "({0:g}, {1:.15g})".format(z, f(z))

z = sort(concatenate((linspace(-10, 20, num=200), samples)), kind='mergesort')
plot(z, f(z), '-')
plot(samples, f(samples), 'o')

ylim(-.1, 1.1)
show()
