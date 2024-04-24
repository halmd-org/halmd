#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2011  Felix Höfling
#
# This file is part of HALMD.
#
# HALMD is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
#

import argparse
import h5py
from numpy import *
from pylab import *

def main():
    # define and parse command line arguments
    parser = argparse.ArgumentParser(prog='plot_h5md.py')
    parser.add_argument('--range', type=int, nargs=2, help='select range of data points')
    parser.add_argument('--dump', metavar='FILENAME', help='dump plot data to filename')
    parser.add_argument('--no-plot', action='store_true', help='do not produce plots, but do the analysis')
    parser.add_argument('--group', help='particle group (default: %(default)s)', default='all')
    parser.add_argument('input', metavar='INPUT', help='H5MD input file with data for state variables')
    args = parser.parse_args()

    # evaluate option --range
    range = args.range or [0, -1]

    # open and read data file
    H5 = h5py.File(args.input, 'r')
    H5obs = H5['observables']
    H5particle = H5['particles'][args.group]
    H5box = H5particle['box']

    # print some details
    print('Particles: {0}'.format(H5obs['particle_number'][()]))
    print('Box size:', ' '.join(['{0:g}'.format(x) for x in diag(H5box['edges'])]))
    print('Density: {0:g}'.format(H5obs['density'][()]))

    # compute and print some averages
    # the simplest selection of a data set looks like this:
    #     temp, temp_err = compute_average(H5obs['temperature'], 'Temperature')
    # full support for slicing (the second pair of brackets) requires conversion to a NumPy array before
    temp, temp_err = compute_average(array(H5obs['temperature/value'])[range[0]:range[1]], 'Temperature')
    pressure, pressure_err = compute_average(array(H5obs['pressure/value'])[range[0]:range[1]], 'Pressure')
    epot, epot_err = compute_average(array(H5obs['potential_energy/value'])[range[0]:range[1]], 'Potential energy')

    # select data for plotting the potential energy as function of time
    x = array(H5obs['potential_energy/time'])[range[0]:range[1]]
    y = array(H5obs['potential_energy/value'])[range[0]:range[1]]

    # append plot data to file
    if args.dump:
        f = open(args.dump, 'a')
        print('# time   E_pot(t)', file=f)
        savetxt(f, array((x, y)).T)
        print('\n', file=f)
        f.close()

    # generate plot
    if not args.no_plot:
        # plot potential energy versus time
        plot(x, y, '-b', label=args.input)

        # plot mean value for comparison
        x = linspace(min(x), max(x), num=50)
        y = zeros_like(x) + epot
        plot(x, y, ':k')

        # add axes labels and finalise plot
        axis('tight')
        xlabel(r'Time $t$')
        ylabel(r'Potential energy $E_\mathrm{pot}$')
        legend(loc='best', frameon=False)
        show()

    H5.close()

def compute_average(data, label, nblocks = 10):
    """ compute and print average of data set
        The error is estimated from grouping the data in shorter blocks """

    # group data in blocks, discard excess data
    data = reshape(data[:(nblocks * (data.shape[0] / nblocks))], (nblocks, -1))

    # compute mean and error
    avg = mean(data)
    err = std(mean(data, axis=1)) / sqrt(nblocks - 1)
    print('{0:s}: {1:.4f} ± {2:.2g}'.format(label, avg, err))

    # return as tuple
    return avg, err

if __name__ == '__main__':
    main()
