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
# Script for plotting 1D density profile from simulation data of wetting process
#

import argparse

import h5py
import matplotlib.pyplot as plt
import numpy as np

from plotting.density_profile import (compute_density_modes_1D, compute_1D_density)


def main():

    # define and parse command line arguments
    parser = argparse.ArgumentParser(prog='plot_density_profile.py')
    parser.add_argument('--range', type=int, nargs=2, help='select range of data points')
    parser.add_argument('--filter-width', type=float, help='define width of the Gaussian filter')
    parser.add_argument('--inputs', nargs='+', help='H5MD input files with data for state variables')
    parser.add_argument('--labels', nargs='+', default=None, help='labels for plots in a figure' )
    parser.add_argument('--dump', metavar='FILENAME', help='dump plot data to filename')
    args = parser.parse_args()

    if args.labels and len(args.labels) != len(args.input):
        raise ValueError("The number of labels is not equal to the number of input files")

    # open and read data file
    range = args.range or [159, 225] # range[0] is the position index of wall

    plt.figure()

    for i, input in enumerate(args.inputs):
        # open and read data file
        H5 = h5py.File(input, 'r')

        box_edges = np.diagonal(np.array(H5['particles/all/box/edges']))
        wavevector = np.array(H5['structure/all/density_mode/wavevector'])
        positions = np.array(H5['particles/all/position/value'])[::5,:,:] # choose positions every 5 time steps

        # print some details
        print('Particles: {0}'.format(H5['observables/particle_number'][()]))
        print('Box size:', end = ' ')
        for x in box_edges:
            print(' {0:g}'.format(x), end = ' ')
        print()
        print('Density: {0:g}'.format(H5['observables/density'][()]))

        # compute density modes
        density_modes = compute_density_modes_1D(wavevector, positions, coord_index=0)

        # compute density profile
        real_grid, density = compute_1D_density(density_modes, box_edges, wavevector, smoothing=True)

        # shift to wall
        real_grid = real_grid - real_grid[range[0]]

        # generate plot
        if args.labels:
            plt.plot(real_grid[range[0]:range[1]], density[range[0]:range[1]], label=args.labels[i])
        else:
            plt.plot(real_grid[range[0]:range[1]], density[range[0]:range[1]])

        H5.close()

    plt.legend()
    plt.ylabel(r"$\rho(x)/\rho_{bulk}$", fontsize=13)
    plt.xlabel(r"$x[\sigma]$", fontsize=13)
    plt.tight_layout()
    plt.savefig(f'{args.dump}.png')
    plt.show()


if __name__ == '__main__':
    main()