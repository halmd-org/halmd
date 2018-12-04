#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Jake Atwell
# Copyright © 2018 Felix Höfling
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

def main(args):
    import h5py
    import numpy as np
    from scipy.optimize import curve_fit

    # open and read H5MD input file
    with h5py.File(args.input, 'r') as H5:
        msd = H5['dynamics/all/mean_square_displacement']

        # get relavant data for MSD calculation
        x = np.array(msd['time']).flatten()
        y = np.array(msd['value']).flatten()
        yerr = np.array(msd['error']).flatten()

    # bring data in chronological order
    idx = np.argsort(x)
    x, y, yerr = x[idx], y[idx], yerr[idx]

    # partition data set for fit using a selection mask
    mask = np.ones_like(x, dtype=bool)

    # deselect according to --range argument
    if args.range:
        mask[:args.range[0]] = False
        mask[args.range[1]:] = False

    # deselect according to --min-time argument, exclude t = 0 data
    min_time = args.min_time or 0
    mask[np.where(x <= min_time)] = False

    # split input data
    fit_x, no_fit_x = x[mask], x[~mask]
    fit_y, no_fit_y = y[mask], y[~mask]
    fit_yerr, no_fit_yerr = yerr[mask], yerr[~mask]

    # fit model
    linear_model = lambda x, D : 2 * args.dimension * D * x

    popt, pcov = curve_fit(linear_model, fit_x, fit_y, sigma=fit_yerr)
    perr = np.sqrt(np.diag(pcov))

    print("D = {0:.5g} ± {1:.2g}".format(popt[0], perr[0]))

    # generate plot
    if not args.no_plot:
        import matplotlib.pyplot as plt

        # time grid for reference lines
        t = np.logspace(np.log10(fit_x[0]) - 1, np.log10(fit_x[-1]) + .3)

        # plot MSD versus time
        if not args.rectify:
            plt.errorbar(fit_x, fit_y, yerr=fit_yerr, fmt='xb', label="data used in fit")
            plt.errorbar(no_fit_x, no_fit_y, yerr=no_fit_yerr, fmt='xg', label="data not used in fit")
            plt.loglog(t, linear_model(t, *popt), '-k')
            ylabel = "Mean-square displacement"
        else:
            plt.errorbar(fit_x, fit_y / fit_x, yerr=fit_yerr / fit_x, fmt='xb', label="data used in fit")
            plt.errorbar(no_fit_x, no_fit_y / no_fit_x, yerr=no_fit_yerr / no_fit_x, fmt='xg', label="data not used in fit")
            plt.semilogx(t, linear_model(t, *popt) / t, '-k')
            ylabel = "MSD(t) / t"

        # add axes labels and finalise plot
        plt.axis('tight')
        plt.xlabel('Time')
        plt.ylabel(ylabel)
        plt.legend(loc='best', frameon=False)
        plt.show()


def parse_args():
    import argparse

    # define and parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--dimension', type=int, default=3, help='space dimension')
    parser.add_argument('--no-plot', action='store_true', help='do not produce plots, but do the analysis')
    parser.add_argument('--rectify', action='store_true', help='rectify plot by showing MSD(t) / t')
    parser.add_argument('--range', type=int, nargs=2, help='range of data points to include in fit')
    parser.add_argument('--min-time', type=float, nargs=1, help='left end of time interval used for fitting')
    parser.add_argument('input', metavar='INPUT', help='H5MD input file with MSD data')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)

