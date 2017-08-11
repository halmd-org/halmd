#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import h5py
from numpy import *
from pylab import *
from scipy.optimize import curve_fit

def main():
    # define and parse command line arguments
    parser = argparse.ArgumentParser(prog='plot_h5md.py')
    parser.add_argument('--range', type=int, nargs=2, help='select range of data points to include in fit')
    parser.add_argument('--start_t', type=float, nargs=1, help='select time value at which fitting should begin')
    parser.add_argument('--dump', metavar='FILENAME', help='dump plot data to filename')
    parser.add_argument('--no-plot', action='store_true', help='do not produce plots, but do the analysis')
    parser.add_argument('--group', help='particle group (default: %(default)s)', default='all')
    parser.add_argument('input', metavar='INPUT', help='H5MD input file with data for state variables')
    args = parser.parse_args()

    # evaluate option --range
    range = args.range or [0, -1] #could also use np.where()
    start_t = args.start_t or 0

    # open and read data file
    H5 = h5py.File(args.input, 'r')

    # obtain relavant data for MSD calculation
    msd = H5['dynamics/all/mean_square_displacement']
    x = array(msd['time']).flatten()
    y = array(msd['value']).flatten()
    dy = array(msd['error']).flatten()

    # organize data into a single array so it can be sorted
    data = array((x, y, dy)).T 
    sort_column_index = 0 # index of colum by which data should be sorted (0-sorting by time, 1-sorting by value)
    x,y,dy = data[argsort(data[:,sort_column_index])].T
    
    full_x = x

    # data selection if --range argument passed
    if args.range: 
        data_divisions = [0, range[0], range[1], -1]

        # data used in the fit
        fit_x = x[data_divisions[1]:data_divisions[2]]
        fit_y = y[data_divisions[1]:data_divisions[2]]
        fit_dy = dy[data_divisions[1]:data_divisions[2]]

        # data not used  in the fit
        no_fit_x = append(x[data_divisions[0]:data_divisions[1]], x[data_divisions[2]:data_divisions[3]])
        no_fit_y = append(y[data_divisions[0]:data_divisions[1]], y[data_divisions[2]:data_divisions[3]])
        no_fit_dy = append(dy[data_divisions[0]:data_divisions[1]], dy[data_divisions[2]:data_divisions[3]])

    # data selection if --start_t argument passed
    if start_t: 
        # data used in the fit
        fit_x = x[where(x > start_t)]
        fit_y = y[where(x > start_t)]
        fit_dy = dy[where(x > start_t)]

        # data not used  in the fit
        no_fit_x = x[where(x <= start_t)]
        no_fit_y = y[where(x <= start_t)]
        no_fit_dy = dy[where(x <= start_t)]


    popt, pcov = curve_fit(lin_fit_no_offset, fit_x, fit_y, sigma=fit_dy)
    perr = sqrt(diag(pcov))

    print "D = {0:.3f} ± {1:.3f}".format(popt[0], perr[0])

    # append plot data to file
    if args.dump:
        f = open(args.dump, 'a')
        print >>f, '# time   mean_square_disp(t)'
        savetxt(f, array((x, y)).T)
        print >>f, '\n'
        f.close()
    
    # generate plot
    if not args.no_plot:
        # plot potential energy versus time
        errorbar(fit_x, fit_y, fmt='xb', yerr=fit_dy, label="data used in fit")
        errorbar(no_fit_x, no_fit_y, fmt='xg', yerr=no_fit_dy, label="data not used in fit")
        plot(full_x, lin_fit_no_offset(full_x, *popt), ':r')
        xscale('log')
        yscale('log')

        # add axes labels and finalise plot
        axis('tight')
        xlabel(r'Time')
        ylabel(r'Mean Square  Displacement')
        legend(loc='best', frameon=False)
        show()



    H5.close()

def lin_fit_no_offset(x, a):
    return a*x

if __name__ == '__main__':
    main()

