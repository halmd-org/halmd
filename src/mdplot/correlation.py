#!/usr/bin/python
#
# mdplot - Molecular Dynamics simulation plotter
#
# Copyright (C) 2008  Peter Colberg
#
# This program is free software: you can redistribute it and/or modify
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

import Gnuplot
import glob
import os, os.path, sys
import subprocess
import tables
import tempfile
import time

import parameter


"""
Plot correlations
"""
def plot(tcf):
    # parameter sets
    sets = [root.parameters for (root, name) in tcf]

    def spawn_gnuplot(title):
        # one gnuplot instance per plot to allow multiple *interactive* windows
        # FIXME write commands to file so we can tinker with it if necessary
        g = Gnuplot.Gnuplot()

        # With the current gnuplot release (4.2.3), the persist option for the
        # wxt terminal is buggy and causes CPU usage to rise to 100% after the
        # main gnuplot process exits. This bug is fixed in gnuplot 4.3 CVS
        # since 2007-04-29.
        g('set terminal wxt persist size 1024,768 enhanced title "ljfluid: %s"' % title)

        g('set key outside vertical center bottom Left reverse box')
        # minimal timestep value of sets
        # g('set xrange [%.3G:]' % max(parameter.values(sets, 'ljfluid/timestep')))
        g('set logscale x')
        g('set format x "10^{%T}"')
        g('set xlabel "{/Symbol Dt}"')
        return g

    # plot commands
    plot_line = '"%s" binary array=inf format="%%float%%float%%*float" using 1:%s notitle with lines lt %d'
    plot_error = '"%s" binary array=inf format="%%float%%float%%float" using 1:%s:%s title "%s" with yerrorbars lt %d'
    # plot titles
    titles = parameter.difference(sets)

    def time_ordered_write(f, dataset):
        data = {}
        for (n, block) in enumerate(dataset):
            for sample in block:
                data[(sample[0], n)] = sample.tostring()

        # sort dataset after time interval
        for k in sorted(data.keys()):
            f.write(data[k])

    # write mean squared displacement data files
    files_MSD = []
    for (root, name) in tcf:
        f = file(name + '_msd.bin', 'wb')
        time_ordered_write(f, root.MSD)
        f.close()
        files_MSD.append(f.name)

    g = spawn_gnuplot('Mean squared displacement')
    g('set logscale y')
    g('set format y "10^{%T}"')
    g('set ylabel "<({/Symbol x}({/Symbol t}+{/Symbol Dt}) - {/Symbol x}({/Symbol t}))^2>_N"')
    plot = []
    for (i, f) in enumerate(files_MSD):
        plot.append(plot_line % (f, '2', i + 1))
        plot.append(plot_error % (f, '2', '3', titles[i], i + 1))
    g('plot ' + ', '.join(plot))

    g = spawn_gnuplot('Diffusion constant')
    g('unset logscale y')
    g('set format y')
    g('set ylabel "<({/Symbol x}({/Symbol t}+{/Symbol Dt}) - {/Symbol x}({/Symbol t}))^2>_N / {/Symbol Dt}"')
    plot = []
    for (i, f) in enumerate(files_MSD):
        plot.append(plot_line % (f, '($2/$1)', i + 1))
        plot.append(plot_error % (f, '($2/$1)', '($3/$1)', titles[i], i + 1))
    g('plot ' + ', '.join(plot))

    # write mean quartic displacement data files
    files_MQD = []
    for (root, name) in tcf:
        f = file(name + '_mqd.bin', 'wb')
        time_ordered_write(f, root.MQD)
        f.close()
        files_MQD.append(f.name)

    g = spawn_gnuplot('Mean quartic displacement')
    g('set logscale y')
    g('set format y "10^{%T}"')
    g('set ylabel "<({/Symbol x}({/Symbol t}+{/Symbol Dt}) - {/Symbol x}({/Symbol t}))^4>_N"')
    plot = []
    for (i, f) in enumerate(files_MQD):
        plot.append(plot_line % (f, '2', i + 1))
        plot.append(plot_error % (f, '2', '3', titles[i], i + 1))
    g('plot ' + ', '.join(plot))

    # write velocity autocorrelation data files
    files_VAC = []
    for (root, name) in tcf:
        f = file(name + '_vac.bin', 'wb')
        time_ordered_write(f, root.VAC)
        f.close()
        files_VAC.append(f.name)

    g = spawn_gnuplot('Velocity autocorrelation')
    g('set logscale y')
    g('set format y "10^{%T}"')
    g('set ylabel "<({/Symbol n}({/Symbol t}+{/Symbol Dt}) {/Symbol n}({/Symbol t}))>_N"')
    plot = []
    for (i, f) in enumerate(files_VAC):
        plot.append(plot_line % (f, '2', i + 1))
        plot.append(plot_error % (f, '2', '3', titles[i], i + 1))
    g('plot ' + ', '.join(plot))


