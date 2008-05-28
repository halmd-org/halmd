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
Plot thermodynamic equilibrium properties
"""
def plot(tep):
    # parameter sets
    sets = [root.parameters for (root, name) in tep]

    def spawn_gnuplot(title):
        # one gnuplot instance per plot to allow multiple *interactive* windows
        g = Gnuplot.Gnuplot()

        # With the current gnuplot release (4.2.3), the persist option for the
        # wxt terminal is buggy and causes CPU usage to rise to 100% after the
        # main gnuplot process exits. This bug is fixed in gnuplot 4.3 CVS
        # since 2007-04-29.
        g('set terminal wxt persist size 1200,900 enhanced title "ljfluid: %s"' % title)

        g('set key outside vertical center bottom Left reverse box')
        g('set xlabel "{/Symbol t}"')
        return g

    # plot command
    plot = '"%s" binary array=inf format="%%float" using (($0 + 1)*%f):%s title "%s" with lines'
    # plot titles
    titles = parameter.difference(sets)

    # mean potential energy per particle
    g = spawn_gnuplot('Mean potential energy per particle')
    g('set ylabel "{/Symbol e}_{pot}({/Symbol t})"')
    plots = []
    for (i, (root, name)) in enumerate(tep):
        f = file(name + '_epot.bin', 'wb')
        f.write(root.EPOT.read().tostring())
        f.close()
        plots.append(plot % (f.name, root.parameters.ljfluid._v_attrs.timestep, '1', titles[i]))
    g('plot ' + ', '.join(plots))

    # mean kinetic energy per particle
    g = spawn_gnuplot('Mean kinetic energy per particle')
    g('set ylabel "{/Symbol e}_{kin}({/Symbol t})"')
    plots = []
    for (i, (root, name)) in enumerate(tep):
        f = file(name + '_ekin.bin', 'wb')
        f.write(root.EKIN.read().tostring())
        f.close()
        plots.append(plot % (f.name, root.parameters.ljfluid._v_attrs.timestep, '1', titles[i]))
    g('plot ' + ', '.join(plots))

    # mean total energy per particle
    g = spawn_gnuplot('Mean total energy per particle')
    g('set ylabel "{/Symbol e}({/Symbol t})"')
    plots = []
    for (i, (root, name)) in enumerate(tep):
        f = file(name + '_etot.bin', 'wb')
        f.write(root.ETOT.read().tostring())
        f.close()
        plots.append(plot % (f.name, root.parameters.ljfluid._v_attrs.timestep, '1', titles[i]))
    g('plot ' + ', '.join(plots))

    # temperature
    g = spawn_gnuplot('Temperature')
    g('set ylabel "T^*({/Symbol t})"')
    plots = []
    for (i, (root, name)) in enumerate(tep):
        f = file(name + '_temp.bin', 'wb')
        f.write(root.TEMP.read().tostring())
        f.close()
        plots.append(plot % (f.name, root.parameters.ljfluid._v_attrs.timestep, '1', titles[i]))
    g('plot ' + ', '.join(plots))

    # pressure
    g = spawn_gnuplot('Pressure')
    g('set ylabel "P^*({/Symbol t})"')
    plots = []
    for (i, (root, name)) in enumerate(tep):
        f = file(name + '_press.bin', 'wb')
        f.write(root.PRESS.read().tostring())
        f.close()
        plots.append(plot % (f.name, root.parameters.ljfluid._v_attrs.timestep, '1', titles[i]))
    g('plot ' + ', '.join(plots))

    # velocity center of mass
    g = spawn_gnuplot('Velocity center of mass')
    g('set ylabel "({/Symbol n}^2({/Symbol t}))^{1/2}"')
    plots = []
    for (i, (root, name)) in enumerate(tep):
        f = file(name + '_vcm.bin', 'wb')
        f.write(root.VCM.read().tostring())
        f.close()
        plots.append(plot % (f.name, root.parameters.ljfluid._v_attrs.timestep, '(sqrt($1*$1+$2*$2))', titles[i]))
    g('plot ' + ', '.join(plots))


