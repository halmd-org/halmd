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
import numpy
import subprocess
import tables
import tempfile
import time


"""
Plot trajectories
"""
def plot(trj):
    for (root, name) in trj:
        render(root, name)


"""
Render trajectory movie
"""
def render(root, basename):
    # remove stray gnuplot output files
    for fn in glob.glob(basename + '_*.png'):
        os.unlink(fn)

    # font requires Debian package ttf-dejavu-core
    if not 'GDFONTPATH' in os.environ:
        os.environ['GDFONTPATH'] = '/usr/share/fonts/truetype/ttf-dejavu'

    g = Gnuplot.Gnuplot()

    output = 'set output "%s"'
    g('set terminal png font "DejaVuSans,12" enhanced size 1080, 1080')
    g('unset key')
    g('unset title')

    # simulation box length
    box = root.parameters.ljfluid._v_attrs.box_length
    # positional coordinates dimension
    dimension = root.parameters.ljfluid._v_attrs.dimension
    # simulation time resolution
    timestep = root.parameters.ljfluid._v_attrs.timestep

    g('set xrange [0:%f]' % box)
    g('set yrange [0:%f]' % box)

    if dimension == 3:
        # FIXME use fraction of cutoff length as particle radius
        plot = 'splot "%s" binary record=inf format="%%%s" using 1:2:3 notitle with points pt 7 ps 2'

        g('set zrange [-%f:%f]' % (box, box))
        g('set grid')
        # draw a complete box around splot
        g('set border 4095')

    else:
        # FIXME use fraction of cutoff length as particle radius
        plot = 'plot "%s" binary array=inf format="%%%s" using 1:2 notitle with points pt 7 ps 2'

        # draw grid for cell-list based MD simulation
        if root.parameters.ljfluid._v_attrs.__contains__("cell_length"):
            g('set format "%.3g"')
            g('set xtics 0, %f' % root.parameters.ljfluid._v_attrs.cell_length)
            g('set x2tics 0, %f' % root.parameters.ljfluid._v_attrs.cell_length)
            g('set ytics 0, %f' % root.parameters.ljfluid._v_attrs.cell_length)
            g('set y2tics 0, %f' % root.parameters.ljfluid._v_attrs.cell_length)
            g('set grid xtics ytics lt -1')
        else:
            g('set x2tics')
            g('set y2tics')

        g('set x2range [0:%f]' % box)
        g('set y2range [0:%f]' % box)

    data_file = None
    frame_fn = basename + '_%06d.png'
    sys.stdout.write('gnuplot: %6sf' % '')
    for (i, sample) in enumerate(root.trajectory.positions):
        # write trajectory sample to temporary binary file
        f = tempfile.NamedTemporaryFile('wb')
        tid = (sample.dtype == numpy.float64) and 'double' or 'float'
        f.write(sample.tostring())
        f.flush()
        # plot trajectory sample
        g(output % (frame_fn % i))
        g('set label 1 front "t^* = %g" font "DejaVuSans,14" at graph 0.85, 0.05' % (i * timestep))
        g(plot % (f.name, tid))

        # FIXME better way to wait for file creation (maybe pyinotify?)
        while not os.path.exists(frame_fn % i):
            time.sleep(1e-2)

        # erase previously printed frame number characters
        sys.stdout.write('\010 \010' * 7)
        sys.stdout.write('%6df' % (i + 1))
        sys.stdout.flush()

        # erase previous data file
        data_file = f

    # wait for gnuplot to finish
    del g
    # remove temporary files
    del data_file

    sys.stdout.write('\n')
    sys.stdout.flush()

    # render movie with ffmpeg2theora
    ffmpeg = subprocess.Popen(['ffmpeg2theora', frame_fn, '-x', '720', '-y', '720', '-S', '0', '-o', '%s.ogg' % basename, '--nosound'], close_fds=True)
    errno = ffmpeg.wait()

    # remove gnuplot output files
    for fn in glob.glob(basename + '_*.png'):
        os.unlink(fn)

    if errno:
        raise Exception('mencoder exited with error code %d' % errno)


