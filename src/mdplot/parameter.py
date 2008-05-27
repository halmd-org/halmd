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


"""
Determine parameter difference of multiple parameter sets
"""
def difference(sets):
    # fully-qualified names of attributes with first value
    visited = {}
    # fully-qualified names of differing attributes
    diff = {}
    # fully-qualified names with value for each parameter set
    attrs = []

    # iterate over parameter sets
    for parameters in sets:
        pairs = []
        # iterate over all nodes in parameters group
        for node in parameters._f_walkNodes():
            # iterate over attribute names in node
            for attr in node._v_attrs._f_list():
                # fully-qualified attribute name
                path = node._v_name + "/" + attr
                value = node._v_attrs.__getattr__(attr)
                pairs.append((path, value))

                if not path in visited:
                    visited[path] = value
                elif not visited[path] == value:
                    diff[path] = True

        attrs.append(pairs)

    # plot labels
    labels = []
    for pairs in attrs:
        l = []
        for (path, value) in pairs:
            if path in diff:
                try:
                    # parameter with doubleing-point value
                    l.append("%s = %.3G" % (label(path), value))
                except TypeError:
                    # parameter with arbitrary value
                    l.append(" = ".join((path, str(value))))

        labels.append(", ".join(l))

    return labels


"""
Returns values of a parameter from multiple parameter sets
"""
def values(sets, path):
    values = []
    nodes = path.split('/')
    attr = nodes.pop()

    for set in sets:
        for node in nodes:
            if not set.__contains__(node):
                return values
            set = set._v_children[node]

        if attr in set._v_attrs:
            values.append(set._v_attrs.__getattr__(attr))

    return values


"""
Translate fully-qualified parameter name to plot label
"""
def label(path):
    # plot labels for fully-qualified parameter names
    labels = {
            'autocorrelation/steps': 'steps',
            'autocorrelation/block_size': 'block-size',
            'autocorrelation/block_shift': 'block-shift',
            'autocorrelation/block_count': 'block-count',
            'autocorrelation/max_samples': 'max-samples',
            'ljfluid/dimension': 'dim',
            'ljfluid/particles': 'N',
            'ljfluid/blocks': 'CUDA blocks',
            'ljfluid/threads': 'CUDA threads',
            'ljfluid/density': '{/Symbol rho}^*',
            'ljfluid/box_length': '{/Symbol L}^*',
            'ljfluid/timestep': '{/Symbol dt}',
            'ljfluid/cutoff_distance': '{/Symbol x}_{cut}',
            'program/name': 'program',
            'program/version': 'version',
    }

    return path in labels and labels[path] or path


