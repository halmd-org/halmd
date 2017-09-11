#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2014 Felix Höfling
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

#
# Particle trajectories can be visualised with VMD (Visual Molecular Dynamics)
# [1] taking advantage of a plugin [2] (under development) for the H5MD output
# files. The plugin is still a bit restrictive, the following patch routine
# tries to meet these restrictions.
#
# [1] http://www.ks.uiuc.edu/Research/vmd
# [2] https://github.com/h5md/VMD-h5mdplugin
#

import argparse
parser = argparse.ArgumentParser(description="patch H5MD output file to work with VMD plugin")
parser.add_argument('input', nargs=1, metavar='INPUT', help='H5MD input file')
parser.add_argument('output', nargs=1, metavar='OUTPUT', help='H5MD output file')
args = parser.parse_args()

# copy input to destination file
from shutil import copyfile
copyfile(args.input[0], args.output[0])

import h5py
import numpy

f = h5py.File(args.output[0], 'r+')

if not "particles" in f:
    SystemExit("No particle trajectory in input file")

# iterate over particle groups
s = None
for p in f["particles"].itervalues():
    # make species time-independent
    if "species" in p.keys():
        species = p["species/value"][0]
        del p["species"]
        p["species"] = species
        s = numpy.append(s, species) if s != None else species

    # add image data and fold back particle positions to periodic box
    if "image" not in p.keys() and "position" in p.keys():
        image = p.create_group("image")
        pos = p["position"]
        image["time"] = pos["time"]     # create hard links
        image["step"] = pos["step"]

        # blow up 2D data to 3 dimensions
        dimension = p["box"].attrs["dimension"]
        length = numpy.ones((3,))
        length[:dimension] = numpy.diagonal(p["box/edges"])

        r = numpy.zeros(pos["value"].shape[:-1] + (3,))
        r[..., :dimension] = pos["value"]

        img = numpy.round(r / length)
        image["value"] = img

        del pos["value"]
        pos["value"] = r - img * length

# VMD data in parameters group
vmd = f.require_group("parameters/vmd_structure")
vmd["indexOfSpecies"] = numpy.unique(s) if type(s) != type(None) else [0,]
nspecies = len(vmd["indexOfSpecies"])
vmd["name"] = numpy.array(("He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"))[:nspecies]

f.close()

