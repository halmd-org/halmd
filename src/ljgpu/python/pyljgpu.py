#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# pyljgpu - wrap ljgpu process as python class
#
# Copyright © 2008-2009  Peter Colberg
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

import subprocess

"""
Wrap ljgpu process as python class
"""
class ljgpu(object):
    def __init__(self, **kwargs):
        self.__opts = {}
        self.__set(**kwargs)

    def set(self, **kwargs):
        self.__set(**kwargs)

    def run(self):
        proc = subprocess.Popen(self.__get())
        try:
            proc.wait()
        finally:
            ret = proc.wait()
            if ret:
                raise SystemExit(ret)

    def __str__(self):
        return ' '.join(self.__get())

    def __get(self):
        args = ['ljgpu']
        for key, val in self.__opts.iteritems():
            # append option name
            opt = '--%s' % key.replace('_', '-')
            if not val is False:
                args.append(opt)
            # append option value
            if hasattr(val, '__iter__'):
                args.append(','.join([str(v) for v in val]))
            elif not isinstance(val, bool):
                args.append(str(val))
        return args

    def __set(self, **kwargs):
        for key, val in kwargs.iteritems():
            self.__opts[key] = val


if __name__ == '__main__':
    import sys
    mdsim = ljgpu(verbose=True, density=0.5, steps=100)
    mdsim.set(binary=(2000,8000), random_seed=42, cell_occupancy=0.2)
    for timestep in (0.001, 0.002, 0.003):
        mdsim.set(timestep=timestep)
        print >> sys.stderr, '>>> %s' % mdsim
        try:
            mdsim.run()
        except SystemExit, e:
            print >> sys.stderr, '>>> ljgpu failed with exit code %s' % e
            raise

