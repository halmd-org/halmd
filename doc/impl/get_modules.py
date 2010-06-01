#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Extract module dependencies and options from source files
#
# Copyright © 2010  Felix Höfling
#
# This file is part of HALMD.
#
# HALMD is free software: you can redistribute it and/or modify
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

import re, os, sys
import scriptutil as su

def get_base(module):
    """Return root of inheritance tree."""
    if module in baseclass:
        return get_base(baseclass[module])
    else:
        return module

def get_result_filename(module, path=''):
    """Return name of result file for given module."""
    base = get_base(module)
    base = base.replace('halmd::', '')
    return os.path.join(path, base.replace('::', '_') + '.txt')

# scan header files for resolve() function
basepath = '../../src'
module_files = su.ffindgrep(basepath,
                            namefs=(lambda fn: fn.endswith('.hpp'),),
                            regexl=(r'static void resolve\(',)).keys()

modules = map(lambda s: s[len(basepath)+1:-4].replace('/', '::'), module_files)
modules.sort()

print """
List of modules
---------------
"""
for s in modules:
    print '    %s' % s
print

baseclass = dict()
typedefs = dict()

# extract base class and typedefs from header file
for m in modules:
    fn = os.path.join(basepath, m.replace('::', '/') + '.hpp')
    if not os.access(fn, os.R_OK):
        print >>sys.stderr, 'Error: cannot read module definition in %s' % fn
        modules.remove(m)
        continue

    fh = open(fn, 'r')
    lines = fh.read().splitlines()
    fh.close()

    expr = re.compile(r':\s+public\s+([\w:]+)')
    for s in lines:
        match = expr.search(s)
        if match:
            baseclass[m] = match.group(1)

    typedefs[m] = dict()
    expr = re.compile(r'typedef\s+(.+)\s+(\w+);')
    for s in lines:
        match = expr.search(s)
        if match:
            typedefs[m][match.group(2)] = match.group(1)

# remove result files from previous runs (dangerous!?)
for m in modules:
    fn = get_result_filename(m)
    if os.access(fn, os.W_OK):
        os.remove(fn)

# determine dependencies from calls to module<...>::required()
for m in modules:
    # check whether cpp file exists and read contents
    fn = os.path.join(basepath, m.replace('::', '/') + '.cpp')
    if not os.access(fn, os.R_OK):
        print >>sys.stderr, 'Note: %s is an abstract base class' % m
        modules.remove(m)
        continue

    fh = open(fn, 'r')
    content = fh.read()
    lines = content.splitlines()
    fh.close()

    # extract dependencies
    dependencies = []
    expr = re.compile(r'module<(.*)>::required')
    for s in lines:
        match = expr.search(s)
        if match:
            dependencies.append(match.group(1))
    dependencies.sort()

    # extract options
    options = []
    '\n'.join(lines)
    match = re.search(r'add_options\(\)\s*\(([^;]+)*\)', '\n'.join(lines))
    if match:
        for opt in re.split(r'(?<=")\)\s*\((?=")', match.group(1)):
            match = re.search(r'"(?P<name>[\w-]+)(?:,(?P<short>\w))?",' +
                              r'.*value<(?P<type>.+)>\(\).*,\s*"(?P<desc>.*)"', opt)
            if match:
                options.append(match.groupdict())

    # collect results in one file per base module
    if dependencies or options:
        fh = open(get_result_filename(m), 'a')
        print >>fh, '**Module %s**\n' % m
#        print >>fh, 'Base module: %s' % get_base(m)

        if dependencies:
            print >>fh, '* Dependencies:\n'
            for type in dependencies:
                if type in typedefs[m]:
                    type = typedefs[m][type]
                type = type.split('<')[0]
                print >>fh, '  - %s' % type
            print >>fh

        if options:
            # create table
            table = []
            for opts in options:
                row = 3 * [None]
                if opts['short'] == None:
                    row[0] = '--%(name)s' % opts
                else:
                    row[0] = '--%(name)s, -%(short)s' % opts
                row[1] = opts['type']
                row[2] = opts['desc']
                table.append(row)
            # determine column widths
            width = 3 * [None]
            rule = '+'
            for i in range(3):
                width[i] = max([len(row[i].decode('utf-8')) for row in table])
                rule += (width[i] + 2) * '-' + '+'
            # print table
            print >>fh, '* Options:\n'
            print >>fh, '  .. table::'
            print >>fh, '\n    ' + rule
            for row in table:
                line = '    '
                for n,s in enumerate(row):
                    line += '| %s ' % s.ljust(width[n])
                line += '|'
                print >>fh, line + '\n    ' + rule
            print >>fh

        fh.close()

