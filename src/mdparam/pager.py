#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# pager.py - Git-style paging if output is longer than screen
#
# Copyright © 2008-2009  Peter Colberg
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
#

import os, sys
import select

"""
The following code is based on pager.c shipped with Git v1.5.6.5.
"""

def run_pager(pager):
    # Work around bug in "less" by not starting it until we have real input
    select.select([0], [], [0])
    os.execlp(pager, pager)
    os.execl('/bin/sh', 'sh', '-c', pager)


def setup_pager():
    if not hasattr(sys.stdout, 'isatty'):
        return

    pager = os.getenv('PAGER')
    if pager is None:
        pager = 'less'
    elif pager == 'cat':
        return

    fd = os.pipe()
    pid = os.fork()
    if pid < 0:
        os.close(fd[0])
        os.close(fd[1])
        return

    # return in the child
    if pid == 0:
        os.dup2(fd[1], 1)
        os.dup2(fd[1], 2)
        os.close(fd[0])
        os.close(fd[1])
        return

    # The original process turns into the PAGER
    os.dup2(fd[0], 0)
    os.close(fd[0])
    os.close(fd[1])

    os.environ['LESS'] = 'FRSX'
    try:
        run_pager(pager)
    except os.OSError:
        print >> sys.stderr, "unable to execute pager '%s'" % pager
        os._exit(255)


def test_pager():
    setup_pager()
    # single screen of output
    from time import sleep
    for line in ('line 1', 'line 2'):
        print line
        sys.stdout.flush()
        sleep(1)
    # multiple screens of output
    s = file(__file__).read()
    sys.stdout.write('%s%s%s%s' % (s, s, s, s))


if __name__ == '__main__':
    test_pager()

