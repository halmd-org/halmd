#!/bin/bash
#
# Copyright © 2011-2017 Felix Höfling
#
# This file is part of HALMD.
#
# HALMD is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General
# Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
#

##
# Extract and print perfomance metrics from HALMD logging output
#
if [ "$1" = "--help" -o $# -eq 0 ]
then
    echo -e "Usage: $0 LOG_FILE\n"
    exit
fi

LOG_FILE=$1

TIMINGS=$(sed -n -e 's/.*MD integration step: \([0-9.]*\) \(..\).*/\1 \2/p' "${LOG_FILE}")
PARTICLES=$(sed -n -e '/number of particles/{s/.*: \([0-9]*\).*/\1/p;q}' "${LOG_FILE}")
echo -e "$TIMINGS" | gawk -v N=$PARTICLES '{
    factor = ($2 == "µs") ? 1e-3 : 1;
    a += $1 * factor; n+=1
}END{\
    a = a/n;
    print N, "particles"; \
    print a, "ms per step"; \
    print 1e6*a/N, "ns per step and particle"; \
    print 1000/a, "steps per second" \
}'

