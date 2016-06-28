#!/bin/bash
#
# Copyright © 2011-2012  Felix Höfling
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

##
# Generate an initial configuration for a benchmarking suite
#

if [ "$1" = "--help" -o $# -eq 0 ]
then
    echo -e "Usage: generate_configuration.sh BENCHMARK_NAME [SUFFIX [HALMD_OPTIONS]]\n"
    exit
fi

SCRIPT_DIR="$(dirname $0)"
BENCHMARK_NAME=$1
SUFFIX=${2:+_$2}
HALMD_OPTIONS=$3

SCRIPT="${SCRIPT_DIR}/${BENCHMARK_NAME}/generate_configuration.lua"
OUTPUT="${BENCHMARK_NAME}/configuration${SUFFIX}"

halmd \
  "${SCRIPT}" \
  --verbose \
  --output "${OUTPUT}" \
  ${HALMD_OPTIONS}
