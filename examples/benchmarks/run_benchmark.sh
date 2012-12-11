#!/bin/bash
#
# Copyright © 2011-2012  Felix Höfling
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

##
# Run a benchmarking suite from an input configuration
#

if [ "$1" = "--help" -o $# -eq 0 ]
then
    echo -e "Usage: run_benchmark.sh BENCHMARK_NAME [COUNT [INPUT_FILE [SUFFIX [DEVICE_NAME [HALMD_OPTIONS]]]]]\n"
    exit
fi

SCRIPT_DIR="$(dirname $0)"
BENCHMARK_NAME=$1
COUNT=${2:-5}
INPUT_FILE=${3:-"${BENCHMARK_NAME}/configuration.h5"}
SUFFIX=${4:+_$4}
DEVICE_NAME=${5:-$(nvidia-smi -a | sed -ne '/Product Name/{s/.*Tesla \([A-Z][0-9]\+\).*/\1/p;q}')}
HALMD_OPTIONS=$6

HALMD_VERSION=$(halmd --version | cut -c 26- | sed -e '1s/.*-g\([a-z0-9]\+\).*/\1/;q')
BENCHMARK_TAG="${DEVICE_NAME}_${HALMD_VERSION}${SUFFIX}"

SCRIPT="${SCRIPT_DIR}/${BENCHMARK_NAME}/run_benchmark.lua"
OUTPUT="${BENCHMARK_NAME}/benchmark_${BENCHMARK_TAG}"

# run benchmark
halmd "${SCRIPT}" \
  --trajectory "${INPUT_FILE}" \
  --output "${OUTPUT}" \
  --count "${COUNT}" \
  --verbose \
  ${HALMD_OPTIONS}

TIMINGS=$(sed -n -e 's/.*MD integration step: \([0-9.]*\).*/\1/p' "${OUTPUT}.log")
PARTICLES=$(sed -n -e 's/.*number of particles: \([0-9]*\).*/\1/p' "${OUTPUT}.log")
echo -e "$TIMINGS" | gawk -v N=$PARTICLES '{a += $1; n+=1}END{\
    a = a/n;
    print N, "particles"; \
    print a, "ms per step"; \
    print 1e6*a/N, "ns per step and particle"; \
    print 1000/a, "steps per second" \
}'
