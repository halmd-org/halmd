#!/bin/bash
#
# Copyright © 2011-2017 Felix Höfling
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
DEVICE_NAME=${5:-$(nvidia-smi -a | sed -ne '/Product Name/{s/.*: [A-Za-z]* \(.*\)/\1/;s/ //g;p;q}')}
HALMD_OPTIONS=$6

HALMD_VERSION=$(halmd --version | cut -c 26- | sed -e '1s/-patch.* \([a-z0-9]\+\)\]/-g\1/;q')
BENCHMARK_TAG="${DEVICE_NAME}_${HALMD_VERSION}${SUFFIX}"

SCRIPT="${SCRIPT_DIR}/${BENCHMARK_NAME}/run_benchmark.lua"
OUTPUT="${BENCHMARK_NAME}/benchmark_${BENCHMARK_TAG}"

# run benchmark
halmd ${HALMD_OPTIONS} "${SCRIPT}" \
  --trajectory "${INPUT_FILE}" \
  --output "${OUTPUT}" \
  --count "${COUNT}" \
  --verbose

# print results
${SCRIPT_DIR}/print_timings.sh "${OUTPUT}.log"

