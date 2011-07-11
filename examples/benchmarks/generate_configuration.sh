#/usr/bin/bash
#
# Copyright © 2011  Felix Höfling
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
# Generate an initial configuration for a benchmarking suite
#

if [ "$1" = "--help" ]
then
    echo -e "Usage: generate_configuration.sh BENCHMARK_NAME DEVICE_NAME SUFFIX\n"
    exit
fi

DEVICE_NAME=`nvidia-smi -a | sed -ne '/Product Name/{s/.*Tesla \([A-Z][0-9]\+\).*/\1/p;q}'`

BENCHMARK_NAME=${1:-"lennard_jones"}
DEVICE_NAME=${2:-`nvidia-smi -a | sed -ne '/Product Name/{s/.*Tesla \([A-Z][0-9]\+\).*/\1/p;q}'`}
BENCHMARK_TAG="${DEVICE_NAME}${3:+_$3}"

INPUT_DIR=$PWD
WORKING_DIR=$PWD/data

CONFIG_DIR=${INPUT_DIR}/${BENCHMARK_NAME}
OUTPUT_DIR=${WORKING_DIR}/${BENCHMARK_NAME}

OUTPUT_PREFIX="${OUTPUT_DIR}/configuration"
halmd \
  --verbose \
  --config "${CONFIG_DIR}/generate_configuration.rc" \
  --output "${OUTPUT_PREFIX}"
