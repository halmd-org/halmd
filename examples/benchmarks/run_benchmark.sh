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
# Run a benchmarking suite
#

if [ "$1" = "--help" ]
then
    echo -e "Usage: run_benchmark.sh BENCHMARK_NAME COUNT DEVICE_NAME SUFFIX\n"
    exit
fi

BENCHMARK_NAME=${1:-"lennard_jones"}
COUNT=${2:-5}
DEVICE_NAME=${3:-$(nvidia-smi -a | sed -ne '/Product Name/{s/.*Tesla \([A-Z][0-9]\+\).*/\1/p;q}')}
SUFFIX=${4:+_$4}
HALMD_VERSION=$(halmd --version | sed -e '1s/.*-g\([a-z0-9]\+\) (.*)$/\1/;q')
BENCHMARK_TAG="${DEVICE_NAME}_${HALMD_VERSION}${SUFFIX}"

INPUT_DIR=$PWD
WORKING_DIR=$PWD/data

CONFIG_DIR=${INPUT_DIR}/${BENCHMARK_NAME}
OUTPUT_DIR=${WORKING_DIR}/${BENCHMARK_NAME}

# run benchmark several times by continuation of the trajectory
PREVIOUS_OUTPUT_PREFIX="${OUTPUT_DIR}/configuration"
for I in $(seq $COUNT)
do
    OUTPUT_PREFIX="${OUTPUT_DIR}/benchmark_${BENCHMARK_TAG}-${I}"
    echo halmd \
      --verbose \
      --config "${CONFIG_DIR}/run_benchmark.rc" \
      --output "${OUTPUT_PREFIX}" \
      trajectory --file "${PREVIOUS_OUTPUT_PREFIX}.trj"

    PREVIOUS_OUTPUT_PREFIX="${OUTPUT_PREFIX}"
done

TIMINGS=`grep -h "MD integration step:" "${OUTPUT_DIR}/benchmark_${BENCHMARK_TAG}"-*.log`
echo -e "$TIMINGS"
echo -e "$TIMINGS" | gawk '{a += $6; n+=1}END{\
    a = a/n;
    print a, "ms per step"; \
    print 1e6*a/64000, "ns per step and particle"; \
    print 1000/a, "steps per second" \
}'
