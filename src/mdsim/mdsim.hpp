/* Event-based Molecular Dynamics simulation of hard spheres
 *
 * Copyright (C) 2008  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MDSIM_MDSIM_HPP
#define MDSIM_MDSIM_HPP

#include "correlation.hpp"
#include "energy.hpp"
#include "hardspheres.hpp"
#include "options.hpp"
#include "perf.hpp"
#include "trajectory.hpp"

namespace mdsim
{

/**
 * Event-based Molecular Dynamics simulation of hard spheres
 */
class mdsim
{
public:
    enum {
	/** HDF5 buffers flush to disk interval in seconds */
	FLUSH_TO_DISK_INTERVAL = 900,
	/** waiting time in seconds before runtime estimate after block completion */
	TIME_ESTIMATE_WAIT_AFTER_BLOCK = 300,
	/** runtime estimate interval in seconds */
	TIME_ESTIMATE_INTERVAL = 1800,
    };

public:
    /** initialize MD simulation program */
    mdsim(options const& opts);
    /** run MD simulation program */
    void operator()();

private:
    /** program options */
    options const& opts;
    /** hard spheres simulation */
    hardspheres fluid;
#ifndef USE_BENCHMARK
    /** block correlations */
    correlation tcf;
    /**  trajectory file writer */
    trajectory<true> traj;
    /** thermodynamic equilibrium properties */
    energy tep;
#endif
    /** performance data */
    perf prf;
};

} // namespace mdsim

#endif /* ! MDSIM_MDSIM_HPP */
