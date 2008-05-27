/* Molecular Dynamics Simulation of a Lennard-Jones fluid
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

#include <exception>
#include <iostream>
#include <vector>
#include "log.hpp"
#include "mdsim.hpp"
#include "options.hpp"
#include "vector2d.hpp"
#include "vector3d.hpp"
#include "version.h"


int main(int argc, char **argv)
{
    mdsim::options opts;

    // parse program options
    try {
	opts.parse(argc, argv);
    }
    catch (mdsim::options::exit_exception const& e) {
	return e.status();
    }

    mdsim::log::init(opts);

    LOG(PROGRAM_NAME " (" PROGRAM_VERSION ")");

    try {
	// initialize molecular dynamics simulation
#ifdef DIM_3D
	mdsim::mdsim<3, vector3d<double> > sim(opts);
#else
	mdsim::mdsim<2, vector2d<double> > sim(opts);
#endif
	// run MD simulation
	sim();
    }
    catch (std::exception const& e) {
	LOG_ERROR(e.what());
	LOG_WARNING(PROGRAM_NAME " aborted");
	return EXIT_FAILURE;
    }

    LOG(PROGRAM_NAME " exit");
    return EXIT_SUCCESS;
}
