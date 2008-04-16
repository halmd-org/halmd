/* Molecular Dynamics simulation program options
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

#include "options.hpp"
#include "version.h"
#include <iostream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;


namespace mdsim {

options::options()
{
    //
    // default program parameters
    //

    npart_ = 10000;
    density_ = 0.05;
    timestep_ = 0.005;
    temp_ = 1.;
    steps_ = 1000;
    avgsteps_ = 100;

    rngseed_ = 123;
}

void options::parse(int argc, char** argv)
{
    po::options_description gen_opts("General options");
    gen_opts.add_options()
	("version,V", "output version and exit")
	("help,h", "display this help and exit")
	;

    po::options_description mdsim_opts("MD simulation parameters");
    mdsim_opts.add_options()
	("particles,N", po::value(&npart_), "number of particles")
	("density,d", po::value(&density_), "particle density")
	("timestep,t", po::value(&timestep_), "simulation timestep")
	("temperature,T", po::value(&temp_), "initial temperature")
	("steps,s", po::value(&steps_), "number of simulation steps")
	("average,S", po::value(&avgsteps_), "number of average accumulation steps")
	;

    po::options_description misc_opts("Miscellaneous options");
    misc_opts.add_options()
	("seed,R", po::value(&rngseed_), "random number generator integer seed")
	;

    po::options_description opts;
    opts.add(gen_opts).add(mdsim_opts).add(misc_opts);

    po::variables_map vm;

    try {
	po::store(po::command_line_parser(argc, argv).options(opts).run(), vm);
    }
    catch (std::exception const& e) {
	std::cerr << PROGRAM_NAME << ": " << e.what() << "\n";
	std::cerr << "Try `" << PROGRAM_NAME << " --help' for more information.\n";
	throw options::exception(EXIT_FAILURE);
    }

    po::notify(vm);

    if (vm.count("help")) {
	std::cout << "Usage: " PROGRAM_NAME " [OPTION]..." << opts;
	throw options::exception(EXIT_SUCCESS);
    }

    if (vm.count("version")) {
	std::cout << PROGRAM_NAME " (" PROGRAM_VERSION ")\n"
	    "\n" PROGRAM_COPYRIGHT "\n" "This is free software. "
	    "You may redistribute copies of it under the terms of\n"
	    "the GNU General Public License "
	    "<http://www.gnu.org/licenses/gpl.html>.\n"
	    "There is NO WARRANTY, to the extent permitted by law.\n";
	throw options::exception(EXIT_SUCCESS);
    }
}

} // namespace mdsim
