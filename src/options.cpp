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


namespace boost { namespace program_options
{

/**
 * program option value validation for double precision floating-point values
 */
template <>
void validate(boost::any& value_store, std::vector<std::string> const& values, double*, long)
{
    std::string const& s = po::validators::get_single_string(values);
    double value;

    try {
	value = lexical_cast<double>(s);
    }
    catch (bad_lexical_cast const& e)
    {
	throw po::invalid_option_value(s);
    }

    // require positive value
    if (value > 0.) {
	value_store = boost::any(value);
    }
    else {
	throw po::invalid_option_value(s);
    }
}

/**
 * program option value validation for unsigned integer values
 */
template <>
void validate(boost::any& value_store, std::vector<std::string> const& values, unsigned int*, long)
{
    std::string const& s = po::validators::get_single_string(values);
    unsigned int value;

    try {
	value = lexical_cast<unsigned int>(s);
    }
    catch (bad_lexical_cast const& e)
    {
	throw po::invalid_option_value(s);
    }

    // require positive value
    if (value > 0) {
	value_store = boost::any(value);
    }
    else {
	throw po::invalid_option_value(s);
    }
}

/**
 * program option value validation for 64-bit unsigned integer values
 */
template <>
void validate(boost::any& value_store, std::vector<std::string> const& values, uint64_t*, long)
{
    std::string const& s = po::validators::get_single_string(values);
    uint64_t value;

    try {
	value = lexical_cast<uint64_t>(s);
    }
    catch (bad_lexical_cast const& e)
    {
	throw po::invalid_option_value(s);
    }

    // require positive value
    if (value > 0) {
	value_store = boost::any(value);
    }
    else {
	throw po::invalid_option_value(s);
    }
}

}} // namespace boost::program_options


namespace mdsim {

/**
 * set default parameters
 */
options::options()
{
    npart_ = 10000;
    density_ = 0.05;
    timestep_ = 0.005;
    temp_ = 1.;
    steps_ = 1000;
    avgsteps_ = 100;

    rngseed_ = 123;
}

/**
 * parse program option values
 */
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
