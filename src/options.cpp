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
#include <fstream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;


namespace boost { namespace program_options
{

/**
 * program option value validation for double precision floating-point values
 */
template <>
void validate(boost::any& value_store, vector<string> const& values, double*, long)
{
    string const& s = validators::get_single_string(values);
    double value;

    try {
	value = lexical_cast<double>(s);
    }
    catch (bad_lexical_cast const& e)
    {
	throw invalid_option_value(s);
    }

    // require positive value
    if (value > 0.) {
	value_store = boost::any(value);
    }
    else {
	throw invalid_option_value(s);
    }
}

/**
 * program option value validation for unsigned integer values
 */
template <>
void validate(boost::any& value_store, vector<string> const& values, unsigned int*, long)
{
    string const& s = validators::get_single_string(values);
    unsigned int value;

    try {
	value = lexical_cast<unsigned int>(s);
    }
    catch (bad_lexical_cast const& e)
    {
	throw invalid_option_value(s);
    }

    // require positive value
    if (value > 0) {
	value_store = boost::any(value);
    }
    else {
	throw invalid_option_value(s);
    }
}

/**
 * program option value validation for 64-bit unsigned integer values
 */
template <>
void validate(boost::any& value_store, vector<string> const& values, uint64_t*, long)
{
    string const& s = validators::get_single_string(values);
    uint64_t value;

    try {
	value = lexical_cast<uint64_t>(s);
    }
    catch (bad_lexical_cast const& e)
    {
	throw invalid_option_value(s);
    }

    // require positive value
    if (value > 0) {
	value_store = boost::any(value);
    }
    else {
	throw invalid_option_value(s);
    }
}

/**
 * Function used to check that 'opt1' and 'opt2' are not specified at the same time.
 */
void conflicting_options(const variables_map& vm, const char* opt1, const char* opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted() && vm.count(opt2) && !vm[opt2].defaulted()) {
	throw logic_error(string("conflicting options '") + opt1 + "' and '" + opt2 + "'");
    }
}

}} // namespace boost::program_options


namespace mdsim {

/**
 * set default parameters
 */
options::options()
{
    // MD simulation parameters
    npart_ = 10000;
    density_ = 0.05;
    timestep_ = 0.005;
    temp_ = 1.;
    steps_ = 1000;
    avgsteps_ = 100;

    // Autocorrelation options
    block_count_ = 6;
    block_size_ = 6;
    block_shift_ = 2;
    max_samples_ = 1000;

    // Other options
    rngseed_ = 42;
    output_file_prefix_ = PROGRAM_NAME;
    quiet_ = false;
}

/**
 * parse program option values
 */
void options::parse(int argc, char** argv)
{
    po::options_description mdsim_opts("MD simulation parameters");
    mdsim_opts.add_options()
	("particles,N", po::value(&npart_), "number of particles")
	("density,d", po::value(&density_), "particle density")
	("timestep,t", po::value(&timestep_), "simulation timestep")
	("temperature,K", po::value(&temp_), "initial temperature")
	("steps,s", po::value(&steps_), "number of simulation steps")
	("average,S", po::value(&avgsteps_), "number of average accumulation steps")
	;

    po::options_description tcf_opts("Autocorrelation options");
    tcf_opts.add_options()
	("block-count", po::value(&block_count_), "block count")
	("block-size", po::value(&block_size_), "block size")
	("block-shift", po::value(&block_shift_), "block shift")
	("max-samples", po::value(&max_samples_), "maximum number of samples per block")
	;

    po::options_description misc_opts("Other options");
    misc_opts.add_options()
	("seed,R", po::value(&rngseed_), "random number generator integer seed")
	("input,i", po::value<vector<string> >(), "parameter input file")
	("output,o", po::value(&output_file_prefix_), "output file prefix")
	("quiet,q", po::bool_switch(&quiet_), "suppress status output")
	("version,V", "output version and exit")
	("help,h", "display this help and exit")
	;

    po::options_description opts;
    opts.add(mdsim_opts).add(tcf_opts).add(misc_opts);

    po::variables_map vm;

    try {
	// parse command line options
	po::store(po::parse_command_line(argc, argv, opts), vm);

	// parse optional parameter input files
	if (vm.count("input")) {
	    vector<string> const& files = vm["input"].as<vector<string> >();

	    for (vector<string>::const_iterator it = files.begin(); it != files.end(); ++it) {
		ifstream ifs(it->c_str());
		if (ifs.fail()) {
		    cerr << PROGRAM_NAME ": could not open parameter input file '" << *it << "'\n";
		    throw options::exit_exception(EXIT_FAILURE);
		}

		po::store(po::parse_config_file(ifs, opts), vm);
	    }
	}
    }
    catch (std::exception const& e) {
	cerr << PROGRAM_NAME ": " << e.what() << "\n";
	cerr << "Try `" PROGRAM_NAME " --help' for more information.\n";
	throw options::exit_exception(EXIT_FAILURE);
    }

    po::notify(vm);

    if (vm.count("help")) {
	cout << "Usage: " PROGRAM_NAME " [OPTION]...\n" << opts << "\n";
	throw options::exit_exception(EXIT_SUCCESS);
    }

    if (vm.count("version")) {
	cout << PROGRAM_NAME " (" PROGRAM_VERSION ")\n"
	    "\n" PROGRAM_COPYRIGHT "\n" "This is free software. "
	    "You may redistribute copies of it under the terms of\n"
	    "the GNU General Public License "
	    "<http://www.gnu.org/licenses/gpl.html>.\n"
	    "There is NO WARRANTY, to the extent permitted by law.\n";
	throw options::exit_exception(EXIT_SUCCESS);
    }
}

} // namespace mdsim
