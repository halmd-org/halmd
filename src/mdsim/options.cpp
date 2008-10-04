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

#include <H5Cpp.h>
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include "H5xx.hpp"
#include "date_time.hpp"
#include "options.hpp"
#include "version.h"
namespace po = boost::program_options;


namespace boost { namespace program_options
{

/**
 * program option value validation for double precision floating-point values
 */
template <>
void validate(boost::any& value_store, std::vector<std::string> const& values, double*, long)
{
    std::string const& s = validators::get_single_string(values);
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
void validate(boost::any& value_store, std::vector<std::string> const& values, unsigned int*, long)
{
    std::string const& s = validators::get_single_string(values);
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
void validate(boost::any& value_store, std::vector<std::string> const& values, uint64_t*, long)
{
    std::string const& s = validators::get_single_string(values);
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
void conflicting_options(const variables_map& vm, char const* opt1, char const* opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted() && vm.count(opt2) && !vm[opt2].defaulted()) {
	throw std::logic_error(std::string("conflicting options '") + opt1 + "' and '" + opt2 + "'");
    }
}

/**
 * Function used to check that 'opt2' is specified if 'opt1' is specified
 */
void dependent_option(const variables_map& vm, char const* opt1, char const* opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted() && !(vm.count(opt2) && !vm[opt2].defaulted())) {
	throw std::logic_error(std::string("option '") + opt1 + "' requires missing option '" + opt2 + "'");
    }
}

/**
 * Accumulating program option value
 */
template <typename T, typename charT = char>
class accumulating_value : public boost::program_options::typed_value<T,charT>
{
    //
    // Originally written by Bryan Green <bgreen@nas.nasa.gov>
    // http://article.gmane.org/gmane.comp.lib.boost.user/29084
    //

public:
    accumulating_value(T* store_to=0) : boost::program_options::typed_value<T,charT>(store_to), origin(0)
    {
	(void) boost::program_options::typed_value<T,charT>::zero_tokens();
    }

    accumulating_value* default_value(const T &v)
    {
	// setting a default value sets the origin to that value
	origin = v;
	(void) boost::program_options::typed_value<T,charT>::default_value(v);
	return this;
    }

    accumulating_value* default_value(const T &v,const std::string& textual)
    {
	// setting a default value sets the origin to that value
	origin = v;
	(void) boost::program_options::typed_value<T,charT>::default_value(v, textual);
	return this;
    }

    void xparse(boost::any& value_store, const std::vector<std::basic_string<charT> >& new_tokens) const
    {
	// if this is the first occurrence of the option, initialize it to the origin
	if (value_store.empty())
	    value_store = boost::any(origin);

	++boost::any_cast<T&>(value_store);
    }

private:
    /** the numeric origin from which to increment upward */
    T origin;
};

template <typename T>
accumulating_value<T>* accum_value(T *v)
{
    accumulating_value<T>* r = new accumulating_value<T>(v);
    return r;
}

template <typename T>
accumulating_value<T>* accum_value()
{
    return accum_value<T>(0);
}

/**
 * returns HDF5 attribute value if attribute exists, or empty value otherwise
 */
template <typename T>
boost::any parse_attribute(H5xx::group const& node, char const* name)
{
    try {
	return boost::any(node[name].as<T>());
    }
    catch (H5::Exception const& e) {
	// discard missing HDF5 attribute for backward compatibility
	return boost::any();
    }
}

/**
 * override variable value if defaulted
 */
void store(boost::any const& value, variable_value& vv) {
    if (vv.defaulted() && !value.empty()) {
	vv = variable_value(value, true);
    }
}

}} // namespace boost::program_options


namespace mdsim {

/**
 * parse program option values
 */
void options::parse(int argc, char** argv)
{
    po::options_description mdsim_opts("MD simulation parameters");
    mdsim_opts.add_options()
	("particles,N", po::value<unsigned int>()->default_value(128), "number of particles")
	("density,d", po::value<double>()->default_value(0.1), "particle density")
	("box-length,L", po::value<double>(), "simulation box length")
	("timestep,r", po::value<double>()->default_value(0.01), "simulation timestep")
	("temperature,K", po::value<double>()->default_value(1.), "initial temperature")
	("rng-seed,R", po::value<unsigned int>()->default_value(42), "random number generator integer seed")
#ifndef USE_BENCHMARK
	("steps,s", po::value<uint64_t>()->default_value(10000), "number of simulation steps")
	("time,t", po::value<double>(), "total simulation time")
	("equilibrate,E", po::value<uint64_t>()->default_value(0), "number of equilibration steps")
#else
	("equilibrate,E", po::value<uint64_t>()->default_value(10000), "number of equilibration steps")
#endif
	("trajectory,I", po::value<std::string>(), "trajectory input file")
	("sample,S", po::value<int64_t>(), "resume from sample in trajectory file")
	;

#ifndef USE_BENCHMARK
    po::options_description tcf_opts("Time correlation function options");
    tcf_opts.add_options()
	("block-size,B", po::value<unsigned int>()->default_value(10), "block size")
	("max-samples,M", po::value<uint64_t>()->default_value(1000), "maximum number of samples per block")
	("q-values", po::value<unsigned int>()->default_value(5), "largest multiple of smallest q-vector for Fourier transformation")
	("dump-trajectories", po::bool_switch(), "dump particle trajectories")
	;
#endif

    po::options_description misc_opts("Other options");
    misc_opts.add_options()
	("output,o", po::value<std::string>()->default_value(PROGRAM_NAME "_%Y%m%d_%H%M%S"), "output file prefix")
	("input,i", po::value<std::vector<std::string> >(), "parameter input file")
	("dry-run,n", po::bool_switch(), "perform a trial run without simulation")
	("daemon", po::bool_switch(), "run program in background")
	("verbose,v", po::accum_value<int>()->default_value(0), "increase verbosity")
	("version,V", "output version and exit")
	("help,h", "display this help and exit")
	;

    po::options_description opts;
    opts.add(mdsim_opts);
#ifndef USE_BENCHMARK
    opts.add(tcf_opts);
#endif
    opts.add(misc_opts);

    try {
	// parse command line options
	po::store(po::parse_command_line(argc, argv, opts), vm);

	// parse optional parameter input files
	if (vm.count("input")) {
	    std::vector<std::string> const& files = vm["input"].as<std::vector<std::string> >();

	    for (std::vector<std::string>::const_iterator it = files.begin(); it != files.end(); ++it) {
		std::ifstream ifs(it->c_str());
		if (ifs.fail()) {
		    std::cerr << PROGRAM_NAME ": could not open parameter input file '" << *it << "'\n";
		    throw options::exit_exception(EXIT_FAILURE);
		}

		po::store(po::parse_config_file(ifs, opts), vm);
	    }
	}
    }
    catch (std::exception const& e) {
	std::cerr << PROGRAM_NAME ": " << e.what() << "\n";
	std::cerr << "Try `" PROGRAM_NAME " --help' for more information.\n";
	throw options::exit_exception(EXIT_FAILURE);
    }

    po::notify(vm);

    try {
	po::conflicting_options(vm, "density", "box-length");
	po::conflicting_options(vm, "steps", "time");
	po::dependent_option(vm, "sample", "trajectory");
    }
    catch (std::exception const& e) {
	std::cerr << PROGRAM_NAME ": " << e.what() << "\n";
	throw options::exit_exception(EXIT_FAILURE);
    }

    // override const operator[] in variables_map
    std::map<std::string, po::variable_value>& vm_ = vm;

    // optionally read parameters from HDF5 input file
    if (vm.count("trajectory")) {
	// store absolute input file path
	boost::filesystem::path path(vm["trajectory"].as<std::string>());
	vm_["trajectory"] = po::variable_value(boost::filesystem::complete(path).string(), false);

	try {
	    H5::H5File file(vm["trajectory"].as<std::string>(), H5F_ACC_RDONLY);
	    H5::Group param(file.openGroup("param")), node;

	    node = param.openGroup("mdsim");
	    po::store(po::parse_attribute<unsigned int>(node, "particles"), vm_["particles"]);
	    po::store(po::parse_attribute<double>(node, "density"), vm_["density"]);
	    po::store(po::parse_attribute<double>(node, "box_length"), vm_["box-length"]);
	    po::store(po::parse_attribute<double>(node, "timestep"), vm_["timestep"]);
	    po::store(po::parse_attribute<double>(node, "temperature"), vm_["temperature"]);

	    node = param.openGroup("correlation");
	    po::store(po::parse_attribute<uint64_t>(node, "steps"), vm_["steps"]);
	    po::store(po::parse_attribute<double>(node, "time"), vm_["time"]);
	    po::store(po::parse_attribute<unsigned int>(node, "block_size"), vm_["block-size"]);
	    po::store(po::parse_attribute<uint64_t>(node, "max_samples"), vm_["max-samples"]);
	    po::store(po::parse_attribute<unsigned int>(node, "q_values"), vm_["q-values"]);
	}
	catch (H5::Exception const& e) {
	    std::cerr << PROGRAM_NAME ": " << "failed to read parameters from HDF5 input file\n";
	    throw options::exit_exception(EXIT_FAILURE);
	}
    }

    // format timestamp in output file prefix
    boost::filesystem::path path(date_time::format(vm["output"].as<std::string>()));
    // store absolute output file path
    vm_["output"] = po::variable_value(boost::filesystem::complete(path).string(), false);

    if (vm.count("help")) {
	std::cout << "Usage: " PROGRAM_NAME " [OPTION]...\n" << opts << "\n";
	throw options::exit_exception(EXIT_SUCCESS);
    }

    if (vm.count("version")) {
	std::cout << PROGRAM_NAME " " PROGRAM_VERSION "\n"
	    "variant: " PROGRAM_VARIANT "\n"
	    "\n" PROGRAM_COPYRIGHT "\n" "This is free software. "
	    "You may redistribute copies of it under the terms of\n"
	    "the GNU General Public License "
	    "<http://www.gnu.org/licenses/gpl.html>.\n"
	    "There is NO WARRANTY, to the extent permitted by law.\n";
	throw options::exit_exception(EXIT_SUCCESS);
    }
}

} // namespace mdsim
