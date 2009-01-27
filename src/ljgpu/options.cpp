/* Molecular Dynamics simulation program options
 *
 * Copyright © 2008-2009  Peter Colberg
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
#include <algorithm>
#include <boost/array.hpp>
#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <fstream>
#include <iostream>
#include <ljgpu/options.hpp>
#include <ljgpu/util/H5xx.hpp>
#include <ljgpu/util/date_time.hpp>
#include <ljgpu/version.h>
#include <map>
namespace po = boost::program_options;

#define foreach BOOST_FOREACH

namespace boost { namespace program_options
{

/**
 * program option value validation for single precision floating-point values
 */
template <>
void validate(boost::any& value_store, std::vector<std::string> const& values, float*, long)
{
    std::string const& s = validators::get_single_string(values);
    float value;

    try {
	value = lexical_cast<float>(s);
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
 * override variable value if defaulted or empty
 */
void store(boost::any const& value, variable_value& vv) {
    if (!value.empty() && (vv.defaulted() || vv.empty())) {
	vv = variable_value(value, true);
    }
}

}} // namespace boost::program_options

namespace std
{

/**
 * extract comma-separated option values into fixed-size array
 */
template <typename T, size_t size>
static istream& operator>>(istream& is, boost::array<T, size>& value)
{
    foreach (T& v, value) {
	string str;
	getline(is, str, ',');
	v = boost::lexical_cast<T>(str);
    }
    return is;
}

template <typename T, size_t size>
static ostream& operator<<(ostream& os, boost::array<T, size> const& value)
{
    foreach (T const& v, value) {
	if (&v != &value.front()) {
	    os << ',';
	}
	os << v;
    }
    return os;
}

/**
 * extract comma-separated option values into variable-size array
 */
template <typename T>
static istream& operator>>(istream& is, boost::multi_array<T, 1>& value)
{
    std::vector<T> v;
    string str;
    while (!is.eof()) {
	getline(is, str, ',');
	v.push_back(boost::lexical_cast<T>(str));
    }
    boost::array<size_t, 1> extents = boost::assign::list_of(v.size());
    value.resize(extents);
    std::copy(v.begin(), v.end(), value.begin());
    return is;
}

template <typename T>
static ostream& operator<<(ostream& os, boost::multi_array<T, 1> const& value)
{
    foreach (T const& v, value) {
	if (&v != &(*value.begin())) {
	    os << ',';
	}
	os << v;
    }
    return os;
}

} // namespace std

namespace ljgpu
{

/**
 * parse program option values
 */
void options::parse(int argc, char** argv)
{
    using namespace std;

    po::options_description desc("Program options");
    desc.add_options()
	("output,o",
	 po::value<string>()->default_value(PROGRAM_NAME "_%Y%m%d_%H%M%S"),
	 "output file prefix")
	("input,I", po::value<vector<string> >(),
	 "parameter input file")
	("processor,P", po::value<vector<int> >(),
	 "CPU core(s)")
	("daemon,b", po::bool_switch(),
	 "run program in background")
	("verbose,v", po::accum_value<int>()->default_value(0),
	 "increase verbosity")
	("version",
	 "output version and exit")
	("help",
	 "display this help and exit")
	("backend",
#ifdef WITH_CUDA
	 po::value<string>()->default_value("gpu_neighbour"),
#else
	 po::value<string>()->default_value("host"),
#endif
	 "MD simulation backend")
	;

    try {
	po::command_line_parser parser(argc, argv);
	po::parsed_options parsed(parser.options(desc).allow_unregistered().run());
	po::store(parsed, vm);
	unparsed = po::collect_unrecognized(parsed.options, po::include_positional);
    }
    catch (exception const& e) {
	cerr << PROGRAM_NAME ": " << e.what() << "\n";
	cerr << "Try `" PROGRAM_NAME " --help' for more information.\n";
	throw options::exit_exception(EXIT_FAILURE);
    }

    if (vm.count("help")) {
	cout << "Usage: " PROGRAM_NAME " [OPTION]...\n" << desc << "\n";
    }

    if (vm.count("version")) {
	cout << PROGRAM_NAME " " PROGRAM_VERSION "\n"
	    "variant: " PROGRAM_VARIANT "\n"
	    "\n" PROGRAM_COPYRIGHT "\n" "This is free software. "
	    "You may redistribute copies of it under the terms of\n"
	    "the GNU General Public License "
	    "<http://www.gnu.org/licenses/gpl.html>.\n"
	    "There is NO WARRANTY, to the extent permitted by law.\n";
	throw options::exit_exception(EXIT_SUCCESS);
    }
}

/**
 * parse backend option values
 */
void options::parse(po::options_description const& opt)
{
    using namespace std;

    try {
	po::command_line_parser parser(unparsed);
	po::store(parser.options(opt).run(), vm);

	// parse optional parameter input files
	if (vm.count("input")) {
	    foreach (string const& fn, vm["input"].as<vector<string> >()) {
		ifstream ifs(fn.c_str());
		if (ifs.fail()) {
		    cerr << PROGRAM_NAME ": could not open parameter input file '" << fn << "'\n";
		    throw options::exit_exception(EXIT_FAILURE);
		}

		po::store(po::parse_config_file(ifs, opt), vm);
	    }
	}
    }
    catch (exception const& e) {
	cerr << PROGRAM_NAME ": " << e.what() << "\n";
	cerr << "Try `" PROGRAM_NAME " --help' for more information.\n";
	throw options::exit_exception(EXIT_FAILURE);
    }

    po::notify(vm);

    if (vm.count("help")) {
	cout << opt << "\n";
	throw options::exit_exception(EXIT_SUCCESS);
    }

    try {
	po::conflicting_options(vm, "density", "box-length");
	po::conflicting_options(vm, "steps", "time");
	po::dependent_option(vm, "trajectory-sample", "trajectory");
	po::dependent_option(vm, "discard-velocities", "trajectory-sample");
	po::dependent_option(vm, "thermostat", "temperature");
	po::conflicting_options(vm, "binary", "particles");
    }
    catch (exception const& e) {
	cerr << PROGRAM_NAME ": " << e.what() << "\n";
	throw options::exit_exception(EXIT_FAILURE);
    }

    // override const operator[] in variables_map
    map<string, po::variable_value>& vm_(vm);

    // optionally read parameters from HDF5 input file
    if (vm.count("trajectory")) {
	// store absolute input file path
	boost::filesystem::path path(vm["trajectory"].as<string>());
	vm_["trajectory"] = po::variable_value(boost::filesystem::complete(path).string(), false);

	try {
	    H5XX_NO_AUTO_PRINT(H5::Exception);
	    H5::H5File file(vm["trajectory"].as<string>(), H5F_ACC_RDONLY);
	    H5::Group param(file.openGroup("param"));

	    H5::Group node(param.openGroup("mdsim"));
	    po::store(po::parse_attribute<unsigned int>(node, "particles"),
		      vm_["particles"]);
	    po::store(po::parse_attribute<boost::array<unsigned int, 2> >(node, "particles"),
		      vm_["binary"]);
	    po::store(po::parse_attribute<float>(node, "density"),
		      vm_["density"]);
	    po::store(po::parse_attribute<float>(node, "box_length"),
		      vm_["box-length"]);
	    po::store(po::parse_attribute<boost::array<float, 3> >(node, "potential_epsilon"),
		      vm_["epsilon"]);
	    po::store(po::parse_attribute<boost::array<float, 3> >(node, "potential_sigma"),
		      vm_["sigma"]);
	    po::store(po::parse_attribute<float>(node, "cutoff_radius"),
		      vm_["cutoff"]);
	    po::store(po::parse_attribute<float>(node, "potential_smoothing"),
		      vm_["smooth"]);
	    po::store(po::parse_attribute<float>(node, "timestep"),
		      vm_["timestep"]);
	    po::store(po::parse_attribute<unsigned int>(node, "threads"),
		      vm_["threads"]);
	    po::store(po::parse_attribute<float>(node, "temperature"),
		      vm_["temperature"]);
	    po::store(po::parse_attribute<float>(node, "cell_occupancy"),
		      vm_["cell-occupancy"]);
	    po::store(po::parse_attribute<float>(node, "neighbour_skin"),
		      vm_["skin"]);
	    po::store(po::parse_attribute<float>(node, "pair_separation"),
		      vm_["pair-separation"]);

	    node = param.openGroup("correlation");
	    po::store(po::parse_attribute<uint64_t>(node, "steps"),
		      vm_["steps"]);
	    po::store(po::parse_attribute<float>(node, "time"),
		      vm_["time"]);
	    po::store(po::parse_attribute<unsigned int>(node, "sample_rate"),
		      vm_["sample-rate"]);
	    po::store(po::parse_attribute<unsigned int>(node, "block_size"),
		      vm_["block-size"]);
	    po::store(po::parse_attribute<uint64_t>(node, "max_samples"),
		      vm_["max-samples"]);
	    po::store(po::parse_attribute<boost::multi_array<float, 1> >(node, "q_values"),
		      vm_["q-values"]);
	    po::store(po::parse_attribute<float>(node, "q_error"),
		      vm_["q-error"]);
	}
	catch (H5::Exception const& e) {
	    cerr << PROGRAM_NAME ": " << "failed to read parameters from HDF5 input file\n";
	    throw options::exit_exception(EXIT_FAILURE);
	}
    }

    // format timestamp in output file prefix
    boost::filesystem::path path(date_time::format(vm["output"].as<string>()));
    // store absolute output file path
    vm_["output"] = po::variable_value(boost::filesystem::complete(path).string(), false);
}

options_description<mdsim_impl>::options_description()
    : po::options_description("MD simulation options")
{
    add_options()
	("particles,N", po::value<unsigned int>()->default_value(1000),
	 "number of particles")
	("binary,M", po::value<boost::array<unsigned int, 2> >(),
	 "binary mixture with A,B particles")
	("dimension", po::value<int>()->default_value(3),
	 "positional coordinates dimension")
	("density,d", po::value<float>()->default_value(0.75),
	 "particle density")
	("box-length,L", po::value<float>(),
	 "simulation box length")
	("timestep,h", po::value<float>()->default_value(0.001),
	 "simulation timestep")
	("random-seed", po::value<unsigned int>(),
	 "random number generator integer seed")
	;

    po::options_description mdsim_desc;
    mdsim_desc.add_options()
	("temperature,K", po::value<float>()->default_value(1.12),
	 "Boltzmann distribution temperature")
	("steps,s", po::value<uint64_t>()->default_value(10000),
	 "number of simulation steps")
	("time,t", po::value<float>(),
	 "total simulation time")
	("trajectory,J", po::value<std::string>(),
	 "trajectory input file")
	("trajectory-sample,S", po::value<int64_t>(),
	 "trajectory sample for initial state")
	("discard-velocities", po::bool_switch(),
	 "reset velocities to Boltzmann distribution")
	;
    add(mdsim_desc);

    po::options_description tcf_desc;
    tcf_desc.add_options()
	("sample-rate", po::value<unsigned int>()->default_value(1),
	 "sample rate for lowest block level")
	("block-size", po::value<unsigned int>()->default_value(10),
	 "block size")
	("max-samples", po::value<uint64_t>()->default_value(10000),
	 "maximum number of samples per block")
	("q-values", po::value<boost::multi_array<float, 1> >(),
	 "wave vector value(s) for correlation functions")
	("q-error", po::value<float>()->default_value(0.001),
	 "relative deviation of averaging wave vectors")
	;
    add(tcf_desc);

    po::options_description desc;
    desc.add_options()
	("disable-correlation", po::bool_switch(),
	 "disable correlation functions")
	("disable-energy", po::bool_switch(),
	 "disable thermal equilibrium properties")
	("enable-trajectory", po::bool_switch(),
	 "dump particle trajectories")
	("dry-run,n", po::bool_switch(),
	 "test parameters")
	;
    add(desc);
}

options_description<ljfluid_impl_base>::options_description()
{
    using namespace boost::assign;
    add_options()
	("cutoff", po::value<float>()->default_value(2.5),
	 "truncate potential at cutoff radius")
	("smooth", po::value<float>(),
	 "C²-potential smoothing factor")
	("thermostat", po::value<float>(),
	 "heat bath collision probability")
	("epsilon", po::value<boost::array<float, 3> >()->default_value(list_of(1.0f)(1.5f)(0.5f)),
	 "potential well depths AA,AB,BB")
	("sigma", po::value<boost::array<float, 3> >()->default_value(list_of(1.0f)(0.8f)(0.88f)),
	 "collision diameters AA,AB,BB")
	;
}

options_description<ljfluid_impl_gpu_base>::options_description()
{
    add_options()
	("tcf-backend", po::value<std::string>()->default_value("gpu"),
	 "compute correlation functions on GPU or host")
	("device,D", po::value<int>()->default_value(0),
	 "CUDA device ordinal")
	("threads,T", po::value<unsigned int>()->default_value(128),
	 "number of CUDA threads per block")
	;
}

options_description<ljfluid_impl_gpu_square>::options_description()
{
}

options_description<ljfluid_impl_gpu_neighbour>::options_description()
{
    add_options()
	("cell-occupancy", po::value<float>()->default_value(0.5),
	 "desired average cell occupancy")
	("skin", po::value<float>()->default_value(0.3),
	 "neighbour list skin")
	;
}

options_description<ljfluid_impl_gpu_cell>::options_description()
{
    add_options()
	("cell-occupancy", po::value<float>()->default_value(0.5),
	 "desired average cell occupancy")
	;
}

options_description<ljfluid_impl_host>::options_description()
{
    add_options()
	("tcf-backend", po::value<std::string>()->default_value("host"),
	 "compute correlation functions on GPU or host")
	("skin", po::value<float>()->default_value(0.3),
	 "neighbour list skin")
	;
}

options_description<hardsphere_impl>::options_description()
{
    add_options()
	("tcf-backend", po::value<std::string>()->default_value("host"),
	 "compute correlation functions on GPU or host")
	("pair-separation,p", po::value<float>()->default_value(0.5),
	 "particle pair separation")
	;
}

} // namespace ljgpu
