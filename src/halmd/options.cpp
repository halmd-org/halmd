/* Molecular Dynamics simulation program options
 *
 * Copyright © 2008-2009  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
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

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/array.hpp>
#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>
#include <fstream>
#include <iostream>
#include <map>

#include <halmd/mdsim/impl.hpp>
#include <halmd/options.hpp>
#include <halmd/util/H5xx.hpp>
#include <halmd/util/date_time.hpp>
#include <halmd/version.h>

namespace po = boost::program_options;

#define foreach BOOST_FOREACH

namespace boost { namespace program_options
{

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

namespace halmd
{

/**
 * parse program option values
 */
void options::parse(int argc, char** argv)
{
    using namespace std;
    using namespace boost::algorithm;

    po::options_description desc("Program options");
    desc.add_options()
        ("output,o",
         po::value<string>()->default_value(to_lower_copy(string(PROGRAM_NAME)) + "_%Y%m%d_%H%M%S"),
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
#ifndef BACKEND_EXECUTABLES
        ("backend",
#ifdef WITH_CUDA
         po::value<string>()->default_value("gpu_neighbour"),
#else
         po::value<string>()->default_value("host"),
#endif
         "MD simulation backend")
#endif /* ! BACKEND_EXECUTABLES */
        ;

    try {
        po::command_line_parser parser(argc, argv);
        po::parsed_options parsed(parser.options(desc).allow_unregistered().run());
        po::store(parsed, vm);
        unparsed_command_line_options = po::collect_unrecognized(parsed.options, po::include_positional);

        // parse optional parameter input files
        if (vm.count("input")) {
            foreach (string const& fn, vm["input"].as<vector<string> >()) {
                ifstream ifs(fn.c_str());
                if (ifs.fail()) {
                    cerr << PROGRAM_NAME ": could not open parameter input file '" << fn << "'\n";
                    throw options::exit_exception(EXIT_FAILURE);
                }

                po::parsed_options parsed(po::parse_config_file(ifs, desc, true));
                po::store(parsed, vm);

                // store unparsed options
                std::vector<po::basic_option<char> > unparsed;
                foreach (po::basic_option<char>& option, parsed.options) {
                    if (option.unregistered) {
                        option.unregistered = false;
                        unparsed.push_back(option);
                    }
                }
                unparsed_config_file_options.push_back(unparsed);
            }
        }
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
        cout << PROGRAM_NAME " (" PROGRAM_DESC ") " PROGRAM_VERSION "\n"
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
        po::command_line_parser parser(unparsed_command_line_options);
        po::store(parser.options(opt).run(), vm);

        // parse unrecognised options from parameter input files
        foreach (std::vector<po::basic_option<char> > const& unparsed, unparsed_config_file_options) {
            po::parsed_options parsed(&opt);
            std::copy(unparsed.begin(), unparsed.end(), std::back_inserter(parsed.options));
            po::store(parsed, vm);
        }

        po::notify(vm);
    }
    catch (exception const& e) {
        cerr << PROGRAM_NAME ": " << e.what() << "\n";
        cerr << "Try `" PROGRAM_NAME " --help' for more information.\n";
        throw options::exit_exception(EXIT_FAILURE);
    }

    if (vm.count("help")) {
        cout << opt << "\n";
        throw options::exit_exception(EXIT_SUCCESS);
    }

    try {
        po::conflicting_options(vm, "density", "box-length");
        po::conflicting_options(vm, "steps", "time");
        po::dependent_option(vm, "trajectory-sample", "trajectory");
        po::dependent_option(vm, "thermostat", "temperature");
        po::conflicting_options(vm, "binary", "particles");
        po::conflicting_options(vm, "disable-correlation", "q-values");
        po::conflicting_options(vm, "disable-correlation", "q-error");
        po::dependent_option(vm, "measure-temperature-after-time", "trajectory-sample");
        po::dependent_option(vm, "measure-temperature-after-time", "energy");
        po::dependent_option(vm, "measure-temperature-after-time", "temperature");
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
            po::store(po::parse_attribute<int>(node, "dimension"),
                      vm_["dimension"]);
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
            // parse scalar cutoff radius for backwards compatibility
            po::store(po::parse_attribute<float>(node, "cutoff_radius"),
                      vm_["cutoff"]);
            po::store(po::parse_attribute<boost::array<float, 3> >(node, "cutoff_radius"),
                      vm_["cutoff"]);
            po::store(po::parse_attribute<float>(node, "potential_smoothing"),
                      vm_["smooth"]);
            po::store(po::parse_attribute<double>(node, "timestep"),
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
            po::store(po::parse_attribute<double>(node, "time"),
                      vm_["time"]);
            po::store(po::parse_attribute<unsigned int>(node, "sample_rate"),
                      vm_["sample-rate"]);
            po::store(po::parse_attribute<unsigned int>(node, "block_size"),
                      vm_["block-size"]);
            po::store(po::parse_attribute<boost::multi_array<uint64_t, 1> >(node, "max_samples"),
                      vm_["max-samples"]);
            po::store(po::parse_attribute<uint64_t>(node, "min_samples"),
                      vm_["min-samples"]);
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

    if (vm.count("energy")) {
        // store absolute input file path
        boost::filesystem::path path(vm["energy"].as<string>());
        vm_["energy"] = po::variable_value(boost::filesystem::complete(path).string(), false);
    }

    // format timestamp in output file prefix
    boost::filesystem::path path(date_time::format(vm["output"].as<string>()));
    // store absolute output file path
    vm_["output"] = po::variable_value(boost::filesystem::complete(path).string(), false);
}

#define IMPL(x) typename mdsim_impl::impl_##x()

template <typename mdsim_impl>
options::description<mdsim_impl>::description() : po::options_description("MD simulation options")
{
    using namespace boost::assign;

    add_options()
        ("particles,N", po::value<unsigned int>()->default_value(1000),
         "number of particles")
        ("dimension", po::value<int>()->default_value(3),
         "positional coordinates dimension")
        ("density,d", po::value<float>()->default_value(0.75),
         "particle density")
        ("box-length,L", po::value<float>(),
         "simulation box length")
        ("timestep,h", po::value<double>()->default_value(0.001),
         "simulation timestep")
        ("random-seed", po::value<unsigned int>(),
         "random number generator integer seed")
        ("random-device", po::value<std::string>()->default_value("/dev/random"),
         "random number generator device")
        ;

    po::options_description mdsim;
    mdsim.add_options()
        ("temperature,K", po::value<float>()->default_value(1.12),
         "Boltzmann distribution temperature")
        ("steps,s", po::value<uint64_t>()->default_value(10000),
         "number of simulation steps")
        ("time,t", po::value<double>(),
         "total simulation time")
        ("trajectory-sample,S", po::value<int64_t>(),
         "trajectory sample for initial state")
        ("trajectory,J", po::value<std::string>(),
         "trajectory input file")
        ;
    add(mdsim);

    boost::array<uint64_t, 2> i = list_of(10000)(10000);
    boost::multi_array<uint64_t, 1> max_samples(boost::extents[i.size()]);
    max_samples.assign(i.begin(), i.end());

    po::options_description tcf;
    tcf.add_options()
        ("sample-rate", po::value<unsigned int>()->default_value(1),
         "sample rate for lowest block level")
        ("block-size", po::value<unsigned int>()->default_value(10),
         "block size")
        ("max-samples", po::value<boost::multi_array<uint64_t, 1> >()->default_value(max_samples),
         "maximum number of samples for lowest blocks")
        ("min-samples", po::value<uint64_t>()->default_value(100),
         "minimum number of trajectory samples")
        ("q-values", po::value<boost::multi_array<float, 1> >(),
         "wave vector values for correlation functions")
        ("minimum-velocity-filter", po::value<boost::multi_array<float, 1> >(),
         "lower boundary for fastest particles velocity")
        ("maximum-velocity-filter", po::value<boost::multi_array<float, 1> >(),
         "upper boundary for slowest particles velocity")
        ("q-error", po::value<float>()->default_value(0.001),
         "relative deviation of averaging wave vectors")
        ;
    add(tcf);

    po::options_description misc;
    misc.add_options()
        ("disable-correlation", po::bool_switch(),
         "disable correlation functions")
        ("disable-energy", po::bool_switch(),
         "disable thermal equilibrium properties")
        ("disable-trajectory", po::bool_switch(),
         "dump only start and end trajectory sample")
        ("dry-run,n", po::bool_switch(),
         "test parameters")
        ;
    add(misc);

    if (IMPL(lennard_jones_potential)) {
        add_options()
            ("cutoff", po::value<boost::array<float, 3> >()->default_value(list_of(2.5f)(2.5f)(2.5f)),
             "truncate potential at cutoff radius")
            ("smooth", po::value<float>(),
             "C²-potential smoothing factor")
            ("binary,M", po::value<boost::array<unsigned int, 2> >(),
             "binary mixture with A,B particles")
            ("epsilon", po::value<boost::array<float, 3> >()->default_value(list_of(1.0f)(1.5f)(0.5f)),
             "potential well depths AA,AB,BB")
            ("sigma", po::value<boost::array<float, 3> >()->default_value(list_of(1.0f)(0.8f)(0.88f)),
             "collision diameters AA,AB,BB")
            ("energy", po::value<std::string>(),
             "energy input file")
            ("measure-temperature-after-time", po::value<double>(),
             "rescale velocities to temperature measured after time")
            ;
    }
    if (IMPL(thermostat)) {
        add_options()
            ("thermostat", po::value<float>(),
             "heat bath collision probability")
            ;
    }
    if (IMPL(gpu)) {
        add_options()
            ("tcf-backend", po::value<std::string>()->default_value("gpu"),
             "correlation functions backend")
#ifndef __DEVICE_EMULATION__
            ("device,D", po::value<boost::multi_array<int, 1> >(),
             "CUDA device(s)")
#endif
            ("threads,T", po::value<unsigned int>()->default_value(128),
             "number of CUDA threads per block")
            ;
    }
    if (IMPL(host)) {
        add_options()
            ("tcf-backend", po::value<std::string>()->default_value("host"),
             "correlation functions backend")
            ;
    }
    if (IMPL(fixed_size_cell_lists)) {
        add_options()
            ("cell-occupancy", po::value<float>()->default_value(0.4),
             "desired average cell occupancy")
            ;
    }
    if (IMPL(neighbour_lists)) {
        add_options()
            ("skin", po::value<float>()->default_value(0.5),
             "neighbour list skin")
            ;
    }
    if (IMPL(hardsphere_potential)) {
        add_options()
            ("pair-separation,p", po::value<float>()->default_value(0.5),
             "particle pair separation")
            ;
    }
}

#undef IMPL

// explicit instantiation
template class options::description<ljfluid_impl_gpu_square>;
template class options::description<ljfluid_impl_gpu_neighbour>;
template class options::description<ljfluid_impl_host>;
template class options::description<hardsphere_impl>;

} // namespace halmd
