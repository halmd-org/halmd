/* Molecular Dynamics simulation of a Lennard-Jones fluid
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

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/proto/args.hpp> //< proto/matches.hpp:95: error: 'N' is not a member of 'boost::proto'
#include <exception>
#include <iostream>

#ifdef WITH_CUDA
# include <cuda_wrapper.hpp>
#endif
#include <H5xx.hpp>

#include <halmd/deprecated/util/exception.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/main.hpp>
#include <halmd/options.hpp>
#include <halmd/script.hpp>
#include <halmd/utility/timer.hpp>
#include <halmd/utility/hostname.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/modules/factory.hpp>
#include <halmd/utility/modules/policy.hpp>
#include <halmd/utility/modules/resolver.hpp>
#include <halmd/utility/modules/writer.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace halmd;
using namespace std;

int main(int argc, char **argv)
{
    static halmd::script script; //< load Lua scripting engine

    options_parser options(script.options());
    try {
        options.parse(argc, argv);
    }
    catch (exit_exception const& e) {
        return e.code();
    }
    po::variables_map vm(options.parsed());

    // FIXME split log_to_console, log_to_file

#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif


#ifndef NDEBUG
    // enable logging as early as possible if debugging
    shared_ptr<logging> logger(new logging(vm));
#endif

    // resolve module dependencies
#ifndef NDEBUG
    if (vm["verbose"].as<int>() >= logging::debug) {
        write_graphviz(vm["output"].as<string>() + "_registry.dot", modules::registry::graph());
    }
#endif
    modules::resolver resolver(modules::registry::graph());
    try {
        resolver.resolve<halmd::main>(vm);
    }
    catch (program_options::error const& e) {
        cerr << PROGRAM_NAME ": " << e.what() << endl;
        return EXIT_FAILURE;
    }
#ifndef NDEBUG
    if (vm["verbose"].as<int>() >= logging::debug) {
        write_graphviz(vm["output"].as<string>() + "_resolver.dot", resolver.graph());
    }
#endif
    modules::policy policy(resolver.graph());
#ifndef NDEBUG
    if (vm["verbose"].as<int>() >= logging::debug) {
        write_graphviz(vm["output"].as<string>() + "_policy.dot", policy.graph());
    }
#endif
    modules::factory factory(policy.graph());

#ifdef NDEBUG
    // enable logging after successful option parsing if not debugging
    shared_ptr<logging> logger(new logging(vm));
#endif

    LOG(PROJECT_NAME " (" PROGRAM_DESC ") " PROGRAM_VERSION);
    LOG("variant: " << PROGRAM_VARIANT);
#ifndef NDEBUG
    LOG_WARNING("built with enabled debugging");
#endif
#ifdef __DEVICE_EMULATION__
    LOG_WARNING("built with device emulation");
#endif

    // print command line
    std::vector<string> cmd(argv, argv + argc);
    LOG("command line: " << boost::algorithm::join(cmd, " "));

    LOG("host name: " << host_name());

    int status_ = halmd::HALMD_EXIT_SUCCESS;
#ifdef NDEBUG
    try {
#endif
        // run MD simulation
        shared_ptr<halmd::main> sampler(modules::fetch<halmd::main>(factory, vm));
        sampler->run();
#ifdef NDEBUG
    }
#ifdef WITH_CUDA
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        LOG_WARNING(PROJECT_NAME " aborted");
        return halmd::HALMD_EXIT_CUDA_ERROR;
    }
#ifndef __DEVICE_EMULATION__
    catch (cuda::driver::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        LOG_WARNING(PROJECT_NAME " aborted");
        return halmd::HALMD_EXIT_CUDA_ERROR;
    }
#endif /* ! __DEVICE_EMULATION__ */
#endif /* WITH_CUDA */
    catch (std::exception const& e) {
        LOG_ERROR(e.what());
        LOG_WARNING(PROJECT_NAME " aborted");
        return halmd::HALMD_EXIT_EXCEPTION;
    }
#endif /* NDEBUG */

    LOG(PROJECT_NAME " exit");
    return status_;
}
