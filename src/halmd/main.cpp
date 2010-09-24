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
#include <exception>
#include <iostream>
#include <libgen.h>
#include <sched.h>

#ifdef WITH_CUDA
# include <cuda_wrapper.hpp>
#endif
#include <H5xx.hpp>

#include <halmd/deprecated/util/exception.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/main.hpp>
#include <halmd/utility/timer.hpp>
#include <halmd/utility/hostname.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/modules/factory.hpp>
#include <halmd/utility/modules/policy.hpp>
#include <halmd/utility/modules/resolver.hpp>
#include <halmd/utility/modules/writer.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace halmd;
using namespace std;

int main(int argc, char **argv)
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif

    // parse program options
    po::options vm;
    po::unparsed_options unparsed;
    try {
        po::parse_options(argc, argv, vm, unparsed);
    }
    catch (po::options_parser_error const& e) {
        return e.status();
    }

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
        resolver.resolve<halmd::main>(vm, unparsed);
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

    LOG(PROGRAM_NAME " (" PROGRAM_DESC ") " PROGRAM_VERSION);
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
    LOG("timer resolution: " << 1.E9 * timer::elapsed_min() << " ns");

    int status_ = halmd::HALMD_EXIT_SUCCESS;
#ifdef NDEBUG
    try {
#endif
        // bind process to CPU core(s)
        if (!vm["processor"].empty()) {
            cpu_set_t mask;
            CPU_ZERO(&mask);
            BOOST_FOREACH(int cpu, vm["processor"].as<std::vector<int> >()) {
                LOG("adding CPU core " << cpu << " to process CPU affinity mask");
                CPU_SET(cpu, &mask);
            }
            if (0 != sched_setaffinity(getpid(), sizeof(cpu_set_t), &mask)) {
                throw std::logic_error("failed to set process CPU affinity mask");
            }
        }

        // run MD simulation
        shared_ptr<halmd::main> script(modules::fetch<halmd::main>(factory, vm));
        script->load_wrapper();
        script->load_library();
        script->run();
#ifdef NDEBUG
    }
#ifdef WITH_CUDA
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        LOG_WARNING(PROGRAM_NAME " aborted");
        return halmd::HALMD_EXIT_CUDA_ERROR;
    }
#ifndef __DEVICE_EMULATION__
    catch (cuda::driver::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        LOG_WARNING(PROGRAM_NAME " aborted");
        return halmd::HALMD_EXIT_CUDA_ERROR;
    }
#endif /* ! __DEVICE_EMULATION__ */
#endif /* WITH_CUDA */
    catch (std::exception const& e) {
        LOG_ERROR(e.what());
        LOG_WARNING(PROGRAM_NAME " aborted");
        return halmd::HALMD_EXIT_EXCEPTION;
    }
#endif /* NDEBUG */

    LOG(PROGRAM_NAME " exit");
    return status_;
}
