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
#include <halmd/core.hpp>
#include <halmd/util/H5xx.hpp>
#include <halmd/util/exception.hpp>
#include <halmd/util/hostname.hpp>
#include <halmd/util/logger.hpp>
#include <halmd/util/timer.hpp>
#include <halmd/utility/module.hpp>
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
    try {
        vm.parse(argc, argv);
    }
    catch (halmd::options::exit_exception const& e) {
        return e.status();
    }

#ifndef NDEBUG
    // enable logging as early as possible if debugging
    halmd::logger::init(vm["output"].as<std::string>() + ".log", vm["verbose"].as<int>());
#endif

    // resolve module dependencies
    try {
        module<core>::resolve(vm);
    }
    catch (std::exception const& e) {
        cerr << PROGRAM_NAME ": " << e.what() << endl;
        return EXIT_FAILURE;
    }

    // parse module options
    try {
        vm.parse(module<>::options());
    }
    catch (halmd::options::exit_exception const& e) {
        return e.status();
    }

#ifdef NDEBUG
    // enable logging after successful option parsing if not debugging
    halmd::logger::init(vm["output"].as<std::string>() + ".log", vm["verbose"].as<int>());
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

    LOG("MD simulation backend: " << vm["backend"].as<string>());
    std::string const hostname(halmd::get_hostname());
    try {
        LOG("host name: " << halmd::get_fqdn_hostname(hostname));
    }
    catch (std::runtime_error const&) {
        LOG("host name: " << hostname);
    }
    LOG("timer resolution: " << 1.E9 * halmd::high_resolution_timer::resolution() << " ns");

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
        module<core>::fetch(vm)->run();
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
