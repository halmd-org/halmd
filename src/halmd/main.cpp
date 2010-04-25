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
#include <halmd/mdlib.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/options.hpp>
#include <halmd/util/H5xx.hpp>
#include <halmd/util/exception.hpp>
#include <halmd/util/hostname.hpp>
#include <halmd/util/logger.hpp>
#include <halmd/util/timer.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace std;

int main(int argc, char **argv)
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif

    // parse program options
    halmd::options opt;
    try {
        opt.parse(argc, argv);
    }
    catch (halmd::options::exit_exception const& e) {
        return e.status();
    }

    // load backend library
    halmd::mdlib mdlib;
#ifndef BACKEND_EXECUTABLES
    try {
        string const backend(opt["backend"].as<string>());
        boost::filesystem::path exe(argv[0]);
        boost::filesystem::path lib("libhalmd_" + backend + ".so");
        mdlib.open((exe.parent_path() / lib));
    }
    catch (std::exception const& e) {
        cerr << PROGRAM_NAME ": " << e.what() << endl;
        return EXIT_FAILURE;
    }
    if (mdlib.version() != PROGRAM_VERSION) {
        cerr << PROGRAM_NAME ": mismatching program and backend version" << endl;
        return EXIT_FAILURE;
    }
#endif

    // parse backend options
    try {
        opt.parse(mdlib.options());
    }
    catch (halmd::options::exit_exception const& e) {
        return e.status();
    }

    halmd::logger::init(opt["output"].as<std::string>() + ".log", opt["verbose"].as<int>());

    LOG(PROGRAM_NAME " (" PROGRAM_DESC ") " PROGRAM_VERSION);
    LOG("variant: " << mdlib.variant());
#ifndef NDEBUG
    LOG_WARNING("built with enabled debugging");
#endif
#ifdef __DEVICE_EMULATION__
    LOG_WARNING("built with device emulation");
#endif

    // print command line
    std::vector<string> cmd(argv, argv + argc);
    LOG("command line: " << boost::algorithm::join(cmd, " "));

    LOG("MD simulation backend: " << mdlib.backend());
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
        if (!opt["processor"].empty()) {
            cpu_set_t mask;
            CPU_ZERO(&mask);
            BOOST_FOREACH(int cpu, opt["processor"].as<std::vector<int> >()) {
                LOG("adding CPU core " << cpu << " to process CPU affinity mask");
                CPU_SET(cpu, &mask);
            }
            if (0 != sched_setaffinity(getpid(), sizeof(cpu_set_t), &mask)) {
                throw std::logic_error("failed to set process CPU affinity mask");
            }
        }

        // run MD simulation
        int const dimension = opt["dimension"].as<int>();
        if (dimension == 3) {
            halmd::mdsim::core<3>::resolve(opt);
            halmd::mdsim::core<3> core(opt);
            core.run();
        }
        else if (dimension == 2) {
            halmd::mdsim::core<2>::resolve(opt);
            halmd::mdsim::core<2> core(opt);
            core.run();
        }
        else {
            throw logic_error("invalid dimension: " + lexical_cast<string>(dimension));
        }
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
