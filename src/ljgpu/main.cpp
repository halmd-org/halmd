/* Molecular Dynamics simulation of a Lennard-Jones fluid
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

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#ifdef WITH_CUDA
# include <cuda_wrapper.hpp>
#endif
#include <exception>
#include <iostream>
#include <libgen.h>
#include <ljgpu/mdlib.hpp>
#include <ljgpu/options.hpp>
#include <ljgpu/util/H5xx.hpp>
#include <ljgpu/util/exception.hpp>
#include <ljgpu/util/log.hpp>
#include <ljgpu/version.h>
#include <sched.h>
using namespace boost;
using namespace std;

int main(int argc, char **argv)
{
#ifdef NDEBUG
    // turns off the automatic error printing from the HDF5 library
    H5::Exception::dontPrint();
#endif

    // parse program options
    ljgpu::options opt;
    try {
        opt.parse(argc, argv);
    }
    catch (ljgpu::options::exit_exception const& e) {
        return e.status();
    }

    // load backend library
    ljgpu::mdlib mdlib;
#ifndef STATIC_BACKEND
    try {
        string const backend(opt["backend"].as<string>());
        boost::filesystem::path exe(argv[0]);
        boost::filesystem::path lib("libljgpu_" + backend + ".so");
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
    catch (ljgpu::options::exit_exception const& e) {
        return e.status();
    }

    ljgpu::log::init(opt["output"].as<std::string>() + ".log", opt["verbose"].as<int>());

    LOG(PROGRAM_NAME " " PROGRAM_VERSION);
    LOG("variant: " << mdlib.variant());
#ifndef NDEBUG
    LOG_WARNING("built with enabled debugging");
#endif
#ifdef __DEVICE_EMULATION__
    LOG_WARNING("built with device emulation");
#endif

    // print command line
    vector<string> cmd(argv, argv + argc);
    LOG("command line: " << boost::algorithm::join(cmd, " "));

    LOG("MD simulation backend: " << mdlib.backend());

    int status_ = ljgpu::LJGPU_EXIT_SUCCESS;
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
        status_ = mdlib.mdsim(opt);
#ifdef NDEBUG
    }
#ifdef WITH_CUDA
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        LOG_WARNING(PROGRAM_NAME " aborted");
        return ljgpu::LJGPU_EXIT_CUDA_ERROR;
    }
#ifndef __DEVICE_EMULATION__
    catch (cuda::driver::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        LOG_WARNING(PROGRAM_NAME " aborted");
        return ljgpu::LJGPU_EXIT_CUDA_ERROR;
    }
#endif /* ! __DEVICE_EMULATION__ */
#endif /* WITH_CUDA */
    catch (std::exception const& e) {
        LOG_ERROR(e.what());
        LOG_WARNING(PROGRAM_NAME " aborted");
        return ljgpu::LJGPU_EXIT_EXCEPTION;
    }
#endif /* NDEBUG */

    LOG(PROGRAM_NAME " exit");
    return status_;
}
