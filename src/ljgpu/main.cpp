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

#include <H5Cpp.h>
#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <cuda_wrapper.hpp>
#include <exception>
#include <iostream>
#include <libgen.h>
#include <ljgpu/mdlib.hpp>
#include <ljgpu/options.hpp>
#include <ljgpu/util/log.hpp>
#include <ljgpu/version.h>
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

    string const backend(opt["backend"].as<string>());
    ljgpu::mdlib mdlib;
    try {
	boost::filesystem::path exe(argv[0]);
	boost::filesystem::path lib("libljgpu_" + backend + ".so");
	mdlib.open((exe.parent_path() / lib));
    }
    catch (std::exception const& e) {
	cerr << e.what() << endl;
	return EXIT_FAILURE;
    }
    if (mdlib.version() != PROGRAM_VERSION) {
	cerr << PROGRAM_NAME << ": mismatching program and backend version" << endl;
	return EXIT_FAILURE;
    }

    // parse backend options
    try {
	opt.parse(mdlib.options());
    }
    catch (ljgpu::options::exit_exception const& e) {
	return e.status();
    }

    ljgpu::log::init(opt["output"].as<std::string>() + ".log", opt["verbose"].as<int>());

    LOG(PROGRAM_NAME " " PROGRAM_VERSION);
    LOG("variant: " << PROGRAM_VARIANT);
#ifndef NDEBUG
    LOG_WARNING("built with enabled debugging");
#endif

    // print command line
    vector<string> cmd(argv, argv + argc);
    LOG("command line: " << boost::algorithm::join(cmd, " "));

    LOG("MD simulation backend: " << backend);

#ifdef NDEBUG
    try {
#endif
	// run MD simulation
	mdlib.mdsim(opt);
#ifdef NDEBUG
    }
    catch (cuda::error const& e) {
	LOG_ERROR("CUDA: " << e.what());
	LOG_WARNING(PROGRAM_NAME " aborted");
	return EXIT_FAILURE;
    }
    catch (std::exception const& e) {
	LOG_ERROR(e.what());
	LOG_WARNING(PROGRAM_NAME " aborted");
	return EXIT_FAILURE;
    }
#endif /* NDEBUG */

    LOG(PROGRAM_NAME " exit");
    return EXIT_SUCCESS;
}
