/* test simulation results for thermodynamic variables in the dilute limit
 *
 * Copyright (C) 2010  Felix Höfling
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE thermodynamics
#include <boost/test/unit_test.hpp>

#include <boost/program_options.hpp>
#include <cuda_wrapper.hpp>
#include <map>
#include <string>
#include <utility>

// #include <halmd/mdlib.hpp>
#include <halmd/core.hpp>
// #include <halmd/mdsim/core.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/box.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/io/logger.hpp>

using namespace halmd;

BOOST_AUTO_TEST_CASE( thermodynamics )
{
    typedef boost::program_options::variable_value variable_value;

    // parse and define program options
    options vm;
/*    try {
        vm.parse(argc, const_cast<char**>(argv));
    }
    catch (halmd::options::exit_exception const& e) {
        BOOST_FAIL("Parsing command line options failed with exit code " << e.status());
    }
    BOOST_TEST_MESSAGE("Still alive");*/

    vm.set("time", variable_value(1e0f, false));
    vm.set("density", variable_value(0.2f, false));
    vm.set("temperature", variable_value(2.0f, false));
    vm.set("dimension", variable_value(3, true));
    vm.set("verbose", variable_value(2, true));

#ifndef NDEBUG
    // enable logging as early as possible if debugging
    io::logger::init(vm);
#endif

    // resolve module dependencies
    try {
        module<halmd::mdsim::host::box<3> >::required(vm);
//         module<halmd::mdsim::host::particle<3, double> >::required(vm);
//         module<core>::required(vm);
    }
    catch (std::exception const& e) {
        BOOST_FAIL(e.what());
    }

    // parse module options
    try {
        vm.parse(module<>::options());
    }
    catch (halmd::options::exit_exception const& e) {
        BOOST_FAIL("exception in options parser: " << e.status());
    }

#ifdef NDEBUG
    // enable logging after successful option parsing if not debugging
    io::logger::init(vm);
#endif

    BOOST_TEST_MESSAGE("Still alive");

    // run MD simulation
    shared_ptr<halmd::core> core(module<halmd::core>::fetch(vm));
    core->run();
}
