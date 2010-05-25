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

#include <halmd/core.hpp>
// #include <halmd/mdsim/core.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/box.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/io/logger.hpp>

using namespace halmd;
using namespace boost::assign;

BOOST_AUTO_TEST_CASE( thermodynamics )
{
    typedef boost::program_options::variable_value variable_value;

    // manually define program option values
    options vm;
    vm["backend"]       = variable_value(std::string(MDSIM_BACKEND), false);
    vm["particles"]     = variable_value(1000u, false);
    vm["steps"]         = variable_value(uint64_t(1000), false);
    vm["timestep"]      = variable_value(0.001, false);
    vm["density"]       = variable_value(0.2f, false);
    vm["temperature"]   = variable_value(2.0f, false);
    vm["dimension"]     = variable_value(3, false);
    vm["verbose"]       = variable_value(2, false);
    vm["epsilon"]       = variable_value(boost::array<float, 3>(list_of(1.0f)(1.5f)(0.5f)), false);
    vm["sigma"]         = variable_value(boost::array<float, 3>(list_of(1.0f)(0.8f)(0.88f)), false);
    vm["cutoff"]        = variable_value(boost::array<float, 3>(list_of(2.5f)(2.5f)(2.5f)), false);
    vm["skin"]          = variable_value(0.5f, false);
    vm["random-seed"]   = variable_value(42u, false);

    BOOST_TEST_MESSAGE("enable logging to console");
    io::logger::init(vm);

    BOOST_TEST_MESSAGE("resolve module dependencies");
    module<core>::required(vm);

    BOOST_TEST_MESSAGE("initialise MD simulation");
    shared_ptr<halmd::core> core(module<halmd::core>::fetch(vm));
    BOOST_TEST_MESSAGE("run MD simulation");
    core->run();
}
