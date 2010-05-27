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
#include <limits>
#include <map>
#include <string>
#include <utility>

#include <halmd/mdsim/core.hpp>
#include <halmd/mdsim/thermodynamics.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/io/logger.hpp>

/** reference values for the state variables of a 3-dim. LJ fluid can be found in
 *  L. Verlet, Phys. Rev. 159, 98 (1967)
 *  and in Hansen & McDonald, Table 4.2
 */

void set_default_options(halmd::options& vm);

const double eps = std::numeric_limits<double>::epsilon();
const float eps_float = std::numeric_limits<float>::epsilon();

const int dim = 3;

BOOST_AUTO_TEST_CASE( thermodynamics )
{
    using namespace halmd;

    typedef boost::program_options::variable_value variable_value;

    // manually define program option values
    options vm;
    set_default_options(vm);

    float density = 0.4;
    float temp = 1.424;

    vm["density"]       = variable_value(density, false);
    vm["temperature"]   = variable_value(temp, false);
    vm["dimension"]     = variable_value(dim, false);

    BOOST_TEST_MESSAGE("use backend " << vm["backend"].as<std::string>());

    BOOST_TEST_MESSAGE("enable logging to console");
    io::logger::init(vm);

    // set up modules
    BOOST_TEST_MESSAGE("resolve module dependencies");
    module<mdsim::core<dim> >::required(vm);
    module<mdsim::thermodynamics<dim> >::required(vm);

    BOOST_TEST_MESSAGE("initialise MD simulation");
    // init core module
    shared_ptr<mdsim::core<dim> > core(module<mdsim::core<dim> >::fetch(vm));
    core->init();

    // do some equilibration
    BOOST_TEST_MESSAGE("equilibrate initial state");
    for (uint64_t i = 0; i < 1000; ++i) {
        core->mdstep();
    }

    // measure thermodynamic properties
    shared_ptr<mdsim::thermodynamics<dim> >
            thermodynamics(module<mdsim::thermodynamics<dim> >::fetch(vm));

    BOOST_CHECK_CLOSE_FRACTION(temp,    (float)thermodynamics->temp(), eps_float);
    BOOST_CHECK_CLOSE_FRACTION(density, (float)thermodynamics->density(), eps_float);
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), eps);
    double en_tot = thermodynamics->en_tot();

    // run simulation
    BOOST_TEST_MESSAGE("run MD simulation");
    for (uint64_t i = 0; i < core->steps(); ++i) {
        core->mdstep();
    }

    double press = thermodynamics->pressure();
    double en_pot = thermodynamics->en_pot();

    BOOST_CHECK_CLOSE_FRACTION(temp,    (float)thermodynamics->temp(), eps_float);
    BOOST_CHECK_CLOSE_FRACTION(density, (float)thermodynamics->density(), eps_float);
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), eps);
    BOOST_CHECK_CLOSE_FRACTION(en_tot, en_pot + dim * temp / 2, eps_float);

    BOOST_TEST_MESSAGE("Density: " << density);
    BOOST_TEST_MESSAGE("Temperature: " << temp);
    BOOST_TEST_MESSAGE("Pressure: " << press);
    BOOST_TEST_MESSAGE("Potential energy: " << en_pot);
    BOOST_TEST_MESSAGE("β P / ρ = " << press / temp / density);
    BOOST_TEST_MESSAGE("β U / N = " << en_pot / temp);
    BOOST_CHECK_CLOSE_FRACTION(press / temp / density, 0.38, 0.1);
    BOOST_CHECK_CLOSE_FRACTION(en_pot / temp, -2.73, 0.1);
}

void set_default_options(halmd::options& vm)
{
    typedef boost::program_options::variable_value variable_value;
    using namespace boost::assign;

    vm["backend"]       = variable_value(std::string(MDSIM_BACKEND), false);
    vm["particles"]     = variable_value(1000u, false);
    vm["steps"]         = variable_value(uint64_t(100), false);
    vm["timestep"]      = variable_value(0.001, false);
    vm["density"]       = variable_value(0.4f, false);
    vm["temperature"]   = variable_value(2.0f, false);
    vm["dimension"]     = variable_value(3, false);
    vm["verbose"]       = variable_value(1, false);
    vm["epsilon"]       = variable_value(boost::array<float, 3>(list_of(1.0f)(1.5f)(0.5f)), false);
    vm["sigma"]         = variable_value(boost::array<float, 3>(list_of(1.0f)(0.8f)(0.88f)), false);
    vm["cutoff"]        = variable_value(boost::array<float, 3>(list_of(2.5f)(2.5f)(2.5f)), false);
    vm["skin"]          = variable_value(0.5f, false);
    vm["random-seed"]   = variable_value(42u, false);
}
