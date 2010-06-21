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

#include <halmd/math/accum.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/mdsim/thermodynamics.hpp>
#include <halmd/mdsim/host/velocity/boltzmann.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/modules/factory.hpp>
#include <halmd/utility/modules/policy.hpp>
#include <halmd/utility/modules/resolver.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/io/logger.hpp>

using namespace halmd;
using namespace std;

/** reference values for the state variables of a 3-dim. LJ fluid can be found in
 *  L. Verlet, Phys. Rev. 159, 98 (1967)
 *  and in Hansen & McDonald, Table 4.2
 *
 *  a more detailed analysis and more accurate values are found in
 *  Johnson, Zollweg, and Gubbins, Mol. Phys. 78, 591 (1993).
 */

void set_default_options(halmd::po::options& vm);

const double eps = numeric_limits<double>::epsilon();
const float eps_float = numeric_limits<float>::epsilon();

const int dim = 3;

/**
 * heat capacity from microcanonical fluctuations of kinetic energy
 * see Lebowitz, Percus, and Verlet, Phys. Rev. 153, 250 (1967) for details
 */
inline double heat_capacity(double en_kin, double variance, unsigned npart)
{
    return 1 / (2./3 - npart * variance / (en_kin * en_kin));
}

/** test Verlet integrator: 'ideal' gas without interactions (setting ε=0) */

BOOST_AUTO_TEST_CASE( ideal_gas )
{
    using namespace boost::assign;

    typedef boost::program_options::variable_value variable_value;

    // manually define program option values
    po::options vm;
    set_default_options(vm);

    float density = 1.;
    float temp = 1.;
    float rc = 0.1;

    // override const operator[] in variables_map
    map<string, variable_value>& vm_(vm);
    vm_["density"]      = variable_value(density, false);
    vm_["temperature"]  = variable_value(temp, false);
    vm_["dimension"]    = variable_value(dim, false);
    vm_["epsilon"]      = variable_value(boost::array<float, 3>(list_of(0.f)(0.f)(0.f)), false);
    vm_["cutoff"]       = variable_value(boost::array<float, 3>(list_of(rc)(rc)(rc)), false);

    // enable logging to console
    io::logger::init(vm);

    BOOST_TEST_MESSAGE("use backend " << vm["backend"].as<string>());

    // set up modules
    BOOST_TEST_MESSAGE("resolve module dependencies");
    po::unparsed_options unparsed;
    modules::resolver resolver(modules::registry::graph());
    resolver.resolve<mdsim::core<dim> >(vm, unparsed);
    resolver.resolve<mdsim::thermodynamics<dim> >(vm, unparsed);
    modules::policy policy(resolver.graph());
    modules::factory factory(policy.graph());

    BOOST_TEST_MESSAGE("initialise MD simulation");
    // init core module
    shared_ptr<mdsim::core<dim> > core(modules::fetch<mdsim::core<dim> >(factory, vm));

    // measure thermodynamic properties
    shared_ptr<mdsim::thermodynamics<dim> >
            thermodynamics(modules::fetch<mdsim::thermodynamics<dim> >(factory, vm));

    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), eps);
    double en_tot = thermodynamics->en_tot();

    // microcanonical simulation run
    BOOST_TEST_MESSAGE("run MD simulation");
    uint64_t steps = 1000;
    for (uint64_t i = 0; i < steps; ++i) {
        core->mdstep();
    }

    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), eps);
    BOOST_CHECK_CLOSE_FRACTION(en_tot, thermodynamics->en_tot(), 10 * eps);

    BOOST_CHECK_CLOSE_FRACTION(temp, (float)thermodynamics->temp(), eps_float);
    BOOST_CHECK_CLOSE_FRACTION(density, (float)thermodynamics->density(), eps_float);
    BOOST_CHECK_CLOSE_FRACTION(thermodynamics->pressure() / temp / density, 1., eps_float);
}

BOOST_AUTO_TEST_CASE( thermodynamics )
{

    typedef boost::program_options::variable_value variable_value;

    // manually define program option values
    po::options vm;
    set_default_options(vm);

    float density = 0.4;
    float temp = 2.0;
    float rc = 3.5;

    using namespace boost::assign;

    // override const operator[] in variables_map
    map<string, variable_value>& vm_(vm);
//     vm_["force"]        = variable_value(string("power-law"), false);
//     vm_["index"]        = variable_value(48, false);
    vm_["density"]      = variable_value(density, false);
    vm_["temperature"]  = variable_value(temp, false);
    vm_["dimension"]    = variable_value(dim, false);
    vm_["particles"]    = variable_value(864u, false);
//     vm_["verbose"]      = variable_value(2, true);
    vm_["cutoff"]       = variable_value(boost::array<float, 3>(list_of(rc)(rc)(rc)), true);

    // enable logging to console
    io::logger::init(vm);

    BOOST_TEST_MESSAGE("use backend " << vm["backend"].as<string>());

    // set up modules
    BOOST_TEST_MESSAGE("resolve module dependencies");
    po::unparsed_options unparsed;
    modules::resolver resolver(modules::registry::graph());
    resolver.resolve<mdsim::core<dim> >(vm, unparsed);
    resolver.resolve<mdsim::thermodynamics<dim> >(vm, unparsed);
    modules::policy policy(resolver.graph());
    modules::factory factory(policy.graph());

    BOOST_TEST_MESSAGE("initialise MD simulation");
    // init core module
    shared_ptr<mdsim::core<dim> > core(modules::fetch<mdsim::core<dim> >(factory, vm));

    // measure thermodynamic properties
    shared_ptr<mdsim::thermodynamics<dim> >
            thermodynamics(modules::fetch<mdsim::thermodynamics<dim> >(factory, vm));

    // poor man's thermostat
    shared_ptr<mdsim::host::velocity::boltzmann<dim, double> >
            boltzmann(modules::fetch<mdsim::host::velocity::boltzmann<dim, double> >(factory, vm));

    // prepare system at given temperature, run for t*=30
    BOOST_TEST_MESSAGE("equilibrate initial state");
    uint64_t steps = static_cast<uint64_t>(round(30 / vm["timestep"].as<double>()));
    for (uint64_t i = 0; i < steps; ++i) {
        core->mdstep();
        if((i+1) % 200 == 0) {
            boltzmann->set();
        }
    }

    // take averages of fluctuating quantities,
    accumulator<double> temp_, press, en_pot;

    // equilibration run, measure temperature in second half
    for (uint64_t i = 0; i < steps; ++i) {
        core->mdstep();
        if(i > steps/2 && i % 10 == 0) {
            temp_ += thermodynamics->temp();
        }
    }
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), steps * eps);

    boltzmann->rescale(sqrt(temp/temp_.mean()));
    double en_tot = thermodynamics->en_tot();

    // microcanonical simulation run
    BOOST_TEST_MESSAGE("run MD simulation");
    steps = 1000;
    for (uint64_t i = 0; i < steps; ++i) {
        core->mdstep();
        if(i % 10 == 0) {
            temp_  += thermodynamics->temp();
            press  += thermodynamics->pressure();
            en_pot += thermodynamics->en_pot();
        }
    }

    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), steps * eps);
    BOOST_CHECK_CLOSE_FRACTION(en_tot, thermodynamics->en_tot(),
                               steps * 1e-10 / fabs(en_tot));

    BOOST_CHECK_CLOSE_FRACTION(temp, temp_.mean(), 5e-3);
    BOOST_CHECK_CLOSE_FRACTION(density, (float)thermodynamics->density(), eps_float);

    double Cv = heat_capacity(temp_.mean(), temp_.var(), vm["particles"].as<unsigned>());

    // corrections for trunctated LJ potential, see e.g. Johnsen et al. (1993)
    double press_corr = 32./9 * M_PI * pow(density, 2) * (pow(rc, -6) - 1.5) * pow(rc, -3);
    double en_corr = 8./9 * M_PI * density * (pow(rc, -6) - 3) * pow(rc, -3);

    BOOST_TEST_MESSAGE("Density: " << density);
    BOOST_TEST_MESSAGE("Temperature: " << temp_.mean() << " ± " << temp_.err());
    BOOST_TEST_MESSAGE("Pressure: " << press.mean() << " ± " << press.err());
    BOOST_TEST_MESSAGE("Pressure (corrected): " << press.mean() + press_corr);
    BOOST_TEST_MESSAGE("Potential energy: " << en_pot.mean() << " ± " << en_pot.err());
    BOOST_TEST_MESSAGE("Potential energy (corrected): " << en_pot.mean() + en_corr);
    BOOST_TEST_MESSAGE("β P / ρ = " << press.mean() / temp_.mean() / density);
    BOOST_TEST_MESSAGE("β U / N = " << en_pot.mean() / temp_.mean());
    BOOST_TEST_MESSAGE("Heat capacity = " << Cv);
    BOOST_CHECK_CLOSE_FRACTION(press.mean() + press_corr, 0.70, 0.01);
    BOOST_CHECK_CLOSE_FRACTION(en_pot.mean() + en_corr, -2.54, 0.01);
}

void set_default_options(halmd::po::options& vm)
{
    typedef boost::program_options::variable_value variable_value;
    using namespace boost::assign;

    // override const operator[] in variables_map
    map<string, variable_value>& vm_(vm);
    vm_["backend"]      = variable_value(string(MDSIM_BACKEND), true);
    vm_["force"]        = variable_value(string("lj"), true);
    vm_["integrator"]   = variable_value(string("verlet"), true);
    vm_["particles"]    = variable_value(1000u, true);
    vm_["timestep"]     = variable_value(0.001, true);
    vm_["smooth"]       = variable_value(0.005f, true);
    vm_["density"]      = variable_value(0.4f, true);
    vm_["temperature"]  = variable_value(2.0f, true);
    vm_["dimension"]    = variable_value(3, true);
    vm_["verbose"]      = variable_value(0, true);
    vm_["epsilon"]      = variable_value(boost::array<float, 3>(list_of(1.0f)(1.5f)(0.5f)), true);
    vm_["sigma"]        = variable_value(boost::array<float, 3>(list_of(1.0f)(0.8f)(0.88f)), true);
    vm_["cutoff"]       = variable_value(boost::array<float, 3>(list_of(2.5f)(2.5f)(2.5f)), true);
    vm_["random-seed"]  = variable_value(42u, true);
}
