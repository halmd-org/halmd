/* test simulation results for thermodynamic variables in the dilute limit
 *
 * Copyright © 2010  Felix Höfling
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
#include <boost/test/parameterized_test.hpp>

#include <boost/program_options.hpp>
#include <limits>
#include <map>
#include <string>
#include <utility>

#include <halmd/numeric/accumulator.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/mdsim/thermodynamics.hpp>
#include <halmd/mdsim/velocity.hpp>
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

/**
 * heat capacity from microcanonical fluctuations of kinetic energy
 * see Lebowitz, Percus, and Verlet, Phys. Rev. 153, 250 (1967) for details
 */
inline double heat_capacity(double en_kin, double variance, unsigned npart)
{
    return 1 / (2./3 - npart * variance / (en_kin * en_kin));
}

/** test Verlet integrator: 'ideal' gas without interactions (setting ε=0) */

template <int dimension>
void ideal_gas(po::options vm)
{
    using namespace boost::assign;

    typedef boost::program_options::variable_value variable_value;

    float density = 1.;
    float temp = 1.;
    float rc = .1;

    // override const operator[] in variables_map
    map<string, variable_value>& vm_(vm);
    vm_["dimension"]    = variable_value(dimension, false);
    vm_["density"]      = variable_value(density, false);
    vm_["temperature"]  = variable_value(temp, false);
    vm_["epsilon"]      = variable_value(boost::array<float, 3>(list_of(0.f)(0.f)(0.f)), false);
    vm_["cutoff"]       = variable_value(boost::array<float, 3>(list_of(rc)(rc)(rc)), false);

    // enable logging to console
    shared_ptr<logging> logger(new logging(vm));

    BOOST_TEST_MESSAGE("using backend '" << vm["backend"].as<string>() << "' in " <<
                       dimension << " dimensions");

    // set up modules
    BOOST_TEST_MESSAGE("resolve module dependencies");
    po::unparsed_options unparsed;
    modules::resolver resolver(modules::registry::graph());
    resolver.resolve<mdsim::core<dimension> >(vm, unparsed);
    resolver.resolve<mdsim::thermodynamics<dimension> >(vm, unparsed);
    modules::policy policy(resolver.graph());
    modules::factory factory(policy.graph());

    // init core module
    BOOST_TEST_MESSAGE("initialise MD simulation");
    shared_ptr<mdsim::core<dimension> > core(modules::fetch<mdsim::core<dimension> >(factory, vm));

    // prepare system with Maxwell-Boltzmann distributed velocities
    BOOST_TEST_MESSAGE("assign positions and velocities");
    core->prepare();

    // measure thermodynamic properties
    shared_ptr<mdsim::thermodynamics<dimension> >
            thermodynamics(modules::fetch<mdsim::thermodynamics<dimension> >(factory, vm));

    double vcm_limit = (vm["backend"].as<string>() == "gpu") ? eps_float : eps;
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    double en_tot = thermodynamics->en_tot();

    // microcanonical simulation run
    BOOST_TEST_MESSAGE("run MD simulation");
    uint64_t steps = 1000;
    for (uint64_t i = 0; i < steps; ++i) {
        core->mdstep();
    }

    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);
    BOOST_CHECK_CLOSE_FRACTION(en_tot, thermodynamics->en_tot(), 10 * eps);

    BOOST_CHECK_CLOSE_FRACTION(temp, (float)thermodynamics->temp(), eps_float);
    BOOST_CHECK_CLOSE_FRACTION(density, (float)thermodynamics->density(), eps_float);
    BOOST_CHECK_CLOSE_FRACTION(thermodynamics->pressure() / temp / density, 1., eps_float);
}

template <int dimension>
void thermodynamics(po::options vm)
{
    typedef boost::program_options::variable_value variable_value;

    float density = 0.3;
    float temp = 3.0;
    float rc = 4.0;

    using namespace boost::assign;

    // override const operator[] in variables_map
    map<string, variable_value>& vm_(vm);
    vm_["dimension"]    = variable_value(dimension, false);
//     vm_["force"]        = variable_value(string("power-law"), false);
//     vm_["index"]        = variable_value(48, false);
    vm_["density"]      = variable_value(density, false);
    vm_["temperature"]  = variable_value(temp, false);
    vm_["particles"]    = variable_value(4000u, false);
    vm_["verbose"]      = variable_value(2, true);
    vm_["cutoff"]       = variable_value(boost::array<float, 3>(list_of(rc)(rc)(rc)), true);

    // enable logging to console
    shared_ptr<logging> logger(new logging(vm));

    BOOST_TEST_MESSAGE("using backend '" << vm["backend"].as<string>() << "' in " <<
                       dimension << " dimensions");

    // set up modules
    BOOST_TEST_MESSAGE("resolve module dependencies");
    po::unparsed_options unparsed;
    modules::resolver resolver(modules::registry::graph());
    resolver.resolve<mdsim::core<dimension> >(vm, unparsed);
    resolver.resolve<mdsim::thermodynamics<dimension> >(vm, unparsed);
    modules::policy policy(resolver.graph());
    modules::factory factory(policy.graph());

    BOOST_TEST_MESSAGE("initialise MD simulation");
    // init core module
    shared_ptr<mdsim::core<dimension> > core(modules::fetch<mdsim::core<dimension> >(factory, vm));

    // measure thermodynamic properties
    shared_ptr<mdsim::thermodynamics<dimension> >
            thermodynamics(modules::fetch<mdsim::thermodynamics<dimension> >(factory, vm));

    // poor man's thermostat
    shared_ptr<mdsim::velocity<dimension> >
            boltzmann(modules::fetch<mdsim::velocity<dimension> >(factory, vm));

    // prepare system at given temperature, run for t*=100
    BOOST_TEST_MESSAGE("equilibrate initial state");
    core->prepare();
    uint64_t steps = static_cast<uint64_t>(round(100 / vm["timestep"].as<double>()));
    for (uint64_t i = 0; i < steps; ++i) {
        core->mdstep();
        if((i+1) % 500 == 0) {
            boltzmann->set();
        }
    }

    // take averages of fluctuating quantities,
    accumulator<double> temp_, press, en_pot;

    // equilibration run, measure temperature in second half
    for (uint64_t i = 0; i < steps; ++i) {
        core->mdstep();
        if(i > steps/2 && i % 10 == 0) {
            temp_(thermodynamics->temp());
        }
    }
    double vcm_limit = (vm["backend"].as<string>() == "gpu") ? eps_float : eps;
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    boltzmann->rescale(sqrt(temp / mean(temp_)));
    double en_tot = thermodynamics->en_tot();

    // microcanonical simulation run
    BOOST_TEST_MESSAGE("run MD simulation");
    steps = 1000;
    for (uint64_t i = 0; i < steps; ++i) {
        core->mdstep();
        if(i % 10 == 0) {
            temp_(thermodynamics->temp());
            press(thermodynamics->pressure());
            en_pot(thermodynamics->en_pot());
        }
    }

    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);
    BOOST_CHECK_CLOSE_FRACTION(en_tot, thermodynamics->en_tot(),
                               steps * 1e-10 / fabs(en_tot));

    BOOST_CHECK_CLOSE_FRACTION(temp, mean(temp_), 5e-3);
    BOOST_CHECK_CLOSE_FRACTION(density, (float)thermodynamics->density(), eps_float);

    if (dimension == 3) {
        double Cv = heat_capacity(mean(temp_), variance(temp_), vm["particles"].as<unsigned>());

        // corrections for trunctated LJ potential, see e.g. Johnsen et al. (1993)
        double press_corr = 32./9 * M_PI * pow(density, 2) * (pow(rc, -6) - 1.5) * pow(rc, -3);
        double en_corr = 8./9 * M_PI * density * (pow(rc, -6) - 3) * pow(rc, -3);

        BOOST_TEST_MESSAGE("Density: " << density);
        BOOST_TEST_MESSAGE("Temperature: " << mean(temp_) << " ± " << error_of_mean(temp_));
        BOOST_TEST_MESSAGE("Pressure: " << mean(press) << " ± " << error_of_mean(press));
        BOOST_TEST_MESSAGE("Pressure (corrected): " << mean(press) + press_corr);
        BOOST_TEST_MESSAGE("Potential energy: " << mean(en_pot) << " ± " << error_of_mean(en_pot));
        BOOST_TEST_MESSAGE("Potential energy (corrected): " << mean(en_pot) + en_corr);
        BOOST_TEST_MESSAGE("β P / ρ = " << mean(press) / mean(temp_) / density);
        BOOST_TEST_MESSAGE("β U / N = " << mean(en_pot) / mean(temp_));
        BOOST_TEST_MESSAGE("Heat capacity = " << Cv);
        BOOST_CHECK_CLOSE_FRACTION(mean(press), 1.023, 0.001);
        BOOST_CHECK_CLOSE_FRACTION(mean(en_pot), -1.673, 0.001);
    }
}

void set_default_options(halmd::po::options& vm)
{
    typedef boost::program_options::variable_value variable_value;
    using namespace boost::assign;

    // override const operator[] in variables_map
    map<string, variable_value>& vm_(vm);
    vm_["force"]        = variable_value(string("lj"), true);
    vm_["integrator"]   = variable_value(string("verlet"), true);
    vm_["particles"]    = variable_value(1000u, true);
    vm_["timestep"]     = variable_value(0.001, true);
    vm_["smooth"]       = variable_value(0.005f, true);
    vm_["density"]      = variable_value(0.4f, true);
    vm_["temperature"]  = variable_value(2.0f, true);
    vm_["verbose"]      = variable_value(1, true);
    vm_["epsilon"]      = variable_value(boost::array<float, 3>(list_of(1.0f)(1.5f)(0.5f)), true);
    vm_["sigma"]        = variable_value(boost::array<float, 3>(list_of(1.0f)(0.8f)(0.88f)), true);
    vm_["cutoff"]       = variable_value(boost::array<float, 3>(list_of(2.5f)(2.5f)(2.5f)), true);
    vm_["random-seed"]  = variable_value(42u, true);
    vm_["output"]       = variable_value(string("halmd_test"), true);
}

int init_unit_test_suite()
{
    typedef boost::program_options::variable_value variable_value;
    using namespace boost::assign;
    using namespace boost::unit_test;
    using namespace boost::unit_test::framework;

    // manually define program option values
    boost::array<po::options, 2> vm;
    for_each(vm.begin(), vm.end(), boost::bind(&set_default_options, _1));

    // parametrize specific program options
    {
        map<string, variable_value>& vm_(vm[0]);
        vm_["backend"]      = variable_value(string("host"), true);
    }
    {
        map<string, variable_value>& vm_(vm[1]);
        vm_["backend"]      = variable_value(string("gpu"), true);
    }

    test_suite* ts1 = BOOST_TEST_SUITE( "host" );

    test_suite* ts11 = BOOST_TEST_SUITE( "2d" );
    ts11->add( BOOST_PARAM_TEST_CASE( &ideal_gas<2>, vm.begin(), vm.begin() + 1 ) );
    ts11->add( BOOST_PARAM_TEST_CASE( &thermodynamics<2>, vm.begin(), vm.begin() + 1 ) );

    test_suite* ts12 = BOOST_TEST_SUITE( "3d" );
    ts12->add( BOOST_PARAM_TEST_CASE( &ideal_gas<3>, vm.begin(), vm.begin() + 1 ) );
    ts12->add( BOOST_PARAM_TEST_CASE( &thermodynamics<3>, vm.begin(), vm.begin() + 1 ) );

    ts1->add( ts11 );
    ts1->add( ts12 );

    test_suite* ts2 = BOOST_TEST_SUITE( "gpu" );

    test_suite* ts21 = BOOST_TEST_SUITE( "2d" );
    ts21->add( BOOST_PARAM_TEST_CASE( &ideal_gas<2>, vm.begin() + 1, vm.end() ) );
    ts21->add( BOOST_PARAM_TEST_CASE( &thermodynamics<2>, vm.begin() + 1, vm.end() ) );

    test_suite* ts22 = BOOST_TEST_SUITE( "3d" );
    ts22->add( BOOST_PARAM_TEST_CASE( &ideal_gas<3>, vm.begin() + 1, vm.end() ) );
    ts22->add( BOOST_PARAM_TEST_CASE( &thermodynamics<3>, vm.begin() + 1, vm.end() ) );

    ts2->add( ts21 );
    ts2->add( ts22 );

    master_test_suite().add( ts1 );
    master_test_suite().add( ts2 );

    return 0;
}

static int _dummy = init_unit_test_suite();
