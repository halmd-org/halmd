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
#include <halmd/observables/thermodynamics.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/utility/module.hpp>
#include <halmd/utility/modules/factory.hpp>
#include <halmd/utility/modules/policy.hpp>
#include <halmd/utility/modules/resolver.hpp>
#include <halmd/options.hpp>
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

void set_default_options(halmd::po::variables_map& vm);

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
void ideal_gas(po::variables_map vm)
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
    modules::resolver resolver(modules::registry::graph());
    resolver.resolve<mdsim::core<dimension> >(vm);
    resolver.resolve<observables::thermodynamics<dimension> >(vm);
    modules::policy policy(resolver.graph());
    modules::factory factory(policy.graph());

    // init core module
    BOOST_TEST_MESSAGE("initialise simulation modules");
    shared_ptr<mdsim::core<dimension> > core(modules::fetch<mdsim::core<dimension> >(factory, vm));

    // prepare system with Maxwell-Boltzmann distributed velocities
    BOOST_TEST_MESSAGE("assign positions and velocities");
    core->prepare();

    // measure thermodynamic properties
    shared_ptr<observables::thermodynamics<dimension> >
            thermodynamics(modules::fetch<observables::thermodynamics<dimension> >(factory, vm));

    double vcm_limit = (vm["backend"].as<string>() == "gpu") ? 0.1 * eps_float : eps;
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    double en_tot = thermodynamics->en_tot();

    // microcanonical simulation run
    BOOST_TEST_MESSAGE("run NVE simulation");
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
void thermodynamics(po::variables_map vm)
{
    typedef boost::program_options::variable_value variable_value;

    float density = 0.3;
    float temp = 3.0;
    float rc = 4.0;
    double timestep = 0.01;    // start with small time step for thermalisation
    unsigned npart = (vm["backend"].as<string>() == "gpu") ? 4000 : 1000;

    using namespace boost::assign;

    // override const operator[] in variables_map
    map<string, variable_value>& vm_(vm);
    vm_["dimension"]    = variable_value(dimension, false);
    vm_["density"]      = variable_value(density, false);
    vm_["temperature"]  = variable_value(temp, false);
    vm_["timestep"]     = variable_value(timestep, false);
    vm_["particles"]    = variable_value(npart, false);
//     vm_["verbose"]      = variable_value(int(logging::info), true);
    vm_["cutoff"]       = variable_value(boost::array<float, 3>(list_of(rc)(rc)(rc)), false);

    // enable logging to console
    shared_ptr<logging> logger(new logging(vm));

    BOOST_TEST_MESSAGE("using backend '" << vm["backend"].as<string>() << "' in " <<
                       dimension << " dimensions");

    // set up modules
    BOOST_TEST_MESSAGE("resolve module dependencies");
    modules::resolver resolver(modules::registry::graph());
    resolver.resolve<mdsim::core<dimension> >(vm);
    resolver.resolve<observables::thermodynamics<dimension> >(vm);
    modules::policy policy(resolver.graph());
    modules::factory factory(policy.graph());

    BOOST_TEST_MESSAGE("initialise simulation modules");
    // init core module and all dependencies
    shared_ptr<mdsim::core<dimension> >
        core(modules::fetch<mdsim::core<dimension> >(factory, vm));

    // keep a handle on particle configuration
    // thus, we can destroy and reload core with different options
    shared_ptr<mdsim::particle<dimension> >
        particle(modules::fetch<mdsim::particle<dimension> >(factory, vm));

    // prepare system at given temperature, run for t*=50, couple every Δt*=0.2
    BOOST_TEST_MESSAGE("thermalise initial state at T=" << temp);
    core->prepare();
    uint64_t steps = static_cast<uint64_t>(ceil(50 / timestep));
    uint64_t period = static_cast<uint64_t>(round(.2 / timestep));
    for (uint64_t i = 0; i < steps; ++i) {
        core->mdstep();
        if((i+1) % period == 0) {
            core->velocity->set();
        }
    }

    // reload core with different timestep, keep particle configuration
    core.reset();
    timestep = 0.001;
    vm_["timestep"] = variable_value(timestep, false);
    core = modules::fetch<mdsim::core<dimension> >(factory, vm);
    BOOST_CHECK_EQUAL(core->integrator->timestep(), timestep);

    // measure thermodynamic properties
    shared_ptr<observables::thermodynamics<dimension> >
            thermodynamics(modules::fetch<observables::thermodynamics<dimension> >(factory, vm));

    // take averages of fluctuating quantities,
    accumulator<double> temp_, press, en_pot;

    // equilibration run, measure temperature in second half
    BOOST_TEST_MESSAGE("equilibrate initial state");
    steps = static_cast<uint64_t>(ceil(30 / timestep));
    period = static_cast<uint64_t>(round(.01 / timestep));
    for (uint64_t i = 0; i < steps; ++i) {
        core->mdstep();
        if(i > steps/2 && i % period == 0) {
            temp_(thermodynamics->temp());
        }
    }
    double vcm_limit = (vm["backend"].as<string>() == "gpu") ? 0.1 * eps_float : 20 * eps;
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    double scale = sqrt(temp / mean(temp_));
    BOOST_TEST_MESSAGE("rescale velocities by factor " << scale);
    core->velocity->rescale(scale);

    double en_tot = thermodynamics->en_tot();
    double max_en_diff = 0; // maximal absolut deviation from initial total energy
    temp_.reset();

    // microcanonical simulation run
    BOOST_TEST_MESSAGE("run NVE simulation");
    for (uint64_t i = 0; i < steps; ++i) {
        core->mdstep();
        if(i % period == 0) {
            temp_(thermodynamics->temp());
            press(thermodynamics->pressure());
            en_pot(thermodynamics->en_pot());
            max_en_diff = max(abs(thermodynamics->en_tot() - en_tot), max_en_diff);
        }
    }
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    // with the first released version of halmd (commit f5283a2),
    // an energy drift of less than 5e-6 ε was obtained over 2e8 MD steps
    // using a smoothed potential (dt*=0.001, h=0.005)
    double en_limit = max(2e-5, steps * 1e-12);
    BOOST_CHECK_SMALL(max_en_diff / fabs(en_tot), en_limit);

    BOOST_CHECK_CLOSE_FRACTION(temp, mean(temp_), 2e-3);
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
        // values from Johnson et al.: P = 1.023, Epot = -1.673  (Npart = 864)
        // values from RFA theory: P = 1.0245, Epot = -1.6717
        BOOST_CHECK_CLOSE_FRACTION(mean(press), 1.023, 3e-3);
        BOOST_CHECK_CLOSE_FRACTION(mean(en_pot), -1.673, 2e-3);
    }
}

void set_default_options(halmd::po::variables_map& vm)
{
    typedef boost::program_options::variable_value variable_value;
    using namespace boost::assign;

    // override const operator[] in variables_map
    map<string, variable_value>& vm_(vm);
    vm_["force"]        = variable_value(string("lj"), true);
    vm_["integrator"]   = variable_value(string("verlet"), true);
    vm_["velocity"]     = variable_value(string("boltzmann"), true);
    vm_["position"]     = variable_value(string("lattice"), true);
    vm_["disable-state-vars"] = variable_value(false, true);
    vm_["particles"]    = variable_value(1000u, true);
    vm_["timestep"]     = variable_value(0.001, true);
    vm_["skin"]         = variable_value(0.5f, true);
    // smoothing modifies the equation of state
//     vm_["smooth"]       = variable_value(0.005f, true);
    vm_["density"]      = variable_value(0.4f, true);
    vm_["temperature"]  = variable_value(2.0f, true);
    vm_["verbose"]      = variable_value(int(logging::warning), true);
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
    boost::array<po::variables_map, 2> vm;
    for_each(vm.begin(), vm.end(), boost::bind(&set_default_options, _1));

    // parametrize specific program options
    {
        map<string, variable_value>& vm_(vm[0]);
        vm_["backend"]      = variable_value(string("host"), true);
    }
#ifdef WITH_CUDA
    {
        map<string, variable_value>& vm_(vm[1]);
        vm_["backend"]      = variable_value(string("gpu"), true);
        vm_["threads"]      = variable_value(128u, true);
        vm_["cell-occupancy"] = variable_value(0.4f, true);
        vm_["random-threads"] = variable_value(32u << DEVICE_SCALE, true);
        vm_["random-blocks"] = variable_value(32u, true);
    }
#endif /* WITH_CUDA */

    test_suite* ts1 = BOOST_TEST_SUITE( "host" );

    test_suite* ts11 = BOOST_TEST_SUITE( "2d" );
    ts11->add( BOOST_PARAM_TEST_CASE( &ideal_gas<2>, vm.begin(), vm.begin() + 1 ) );
    ts11->add( BOOST_PARAM_TEST_CASE( &thermodynamics<2>, vm.begin(), vm.begin() + 1 ) );

    test_suite* ts12 = BOOST_TEST_SUITE( "3d" );
    ts12->add( BOOST_PARAM_TEST_CASE( &ideal_gas<3>, vm.begin(), vm.begin() + 1 ) );
    ts12->add( BOOST_PARAM_TEST_CASE( &thermodynamics<3>, vm.begin(), vm.begin() + 1 ) );

    ts1->add( ts11 );
    ts1->add( ts12 );

#ifdef WITH_CUDA
    test_suite* ts2 = BOOST_TEST_SUITE( "gpu" );

    test_suite* ts21 = BOOST_TEST_SUITE( "2d" );
    ts21->add( BOOST_PARAM_TEST_CASE( &ideal_gas<2>, vm.begin() + 1, vm.end() ) );
    ts21->add( BOOST_PARAM_TEST_CASE( &thermodynamics<2>, vm.begin() + 1, vm.end() ) );

    test_suite* ts22 = BOOST_TEST_SUITE( "3d" );
    ts22->add( BOOST_PARAM_TEST_CASE( &ideal_gas<3>, vm.begin() + 1, vm.end() ) );
    ts22->add( BOOST_PARAM_TEST_CASE( &thermodynamics<3>, vm.begin() + 1, vm.end() ) );

    ts2->add( ts21 );
    ts2->add( ts22 );
#endif /* WITH_CUDA */

    master_test_suite().add( ts1 );
#ifdef WITH_CUDA
    master_test_suite().add( ts2 );
#endif /* WITH_CUDA */

    return 0;
}

static int _dummy = init_unit_test_suite();
