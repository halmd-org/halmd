/*
 * Copyright © 2010  Felix Höfling and Peter Colberg
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

#include <limits>
#include <map>
#include <string>
#include <utility>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <test/modules.hpp>

using namespace boost;
using namespace halmd;
using namespace halmd::test;
using namespace std;

/**
 * test simulation results for thermodynamic variables in the dilute limit
 *
 * reference values for the state variables of a 3-dim. LJ fluid can be found in
 * L. Verlet, Phys. Rev. 159, 98 (1967)
 * and in Hansen & McDonald, Table 4.2
 *
 * a more detailed analysis and more accurate values are found in
 * Johnson, Zollweg, and Gubbins, Mol. Phys. 78, 591 (1993).
 */

const double eps = numeric_limits<double>::epsilon();
const float eps_float = numeric_limits<float>::epsilon();

/**
 * heat capacity from microcanonical fluctuations of kinetic energy
 * see Lebowitz, Percus, and Verlet, Phys. Rev. 153, 250 (1967) for details
 */
inline double heat_capacity_nve(double en_kin, double variance, unsigned npart)
{
    return 1 / (2./3 - npart * variance / (en_kin * en_kin));
}

/** test Verlet integrator: 'ideal' gas without interactions (setting ε=0) */

template <int dimension>
void ideal_gas(string const& backend)
{
    using namespace boost::assign;

    float density = 1.;
    float temp = 1.;
    float rc = .1;
    double timestep = 0.001;
    unsigned int npart = 1000;
    unsigned int random_seed = 42;
    fixed_vector<double, dimension> box_ratios;
    box_ratios[0] = 1;
    box_ratios[1] = 2;
    if (dimension == 3) {
        box_ratios[2] = 1.01;
    }

    // enable logging to console
    shared_ptr<logger> log(new logger);
    log->log_to_console(
#ifdef NDEBUG
        logger::warning
#else
        logger::debug
#endif
    );

    BOOST_TEST_MESSAGE("using backend '" << backend << "' in " <<
                       dimension << " dimensions");

#ifdef WITH_CUDA
    shared_ptr<utility::gpu::device> device = make_device(backend);
#endif /* WITH_CUDA */
    shared_ptr<halmd::random::random> random = make_random(
        backend
      , random_seed
    );

    // init core module
    BOOST_TEST_MESSAGE("initialise simulation modules");
    shared_ptr<mdsim::core<dimension> > core(new mdsim::core<dimension>);
    core->particle = make_particle<dimension>(
        backend
      , npart
    );
    core->box = make_shared<mdsim::box<dimension> >(
        core->particle
      , density
      , box_ratios
    );
    core->integrator = make_verlet_integrator<dimension>(
        backend
      , core->particle
      , core->box
      , timestep
    );
    core->force = make_lennard_jones_force<dimension>(
        backend
      , core->particle
      , core->box
      , list_of(rc)(rc)(rc)      /* cutoff */
      , list_of(0.f)(0.f)(0.f)   /* epsilon */
      , list_of(1.f)(0.f)(0.f)   /* sigma */
    );
    core->neighbour = make_neighbour(
        backend
      , core->particle
      , core->box
      , core->force
    );
    core->position = make_lattice(
        backend
      , core->particle
      , core->box
      , random
    );
    core->velocity = make_boltzmann(
        backend
      , core->particle
      , random
      , temp
    );

    // prepare system with Maxwell-Boltzmann distributed velocities
    BOOST_TEST_MESSAGE("assign positions and velocities");
    core->force->aux_enable();              // enable computation of potential energy
    core->prepare();

    // measure thermodynamic properties
    shared_ptr<observables::thermodynamics<dimension> > thermodynamics = make_thermodynamics(
        backend
      , core->particle
      , core->box
      , core->force
    );

    const double vcm_limit = (backend == "gpu") ? 0.1 * eps_float : eps;
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    double en_tot = thermodynamics->en_tot();

    // microcanonical simulation run
    BOOST_TEST_MESSAGE("run NVE simulation");
    uint64_t steps = 1000;
    core->force->aux_disable();             // disable auxiliary variables
    for (uint64_t i = 0; i < steps; ++i) {
        // last step: evaluate auxiliary variables (potential energy, virial, ...)
        if (i == steps - 1) {
            core->force->aux_enable();
        }
        core->mdstep();
    }

    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);
    BOOST_CHECK_CLOSE_FRACTION(en_tot, thermodynamics->en_tot(), 10 * eps);

    BOOST_CHECK_CLOSE_FRACTION(temp, (float)thermodynamics->temp(), eps_float);
    BOOST_CHECK_CLOSE_FRACTION(density, (float)thermodynamics->density(), eps_float);
    BOOST_CHECK_CLOSE_FRACTION(thermodynamics->pressure() / temp / density, 1., eps_float);
}

template <int dimension>
void thermodynamics(string const& backend)
{
    float density = 0.3;
    float temp = 3.0;
    float rc = 4.0;
    double timestep = 0.001;
    unsigned npart = (backend == "gpu") ? 4000 : 1500;
    unsigned int random_seed = 42;
    fixed_vector<double, dimension> box_ratios;
    box_ratios[0] = 1;
    box_ratios[1] = 2;
    if (dimension == 3) {
        box_ratios[2] = 1.01;
    }

    using namespace boost::assign;

    // enable logging to console
    shared_ptr<logger> log(new logger);
    log->log_to_console(
#ifdef NDEBUG
        logger::warning
#else
        logger::debug
#endif
    );

    BOOST_TEST_MESSAGE("using backend '" << backend << "' in " <<
                       dimension << " dimensions");

#ifdef WITH_CUDA
    shared_ptr<utility::gpu::device> device = make_device(backend);
#endif /* WITH_CUDA */
    shared_ptr<halmd::random::random> random = make_random(
        backend
      , random_seed
    );

    BOOST_TEST_MESSAGE("initialise simulation modules");
    // init core module and all dependencies
    shared_ptr<mdsim::core<dimension> > core(new mdsim::core<dimension>);

    core->particle = make_particle<dimension>(backend, npart);

    core->box = make_box<dimension>(core->particle, density, box_ratios);

    core->integrator = make_verlet_nvt_andersen_integrator<dimension>(
        backend, core->particle, core->box, random
      , 0.01 /* time step */, temp, 5. /* collision rate */
    );

    core->force = make_lennard_jones_force<dimension>(
        backend, core->particle, core->box
      , list_of(rc)(rc)(rc)     /* cutoff */
      , list_of(1.f)(0.f)(0.f)  /* epsilon */
      , list_of(1.f)(0.f)(0.f)  /* sigma */
    );

    core->neighbour = make_neighbour(backend, core->particle, core->box, core->force);

    core->position = make_lattice(backend, core->particle, core->box, random);

    core->velocity = make_boltzmann(backend, core->particle, random, temp);

    // relax configuration and thermalise at given temperature, run for t*=50
    BOOST_TEST_MESSAGE("thermalise initial state at T=" << temp);
    core->prepare();
    uint64_t steps = static_cast<uint64_t>(ceil(50 / core->integrator->timestep()));
    for (uint64_t i = 0; i < steps; ++i) {
        core->mdstep();
    }

    // set different timestep and choose NVE integrator
    core->integrator = make_verlet_integrator<dimension>(
        backend, core->particle, core->box, timestep
    );
    BOOST_CHECK_CLOSE_FRACTION(core->integrator->timestep(), timestep, eps_float);

    // measure thermodynamic properties
    shared_ptr<observables::thermodynamics<dimension> > thermodynamics =
        make_thermodynamics(backend, core->particle, core->box, core->force);

    // stochastic thermostat => centre particle velocities around zero
    core->velocity->shift(-thermodynamics->v_cm());

    const double vcm_limit = (backend == "gpu") ? 0.1 * eps_float : 20 * eps;
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    // take averages of fluctuating quantities,
    accumulator<double> temp_, press, en_pot;

    // equilibration run, measure temperature in second half
    BOOST_TEST_MESSAGE("equilibrate initial state");
    steps = static_cast<uint64_t>(ceil(30 / timestep));
    uint64_t period = static_cast<uint64_t>(round(.01 / timestep));
    core->force->aux_disable();                     // disable auxiliary variables
    for (uint64_t i = 0; i < steps; ++i) {
        if (i == steps - 1) {
            core->force->aux_enable();              // enable auxiliary variables in last step
        }
        core->mdstep();
        if(i > steps/2 && i % period == 0) {
            temp_(thermodynamics->temp());
        }
    }
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    // rescale velocities to match the exact temperature
    double v_scale = sqrt(temp / mean(temp_));
    BOOST_TEST_MESSAGE("rescale velocities by factor " << v_scale);
    core->velocity->rescale(v_scale);

    double en_tot = thermodynamics->en_tot();
    double max_en_diff = 0; // maximal absolut deviation from initial total energy
    temp_.reset();

    // microcanonical simulation run
    BOOST_TEST_MESSAGE("run NVE simulation");
    core->force->aux_disable();
    for (uint64_t i = 0; i < steps; ++i) {
        // turn on evaluation of potential energy, virial, etc.
        if(i % period == 0) {
            core->force->aux_enable();
        }

        // perform MD step
        core->mdstep();

        // measurement
        if(i % period == 0) {
            temp_(thermodynamics->temp());
            press(thermodynamics->pressure());
            en_pot(thermodynamics->en_pot());
            max_en_diff = max(abs(thermodynamics->en_tot() - en_tot), max_en_diff);
            core->force->aux_disable();
        }
    }
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    // with the first released version of halmd (commit f5283a2),
    // an energy drift of less than 5e-6 ε was obtained over 2e8 MD steps
    // using a smoothed potential (dt*=0.001, h=0.005)
    const double en_limit = max(2e-5, steps * 1e-12);
    BOOST_CHECK_SMALL(max_en_diff / fabs(en_tot), en_limit);

    BOOST_CHECK_CLOSE_FRACTION(temp, mean(temp_), 2e-3);
    BOOST_CHECK_CLOSE_FRACTION(density, (float)thermodynamics->density(), eps_float);

    if (dimension == 3) {
        double cV = heat_capacity_nve(mean(temp_), variance(temp_), npart);

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
        BOOST_TEST_MESSAGE("Heat capacity c_V = " << cV);
        // values from Johnson et al.: P = 1.023, Epot = -1.673  (Npart = 864)
        // values from RFA theory: P = 1.0245, Epot = -1.6717
        BOOST_CHECK_CLOSE_FRACTION(mean(press), 1.023, 3e-3);
        BOOST_CHECK_CLOSE_FRACTION(mean(en_pot), -1.673, 2e-3);
    }
}

static void __attribute__((constructor)) init_unit_test_suite()
{
    using namespace boost::assign;
    using namespace boost::unit_test;
    using namespace boost::unit_test::framework;

    // parametrize specific program options
    vector<string> backend = list_of
        ("host")
#ifdef WITH_CUDA
        ("gpu")
#endif /* WITH_CUDA */
        ;

    test_suite* ts1 = BOOST_TEST_SUITE( "host" );

    test_suite* ts11 = BOOST_TEST_SUITE( "2d" );
    ts11->add( BOOST_PARAM_TEST_CASE( &ideal_gas<2>, backend.begin(), backend.begin() + 1 ) );
    ts11->add( BOOST_PARAM_TEST_CASE( &thermodynamics<2>, backend.begin(), backend.begin() + 1 ) );

    test_suite* ts12 = BOOST_TEST_SUITE( "3d" );
    ts12->add( BOOST_PARAM_TEST_CASE( &ideal_gas<3>, backend.begin(), backend.begin() + 1 ) );
    ts12->add( BOOST_PARAM_TEST_CASE( &thermodynamics<3>, backend.begin(), backend.begin() + 1 ) );

    ts1->add( ts11 );
    ts1->add( ts12 );

#ifdef WITH_CUDA
    test_suite* ts2 = BOOST_TEST_SUITE( "gpu" );

    test_suite* ts21 = BOOST_TEST_SUITE( "2d" );
    ts21->add( BOOST_PARAM_TEST_CASE( &ideal_gas<2>, backend.begin() + 1, backend.end() ) );
    ts21->add( BOOST_PARAM_TEST_CASE( &thermodynamics<2>, backend.begin() + 1, backend.end() ) );

    test_suite* ts22 = BOOST_TEST_SUITE( "3d" );
    ts22->add( BOOST_PARAM_TEST_CASE( &ideal_gas<3>, backend.begin() + 1, backend.end() ) );
    ts22->add( BOOST_PARAM_TEST_CASE( &thermodynamics<3>, backend.begin() + 1, backend.end() ) );

    ts2->add( ts21 );
    ts2->add( ts22 );
#endif /* WITH_CUDA */

    master_test_suite().add( ts1 );
#ifdef WITH_CUDA
    master_test_suite().add( ts2 );
#endif /* WITH_CUDA */
}
