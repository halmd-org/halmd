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

#include <boost/program_options.hpp>
#include <limits>
#include <map>
#include <string>
#include <utility>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/mdsim/host/forces/lennard_jones.hpp>
#include <halmd/mdsim/host/integrators/verlet.hpp>
#include <halmd/mdsim/host/neighbour.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/position/lattice.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/observables/thermodynamics.hpp>
#include <halmd/options.hpp>
#include <halmd/random/host/random.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/forces/lennard_jones.hpp>
# include <halmd/mdsim/gpu/integrators/verlet.hpp>
# include <halmd/mdsim/gpu/neighbour.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/position/lattice.hpp>
# include <halmd/mdsim/gpu/velocities/boltzmann.hpp>
# include <halmd/observables/gpu/thermodynamics.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif

using namespace boost;
using namespace halmd;
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

#ifdef WITH_CUDA
shared_ptr<utility::gpu::device> make_device(
    string const& backend
  , vector<int> devices
  , unsigned int threads
)
{
    if (backend == "gpu") {
        static weak_ptr<utility::gpu::device> device;
        shared_ptr<utility::gpu::device> device_(device.lock());
        if (!device_) {
            device_ = make_shared<utility::gpu::device>(
                devices
              , threads
            );
            device = device_;
        }
        return device_;
    }
    if (backend == "host") {
        return shared_ptr<utility::gpu::device>(); // null pointer
    }
    throw runtime_error("unknown backend: " + backend);
}
#endif /* WITH_CUDA */

template <int dimension>
shared_ptr<mdsim::particle<dimension> > make_particle(
    string const& backend
  , unsigned int npart
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return make_shared<mdsim::gpu::particle<dimension, float> >(
            make_device(
                backend
              , utility::gpu::device::default_devices()
              , utility::gpu::device::default_threads()
            )
          , vector<unsigned int>(1, npart)
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return make_shared<mdsim::host::particle<dimension, double> >(
            vector<unsigned int>(1, npart)
        );
    }
    throw runtime_error("unknown backend: " + backend);
}

template <int dimension>
shared_ptr<mdsim::integrator<dimension> > make_verlet_integrator(
    string const& backend
  , shared_ptr<mdsim::particle<dimension> > particle
  , shared_ptr<mdsim::box<dimension> > box
  , double timestep
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return make_shared<mdsim::gpu::integrators::verlet<dimension, float> >(
            dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
          , timestep
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return make_shared<mdsim::host::integrators::verlet<dimension, double> >(
            dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
          , timestep
        );
    }
    throw runtime_error("unknown backend: " + backend);
}

template <int dimension>
shared_ptr<mdsim::force<dimension> > make_lennard_jones_force(
    string const& backend
  , shared_ptr<mdsim::particle<dimension> > particle
  , shared_ptr<mdsim::box<dimension> > box
  , array<float, 3> cutoff
  , array<float, 3> epsilon
  , array<float, 3> sigma
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        typedef mdsim::gpu::forces::lennard_jones<float> potential_type;
        typedef mdsim::gpu::forces::pair_trunc<dimension, float, potential_type> force_type;
        return make_shared<force_type>(
            make_shared<potential_type>(particle->ntype, cutoff, epsilon, sigma)
          , dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        typedef mdsim::host::forces::lennard_jones<double> potential_type;
        typedef mdsim::host::forces::pair_trunc<dimension, double, potential_type> force_type;
        return make_shared<force_type>(
            make_shared<potential_type>(particle->ntype, cutoff, epsilon, sigma)
          , dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
        );
    }
    throw runtime_error("unknown backend: " + backend);
}

template <int dimension>
shared_ptr<mdsim::neighbour<dimension> > make_neighbour(
    string const& backend
  , shared_ptr<mdsim::particle<dimension> > particle
  , shared_ptr<mdsim::box<dimension> > box
  , shared_ptr<mdsim::force<dimension> > force
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return make_shared<mdsim::gpu::neighbour<dimension, float> >(
            dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
          , dynamic_pointer_cast<mdsim::gpu::force<dimension, float> >(force)
          , mdsim::host::neighbour<dimension, double>::default_skin
          , mdsim::gpu::neighbour<dimension, float>::default_cell_occupancy
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return make_shared<mdsim::host::neighbour<dimension, double> >(
            dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
          , dynamic_pointer_cast<mdsim::host::force<dimension, double> >(force)
          , mdsim::host::neighbour<dimension, double>::default_skin
        );
    }
    throw runtime_error("unknown backend: " + backend);
}

shared_ptr<halmd::random::random> make_random(
    string const& backend
  , unsigned int seed
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
        return make_shared<random_type>(
            make_device(
                backend
              , utility::gpu::device::default_devices()
              , utility::gpu::device::default_threads()
            )
          , seed
          , random_type::default_blocks()
          , random_type::default_threads()
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return make_shared<halmd::random::host::random>(
            seed
        );
    }
    throw runtime_error("unknown backend: " + backend);
}

template <int dimension>
shared_ptr<mdsim::position<dimension> > make_lattice(
    string const& backend
  , shared_ptr<mdsim::particle<dimension> > particle
  , shared_ptr<mdsim::box<dimension> > box
  , shared_ptr<halmd::random::random> random
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return make_shared<mdsim::gpu::position::lattice<dimension, float, halmd::random::gpu::rand48> >(
            dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
          , dynamic_pointer_cast<halmd::random::gpu::random<halmd::random::gpu::rand48> >(random)
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return make_shared<mdsim::host::position::lattice<dimension, double> >(
            dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
          , dynamic_pointer_cast<halmd::random::host::random>(random)
        );
    }
    throw runtime_error("unknown backend: " + backend);
}

template <int dimension>
shared_ptr<mdsim::velocity<dimension> > make_boltzmann(
    string const& backend
  , shared_ptr<mdsim::particle<dimension> > particle
  , shared_ptr<halmd::random::random> random
  , double temperature
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return make_shared<mdsim::gpu::velocities::boltzmann<dimension, float, halmd::random::gpu::rand48> >(
            dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , dynamic_pointer_cast<halmd::random::gpu::random<halmd::random::gpu::rand48> >(random)
          , temperature
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return make_shared<mdsim::host::velocities::boltzmann<dimension, double> >(
            dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , dynamic_pointer_cast<halmd::random::host::random>(random)
          , temperature
        );
    }
    throw runtime_error("unknown backend: " + backend);
}

template <int dimension>
shared_ptr<observables::thermodynamics<dimension> > make_thermodynamics(
    string const& backend
  , shared_ptr<mdsim::particle<dimension> > particle
  , shared_ptr<mdsim::box<dimension> > box
  , shared_ptr<mdsim::force<dimension> > force
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return make_shared<observables::gpu::thermodynamics<dimension, float> >(
            dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
          , dynamic_pointer_cast<mdsim::gpu::force<dimension, float> >(force)
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return make_shared<observables::host::thermodynamics<dimension, double> >(
            dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
          , dynamic_pointer_cast<mdsim::host::force<dimension, double> >(force)
        );
    }
    throw runtime_error("unknown backend: " + backend);
}

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
    shared_ptr<utility::gpu::device> device = make_device(
        backend
      , utility::gpu::device::default_devices()
      , utility::gpu::device::default_threads()
    );
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
      , list_of(rc)(rc)(rc) /* cutoff */
      , list_of(0.f)(0.f)(0.f) /* epsilon */
      , mdsim::host::forces::lennard_jones<double>::default_sigma()
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
    core->prepare();

    // measure thermodynamic properties
    shared_ptr<observables::thermodynamics<dimension> > thermodynamics = make_thermodynamics(
        backend
      , core->particle
      , core->box
      , core->force
    );

    double vcm_limit = (backend == "gpu") ? 0.1 * eps_float : eps;
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
void thermodynamics(string const& backend)
{
    float density = 0.3;
    float temp = 3.0;
    float rc = 4.0;
    double timestep = 0.01;    // start with small time step for thermalisation
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
    shared_ptr<utility::gpu::device> device = make_device(
        backend
      , utility::gpu::device::default_devices()
      , utility::gpu::device::default_threads()
    );
#endif /* WITH_CUDA */
    shared_ptr<halmd::random::random> random = make_random(
        backend
      , random_seed
    );

    BOOST_TEST_MESSAGE("initialise simulation modules");
    // init core module and all dependencies
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
      , list_of(rc)(rc)(rc) /* cutoff */
      , mdsim::host::forces::lennard_jones<double>::default_epsilon()
      , mdsim::host::forces::lennard_jones<double>::default_sigma()
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

    // set different timestep
    timestep = 0.001;
    core->integrator->timestep(timestep);
    BOOST_CHECK_CLOSE_FRACTION(core->integrator->timestep(), timestep, eps_float);

    // measure thermodynamic properties
    shared_ptr<observables::thermodynamics<dimension> > thermodynamics = make_thermodynamics(
        backend
      , core->particle
      , core->box
      , core->force
    );

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
    double vcm_limit = (backend == "gpu") ? 0.1 * eps_float : 20 * eps;
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
        double Cv = heat_capacity(mean(temp_), variance(temp_), npart);

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

static void __attribute__((constructor)) init_unit_test_suite()
{
    typedef boost::program_options::variable_value variable_value;
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
