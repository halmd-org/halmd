/*
 * Copyright © 2010-2011  Felix Höfling and Peter Colberg
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

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE verlet
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <limits>
#include <numeric>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/mdsim/host/integrators/verlet.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/random/host/random.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/integrators/verlet.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/positions/lattice.hpp>
# include <halmd/mdsim/gpu/velocities/boltzmann.hpp>
# include <halmd/observables/gpu/thermodynamics.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif
#include <test/tools/ctest.hpp>

using namespace boost;
using namespace boost::assign; // list_of
using namespace halmd;
using namespace std;

const double eps = numeric_limits<double>::epsilon();
const float eps_float = numeric_limits<float>::epsilon();

/** test Verlet integrator: 'ideal' gas without interactions (setting ε=0) */

template <typename modules_type>
struct ideal_gas
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::integrator_type integrator_type;
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::position_type position_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::thermodynamics_type thermodynamics_type;
    typedef typename modules_type::velocity_type velocity_type;
    static bool const gpu = modules_type::gpu;

    typedef mdsim::clock clock_type;
    typedef typename clock_type::time_type time_type;
    typedef typename clock_type::step_type step_type;
    typedef mdsim::core core_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    static unsigned int const dimension = vector_type::static_size;

    float density;
    float temp;
    double timestep;
    unsigned int npart;
    fixed_vector<double, dimension> box_ratios;
    fixed_vector<double, dimension> slab;

    boost::shared_ptr<box_type> box;
    boost::shared_ptr<clock_type> clock;
    boost::shared_ptr<core_type> core;
    boost::shared_ptr<integrator_type> integrator;
    boost::shared_ptr<particle_type> particle;
    boost::shared_ptr<position_type> position;
    boost::shared_ptr<random_type> random;
    boost::shared_ptr<thermodynamics_type> thermodynamics;
    boost::shared_ptr<velocity_type> velocity;

    void test();
    ideal_gas();
    void connect();
};

template <typename modules_type>
void ideal_gas<modules_type>::test()
{
    // prepare system with Maxwell-Boltzmann distributed velocities
    BOOST_TEST_MESSAGE("assign positions and velocities");
    particle->aux_enable();              // enable computation of potential energy
    core->setup();

    const double vcm_tolerance = gpu ? 0.1 * eps_float : eps;
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_tolerance);

    double en_tot = thermodynamics->en_tot();

    // microcanonical simulation run
    BOOST_TEST_MESSAGE("run NVE simulation");
    clock->set_timestep(integrator->timestep());
    step_type steps = 1000;
    for (step_type i = 0; i < steps; ++i) {
        // last step: evaluate auxiliary variables (potential energy, virial, ...)
        if (i == steps - 1) {
            particle->aux_enable();
        }
        clock->advance();
        core->mdstep();
    }

    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_tolerance);
    BOOST_CHECK_CLOSE_FRACTION(en_tot, thermodynamics->en_tot(), 10 * eps);

    BOOST_CHECK_CLOSE_FRACTION(temp, (float)thermodynamics->temp(), eps_float);
    BOOST_CHECK_CLOSE_FRACTION(density, (float)thermodynamics->density(), eps_float);
    BOOST_CHECK_CLOSE_FRACTION(thermodynamics->pressure() / temp / density, 1., eps_float);
}

template <typename modules_type>
ideal_gas<modules_type>::ideal_gas()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    density = 1;
    temp = 1;
    time_type timestep = 0.001;
    npart = 1000;
    box_ratios = (dimension == 3) ? list_of(1.)(2.)(1.01) : list_of(1.)(2.);
    double det = accumulate(box_ratios.begin(), box_ratios.end(), 1., multiplies<double>());
    double volume = npart / density;
    double edge_length = pow(volume / det, 1. / dimension);
    typename box_type::vector_type box_length = edge_length * box_ratios;
    slab = 1;

    // create modules
    particle = boost::make_shared<particle_type>(npart);
    box = boost::make_shared<box_type>(box_length);
    random = boost::make_shared<random_type>();
    position = boost::make_shared<position_type>(particle, box, slab);
    velocity = boost::make_shared<velocity_type>(particle, random, temp);
    integrator = boost::make_shared<integrator_type>(particle, box, timestep);
    clock = boost::make_shared<clock_type>();
    thermodynamics = boost::make_shared<thermodynamics_type>(particle, box, clock);

    // create core and connect module slots to core signals
    this->connect();
}

template <typename modules_type>
void ideal_gas<modules_type>::connect()
{
    core = boost::make_shared<core_type>();
    // system preparation
    core->on_prepend_setup( bind(&particle_type::set, particle) );
    core->on_setup( bind(&position_type::set, position) );
    core->on_setup( bind(&velocity_type::set, velocity) );
    core->on_append_setup( bind(&particle_type::prepare, particle) );
    // integration step
    core->on_integrate( bind(&integrator_type::integrate, integrator) );
    core->on_finalize( bind(&integrator_type::finalize, integrator) );
}

template <int dimension, typename float_type>
struct host_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::integrators::verlet<dimension, float_type> integrator_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::positions::lattice<dimension, float_type> position_type;
    typedef halmd::random::host::random random_type;
    typedef mdsim::host::velocities::boltzmann<dimension, float_type> velocity_type;
    typedef observables::host::thermodynamics<dimension, float_type> thermodynamics_type;
    static bool const gpu = false;
};

BOOST_AUTO_TEST_CASE( ideal_gas_host_2d ) {
    ideal_gas<host_modules<2, double> >().test();
}
BOOST_AUTO_TEST_CASE( ideal_gas_host_3d ) {
    ideal_gas<host_modules<3, double> >().test();
}

#ifdef HALMD_WITH_GPU
template <int dimension, typename float_type>
struct gpu_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::integrators::verlet<dimension, float_type> integrator_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::positions::lattice<dimension, float_type> position_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    typedef observables::gpu::thermodynamics<dimension, float_type> thermodynamics_type;
    typedef mdsim::gpu::velocities::boltzmann<dimension, float_type, halmd::random::gpu::rand48> velocity_type;
    static bool const gpu = true;
};

BOOST_FIXTURE_TEST_CASE( ideal_gas_gpu_2d, device ) {
    ideal_gas<gpu_modules<2, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( ideal_gas_gpu_3d, device ) {
    ideal_gas<gpu_modules<3, float> >().test();
}
#endif // HALMD_WITH_GPU
