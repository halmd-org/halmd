/*
 * Copyright © 2011  Felix Höfling and Peter Colberg
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

#define BOOST_TEST_MODULE boltzmann
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/make_shared.hpp>
#include <limits>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/host/particle_group.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/random/host/random.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/particle_group.hpp>
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

/**
 * test initialisation of particle velocities: boltzmann module
 */

template <typename modules_type>
struct boltzmann
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::particle_group_type particle_group_type;
    typedef typename particle_group_type::particle_type particle_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::thermodynamics_type thermodynamics_type;
    typedef typename modules_type::velocity_type velocity_type;
    typedef mdsim::clock clock_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    static unsigned int const dimension = vector_type::static_size;
    static bool const gpu = modules_type::gpu;

    unsigned npart;
    double temp;
    double density;

    shared_ptr<box_type> box;
    shared_ptr<clock_type> clock;
    shared_ptr<particle_type> particle;
    shared_ptr<random_type> random;
    shared_ptr<thermodynamics_type> thermodynamics;
    shared_ptr<velocity_type> velocity;

    void test();
    boltzmann();
};

template <typename modules_type>
void boltzmann<modules_type>::test()
{
    // generate velocity distribution
    BOOST_TEST_MESSAGE("set particle tags");
    particle->set();
    BOOST_TEST_MESSAGE("generate Maxwell-Boltzmann distribution");
    velocity->set();

    //
    // test velocity distribution of final state
    //
    // centre-of-mass velocity ⇒ mean of velocity distribution
    // each particle is an independent "measurement"
    // limit is 3σ, σ = √(<v_x²> / (N - 1)) where <v_x²> = k T
    double vcm_limit = 3 * sqrt(temp / (npart - 1));
    BOOST_TEST_MESSAGE("Absolute limit on instantaneous centre-of-mass velocity: " << vcm_limit);
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);  //< norm_inf tests the max. value

    // temperature ⇒ variance of velocity distribution
    // we have only one measurement of the variance
    // limit is 3σ, σ = √<ΔT²> where <ΔT²> / T² = 2 / (dimension × N)
    double rel_temp_limit = 3 * sqrt(2. / (dimension * npart)) / temp;
    BOOST_TEST_MESSAGE("Relative limit on instantaneous temperature: " << rel_temp_limit);
    BOOST_CHECK_CLOSE_FRACTION(thermodynamics->temp(), temp, rel_temp_limit);

    //
    // test shifting and rescaling
    //
    // multiplication of the velocities by a constant factor
    double scale = 1.5;
    velocity->rescale(scale);
    thermodynamics->clear_cache(); //< reset caches after rescaling the velocities
    BOOST_CHECK_CLOSE_FRACTION(thermodynamics->temp(), scale * scale * temp, rel_temp_limit);

    // shift mean velocity to zero
    fixed_vector<double, dimension> v_cm = thermodynamics->v_cm();
    velocity->shift(-v_cm);
    thermodynamics->clear_cache(); //< reset caches after shifting the velocities
    vcm_limit = gpu ? 0.1 * eps_float : 2 * eps;
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    // first shift, then rescale in one step
    velocity->shift_rescale(v_cm, 1 / scale);
    thermodynamics->clear_cache(); //< reset caches after modifying the velocities
    BOOST_CHECK_CLOSE_FRACTION(thermodynamics->temp(), temp, rel_temp_limit);
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm() - v_cm), vcm_limit);
}

template <typename modules_type>
boltzmann<modules_type>::boltzmann()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    npart = gpu ? 10000 : 3000;
    temp = 2.0;
    density = 0.3;
    typename box_type::vector_type box_length = pow(npart / density, 1. / dimension);

    particle = make_shared<particle_type>(npart);
    box = make_shared<box_type>(box_length);
    random = make_shared<random_type>();
    velocity = make_shared<velocity_type>(particle, random, temp);
    clock = make_shared<clock_type>(0); // bogus time-step
    thermodynamics = make_shared<thermodynamics_type>(make_shared<particle_group_type>(particle), box, clock);
}

template <int dimension, typename float_type>
struct host_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::particle_group_all<dimension, float_type> particle_group_type;
    typedef halmd::random::host::random random_type;
    typedef mdsim::host::velocities::boltzmann<dimension, float_type> velocity_type;
    typedef observables::host::thermodynamics<dimension, float_type> thermodynamics_type;
    static bool const gpu = false;
};

BOOST_AUTO_TEST_CASE( boltzmann_host_2d ) {
    boltzmann<host_modules<2, double> >().test();
}
BOOST_AUTO_TEST_CASE( boltzmann_host_3d ) {
    boltzmann<host_modules<3, double> >().test();
}

#ifdef WITH_CUDA
template <int dimension, typename float_type>
struct gpu_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle_group_all<dimension, float_type> particle_group_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    typedef observables::gpu::thermodynamics<dimension, float_type> thermodynamics_type;
    typedef mdsim::gpu::velocities::boltzmann<dimension, float_type, halmd::random::gpu::rand48> velocity_type;
    static bool const gpu = true;
};

BOOST_FIXTURE_TEST_CASE( boltzmann_gpu_2d, device ) {
    boltzmann<gpu_modules<2, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( boltzmann_gpu_3d, device ) {
    boltzmann<gpu_modules<3, float> >().test();
}
#endif // WITH_CUDA
