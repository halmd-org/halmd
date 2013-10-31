/*
 * Copyright © 2010-2011  Felix Höfling
 * Copyright © 2013       Nicolas Höft
 * Copyright © 2010-2011  Peter Colberg
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
#include <boost/numeric/ublas/banded.hpp>
#include <limits>
#include <numeric>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/integrators/verlet.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/all.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/random/host/random.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/integrators/verlet.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_groups/all.hpp>
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
    typedef typename modules_type::particle_group_type particle_group_type;
    typedef typename modules_type::position_type position_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::thermodynamics_type thermodynamics_type;
    typedef typename modules_type::velocity_type velocity_type;
    static bool const gpu = modules_type::gpu;

    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    static unsigned int const dimension = vector_type::static_size;

    float density;
    float temp;
    double timestep;
    unsigned int npart;
    fixed_vector<double, dimension> box_ratios;
    fixed_vector<double, dimension> slab;

    std::shared_ptr<box_type> box;
    std::shared_ptr<integrator_type> integrator;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<position_type> position;
    std::shared_ptr<random_type> random;
    std::shared_ptr<thermodynamics_type> thermodynamics;
    std::shared_ptr<velocity_type> velocity;

    void test();
    ideal_gas();
};

template <typename modules_type>
void ideal_gas<modules_type>::test()
{
    // prepare system with Maxwell-Boltzmann distributed velocities
    BOOST_TEST_MESSAGE("assign positions and velocities");
    position->set();
    velocity->set();

    const double vcm_tolerance = gpu ? 0.1 * eps_float : eps;
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_tolerance);

    double en_kin = thermodynamics->en_kin();

    // microcanonical simulation run
    BOOST_TEST_MESSAGE("run NVE simulation");
    unsigned int constexpr steps = 1000;
    for (unsigned int i = 0; i < steps; ++i) {
        integrator->integrate();
        if (i == steps -1) {
            particle->aux_enable();
        }
        integrator->finalize();
    }

    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_tolerance);
    BOOST_CHECK_CLOSE_FRACTION(en_kin, thermodynamics->en_kin(), 10 * eps);

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
    double timestep = 0.001;
    npart = 1000;
    box_ratios = (dimension == 3) ? list_of(1.)(2.)(1.01) : list_of(1.)(2.);
    double det = accumulate(box_ratios.begin(), box_ratios.end(), 1., multiplies<double>());
    double volume = npart / density;
    double edge_length = pow(volume / det, 1. / dimension);
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = edge_length * box_ratios[i];
    }
    slab = 1;

    // create modules
    particle = std::make_shared<particle_type>(npart, 1);
    box = std::make_shared<box_type>(edges);
    random = std::make_shared<random_type>();
    position = std::make_shared<position_type>(particle, box, slab);
    velocity = std::make_shared<velocity_type>(particle, random, temp);
    integrator = std::make_shared<integrator_type>(particle, box, timestep);
    std::shared_ptr<particle_group_type> group = std::make_shared<particle_group_type>(particle);
    thermodynamics = std::make_shared<thermodynamics_type>(particle, group, box);
}

/**
 * Construct host array from host particle.
 */
template <typename array_type, typename particle_type>
inline typename std::enable_if<
    std::is_convertible<
        typename std::iterator_traits<typename array_type::iterator>::iterator_category
      , std::random_access_iterator_tag
    >::value
  , void>::type
make_array_from_particle(array_type& array, particle_type const& particle)
{
    array_type output(particle.nparticle());
    std::fill(output.begin(), output.end(), 0);
    array = std::move(output);
}

/**
 * Construct stress pot tensor array from particle.
 */
template <typename particle_type>
inline typename std::enable_if<
    std::is_convertible<
        typename std::iterator_traits<typename particle_type::stress_pot_array_type::iterator>::iterator_category
      , std::random_access_iterator_tag
    >::value
  , void>::type
make_stress_pot_from_particle(
    typename particle_type::stress_pot_array_type& array,
    particle_type const& particle
)
{
   make_array_from_particle(array, particle);
}

#ifdef HALMD_WITH_GPU
/**
 * Construct GPU array from GPU particle.
 */
template <typename array_type, typename particle_type>
inline typename std::enable_if<
    std::is_convertible<
        typename std::iterator_traits<typename array_type::iterator>::iterator_category
      , cuda::device_random_access_iterator_tag
    >::value
  , void>::type
make_array_from_particle(array_type& array, particle_type const& particle)
{
    array_type g_output(particle.nparticle());
    g_output.reserve(particle.dim.threads());
    cuda::memset(g_output.begin(), g_output.end(), 0);
    array = std::move(g_output);
}

/**
 * Construct stress tensor GPU array from GPU particle.
 */
template <typename particle_type>
inline typename std::enable_if<
    std::is_convertible<
        typename std::iterator_traits<typename particle_type::stress_pot_array_type::iterator>::iterator_category
      , cuda::device_random_access_iterator_tag
    >::value
  , void>::type
make_stress_pot_from_particle(
    typename particle_type::stress_pot_array_type& array,
    particle_type const& particle
)
{
    int constexpr stress_pot_size = particle_type::stress_pot_type::static_size;
    typename particle_type::stress_pot_array_type g_output(particle.nparticle());
    g_output.reserve(particle.dim.threads() * stress_pot_size);
    cuda::memset(g_output.begin(), g_output.begin() + g_output.capacity(), 0);
    array = std::move(g_output);
}
#endif

template <int dimension, typename float_type>
struct host_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::integrators::verlet<dimension, float_type> integrator_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::particle_groups::all<particle_type> particle_group_type;
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
    typedef mdsim::gpu::particle_groups::all<particle_type> particle_group_type;
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
