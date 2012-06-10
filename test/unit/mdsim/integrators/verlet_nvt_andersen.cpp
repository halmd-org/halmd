/*
 * Copyright © 2010-2012  Felix Höfling
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

#define BOOST_TEST_MODULE verlet_nvt_andersen
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <cmath>
#include <numeric>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/integrators/verlet_nvt_andersen.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/all.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/random/host/random.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/integrators/verlet_nvt_andersen.hpp>
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

/**
 * test NVT Verlet integrator with stochastic Andersen thermostat
 */
template <typename modules_type>
struct verlet_nvt_andersen
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::force_type force_type;
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

    double timestep;
    float density;
    float temp;
    double coll_rate;
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
    verlet_nvt_andersen();
    void connect();
};

template <typename modules_type>
void verlet_nvt_andersen<modules_type>::test()
{
    // run for Δt*=500
    unsigned int steps = static_cast<unsigned int>(ceil(500 / timestep));
    // ensure that sampling period is sufficiently large such that
    // the samples can be considered independent
    unsigned int period = static_cast<unsigned int>(round(3. / (coll_rate * timestep)));
    accumulator<double> temp_;
    boost::array<accumulator<double>, dimension> v_cm;   //< accumulate velocity component-wise

    position->set();
    velocity->set();
    BOOST_TEST_MESSAGE("run NVT integrator over " << steps << " steps");
    for (unsigned int i = 0; i < steps; ++i) {
        integrator->integrate();
        integrator->finalize();
        if(i % period == 0) {
            temp_(thermodynamics->temp());
            fixed_vector<double, dimension> v(thermodynamics->v_cm());
            for (unsigned int i = 0; i < dimension; ++i) {
                v_cm[i](v[i]);
            }
        }
    }

    //
    // test velocity distribution of final state
    //
    // centre-of-mass velocity ⇒ mean of velocity distribution
    // each particle is an independent "measurement",
    // tolerance is 4.5σ, σ = √(<v_x²> / (N - 1)) where <v_x²> = k T,
    // with this choice, a single test passes with 99.999% probability
    double vcm_tolerance = 4.5 * sqrt(temp / (npart - 1));
    BOOST_TEST_MESSAGE("Absolute tolerance on instantaneous centre-of-mass velocity: " << vcm_tolerance);
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_tolerance);  //< norm_inf tests the max. value

    // temperature ⇒ variance of velocity distribution
    // we have only one measurement of the variance,
    // tolerance is 4.5σ, σ = √<ΔT²> where <ΔT²> / T² = 2 / (dimension × N)
    double rel_temp_tolerance = 4.5 * sqrt(2. / (dimension * npart)) / temp;
    BOOST_TEST_MESSAGE("Relative tolerance on instantaneous temperature: " << rel_temp_tolerance);
    BOOST_CHECK_CLOSE_FRACTION(thermodynamics->temp(), temp, rel_temp_tolerance);

    //
    // test velocity distribution averaged over the whole simulation run
    //
    // centre-of-mass velocity ⇒ mean of velocity distribution
    // #measurements = #particles × #samples,
    // tolerance is 4.5σ, σ = √(<v_x²> / (N × C - 1)) where <v_x²> = k T
    vcm_tolerance = 4.5 * sqrt(temp / (npart * count(v_cm[0]) - 1));
    BOOST_TEST_MESSAGE("Absolute tolerance on centre-of-mass velocity: " << vcm_tolerance);
    for (unsigned int i = 0; i < dimension; ++i) {
        BOOST_CHECK_SMALL(mean(v_cm[i]), vcm_tolerance);
        BOOST_CHECK_SMALL(error_of_mean(v_cm[i]), vcm_tolerance);
    }

    // mean temperature ⇒ variance of velocity distribution
    // each sample should constitute an independent measurement,
    // tolerance is 4.5σ, σ = √(<ΔT²> / (C - 1)) where <ΔT²> / T² = 2 / (dimension × N)
    rel_temp_tolerance = 4.5 * sqrt(2. / (dimension * npart * (count(temp_) - 1))) / temp;
    BOOST_TEST_MESSAGE("Relative tolerance on temperature: " << rel_temp_tolerance);
    BOOST_CHECK_CLOSE_FRACTION(mean(temp_), temp, rel_temp_tolerance);

    // specific heat per particle ⇒ temperature fluctuations
    // c_V = k × (dimension × N / 2)² <ΔT²> / T² / N = k × dimension / 2
    // where we have used <ΔT²> / T² = 2 / (dimension × N),
    // tolerance is 4.5σ, with the approximation
    // σ² = Var[ΔE² / (k T²)] / C → (dimension / 2) × (dimension + 6 / N) / C
    // (one measurement only from the average over C samples)
    double cv = pow(.5 * dimension, 2.) * npart * variance(temp_);
    double cv_variance=  (.5 * dimension) * (dimension + 6. / npart) / count(temp_);
    double rel_cv_tolerance = 4.5 * sqrt(cv_variance) / (.5 * dimension);
    BOOST_TEST_MESSAGE("Relative tolerance on specific heat: " << rel_cv_tolerance);
    BOOST_CHECK_CLOSE_FRACTION(cv, .5 * dimension, rel_cv_tolerance);
}

template <typename modules_type>
verlet_nvt_andersen<modules_type>::verlet_nvt_andersen()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    density = 0.3;
    timestep = 0.01;
    temp = 1;
    coll_rate = 10;
    npart = gpu ? 5000 : 1500;
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
    std::shared_ptr<force_type> force = std::make_shared<force_type>(*particle);
    integrator = std::make_shared<integrator_type>(particle, force, box, random, timestep, temp, coll_rate);
    std::shared_ptr<particle_group_type> group = std::make_shared<particle_group_type>(particle);
    thermodynamics = std::make_shared<thermodynamics_type>(particle, force, group, box);
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
#endif

/**
 * Zero force.
 */
template <typename force_type>
class zero_force
  : public force_type
{
public:
    typedef typename force_type::net_force_array_type net_force_array_type;
    typedef typename force_type::en_pot_array_type en_pot_array_type;
    typedef typename force_type::stress_pot_array_type stress_pot_array_type;
    typedef typename force_type::hypervirial_array_type hypervirial_array_type;

    template <typename particle_type>
    zero_force(particle_type const& particle)
    {
        halmd::cache_proxy<net_force_array_type> net_force = net_force_;
        halmd::cache_proxy<stress_pot_array_type> stress_pot = stress_pot_;
        make_array_from_particle(*net_force, particle);
        make_array_from_particle(*stress_pot, particle);
    }

    virtual halmd::cache<net_force_array_type> const& net_force()
    {
        return net_force_;
    }

    virtual halmd::cache<en_pot_array_type> const& en_pot()
    {
        throw std::runtime_error("not implemented");
    }

    virtual halmd::cache<stress_pot_array_type> const& stress_pot()
    {
        return stress_pot_;
    }

    virtual halmd::cache<hypervirial_array_type> const& hypervirial()
    {
        throw std::runtime_error("not implemented");
    }

private:
    halmd::cache<net_force_array_type> net_force_;
    halmd::cache<stress_pot_array_type> stress_pot_;
};

template <int dimension, typename float_type>
struct host_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef zero_force<mdsim::host::force<dimension, float_type>> force_type;
    typedef mdsim::host::integrators::verlet_nvt_andersen<dimension, float_type> integrator_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::particle_groups::all<particle_type> particle_group_type;
    typedef mdsim::host::positions::lattice<dimension, float_type> position_type;
    typedef halmd::random::host::random random_type;
    typedef mdsim::host::velocities::boltzmann<dimension, float_type> velocity_type;
    typedef observables::host::thermodynamics<dimension, float_type> thermodynamics_type;
    static bool const gpu = false;
};

BOOST_AUTO_TEST_CASE( verlet_nvt_andersen_host_2d ) {
    verlet_nvt_andersen<host_modules<2, double> >().test();
}
BOOST_AUTO_TEST_CASE( verlet_nvt_andersen_host_3d ) {
    verlet_nvt_andersen<host_modules<3, double> >().test();
}

#ifdef HALMD_WITH_GPU
template <int dimension, typename float_type>
struct gpu_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef zero_force<mdsim::gpu::force<dimension, float_type>> force_type;
    typedef mdsim::gpu::integrators::verlet_nvt_andersen<dimension, float_type, halmd::random::gpu::rand48> integrator_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::particle_groups::all<particle_type> particle_group_type;
    typedef mdsim::gpu::positions::lattice<dimension, float_type> position_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    typedef observables::gpu::thermodynamics<dimension, float_type> thermodynamics_type;
    typedef mdsim::gpu::velocities::boltzmann<dimension, float_type, halmd::random::gpu::rand48> velocity_type;
    static bool const gpu = true;
};

BOOST_FIXTURE_TEST_CASE( verlet_nvt_andersen_gpu_2d, device ) {
    verlet_nvt_andersen<gpu_modules<2, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( verlet_nvt_andersen_gpu_3d, device ) {
    verlet_nvt_andersen<gpu_modules<3, float> >().test();
}
#endif // HALMD_WITH_GPU
