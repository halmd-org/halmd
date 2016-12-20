/*
 * Copyright © 2011-2012 Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE phase_space
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <boost/numeric/ublas/banded.hpp>
#include <limits>
#include <numeric>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/all.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/phase_space.hpp>
#include <halmd/observables/host/samples/sample.hpp>
#include <halmd/random/host/random.hpp>
#ifdef HALMD_WITH_GPU
# include <cuda_wrapper/cuda_wrapper.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_groups/all.hpp>
# include <halmd/observables/gpu/phase_space.hpp>
# include <halmd/observables/gpu/samples/sample.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif
#include <test/tools/ctest.hpp>

/**
 * test acquisition of phase space samples
 */

template<typename modules_type>
struct host_samples {
    typedef typename modules_type::input_position_sample_type position_sample_type;
    typedef typename modules_type::input_velocity_sample_type velocity_sample_type;
    typedef typename modules_type::input_species_sample_type species_sample_type;

    host_samples(typename modules_type::phase_space_type&& phase_space) {
        position = phase_space.template acquire<position_sample_type>("position");
        velocity = phase_space.template acquire<velocity_sample_type>("velocity");
        species = phase_space.template acquire<species_sample_type>("species");
    }

    std::shared_ptr<position_sample_type const> position;
    std::shared_ptr<velocity_sample_type const> velocity;
    std::shared_ptr<species_sample_type const> species;
};

#ifdef HALMD_WITH_GPU
/**
 * copy GPU sample to host sample
 */
template<typename modules_type>
struct gpu_samples {
    typedef typename modules_type::input_position_sample_type position_sample_type;
    typedef typename modules_type::input_velocity_sample_type velocity_sample_type;
    typedef typename modules_type::input_species_sample_type species_sample_type;

    typedef typename halmd::observables::gpu::samples::sample<modules_type::dimension, float4> gpu_sample_type;

    gpu_samples(typename modules_type::phase_space_type&& phase_space) {
        auto g_position = phase_space.template acquire<gpu_sample_type>("g_position");
        auto g_velocity = phase_space.template acquire<gpu_sample_type>("g_velocity");

        cuda::host::vector<float4> h_buf(g_position->data().size());
        position = std::make_shared<position_sample_type>(g_position->data().size());
        species = std::make_shared<species_sample_type>(g_position->data().size());
        velocity = std::make_shared<velocity_sample_type>(g_velocity->data().size());

        cuda::copy(g_position->data().begin(), g_position->data().end(), h_buf.begin());
        cuda::thread::synchronize();
        for(size_t i = 0; i < h_buf.size(); ++i) {
            tie(position->data()[i], species->data()[i]) <<= h_buf[i];
        }

        cuda::copy(g_velocity->data().begin(), g_velocity->data().end(), h_buf.begin());
        cuda::thread::synchronize();
        std::copy(h_buf.begin(), h_buf.end(), velocity->data().begin());
    }

    std::shared_ptr<position_sample_type> position;
    std::shared_ptr<velocity_sample_type> velocity;
    std::shared_ptr<species_sample_type> species;
};
#endif

/**
 * randomly shuffle particle arrays
 */
#ifdef HALMD_WITH_GPU
template <typename particle_type, typename rng_type>
void shuffle(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<halmd::random::gpu::random<rng_type> > random
)
{
    // generate random permutation
    cuda::vector<unsigned int> g_index(particle->nparticle());
    cuda::host::vector<unsigned int> h_index(g_index.size());
    std::iota(h_index.begin(), h_index.end(), 0);
    cuda::copy(h_index.begin(), h_index.end(), g_index.begin());
    random->shuffle(g_index);

    // shuffle particles
    particle->rearrange(g_index);
}
#endif

template <typename particle_type>
void shuffle(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<halmd::random::host::random> random
)
{
    // generate random permutation
    std::vector<unsigned int> index(particle->nparticle());
    std::iota(index.begin(), index.end(), 0);
    random->shuffle(index.begin(), index.end());

    // shuffle particles
    particle->rearrange(index);
}

template <typename modules_type>
struct phase_space
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::particle_group_type particle_group_type;
    typedef typename modules_type::phase_space_type phase_space_type;
    typedef typename modules_type::random_type random_type;
    static bool const gpu = modules_type::gpu;

    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    enum { dimension = vector_type::static_size };

    typedef typename modules_type::input_position_sample_type input_position_sample_type;
    typedef typename modules_type::input_velocity_sample_type input_velocity_sample_type;
    typedef typename modules_type::input_species_sample_type input_species_sample_type;
    typedef typename modules_type::input_mass_sample_type input_mass_sample_type;

    std::vector<unsigned int> npart;

    std::shared_ptr<box_type> box;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<input_position_sample_type> input_position_sample;
    std::shared_ptr<input_velocity_sample_type> input_velocity_sample;
    std::shared_ptr<input_species_sample_type> input_species_sample;
    std::shared_ptr<input_mass_sample_type> input_mass_sample;
    std::shared_ptr<random_type> random;

    void test();

    phase_space();
};

template <typename modules_type>
void phase_space<modules_type>::test()
{
    float_type const epsilon = std::numeric_limits<float_type>::epsilon();

    auto& input_position = input_position_sample->data();
    auto& input_velocity = input_velocity_sample->data();
    auto& input_species = input_species_sample->data();

    // prepare input sample
    BOOST_CHECK_EQUAL(input_position.size(), accumulate(npart.begin(), npart.end(), 0u));
    BOOST_CHECK_EQUAL(input_velocity.size(), accumulate(npart.begin(), npart.end(), 0u));
    for (unsigned int i = 0, n = 0; i < npart.size(); ++i) { // iterate over particle species
        for (unsigned int j = 0; j < npart[i]; ++n, ++j) { // iterate over particles
            vector_type& r = input_position[n];
            vector_type& v = input_velocity[n];
            unsigned int& type = input_species[n];
            r[0] = float_type(j) + float_type(1) / (i + 1); //< a large, non-integer value
            r[1] = 0;
            r[dimension - 1] = - static_cast<float_type>(j);
            v[0] = static_cast<float_type>(i);
            v[1] = 0;
            v[dimension - 1] = float_type(1) / (j + 1);
            type = i;
        }
    }

    // copy input sample to particle
    std::shared_ptr<particle_group_type> particle_group = std::make_shared<particle_group_type>(particle);
    {
        auto phase_space = phase_space_type(particle, particle_group, box);
        phase_space.set("position", input_position_sample);
        phase_space.set("velocity", input_velocity_sample);
        phase_space.set("species", input_species_sample);
        phase_space.set("mass", input_mass_sample);
    }

    // randomly permute particles in memory, do it three times since permutations are
    // not commutative
    shuffle(particle, random);
    shuffle(particle, random);
    shuffle(particle, random);

    // compare output and input, copy GPU sample to host before
    typename modules_type::samples_type result(phase_space_type(particle, particle_group, box));
    auto const& result_position = result.position->data();
    auto const& result_velocity = result.velocity->data();
    auto const& result_species = result.species->data();

    BOOST_CHECK_EQUAL(result_position.size(), accumulate(npart.begin(), npart.end(), 0u));
    for (unsigned int i = 0, n = 0; i < npart.size(); ++i) { // iterate over particle species
        for (unsigned int j = 0; j < npart[i]; ++n, ++j) { // iterate over particles
            // compare positions with a tolerance due to mapping to and from the periodic box
            for (unsigned int k = 0; k < dimension; ++k) {
                BOOST_CHECK_CLOSE_FRACTION(result_position[n][k], input_position[n][k], 10 * epsilon);
            }
        }
    }
    // compare velocities directly as they should not have been modified
    BOOST_CHECK_EQUAL_COLLECTIONS(
        result_velocity.begin(), result_velocity.end()
      , input_velocity.begin(), input_velocity.end()
    );
    // compare particle species
    BOOST_CHECK_EQUAL_COLLECTIONS(
        result_species.begin(), result_species.end()
      , input_species.begin(), input_species.end()
    );
}

template <typename modules_type>
phase_space<modules_type>::phase_space()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // choose a value smaller than warp size and some limiting values
    npart.push_back(1024);
    npart.push_back(512);
    npart.push_back(30);
    npart.push_back(1);

    // choose a box length with is not an exactly representable as a
    // floating-point number and which is small enough to have some overflow
    // from the periodic box. In addition, some of the coordinates should sit
    // precisely at the edge.
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = 40./3;
    }

    // create modules
    particle = std::make_shared<particle_type>(accumulate(npart.begin(), npart.end(), 0), npart.size());
    box = std::make_shared<box_type>(edges);
    input_position_sample = std::make_shared<input_position_sample_type>(particle->nparticle());
    input_velocity_sample = std::make_shared<input_velocity_sample_type>(particle->nparticle());
    input_species_sample = std::make_shared<input_species_sample_type>(particle->nparticle());
    input_mass_sample = std::make_shared<input_mass_sample_type>(particle->nparticle());

    random = std::make_shared<random_type>();
}

template <int dimension_, typename float_type>
struct host_modules
{
    static constexpr int dimension = dimension_;
    typedef halmd::mdsim::box<dimension> box_type;
    typedef halmd::mdsim::host::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::host::particle_groups::all<particle_type> particle_group_type;
    typedef halmd::observables::host::phase_space<dimension, float_type> phase_space_type;
    typedef typename halmd::observables::host::samples::sample<dimension, float_type> input_position_sample_type;
    typedef typename halmd::observables::host::samples::sample<dimension, float_type> input_velocity_sample_type;
    typedef typename halmd::observables::host::samples::sample<1, unsigned int> input_species_sample_type;
    typedef typename halmd::observables::host::samples::sample<1, float_type> input_mass_sample_type;
    typedef host_samples<host_modules> samples_type;
    typedef halmd::random::host::random random_type;
    static bool const gpu = false;
};

#ifndef USE_HOST_SINGLE_PRECISION
BOOST_AUTO_TEST_CASE( phase_space_host_2d ) {
    phase_space<host_modules<2, double> >().test();
}
BOOST_AUTO_TEST_CASE( phase_space_host_3d ) {
    phase_space<host_modules<3, double> >().test();
}
#else
BOOST_AUTO_TEST_CASE( phase_space_host_2d ) {
    phase_space<host_modules<2, float> >().test();
}
BOOST_AUTO_TEST_CASE( phase_space_host_3d ) {
    phase_space<host_modules<3, float> >().test();
}
#endif

#ifdef HALMD_WITH_GPU

template <int dimension_, typename float_type>
struct gpu_host_modules
{
    static constexpr int dimension = dimension_;
    typedef halmd::mdsim::box<dimension> box_type;
    typedef halmd::mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::gpu::particle_groups::all<particle_type> particle_group_type;
    typedef halmd::observables::gpu::phase_space<dimension, float_type> phase_space_type;
    typedef typename halmd::observables::host::samples::sample<dimension, float_type> input_position_sample_type;
    typedef typename halmd::observables::host::samples::sample<dimension, float_type> input_velocity_sample_type;
    typedef typename halmd::observables::host::samples::sample<1, unsigned int> input_species_sample_type;
    typedef typename halmd::observables::host::samples::sample<1, float_type> input_mass_sample_type;
    typedef host_samples<gpu_host_modules> samples_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    static bool const gpu = true;
};

template <int dimension_, typename float_type>
struct gpu_gpu_modules
{
    static constexpr int dimension = dimension_;
    typedef halmd::mdsim::box<dimension> box_type;
    typedef halmd::mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::gpu::particle_groups::all<particle_type> particle_group_type;
    typedef halmd::observables::gpu::phase_space<dimension, float_type> phase_space_type;
    typedef typename halmd::observables::host::samples::sample<dimension, float_type> input_position_sample_type;
    typedef typename halmd::observables::host::samples::sample<dimension, float_type> input_velocity_sample_type;
    typedef typename halmd::observables::host::samples::sample<1, unsigned int> input_species_sample_type;
    typedef typename halmd::observables::host::samples::sample<1, float_type> input_mass_sample_type;
    typedef gpu_samples<gpu_gpu_modules> samples_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    static bool const gpu = true;
};

BOOST_FIXTURE_TEST_CASE( phase_space_gpu_host_2d, halmd::device ) {
    phase_space<gpu_host_modules<2, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( phase_space_gpu_host_3d, halmd::device ) {
    phase_space<gpu_host_modules<3, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( phase_space_gpu_gpu_2d, halmd::device ) {
    phase_space<gpu_gpu_modules<2, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( phase_space_gpu_gpu_3d, halmd::device ) {
    phase_space<gpu_gpu_modules<3, float> >().test();
}
#endif // HALMD_WITH_GPU
