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
#include <halmd/mdsim/clock.hpp>
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
# include <halmd/observables/gpu/samples/phase_space.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif
#include <test/tools/ctest.hpp>

/**
 * test acquisition of phase space samples
 */

/**
 * helper structure to capture all host phase space samples
 */
template<int dimension, typename float_type>
struct host_samples
{
    typedef halmd::observables::host::samples::sample<dimension, float_type> position_sample_type;
    typedef halmd::observables::host::samples::sample<dimension, float_type> velocity_sample_type;
    typedef halmd::observables::host::samples::sample<1, unsigned int> species_sample_type;
    typedef halmd::observables::host::samples::sample<1, float_type> mass_sample_type;

    std::shared_ptr<position_sample_type const> position_sample;
    std::shared_ptr<velocity_sample_type const> velocity_sample;
    std::shared_ptr<species_sample_type const> species_sample;
    std::shared_ptr<mass_sample_type const> mass_sample;
};

/**
 * do nothing if input is host sample
 */
template<typename phase_space_type, typename modules_type>
struct copy_samples
{
    static std::shared_ptr<typename modules_type::input_sample_type const>
    copy(typename modules_type::output_phase_space_type&& phase_space)
    {
        auto sample = std::make_shared<typename modules_type::output_sample_type>();
        sample->position_sample = phase_space.acquire_position();
        sample->velocity_sample = phase_space.acquire_velocity();
        sample->species_sample = phase_space.acquire_species();
        sample->mass_sample = phase_space.acquire_mass();
        return sample;
    }
};

#ifdef HALMD_WITH_GPU
/**
 * copy GPU sample to host sample
 */
template<int dimension, typename float_type, typename modules_type>
struct copy_samples<halmd::observables::gpu::phase_space<halmd::observables::gpu::samples::phase_space<dimension, float_type>>, modules_type>
{
    static std::shared_ptr<typename modules_type::input_sample_type const>
    copy(typename modules_type::output_phase_space_type&& phase_space)
    {
        auto sample = phase_space.acquire();
        BOOST_CHECK(sample->step() == 1);

        typedef host_samples<dimension, float_type> host_sample_type;
        typedef typename host_sample_type::position_sample_type host_position_sample_type;
        typedef typename host_sample_type::velocity_sample_type host_velocity_sample_type;
        typedef typename host_sample_type::species_sample_type host_species_sample_type;
        typedef typename host_sample_type::mass_sample_type host_mass_sample_type;
        typedef halmd::observables::gpu::samples::phase_space<dimension, float_type> gpu_sample_type;
        typedef typename gpu_sample_type::position_array_type::value_type gpu_vector_type;

        // allocate memory
        std::shared_ptr<host_position_sample_type> result_position = std::make_shared<host_position_sample_type>(sample->position().size());
        std::shared_ptr<host_velocity_sample_type> result_velocity = std::make_shared<host_velocity_sample_type>(sample->position().size());
        std::shared_ptr<host_species_sample_type> result_species = std::make_shared<host_species_sample_type>(sample->position().size());
        std::shared_ptr<host_mass_sample_type> result_mass = std::make_shared<host_mass_sample_type>(sample->position().size());
        std::shared_ptr<host_sample_type> result = std::make_shared<host_sample_type>();
        result->position_sample = result_position;
        result->velocity_sample = result_velocity;
        result->species_sample = result_species;
        result->mass_sample = result_mass;
        cuda::host::vector<gpu_vector_type> h_buf(sample->position().size());

        // copy from GPU to host via page-locked memory

        // positions and types
        cuda::copy(sample->position(), h_buf);
        cuda::thread::synchronize();
        for (size_t i = 0; i < h_buf.size(); ++i) {
            tie(result_position->data()[i], result_species->data()[i]) <<= h_buf[i];
        }

        // velocities
        cuda::copy(sample->velocity(), h_buf);
        cuda::thread::synchronize();
        std::copy(h_buf.begin(), h_buf.end(), result_velocity->data().begin());

        return result;
    }
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
    typedef typename modules_type::input_phase_space_type input_phase_space_type;
    typedef typename modules_type::input_sample_type input_sample_type;
    typedef typename input_sample_type::position_sample_type input_position_sample_type;
    typedef typename input_sample_type::velocity_sample_type input_velocity_sample_type;
    typedef typename input_sample_type::species_sample_type input_species_sample_type;
    typedef typename input_sample_type::mass_sample_type input_mass_sample_type;
    typedef typename modules_type::output_phase_space_type output_phase_space_type;
    typedef typename modules_type::output_sample_type output_sample_type;
    typedef typename modules_type::random_type random_type;
    static bool const gpu = modules_type::gpu;

    typedef halmd::mdsim::clock clock_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    enum { dimension = vector_type::static_size };

    std::vector<unsigned int> npart;

    std::shared_ptr<box_type> box;
    std::shared_ptr<clock_type> clock;
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

    typename input_position_sample_type::array_type& input_position = input_position_sample->data();
    typename input_velocity_sample_type::array_type& input_velocity = input_velocity_sample->data();
    typename input_species_sample_type::array_type& input_species = input_species_sample->data();

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
        auto phase_space = input_phase_space_type(particle, particle_group, box, clock);
        phase_space.set_position(input_position_sample);
        phase_space.set_velocity(input_velocity_sample);
        phase_space.set_species(input_species_sample);
        phase_space.set_mass(input_mass_sample);
    }

    // randomly permute particles in memory, do it three times since permutations are
    // not commutative
    shuffle(particle, random);
    shuffle(particle, random);
    shuffle(particle, random);

    // acquire sample from particle, construct temporary sampler module
    clock->set_timestep(0); // bogus time-step
    clock->advance();

    // compare output and input, copy GPU sample to host before
    std::shared_ptr<input_sample_type const> result = copy_samples<output_phase_space_type, modules_type>::copy
                    (output_phase_space_type(particle, particle_group, box, clock));
    typename input_position_sample_type::array_type const& result_position = result->position_sample->data();
    typename input_velocity_sample_type::array_type const& result_velocity = result->velocity_sample->data();
    typename input_species_sample_type::array_type const& result_species = result->species_sample->data();

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
    clock = std::make_shared<clock_type>();
    random = std::make_shared<random_type>();
}

template <int dimension, typename float_type>
struct host_modules
{
    typedef halmd::mdsim::box<dimension> box_type;
    typedef halmd::mdsim::host::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::host::particle_groups::all<particle_type> particle_group_type;
    typedef halmd::observables::host::phase_space<dimension, float_type> input_phase_space_type;
    typedef host_samples<dimension, float_type> input_sample_type;
    typedef input_phase_space_type output_phase_space_type;
    typedef input_sample_type output_sample_type;
    typedef halmd::random::host::random random_type;
    static bool const gpu = false;
};

BOOST_AUTO_TEST_CASE( phase_space_host_2d ) {
    phase_space<host_modules<2, double> >().test();
}
BOOST_AUTO_TEST_CASE( phase_space_host_3d ) {
    phase_space<host_modules<3, double> >().test();
}

#ifdef HALMD_WITH_GPU
template <int dimension, typename float_type>
struct gpu_host_modules
{
    typedef halmd::mdsim::box<dimension> box_type;
    typedef halmd::mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::gpu::particle_groups::all<particle_type> particle_group_type;
    typedef halmd::observables::gpu::phase_space<halmd::observables::gpu::host_sample<dimension, float_type> > input_phase_space_type;
    typedef host_samples<dimension, float_type> input_sample_type;
    typedef input_phase_space_type output_phase_space_type;
    typedef input_sample_type output_sample_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    static bool const gpu = true;
};

template <int dimension, typename float_type>
struct gpu_gpu_modules
{
    typedef halmd::mdsim::box<dimension> box_type;
    typedef halmd::mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef halmd::mdsim::gpu::particle_groups::all<particle_type> particle_group_type;
    typedef halmd::observables::gpu::phase_space<halmd::observables::gpu::host_sample<dimension, float_type> > input_phase_space_type;
    typedef host_samples<dimension, float_type> input_sample_type;
    typedef halmd::observables::gpu::phase_space<halmd::observables::gpu::samples::phase_space<dimension, float_type> > output_phase_space_type;
    typedef typename output_phase_space_type::sample_type output_sample_type;
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
