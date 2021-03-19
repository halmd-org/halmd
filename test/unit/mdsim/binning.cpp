/*
 * Copyright © 2014 Felix Höfling
 * Copyright © 2020 Jaslo Ziska
 * Copyright © 2012 Peter Colberg
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

#define BOOST_TEST_MODULE binning
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <algorithm>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <cmath>

#include <halmd/mdsim/host/binning.hpp>
#include <halmd/mdsim/positions/lattice_primitive.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <test/tools/ctest.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/binning.hpp>
# include <test/tools/cuda.hpp>
#endif

/**
 * Partially compress system by applying sine shift to positions.
 *
 * This functor transforms the components of a particle's position according to
 *
 * f(r_i) = r_i + s \lambda_i \sin(r_i \frac{2\pi}{L_i})
 *
 * where r_i is the i-th component of the position vector.
 *
 * For each component i, the value \lambda_i is a constant at which the
 * transformation function f(r_i) has an inflexion point at the middle
 * of the box, i.e. f'(L_i / 2) = f''(L_i / 2) = 0.
 *
 * For scalar values 0 ≤ s ≤ 1, the transformation is monotonic, i.e.
 * the relative ordering of particles in space is conserved.  A value of
 * s = 0 yields a uniform density across the box, while s = 1 yields a
 * non-uniform density, with maximum density at the centre of the box.
 *
 * With a sine shift, the particle density is increased in the middle
 * of the box, to model non-uniform density distributions, while the
 * particles still occupy the entire volume of the box.
 */
template <typename particle_type, typename box_type>
static void
transform_density(particle_type& particle, box_type const& box, float scale)
{
    typedef typename particle_type::position_type position_type;

    BOOST_TEST_MESSAGE( "compress using sine shift with factor " << scale );

    // get particle positions
    std::vector<position_type> position;
    position.reserve(particle.nparticle());
    get_position(particle, std::back_inserter(position));

    // wave number per dimension
    position_type k = element_div(position_type(2 * M_PI), position_type(box.length()));
    // factors where f(r_i) have inflexion point at middle of box
    position_type lambda = element_div(position_type(scale), k);

    // set particle positions by applying density transform
    set_position(
        particle
      , boost::make_transform_iterator(position.begin(), [&](position_type const& r) {
            return r + element_prod(lambda, sin(element_prod(k, r)));
        })
    );
}

/**
 * Test particle binning.
 */
template <typename binning_type, typename particle_type, typename box_type>
static void
test_binning(binning_type& binning, particle_type const& particle, box_type const& box)
{
    typedef typename particle_type::vector_type vector_type;
    typedef typename binning_type::cell_size_type cell_size_type;

    // get particle positions
    std::vector<vector_type> position;
    position.reserve(particle.nparticle());
    get_position(particle, back_inserter(position));

    // number of cells per dimension
    cell_size_type shape = binning.ncell();
    // total number of cells
    unsigned int ncell = std::accumulate(shape.begin(), shape.end(), 1u, std::multiplies<unsigned int>());

    // allocate output array for indices of binned particles
    std::vector<unsigned int> particle_index;
    particle_index.reserve(particle.nparticle());

    // setup profiling timer and connect to profiler
    // (we have no access to the private struct binning.runtime)
    halmd::utility::profiler profiler;
    std::shared_ptr<halmd::accumulator<double>> runtime = std::make_shared<halmd::accumulator<double>>();
    profiler.on_profile(runtime, "update cell list");

    // number of cells per unit length
    vector_type unit_ncell = element_div(vector_type(shape), vector_type(box.length()));

    // Given a multi-dimensional cell index and a reference to a count
    // variable initialised to zero, this function returns a function
    // output iterator that accepts a particle index.
    //
    // The particle index is used to look up the particle's position,
    // derive and validate the cell index, append the cell index to
    // an array with particle indices, and increase the particle count
    // of the cell.
    //
    auto make_cell_iterator = [&](cell_size_type const& cell, unsigned int& count) {
        return boost::make_function_output_iterator(
            [&](unsigned int index) {
                assert( index < particle.nparticle() );
                BOOST_CHECK_EQUAL( cell, floor(element_prod(position[index], unit_ncell)) );
                particle_index.push_back(index);
                ++count;
            }
        );
    };

    // output array for counting the particles per cell
    std::vector<unsigned int> cell_count;

    BOOST_TEST_MESSAGE( "bin particles into " << ncell << " cells" );

    {
        halmd::utility::profiler::scoped_timer_type timer(*runtime);

        // get_cell() triggers the update of the cell list if out of date
        //
        // check that particles are in the correct cell, and
        // output particle indices and cell counts to arrays
        get_cell(
            binning
          , [&](cell_size_type const& cell) {
                cell_count.push_back(0);
                return make_cell_iterator(cell, cell_count.back());
            }
        );

    }

    // output profiling timers
    profiler.profile();

    // print statistics on cell occupation, which depends on density distribution
    halmd::accumulator<double> acc;
    unsigned int maximum = cell_count.front();
    unsigned int minimum = cell_count.front();

    for (unsigned int count : cell_count) {
        acc(count);
        maximum = std::max(maximum, count);
        minimum = std::min(minimum, count);
    }

    BOOST_TEST_MESSAGE( "average number of particles per cell " << mean(acc) << " (σ = " << sigma(acc) << ")" );
    BOOST_TEST_MESSAGE( "minimum number of particles per cell " << minimum );
    BOOST_TEST_MESSAGE( "maximum number of particles per cell " << maximum );

    // check that all particles were binned
    sort(particle_index.begin(), particle_index.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(
        particle_index.begin()
      , particle_index.end()
      , boost::counting_iterator<unsigned int>(0)
      , boost::counting_iterator<unsigned int>(particle.nparticle())
    );
}

/**
 * Make binning unit test.
 *
 * @param shape number of lattice unit cells per dimension
 * @param length lower bound for edge length of cells
 * @param scale scaling parameter for sinoidal transform
 *
 * This fixture creates a lattice of the given shape, a simulation domain
 * with edge lengths equal to the extents of the lattice with unit lattice
 * constant, a system of particles of number of lattice points, a binning
 * instance with given lower bound for edge length of cells, and places the
 * particle on the lattice.
 */
template <typename binning_type>
static void
test_non_uniform_density(typename binning_type::cell_size_type const& shape, float length, float scale)
{
    typedef typename binning_type::particle_type particle_type;
    typedef typename binning_type::matrix_type matrix_type;
    typedef typename binning_type::box_type box_type;
    typedef typename binning_type::cell_size_type shape_type;
    typedef typename particle_type::position_type vector_type;
    typedef typename shape_type::value_type size_type;

    BOOST_TEST_MESSAGE( "number of lattice unit cells " << shape );
    BOOST_TEST_MESSAGE( "lower bound for edge length of cells " << length );

    // convert box edge lengths to edge vectors
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(shape.size());
    for (unsigned int i = 0; i < shape.size(); ++i) {
        edges(i, i) = shape[i];
    }

    // create close-packed lattice of given shape
    halmd::close_packed_lattice<vector_type, shape_type> lattice(shape);
    // create simulation domain
    std::shared_ptr<box_type> box(new box_type(edges));
    // create system of particles of number of lattice points
    std::shared_ptr<particle_type> particle(new particle_type(lattice.size(), 1));
    // create particle binning
    binning_type binning(particle, box, matrix_type(1, 1, length), 0);

    BOOST_TEST_MESSAGE( "number density " << particle->nparticle() / box->volume() );

    // place particles on lattice
    set_position(
        *particle
      , boost::make_transform_iterator(boost::make_counting_iterator(size_type(0)), lattice)
    );
    // transform density using sine shift
    transform_density(*particle, *box, scale);

    // bin particles and test output
    test_binning(binning, *particle, *box);
}

/**
 * Data-driven test case registration.
 */
using namespace boost::unit_test;

unsigned int const DATA_ARRAY_UNIT[] = {1, 2, 4, 8};
float const DATA_ARRAY_COMPRESSION[] = {0., 0.25, 0.5, 0.75, 1.};
auto dataset = data::make(DATA_ARRAY_UNIT) * data::make(DATA_ARRAY_COMPRESSION);
float cell_length = 2;

BOOST_AUTO_TEST_SUITE( host )
    BOOST_DATA_TEST_CASE( two, dataset, unit, compression ) {
#ifdef USE_HOST_SINGLE_PRECISION
        typedef halmd::mdsim::host::binning<2, float> binning_type;
#else
        typedef halmd::mdsim::host::binning<2, double> binning_type;
#endif
        test_non_uniform_density<binning_type>(
            {2 * unit, 3 * unit} // non-square box with coprime edge lengths
          , cell_length
          , compression
        );
    }
    BOOST_DATA_TEST_CASE( three, dataset, unit, compression ) {
#ifdef USE_HOST_SINGLE_PRECISION
        typedef halmd::mdsim::host::binning<3, float> binning_type;
#else
        typedef halmd::mdsim::host::binning<3, double> binning_type;
#endif
        test_non_uniform_density<binning_type>(
            {2 * unit, 5 * unit, 3 * unit} // non-cubic box with coprime edge lengths
          , cell_length
          , compression
        );
    }
BOOST_AUTO_TEST_SUITE_END()

#ifdef HALMD_WITH_GPU
BOOST_FIXTURE_TEST_SUITE( gpu, set_cuda_device )
    BOOST_AUTO_TEST_SUITE( two )
# ifdef USE_GPU_SINGLE_PRECISION
        BOOST_DATA_TEST_CASE( type_float, dataset, unit, compression ) {
            typedef halmd::mdsim::gpu::binning<2, float> binning_type;
            test_non_uniform_density<binning_type>(
                {2 * unit, 3 * unit} // non-square box with coprime edge lengths
              , cell_length
              , compression
            );
        }
# endif
# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
        BOOST_DATA_TEST_CASE( type_dsfloat, dataset, unit, compression) {
            typedef halmd::mdsim::gpu::binning<2, halmd::dsfloat> binning_type;
            test_non_uniform_density<binning_type>(
                {2 * unit, 3 * unit} // non-square box with coprime edge lengths
              , cell_length
              , compression
            );
        }
# endif
    BOOST_AUTO_TEST_SUITE_END()
    BOOST_AUTO_TEST_SUITE( three )
# ifdef USE_GPU_SINGLE_PRECISION
        BOOST_DATA_TEST_CASE( type_float, dataset, unit, compression ) {
            typedef halmd::mdsim::gpu::binning<3, float> binning_type;
            test_non_uniform_density<binning_type>(
                {2 * unit, 5 * unit, 3 * unit} // non-cubic box with coprime edge lengths
              , cell_length
              , compression
            );
        }
# endif
# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
        BOOST_DATA_TEST_CASE( type_dsfloat, dataset, unit, compression) {
            typedef halmd::mdsim::gpu::binning<3, halmd::dsfloat> binning_type;
            test_non_uniform_density<binning_type>(
                {2 * unit, 5 * unit, 3 * unit} // non-cubic box with coprime edge lengths
              , cell_length
              , compression
            );
        }
# endif
    BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif
