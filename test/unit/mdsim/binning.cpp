/*
 * Copyright © 2012  Peter Colberg
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

#define BOOST_TEST_MODULE binning
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <boost/assign/list_of.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <cmath>

#include <halmd/config.hpp>
#include <halmd/mdsim/host/binning.hpp>
#include <halmd/mdsim/positions/lattice_primitive.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/init.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/binning.hpp>
# include <test/tools/cuda.hpp>
#endif

/**
 * Functor to convert particle index to cell index.
 */
template <typename position_array_type, typename shape_type>
class particle_to_cell
{
private:
    typedef typename position_array_type::value_type vector_type;

public:
    typedef shape_type result_type;
    typedef typename position_array_type::size_type argument_type;

    /**
     * Construct functor to convert particle index to cell index.
     *
     * @param position array with particle positions
     * @param unit_ncell number of cells per unit length for each dimension
     */
    particle_to_cell(
        position_array_type const& position
      , vector_type const& unit_ncell
    )
      : position_(position)
      , unit_ncell_(unit_ncell) {}

    /**
     * Returns cell index of the particle with given index.
     */
    result_type operator()(argument_type n) const
    {
        return result_type(floor(element_prod(position_[n], unit_ncell_)));
    }

private:
    /** array with particle positions */
    position_array_type const& position_;
    /** number of cells per unit length for each dimension */
    vector_type unit_ncell_;
};

/**
 * Functor to accumulate particles of a cell.
 *
 * This functor is constructed for a given cell, and called for each
 * particle in the cell. The functor validates that the particle is
 * in the correct cell, appends the particle index to a given output
 * array, and increments the particle count for this cell.
 */
template <typename particle_to_cell_type, typename particle_iterator, typename cell_count_type>
class cell_accumulator
{
private:
    typedef typename particle_to_cell_type::result_type shape_type;
    typedef typename particle_to_cell_type::argument_type argument_type;

public:
    /**
     * Construct cell accumulator.
     *
     * @param index multi-dimensional cell index
     * @param particle_to_cell functor to convert particle index to cell index
     * @param particle output iterator for particle indices
     * @param count reference to cell count variable
     */
    cell_accumulator(
        particle_to_cell_type const& particle_to_cell
      , shape_type const& index
      , particle_iterator const& particle
      , cell_count_type& count
    )
      : particle_to_cell_(particle_to_cell)
      , index_(index)
      , particle_(particle)
      , count_(count) {}

    /**
     * Accumulate particle.
     *
     * @param index particle index
     */
    void operator()(argument_type n)
    {
        BOOST_CHECK_EQUAL( index_, particle_to_cell_(n) );
        *particle_++ = n;
        ++count_;
    }

private:
    /** functor to convert particle index to cell index */
    particle_to_cell_type const& particle_to_cell_;
    /** multi-dimensional cell index */
    shape_type index_;
    /** output iterator for particle indices */
    particle_iterator particle_;
    /** reference to cell count variable */
    cell_count_type& count_;
};

/**
 * Functor to visit particles of a cell.
 *
 * This functor is invoked by binning::get_cell() for each cell.
 * Given a multi-dimensional index, the functor creates a cell
 * accumulator for that cell, which validates the particles
 * in the cell.
 */
template <typename particle_to_cell_type, typename particle_iterator, typename cell_count_array_type>
class cell_index_iterator
{
private:
    typedef typename particle_to_cell_type::result_type shape_type;
    typedef typename cell_count_array_type::value_type cell_count_type;
    typedef cell_accumulator<particle_to_cell_type, particle_iterator, cell_count_type> cell_accumulator_type;

public:
    typedef boost::function_output_iterator<cell_accumulator_type> result_type;

    /**
     * Construct cell index iterator.
     *
     * @param particle output iterator for particle indices
     * @param particle_to_cell functor to convert particle index to cell index
     * @param count cell count array
     */
    cell_index_iterator(
        particle_to_cell_type const& particle_to_cell
      , particle_iterator const& particle
      , cell_count_array_type& count
    )
      : particle_to_cell_(particle_to_cell)
      , particle_(particle)
      , count_(count) {}

    result_type operator()(shape_type const& index)
    {
        count_.push_back(0);
        cell_accumulator_type acc(particle_to_cell_, index, particle_, count_.back());
        return boost::make_function_output_iterator(acc);
    }

private:
    /** functor to convert particle index to cell index */
    particle_to_cell_type const& particle_to_cell_;
    /** output iterator for particle indices */
    particle_iterator particle_;
    /** reference to cell count array */
    cell_count_array_type& count_;
};

template <typename particle_to_cell_type, typename particle_iterator, typename cell_count_array_type>
inline cell_index_iterator<particle_to_cell_type, particle_iterator, cell_count_array_type>
make_cell_index_iterator(
    particle_to_cell_type const& particle_to_cell
  , particle_iterator const& particle
  , cell_count_array_type& count
)
{
    return cell_index_iterator<particle_to_cell_type, particle_iterator, cell_count_array_type>(
        particle_to_cell
      , particle
      , count
    );
}

/**
 * Fixture for binning unit test.
 *
 * This fixture creates a lattice primitive of the given shape, a
 * simulation domain with edge lengths equal to the extents of the
 * lattice with unit lattice constant, a system of particles of number
 * of lattice points, a binning instance with given lower bound for
 * edge length of cells, and places the particle on the lattice.
 */
template <typename binning_type>
class binning_fixture
{
private:
    typedef typename binning_type::matrix_type matrix_type;

public:
    typedef typename binning_type::particle_type particle_type;
    typedef typename binning_type::box_type box_type;
    typedef typename binning_type::cell_size_type shape_type;
    typedef typename particle_type::position_type vector_type;
    typedef typename vector_type::value_type float_type;
    typedef typename shape_type::value_type size_type;
    typedef halmd::close_packed_lattice<vector_type, shape_type> lattice_type;

    /**
     * Construct fixture for binning unit test.
     *
     * @param shape number of lattice unit cells per dimension
     * @param length lower bound for edge length of cells
     */
    binning_fixture(shape_type const& shape, float_type length)
      : lattice(shape)
      , box(boost::make_shared<box_type>(typename box_type::vector_type(lattice.shape())))
      , particle(boost::make_shared<particle_type>(lattice.size()))
      , binning(boost::make_shared<binning_type>(particle, box, matrix_type(1, 1, length), 0))
    {
        particle->set_position(
            boost::make_transform_iterator(boost::make_counting_iterator(size_type(0)), lattice)
          , boost::make_transform_iterator(boost::make_counting_iterator(lattice.size()), lattice)
        );
    }

    /** close-packed lattice primitive */
    lattice_type const lattice;
    /** simulation domain */
    boost::shared_ptr<box_type> const box;
    /** particle arrays */
    boost::shared_ptr<particle_type> const particle;
    /** particle binning */
    boost::shared_ptr<binning_type> const binning;
};

/**
 * Functor to partially compress system by applying sinoidal shift.
 *
 * This functor transforms the components of a particle's position according to
 *
 * f(\xi_i) = \xi_i + c \lambda_i \sin(\xi_i \frac{2\pi}{L_i})
 *
 * where \xi is a placeholder for the i-th component of the position vector.
 *
 * For each component i, the value \lambda_i is a constant at which the
 * transformation function f(\xi_i) has an inflexion point at the middle
 * of the box, i.e. f'(L_i / 2) = f''(L_i / 2) = 0.
 *
 * For scalar values 0 ≤ c ≤ 1, the transformation is monotonic, i.e.
 * the relative ordering of particles in space is conserved.  A value of
 * c = 0 yields a uniform density across the box, while c = 1 yields a
 * non-uniform density, with maximum density at the centre of the box.
 *
 * With a sinoidal shift, the particle density is increased in the middle
 * of the box, to model non-uniform density distributions, while the
 * particles still occupy the entire volume of the box.
 */
template <typename Vector>
class density_transform_sinoidal
{
public:
    typedef Vector result_type;
    typedef typename result_type::value_type float_type;

    /**
     * Construct sinoidal compression functor.
     *
     * @param compression compression factor c, with 0 ≤ c ≤ 1
     * @param box_length edge lengths of the simulation domain
     */
    density_transform_sinoidal(float_type compression, result_type const& box_length)
      : box_length_(box_length)
      , lambda_(compression * box_length_ / (2 * M_PI)) {}

    result_type operator()(result_type const& r) const
    {
        return r + element_prod(lambda_, sin((2 * M_PI) * element_div(r, box_length_)));
    }

private:
    /** edge lengths of the simulation domain */
    result_type box_length_;
    /** values where f(\xi_i) have inflexion point at middle of box for c = 1 */
    result_type lambda_;
};

template <typename Fixture, typename Transform>
class binning_with_density_transform
  : private Fixture
{
private:
    typedef typename Fixture::float_type float_type;
    typedef typename Fixture::vector_type vector_type;
    typedef typename Fixture::shape_type shape_type;
    typedef std::vector<vector_type> position_array_type;
    typedef std::vector<unsigned int> particle_index_array_type;
    typedef std::vector<unsigned int> cell_count_array_type;
    typedef halmd::accumulator<double> accumulator_type;

public:
    binning_with_density_transform(Fixture const& f, Transform const& trans)
      : Fixture(f)
      , trans_(trans) {}

    void operator()()
    {
        BOOST_TEST_MESSAGE( "number density " << Fixture::particle->nparticle() / Fixture::box->volume() );

        // apply density transformation to particle positions
        position_array_type position;
        position.reserve(Fixture::particle->nparticle());
        Fixture::particle->get_position(back_inserter(position));
        std::transform(
            position.begin()
          , position.end()
          , position.begin()
          , trans_
        );
        Fixture::particle->set_position(
            position.begin()
          , position.end()
        );

        // number of cells per dimension
        shape_type shape = Fixture::binning->ncell();
        // total number of cells
        unsigned int ncell = std::accumulate(shape.begin(), shape.end(), 1u, std::multiplies<unsigned int>());

        BOOST_TEST_MESSAGE( "bin particles into " << ncell << " cells" );
        Fixture::binning->update();

        // allocate output array for indices of binned particles
        particle_index_array_type particle_index;
        particle_index.reserve(Fixture::particle->nparticle());
        // allocate output array for particles per cell counts
        cell_count_array_type cell_count;
        cell_count.reserve(ncell);

        // number of cells per unit length
        vector_type unit_ncell = element_div(vector_type(shape), vector_type(Fixture::box->length()));

        // check that particles are in the correct cell, and
        // output particle indices and cell counts to arrays
        Fixture::binning->get_cell(
            make_cell_index_iterator(
                particle_to_cell<position_array_type, shape_type>(position, unit_ncell)
              , back_inserter(particle_index)
              , cell_count
            )
        );

        // print statistics on cell occupation, which depends on density distribution
        accumulator_type acc = std::for_each(cell_count.begin(), cell_count.end(), accumulator_type());
        BOOST_CHECK_EQUAL( count(acc), ncell );
        BOOST_TEST_MESSAGE( "average number of particles per cell " << mean(acc) << " (σ = " << sigma(acc) << ")" );
        BOOST_TEST_MESSAGE( "minimum number of particles per cell " << *min_element(cell_count.begin(), cell_count.end()) );
        BOOST_TEST_MESSAGE( "maximum number of particles per cell " << *max_element(cell_count.begin(), cell_count.end()) );

        // check that all particles were binned
        sort(particle_index.begin(), particle_index.end());
        BOOST_CHECK_EQUAL_COLLECTIONS(
            particle_index.begin()
          , particle_index.end()
          , boost::counting_iterator<unsigned int>(0)
          , boost::counting_iterator<unsigned int>(Fixture::particle->nparticle())
        );
    }

private:
    /** unary particle density transform */
    Transform trans_;
};

template <typename binning_type>
class test_binning_non_uniform_density
{
public:
    typedef binning_fixture<binning_type> fixture_type;
    typedef typename fixture_type::shape_type shape_type;
    typedef typename fixture_type::vector_type vector_type;
    typedef density_transform_sinoidal<vector_type> transform_type;
    typedef binning_with_density_transform<fixture_type, transform_type> test_type;

    /**
     * Construct unit test for binning with non-uniform particle density.
     *
     * @param shape number of lattice unit cells per dimension
     * @param length lower bound for edge length of cells
     * @param compression scaling parameter for sinoidal transform
     */
    test_binning_non_uniform_density(
        shape_type const& shape
      , float length
      , float compression
    )
      : shape_(shape)
      , length_(length)
      , compression_(compression) {}

    /**
     * Run unit test.
     */
    void operator()() const
    {
        BOOST_TEST_MESSAGE( "number of lattice unit cells " << shape_ );
        BOOST_TEST_MESSAGE( "lower bound for edge length of cells " << length_ );
        BOOST_TEST_MESSAGE( "compress using sine shift with factor " << compression_ );
        fixture_type fixture(shape_, length_);
        transform_type transform(compression_, vector_type(fixture.box->length()));
        test_type(fixture, transform)();
    }

private:
    /** number of lattice unit cells per dimension */
    shape_type shape_;
    /** lower bound for edge length of cells */
    float length_;
    /** compression scaling parameter for sinoidal transform */
    float compression_;
};

#ifdef HALMD_WITH_GPU
template <typename test_type>
class test_with_gpu
{
public:
    test_with_gpu(test_type const& test) : test_(test) {}

    void operator()() const
    {
        set_cuda_device device;
        test_();
    }

private:
    test_type test_;
};

template <typename test_type>
inline test_with_gpu<test_type> make_test_with_gpu(test_type const& test)
{
    return test_with_gpu<test_type>(test);
}
#endif /* HALMD_WITH_GPU */

/**
 * Manual test case registration.
 */
HALMD_TEST_INIT( binning )
{
    using namespace boost::assign;
    using namespace boost::unit_test;

    test_suite* ts_host = BOOST_TEST_SUITE( "host" );
    framework::master_test_suite().add(ts_host);

    test_suite* ts_host_two = BOOST_TEST_SUITE( "two" );
    ts_host->add(ts_host_two);

    test_suite* ts_host_three = BOOST_TEST_SUITE( "three" );
    ts_host->add(ts_host_three);

#ifdef HALMD_WITH_GPU
    test_suite* ts_gpu = BOOST_TEST_SUITE( "gpu" );
    framework::master_test_suite().add(ts_gpu);

    test_suite* ts_gpu_two = BOOST_TEST_SUITE( "two" );
    ts_gpu->add(ts_gpu_two);

    test_suite* ts_gpu_three = BOOST_TEST_SUITE( "three" );
    ts_gpu->add(ts_gpu_three);
#endif

    // lower bound on edge length of cells
    float cell_length = 2;

    for (unsigned int unit = 1; unit <= 8; unit *= 2) {
        for (float compression = 0; compression <= 1; compression += 0.2) {
            {
#ifdef USE_HOST_SINGLE_PRECISION
                typedef test_binning_non_uniform_density<halmd::mdsim::host::binning<2, float> > test_type;
#else
                typedef test_binning_non_uniform_density<halmd::mdsim::host::binning<2, double> > test_type;
#endif
                // non-square box with coprime edge lengths
                test_type::shape_type shape = list_of(2 * unit)(3 * unit);

                callback0<> non_uniform_density = test_type(shape, cell_length, compression);
                ts_host_two->add(BOOST_TEST_CASE( non_uniform_density ));
            }
            {
#ifdef USE_HOST_SINGLE_PRECISION
                typedef test_binning_non_uniform_density<halmd::mdsim::host::binning<3, float> > test_type;
#else
                typedef test_binning_non_uniform_density<halmd::mdsim::host::binning<3, double> > test_type;
#endif
                // non-cubic box with coprime edge lengths
                test_type::shape_type shape = list_of(2 * unit)(5 * unit)(3 * unit);

                callback0<> non_uniform_density = test_type(shape, cell_length, compression);
                ts_host_three->add(BOOST_TEST_CASE( non_uniform_density ));
            }
#ifdef HALMD_WITH_GPU
            {
                typedef test_binning_non_uniform_density<halmd::mdsim::gpu::binning<2, float> > test_type;
                // non-square box with coprime edge lengths
                test_type::shape_type shape = list_of(2 * unit)(3 * unit);

                callback0<> non_uniform_density = make_test_with_gpu(test_type(shape, cell_length, compression));
                ts_gpu_two->add(BOOST_TEST_CASE( non_uniform_density ));
            }
            {
                typedef test_binning_non_uniform_density<halmd::mdsim::gpu::binning<3, float> > test_type;
                // non-cubic box with coprime edge lengths
                test_type::shape_type shape = list_of(2 * unit)(5 * unit)(3 * unit);

                callback0<> non_uniform_density = make_test_with_gpu(test_type(shape, cell_length, compression));
                ts_gpu_three->add(BOOST_TEST_CASE( non_uniform_density ));
            }
#endif
        }
    }
}
