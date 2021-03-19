/*
 * Copyright © 2012 Peter Colberg
 * Copyright © 2020 Jaslo Ziska
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

#define BOOST_TEST_MODULE particle
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/positions/lattice_primitive.hpp>
#include <test/tools/constant_iterator.hpp>
#include <test/tools/ctest.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/particle.hpp>
# include <test/tools/cuda.hpp>
#endif

#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <algorithm>
#include <cmath>

/**
 * Primitive lattice with equal number of lattice points per dimension.
 *
 * This lattice functor generates two, three or four-dimensional vectors.
 */
template <typename vector_type>
class equilateral_lattice
  : public halmd::primitive_lattice<vector_type, halmd::fixed_vector<size_t, vector_type::static_size> >
{
public:
    equilateral_lattice(size_t nparticle) : lattice_type(make_lattice(nparticle)) {}

private:
    typedef halmd::primitive_lattice<vector_type, halmd::fixed_vector<size_t, vector_type::static_size> > lattice_type;

    static lattice_type make_lattice(unsigned int nparticle)
    {
        return lattice_type(std::ceil(std::pow(nparticle, 1. / vector_type::static_size)));
    }
};

/**
 * Make lattice iterator given a lattice primitive and particle index.
 */
template <typename lattice_type>
inline boost::transform_iterator<lattice_type, boost::counting_iterator<size_t> >
make_lattice_iterator(lattice_type const& lattice, size_t n)
{
    return boost::make_transform_iterator(boost::make_counting_iterator(n), lattice);
}


/**
 * Construct host array from iterator range.
 */
template <typename array_type, typename iterator_type>
inline typename std::enable_if<
    std::is_convertible<
        typename std::iterator_traits<typename array_type::iterator>::iterator_category
      , std::random_access_iterator_tag
    >::value
  , void>::type
make_array_from_iterator_range(array_type& array, iterator_type const& first, iterator_type const& last)
{
    array_type output(last - first);
    std::copy(first, last, output.begin());
    array = std::move(output);
}

/**
 * Construct host stress pot array from iterator range.
 */
template <typename particle_type, typename iterator_type>
inline typename std::enable_if<
    std::is_convertible<
        typename std::iterator_traits<typename particle_type::stress_pot_array_type::iterator>::iterator_category
      , std::random_access_iterator_tag
    >::value
  , void>::type
make_stress_array_from_iterator_range(
    typename particle_type::stress_pot_array_type& array
  , iterator_type const& first, iterator_type const& last
)
{
    make_array_from_iterator_range(array, first, last);
}

#ifdef HALMD_WITH_GPU
/**
 * Construct GPU array from iterator range.
 */
template <typename array_type, typename iterator_type>
inline typename std::enable_if<
    std::is_convertible<
        typename std::iterator_traits<typename array_type::iterator>::iterator_category
      , cuda::device_random_access_iterator_tag
    >::value
  , void>::type
make_array_from_iterator_range(array_type& array, iterator_type const& first, iterator_type const& last)
{
    array_type g_output(last - first);
    cuda::memory::host::vector<typename array_type::value_type> h_output(g_output.size());
    std::copy(first, last, h_output.begin());
    cuda::copy(h_output.begin(), h_output.end(), g_output.begin());
    array = std::move(g_output);
}

/**
 * Construct GPU stress tensor array from iterator range.
 */
template <typename particle_type, typename iterator_type>
inline typename std::enable_if<
    std::is_convertible<
        typename std::iterator_traits<typename particle_type::stress_pot_array_type::iterator>::iterator_category
      , cuda::device_random_access_iterator_tag
    >::value
  , void>::type
make_stress_array_from_iterator_range(
    typename particle_type::stress_pot_array_type& array
  , iterator_type const& first, iterator_type const& last
)
{
    int constexpr stress_pot_size = particle_type::stress_pot_type::static_size;
    typename particle_type::stress_pot_array_type g_stress_pot(last - first);
    g_stress_pot.reserve(g_stress_pot.size() * stress_pot_size);
    cuda::memory::host::vector<typename particle_type::stress_pot_array_type::value_type> h_stress_pot(g_stress_pot.size());
    h_stress_pot.reserve(g_stress_pot.capacity());

    // convert stress tensor to column-major memory layout
    unsigned int stride = h_stress_pot.capacity() / stress_pot_size;
    for (iterator_type it = first; it != last; ++it) {
        halmd::mdsim::write_stress_tensor(&h_stress_pot[it - first], *it, stride);
    }
    // copy from memory from host to GPU
    cuda::copy(h_stress_pot.begin(), h_stress_pot.begin() + h_stress_pot.capacity(), g_stress_pot.begin());
    array = std::move(g_stress_pot);
}
#endif // HALMD_WITH_GPU

template <typename particle_type, typename iterator_type>
void set_force(particle_type& particle, iterator_type const& first, iterator_type const& last)
{
    auto force = halmd::make_cache_mutable(particle.mutable_force());
    make_array_from_iterator_range(*force, first, last);
}

template <typename particle_type, typename iterator_type>
void set_potential_energy(particle_type& particle, iterator_type const& first, iterator_type const& last)
{
    auto en_pot = halmd::make_cache_mutable(particle.mutable_potential_energy());
    make_array_from_iterator_range(*en_pot, first, last);
}

template <typename particle_type, typename iterator_type>
void set_stress_pot(particle_type& particle, iterator_type const& first, iterator_type const& last)
{
    auto stress_pot = halmd::make_cache_mutable(particle.mutable_stress_pot());
    make_stress_array_from_iterator_range<particle_type>(*stress_pot, first, last);
}


/**
 * Test initialisation, getter and setter of particle positions.
 */
template <typename particle_type>
static void test_position(particle_type& particle)
{
    typedef typename particle_type::position_type position_type;
    typedef typename particle_type::species_type species_type;
    particle_type const& const_particle = particle;

    // set species to ascending sequence of integers starting at 1 ≠ 0
    BOOST_CHECK(
        set_species(particle, boost::counting_iterator<species_type>(1))
            == boost::counting_iterator<species_type>(particle.nparticle() + 1)
    );

    // check that positions are initialised to zero
    std::vector<position_type> position(particle.nparticle());
    BOOST_CHECK(
        get_position(const_particle, position.begin()) == position.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        position.begin()
      , position.end()
      , make_constant_iterator(position_type(0), 0)
      , make_constant_iterator(position_type(0), particle.nparticle())
    );

    // assign square/cubic lattice vectors
    equilateral_lattice<position_type> lattice(particle.nparticle());
    BOOST_CHECK(
        set_position(particle, make_lattice_iterator(lattice, 0))
            == make_lattice_iterator(lattice, particle.nparticle())
    );
    BOOST_CHECK(
        get_position(const_particle, position.begin()) == position.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        position.begin()
      , position.end()
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, particle.nparticle())
    );

    // check that particle species are preserved, since positions
    // and species are stored in the same array in gpu::particle
    std::vector<species_type> species(particle.nparticle());
    BOOST_CHECK(
        get_species(const_particle, species.begin()) == species.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        species.begin()
      , species.end()
      , boost::counting_iterator<species_type>(1)
      , boost::counting_iterator<species_type>(particle.nparticle() + 1)
    );
}

/**
 * Test initialisation, getter and setter of particle images.
 */
template <typename particle_type>
static void test_image(particle_type& particle)
{
    typedef typename particle_type::image_type image_type;
    particle_type const& const_particle = particle;

    // check that images are initialised to zero
    std::vector<image_type> image(particle.nparticle());
    BOOST_CHECK(
        get_image(const_particle, image.begin()) == image.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        image.begin()
      , image.end()
      , make_constant_iterator(image_type(0), 0)
      , make_constant_iterator(image_type(0), particle.nparticle())
    );

    // assign square/cubic lattice vectors
    equilateral_lattice<image_type> lattice(particle.nparticle());
    BOOST_CHECK(
        set_image(particle, make_lattice_iterator(lattice, 0))
            == make_lattice_iterator(lattice, particle.nparticle())
    );
    BOOST_CHECK(
        get_image(const_particle, image.begin()) == image.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        image.begin()
      , image.end()
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of particle velocities.
 */
template <typename particle_type>
static void test_velocity(particle_type& particle)
{
    typedef typename particle_type::velocity_type velocity_type;
    typedef typename particle_type::mass_type mass_type;
    particle_type const& const_particle = particle;

    // set masses to ascending sequence of integers starting at 2 ≠ 1
    BOOST_CHECK(
        set_mass(particle, boost::counting_iterator<mass_type>(2))
            == boost::counting_iterator<mass_type>(particle.nparticle() + 2)
    );

    // check that velocities are initialised to zero
    std::vector<velocity_type> velocity(particle.nparticle());
    BOOST_CHECK(
        get_velocity(const_particle, velocity.begin()) == velocity.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        velocity.begin()
      , velocity.end()
      , make_constant_iterator(velocity_type(0), 0)
      , make_constant_iterator(velocity_type(0), particle.nparticle())
    );

    // assign square/cubic lattice vectors
    equilateral_lattice<velocity_type> lattice(particle.nparticle());
    BOOST_CHECK(
        set_velocity(particle, make_lattice_iterator(lattice, 0))
            == make_lattice_iterator(lattice, particle.nparticle())
    );
    BOOST_CHECK(
        get_velocity(const_particle, velocity.begin()) == velocity.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        velocity.begin()
      , velocity.end()
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, particle.nparticle())
    );

    // check that particle masses are preserved, since velocities
    // and masses are stored in the same array in gpu::particle
    std::vector<mass_type> mass(particle.nparticle());
    BOOST_CHECK(
        get_mass(const_particle, mass.begin()) == mass.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        mass.begin()
      , mass.end()
      , boost::counting_iterator<mass_type>(2)
      , boost::counting_iterator<mass_type>(particle.nparticle() + 2)
    );
}

/**
 * Test initialisation, getter and setter of particle ids.
 */
template <typename particle_type>
static void test_id(particle_type& particle)
{
    typedef typename particle_type::id_type id_type;
    particle_type const& const_particle = particle;

    // check that ids default to ascending sequence of integers
    std::vector<id_type> id(particle.nparticle());
    BOOST_CHECK(
        get_id(const_particle, id.begin()) == id.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        id.begin()
      , id.end()
      , boost::counting_iterator<id_type>(0)
      , boost::counting_iterator<id_type>(particle.nparticle())
    );

    // reverse order of particle ids
    BOOST_CHECK(
        set_id(particle, id.rbegin()) == id.rend()
    );
    BOOST_CHECK(
        get_id(const_particle, id.begin()) == id.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        id.rbegin()
      , id.rend()
      , boost::counting_iterator<id_type>(0)
      , boost::counting_iterator<id_type>(particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of particle reverse ids.
 */
template <typename particle_type>
static void test_reverse_id(particle_type& particle)
{
    typedef typename particle_type::reverse_id_type reverse_id_type;
    particle_type const& const_particle = particle;

    // check that reverse ids default to ascending sequence of integers
    std::vector<reverse_id_type> reverse_id(particle.nparticle());
    BOOST_CHECK(
        get_reverse_id(const_particle, reverse_id.begin()) == reverse_id.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        reverse_id.begin()
      , reverse_id.end()
      , boost::counting_iterator<reverse_id_type>(0)
      , boost::counting_iterator<reverse_id_type>(particle.nparticle())
    );

    // reverse order of reverse particle ids
    BOOST_CHECK(
        set_reverse_id(particle, reverse_id.rbegin()) == reverse_id.rend()
    );
    BOOST_CHECK(
        get_reverse_id(const_particle, reverse_id.begin()) == reverse_id.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        reverse_id.rbegin()
      , reverse_id.rend()
      , boost::counting_iterator<reverse_id_type>(0)
      , boost::counting_iterator<reverse_id_type>(particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of particle species.
 */
template <typename particle_type>
static void test_species(particle_type& particle)
{
    typedef typename particle_type::species_type species_type;
    typedef typename particle_type::position_type position_type;
    particle_type const& const_particle = particle;

    // check default of one species
    BOOST_CHECK_EQUAL( particle.nspecies(), 1u );

    // assign square/cubic lattice vectors
    equilateral_lattice<position_type> lattice(particle.nparticle());
    BOOST_CHECK(
        set_position(particle, make_lattice_iterator(lattice, 0))
            == make_lattice_iterator(lattice, particle.nparticle())
    );

    // check that species are initialised to zero
    std::vector<species_type> species(particle.nparticle());
    BOOST_CHECK(
        get_species(const_particle, species.begin()) == species.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        species.begin()
      , species.end()
      , make_constant_iterator(species_type(0), 0)
      , make_constant_iterator(species_type(0), particle.nparticle())
    );

    // set species to ascending sequence of integers starting at 1 ≠ 0
    BOOST_CHECK(
        set_species(particle, boost::counting_iterator<species_type>(1))
            == boost::counting_iterator<species_type>(particle.nparticle() + 1)
    );
    BOOST_CHECK(
        get_species(const_particle, species.begin()) == species.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        species.begin()
      , species.end()
      , boost::counting_iterator<species_type>(1)
      , boost::counting_iterator<species_type>(particle.nparticle() + 1)
    );

    // check that particle positions are preserved, since positions
    // and species are stored in the same array in gpu::particle
    std::vector<position_type> position(particle.nparticle());
    BOOST_CHECK(
        get_position(const_particle, position.begin()) == position.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        position.begin()
      , position.end()
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of particle masses.
 */
template <typename particle_type>
static void test_mass(particle_type& particle)
{
    typedef typename particle_type::mass_type mass_type;
    typedef typename particle_type::velocity_type velocity_type;
    particle_type const& const_particle = particle;

    // assign square/cubic lattice vectors
    equilateral_lattice<velocity_type> lattice(particle.nparticle());
    BOOST_CHECK(
        set_velocity(particle, make_lattice_iterator(lattice, 0))
            == make_lattice_iterator(lattice, particle.nparticle())
    );

    // check that masses are initialised to unit mass
    std::vector<mass_type> mass(particle.nparticle());
    BOOST_CHECK(
        get_mass(const_particle, mass.begin()) == mass.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        mass.begin()
      , mass.end()
      , make_constant_iterator(mass_type(1), 0)
      , make_constant_iterator(mass_type(1), particle.nparticle())
    );

    // set masses to ascending sequence of integers starting at 2 ≠ 1
    BOOST_CHECK(
        set_mass(particle, boost::counting_iterator<mass_type>(2))
            == boost::counting_iterator<mass_type>(particle.nparticle() + 2)
    );
    BOOST_CHECK(
        get_mass(const_particle, mass.begin()) == mass.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        mass.begin()
      , mass.end()
      , boost::counting_iterator<mass_type>(2)
      , boost::counting_iterator<mass_type>(particle.nparticle() + 2)
    );

    // check that particle velocities are preserved, since velocities
    // and masses are stored in the same array in gpu::particle
    std::vector<velocity_type> velocity(particle.nparticle());
    BOOST_CHECK(
        get_velocity(const_particle, velocity.begin()) == velocity.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        velocity.begin()
      , velocity.end()
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of particle forces.
 */
template <typename particle_type>
static void test_force(particle_type& particle)
{
    typedef typename particle_type::force_type force_type;

    // initialise forces to zero
    set_force(
        particle
      , make_constant_iterator(force_type(0), 0)
      , make_constant_iterator(force_type(0), particle.nparticle())
    );
    std::vector<force_type> force(particle.nparticle());
    BOOST_CHECK( get_force(
        particle
      , force.begin()) == force.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        force.begin()
      , force.end()
      , make_constant_iterator(force_type(0), 0)
      , make_constant_iterator(force_type(0), particle.nparticle())
    );

    // assign square/cubic lattice vectors
    equilateral_lattice<force_type> lattice(particle.nparticle());
    set_force(
        particle
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, particle.nparticle())
    );
    BOOST_CHECK( get_force(
        particle
      , force.begin()) == force.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        force.begin()
      , force.end()
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of potential energy per particle.
 */
template <typename particle_type>
static void test_en_pot(particle_type& particle)
{
    typedef typename particle_type::en_pot_type en_pot_type;

    // initialise potential energies to zero
    set_potential_energy(
        particle
      , make_constant_iterator(en_pot_type(0), 0)
      , make_constant_iterator(en_pot_type(0), particle.nparticle())
    );
    std::vector<en_pot_type> en_pot(particle.nparticle());
    BOOST_CHECK( get_potential_energy(
        particle
      , en_pot.begin()) == en_pot.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        en_pot.begin()
      , en_pot.end()
      , make_constant_iterator(en_pot_type(0), 0)
      , make_constant_iterator(en_pot_type(0), particle.nparticle())
    );

    // set potential energies to ascending sequence of integers starting at 1 ≠ 0
    set_potential_energy(
        particle
      , boost::counting_iterator<en_pot_type>(1)
      , boost::counting_iterator<en_pot_type>(particle.nparticle() + 1)
    );
    BOOST_CHECK( get_potential_energy(
        particle
      , en_pot.begin()) == en_pot.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        en_pot.begin()
      , en_pot.end()
      , boost::counting_iterator<en_pot_type>(1)
      , boost::counting_iterator<en_pot_type>(particle.nparticle() + 1)
    );
}


/**
 * Test initialisation, getter and setter of potential part of stress tensor per particle.
 */
template <typename particle_type>
static void test_stress_pot(particle_type& particle)
{
    typedef typename particle_type::stress_pot_type stress_pot_type;

    // initialise stress tensors to zero
    set_stress_pot(
        particle
      , make_constant_iterator(stress_pot_type(0), 0)
      , make_constant_iterator(stress_pot_type(0), particle.nparticle())
    );
    std::vector<stress_pot_type> stress_pot(particle.nparticle());
    BOOST_CHECK( get_stress_pot(
        particle
      , stress_pot.begin()) == stress_pot.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        stress_pot.begin()
      , stress_pot.end()
      , make_constant_iterator(stress_pot_type(0), 0)
      , make_constant_iterator(stress_pot_type(0), particle.nparticle())
    );

    // assign square/four-dimensional cubic lattice vectors
    equilateral_lattice<stress_pot_type> lattice(particle.nparticle());
    set_stress_pot(
        particle
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, particle.nparticle())
    );
    BOOST_CHECK( get_stress_pot(
        particle
      , stress_pot.begin()) == stress_pot.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        stress_pot.begin()
      , stress_pot.end()
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, particle.nparticle())
    );
}

/**
 * BOOST_AUTO_TEST_SUITE only allows test cases to be registered inside it, no function calls.
 * For this reason the old test_suite_{host,gpu} function had to be replaced with these macros.
 * This way no function is called and BOOST_DATA_TEST_CASE is directly registered
 * inside the test suite.
 */
#define TEST_SUITE_HOST(particle_type, dataset, nspecies)   \
    BOOST_DATA_TEST_CASE( position, dataset, nparticle ) {  \
        particle_type particle(nparticle, nspecies);        \
        test_position(particle);                            \
    }                                                       \
    BOOST_DATA_TEST_CASE( image, dataset, nparticle ) {     \
        particle_type particle(nparticle, nspecies);        \
        test_image(particle);                               \
    }                                                       \
    BOOST_DATA_TEST_CASE( velocity, dataset, nparticle ) {  \
        particle_type particle(nparticle, nspecies);        \
        test_velocity(particle);                            \
    }                                                       \
    BOOST_DATA_TEST_CASE( id, dataset, nparticle )  {       \
        particle_type particle(nparticle, nspecies);        \
        test_id(particle);                                  \
    }                                                       \
    BOOST_DATA_TEST_CASE( reverse_id, dataset, nparticle ) {\
        particle_type particle(nparticle, nspecies);        \
        test_reverse_id(particle);                          \
    }                                                       \
    BOOST_DATA_TEST_CASE( species, dataset, nparticle ) {   \
        particle_type particle(nparticle, nspecies);        \
        test_species(particle);                             \
    }                                                       \
    BOOST_DATA_TEST_CASE( mass, dataset, nparticle ) {      \
        particle_type particle(nparticle, nspecies);        \
        test_mass(particle);                                \
    }                                                       \
    BOOST_DATA_TEST_CASE( force, dataset, nparticle ) {     \
        particle_type particle(nparticle, nspecies);        \
        test_force(particle);                               \
    }                                                       \
    BOOST_DATA_TEST_CASE( en_pot, dataset, nparticle ) {    \
        particle_type particle(nparticle, nspecies);        \
        test_en_pot(particle);                              \
    }                                                       \
    BOOST_DATA_TEST_CASE( stress_pot, dataset, nparticle ) {\
        particle_type particle(nparticle, nspecies);        \
        test_stress_pot(particle);                          \
    }

#ifdef HALMD_WITH_GPU
# define TEST_SUITE_GPU(particle_type, dataset, nspecies)   \
    BOOST_DATA_TEST_CASE( position, dataset, nparticle ) {  \
        particle_type particle(nparticle, nspecies);        \
        test_position(particle);                            \
    }                                                       \
    BOOST_DATA_TEST_CASE( image, dataset, nparticle ) {     \
        particle_type particle(nparticle, nspecies);        \
        test_image(particle);                               \
    }                                                       \
    BOOST_DATA_TEST_CASE( velocity, dataset, nparticle ) {  \
        particle_type particle(nparticle, nspecies);        \
        test_velocity(particle);                            \
    }                                                       \
    BOOST_DATA_TEST_CASE( id, dataset, nparticle ) {        \
        particle_type particle(nparticle, nspecies);        \
        test_id(particle);                                  \
    }                                                       \
    BOOST_DATA_TEST_CASE( reverse_id, dataset, nparticle ) {\
        particle_type particle(nparticle, nspecies);        \
        test_reverse_id(particle);                          \
    }                                                       \
    BOOST_DATA_TEST_CASE( species, dataset, nparticle ) {   \
        particle_type particle(nparticle, nspecies);        \
        test_species(particle);                             \
    }                                                       \
    BOOST_DATA_TEST_CASE( mass, dataset, nparticle ) {      \
        particle_type particle(nparticle, nspecies);        \
        test_mass(particle);                                \
    }                                                       \
    BOOST_DATA_TEST_CASE( force, dataset, nparticle ) {     \
        particle_type particle(nparticle, nspecies);        \
        test_force(particle);                               \
    }                                                       \
    BOOST_DATA_TEST_CASE( en_pot, dataset, nparticle ) {    \
        particle_type particle(nparticle, nspecies);        \
        test_en_pot(particle);                              \
    }                                                       \
    BOOST_DATA_TEST_CASE( stress_pot, dataset, nparticle ) {\
        particle_type particle(nparticle, nspecies);        \
        test_stress_pot(particle);                          \
    }
#endif

/**
 * Data-driven test case registration.
 */
using namespace boost::unit_test;

unsigned int const DATA_ARRAY[] = {109, 4789, 42589};
auto dataset = data::make(DATA_ARRAY);
unsigned int const nspecies = 1;

BOOST_AUTO_TEST_SUITE( host )
    BOOST_AUTO_TEST_SUITE( two )
#ifdef USE_HOST_SINGLE_PRECISION
        typedef halmd::mdsim::host::particle<2, float> particle_type;
#else
        typedef halmd::mdsim::host::particle<2, double> particle_type;
#endif
        TEST_SUITE_HOST(particle_type, dataset, nspecies)
    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE( three )
#ifdef USE_HOST_SINGLE_PRECISION
        typedef halmd::mdsim::host::particle<3, float> particle_type;
#else
        typedef halmd::mdsim::host::particle<3, double> particle_type;
#endif
        TEST_SUITE_HOST(particle_type, dataset, nspecies)
    BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#ifdef HALMD_WITH_GPU
BOOST_FIXTURE_TEST_SUITE( gpu, set_cuda_device )
    BOOST_AUTO_TEST_SUITE( two )
# ifdef USE_GPU_SINGLE_PRECISION
        BOOST_AUTO_TEST_SUITE( type_float )
            typedef halmd::mdsim::gpu::particle<2, float> particle_type;
        TEST_SUITE_GPU(particle_type, dataset, nspecies)
        BOOST_AUTO_TEST_SUITE_END()
# endif
# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
        BOOST_AUTO_TEST_SUITE( type_dsfloat )
            typedef halmd::mdsim::gpu::particle<2, halmd::dsfloat> particle_type;
        TEST_SUITE_GPU(particle_type, dataset, nspecies)
        BOOST_AUTO_TEST_SUITE_END()
# endif
    BOOST_AUTO_TEST_SUITE_END()

    BOOST_AUTO_TEST_SUITE( three )
# ifdef USE_GPU_SINGLE_PRECISION
        BOOST_AUTO_TEST_SUITE( type_float )
            typedef halmd::mdsim::gpu::particle<3, float> particle_type;
        TEST_SUITE_GPU(particle_type, dataset, nspecies)
        BOOST_AUTO_TEST_SUITE_END()
# endif
# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
        BOOST_AUTO_TEST_SUITE( type_dsfloat )
            typedef halmd::mdsim::gpu::particle<3, halmd::dsfloat> particle_type;
        TEST_SUITE_GPU(particle_type, dataset, nspecies)
        BOOST_AUTO_TEST_SUITE_END()
# endif
    BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif
