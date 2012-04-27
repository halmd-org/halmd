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

#define BOOST_TEST_MODULE particle
#include <boost/test/unit_test.hpp>

#include <algorithm> // std::fill
#include <boost/array.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <cmath> // std::ceil, std::pow
#include <iterator> // std::back_inserter

#include <halmd/config.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <test/tools/constant_iterator.hpp>
#include <test/tools/ctest.hpp>
#include <test/unit/mdsim/positions/lattice_iterator.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/particle.hpp>
# include <test/tools/cuda.hpp>
#endif

using namespace boost;
using namespace std;

/**
 * Test initialisation, getter and setter of particle positions.
 */
template <typename particle_type>
void particle_position(particle_type& particle)
{
    typedef typename particle_type::position_type position_type;
    typedef typename particle_type::species_type species_type;

    // set species to ascending sequence of integers starting at 1 ≠ 0
    vector<species_type> species;
    species.reserve(particle.nparticle());
    copy(
        counting_iterator<species_type>(1)
      , counting_iterator<species_type>(particle.nparticle() + 1)
      , back_inserter(species)
    );
    particle.set_species(species);

    // check that species are initialised to zero
    vector<position_type> position;
    particle.get_position(position);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        position.begin()
      , position.end()
      , constant_iterator<position_type>(0, 0)
      , constant_iterator<position_type>(0, particle.nparticle())
    );

    // check that setter method validates input array size
    BOOST_CHECK_THROW(
        particle.set_position(vector<position_type>(position.begin(), position.end() - 1))
      , invalid_argument
    );

    // assign square/cubic lattice vectors
    copy(
         lattice_iterator<position_type>(particle.nparticle(), 0)
       , lattice_iterator<position_type>(particle.nparticle(), particle.nparticle())
       , position.begin()
    );
    particle.set_position(position);

    // zero array, and resize array to mismatching, non-zero size
    // to ensure that getter method properly resizes input array
    fill(position.begin(), position.end(), 0);
    position.resize(particle.nparticle() / 2);

    particle.get_position(position);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        position.begin()
      , position.end()
      , lattice_iterator<position_type>(particle.nparticle(), 0)
      , lattice_iterator<position_type>(particle.nparticle(), particle.nparticle())
    );

    // check that particle species are preserved, since positions
    // and species are stored in the same array in gpu::particle
    species.clear();
    particle.get_species(species);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        species.begin()
      , species.end()
      , counting_iterator<species_type>(1)
      , counting_iterator<species_type>(particle.nparticle() + 1)
    );
}

/**
 * Test initialisation, getter and setter of particle images.
 */
template <typename particle_type>
void particle_image(particle_type& particle)
{
    typedef typename particle_type::image_type image_type;

    // check that images are initialised to zero
    vector<image_type> image;
    particle.get_image(image);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        image.begin()
      , image.end()
      , constant_iterator<image_type>(0, 0)
      , constant_iterator<image_type>(0, particle.nparticle())
    );

    // check that setter method validates input array size
    BOOST_CHECK_THROW(
        particle.set_image(vector<image_type>(image.begin(), image.end() - 1))
      , invalid_argument
    );

    // assign square/cubic lattice vectors
    copy(
         lattice_iterator<image_type>(particle.nparticle(), 0)
       , lattice_iterator<image_type>(particle.nparticle(), particle.nparticle())
       , image.begin()
    );
    particle.set_image(image);

    // zero array, and resize array to mismatching, non-zero size
    // to ensure that getter method properly resizes input array
    fill(image.begin(), image.end(), 0);
    image.resize(particle.nparticle() / 2);

    particle.get_image(image);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        image.begin()
      , image.end()
      , lattice_iterator<image_type>(particle.nparticle(), 0)
      , lattice_iterator<image_type>(particle.nparticle(), particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of particle velocities.
 */
template <typename particle_type>
void particle_velocity(particle_type& particle)
{
    typedef typename particle_type::velocity_type velocity_type;
    typedef typename particle_type::mass_type mass_type;

    // set masses to ascending sequence of integers starting at 2 ≠ 1
    vector<mass_type> mass;
    mass.reserve(particle.nparticle());
    copy(
        counting_iterator<mass_type>(2)
      , counting_iterator<mass_type>(particle.nparticle() + 2)
      , back_inserter(mass)
    );
    particle.set_mass(mass);

    // check that velocities are initialised to zero
    vector<velocity_type> velocity;
    particle.get_velocity(velocity);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        velocity.begin()
      , velocity.end()
      , constant_iterator<velocity_type>(0, 0)
      , constant_iterator<velocity_type>(0, particle.nparticle())
    );

    // check that setter method validates input array size
    BOOST_CHECK_THROW(
        particle.set_velocity(vector<velocity_type>(velocity.begin(), velocity.end() - 1))
      , invalid_argument
    );

    // assign square/cubic lattice vectors
    copy(
         lattice_iterator<velocity_type>(particle.nparticle(), 0)
       , lattice_iterator<velocity_type>(particle.nparticle(), particle.nparticle())
       , velocity.begin()
    );
    particle.set_velocity(velocity);

    // zero array, and resize array to mismatching, non-zero size
    // to ensure that getter method properly resizes input array
    fill(velocity.begin(), velocity.end(), 0);
    velocity.resize(particle.nparticle() / 2);

    particle.get_velocity(velocity);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        velocity.begin()
      , velocity.end()
      , lattice_iterator<velocity_type>(particle.nparticle(), 0)
      , lattice_iterator<velocity_type>(particle.nparticle(), particle.nparticle())
    );

    // check that particle masses are preserved, since velocities
    // and masses are stored in the same array in gpu::particle
    mass.clear();
    particle.get_mass(mass);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        mass.begin()
      , mass.end()
      , counting_iterator<mass_type>(2)
      , counting_iterator<mass_type>(particle.nparticle() + 2)
    );
}

/**
 * Test initialisation, getter and setter of particle tags.
 */
template <typename particle_type>
void particle_tag(particle_type& particle)
{
    typedef typename particle_type::tag_type tag_type;

    // check that tags default to ascending sequence of integers
    vector<tag_type> tag;
    particle.get_tag(tag);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        tag.begin()
      , tag.end()
      , counting_iterator<tag_type>(0)
      , counting_iterator<tag_type>(particle.nparticle())
    );

    // check that setter method validates input array size
    BOOST_CHECK_THROW(
        particle.set_tag(vector<tag_type>(tag.begin(), tag.end() - 1))
      , invalid_argument
    );

    // reverse order of tag array copy and set particle tags
    reverse(tag.begin(), tag.end());
    particle.set_tag(tag);

    // zero array, and resize array to mismatching, non-zero size
    // to ensure that getter method properly resizes input array
    fill(tag.begin(), tag.end(), 0);
    tag.resize(particle.nparticle() / 2);

    particle.get_tag(tag);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        tag.rbegin()
      , tag.rend()
      , counting_iterator<tag_type>(0)
      , counting_iterator<tag_type>(particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of particle reverse tags.
 */
template <typename particle_type>
void particle_reverse_tag(particle_type& particle)
{
    typedef typename particle_type::reverse_tag_type reverse_tag_type;

    // check that reverse tags default to ascending sequence of integers
    vector<reverse_tag_type> reverse_tag;
    particle.get_reverse_tag(reverse_tag);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        reverse_tag.begin()
      , reverse_tag.end()
      , counting_iterator<reverse_tag_type>(0)
      , counting_iterator<reverse_tag_type>(particle.nparticle())
    );

    // check that setter method validates input array size
    BOOST_CHECK_THROW(
        particle.set_reverse_tag(vector<reverse_tag_type>(reverse_tag.begin(), reverse_tag.end() - 1))
      , invalid_argument
    );

    // reverse order of reverse tag array copy and set particle tags
    reverse(reverse_tag.begin(), reverse_tag.end());
    particle.set_reverse_tag(reverse_tag);

    // zero array, and resize array to mismatching, non-zero size
    // to ensure that getter method properly resizes input array
    fill(reverse_tag.begin(), reverse_tag.end(), 0);
    reverse_tag.resize(particle.nparticle() / 2);

    particle.get_reverse_tag(reverse_tag);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        reverse_tag.rbegin()
      , reverse_tag.rend()
      , counting_iterator<reverse_tag_type>(0)
      , counting_iterator<reverse_tag_type>(particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of particle species.
 */
template <typename particle_type>
void particle_species(particle_type& particle)
{
    typedef typename particle_type::species_type species_type;
    typedef typename particle_type::position_type position_type;

    // check default of one species
    BOOST_CHECK_EQUAL( particle.nspecies(), 1 );

    // assign square/cubic lattice vectors
    vector<position_type> position;
    position.reserve(particle.nparticle());
    copy(
         lattice_iterator<position_type>(particle.nparticle(), 0)
       , lattice_iterator<position_type>(particle.nparticle(), particle.nparticle())
       , back_inserter(position)
    );
    particle.set_position(position);

    // check that species are initialised to zero
    vector<species_type> species;
    particle.get_species(species);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        species.begin()
      , species.end()
      , constant_iterator<species_type>(0, 0)
      , constant_iterator<species_type>(0, particle.nparticle())
    );

    // check that setter method validates input array size
    BOOST_CHECK_THROW(
        particle.set_species(vector<species_type>(species.begin(), species.end() - 1))
      , invalid_argument
    );

    // set species to ascending sequence of integers starting at 1 ≠ 0
    copy(
        counting_iterator<species_type>(1)
      , counting_iterator<species_type>(particle.nparticle() + 1)
      , species.begin()
    );
    particle.set_species(species);

    // zero array, and resize array to mismatching, non-zero size
    // to ensure that getter method properly resizes input array
    fill(species.begin(), species.end(), 0);
    species.resize(particle.nparticle() / 2);

    particle.get_species(species);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        species.begin()
      , species.end()
      , counting_iterator<species_type>(1)
      , counting_iterator<species_type>(particle.nparticle() + 1)
    );

    // check that particle positions are preserved, since positions
    // and species are stored in the same array in gpu::particle
    position.clear();
    particle.get_position(position);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        position.begin()
      , position.end()
      , lattice_iterator<position_type>(particle.nparticle(), 0)
      , lattice_iterator<position_type>(particle.nparticle(), particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of particle masses.
 */
template <typename particle_type>
void particle_mass(particle_type& particle)
{
    typedef typename particle_type::mass_type mass_type;
    typedef typename particle_type::velocity_type velocity_type;

    // assign square/cubic lattice vectors
    vector<velocity_type> velocity;
    velocity.reserve(particle.nparticle());
    copy(
         lattice_iterator<velocity_type>(particle.nparticle(), 0)
       , lattice_iterator<velocity_type>(particle.nparticle(), particle.nparticle())
       , back_inserter(velocity)
    );
    particle.set_velocity(velocity);

    // check that masses are initialised to unit mass
    vector<mass_type> mass;
    particle.get_mass(mass);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        mass.begin()
      , mass.end()
      , constant_iterator<mass_type>(1, 0)
      , constant_iterator<mass_type>(1, particle.nparticle())
    );

    // check that setter method validates input array size
    BOOST_CHECK_THROW(
        particle.set_mass(vector<mass_type>(mass.begin(), mass.end() - 1))
      , invalid_argument
    );

    // set masses to ascending sequence of integers starting at 2 ≠ 1
    copy(
        counting_iterator<mass_type>(2)
      , counting_iterator<mass_type>(particle.nparticle() + 2)
      , mass.begin()
    );
    particle.set_mass(mass);

    // zero array, and resize array to mismatching, non-zero size
    // to ensure that getter method properly resizes input array
    fill(mass.begin(), mass.end(), 0);
    mass.resize(particle.nparticle() / 2);

    particle.get_mass(mass);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        mass.begin()
      , mass.end()
      , counting_iterator<mass_type>(2)
      , counting_iterator<mass_type>(particle.nparticle() + 2)
    );

    // check that particle velocities are preserved, since velocities
    // and masses are stored in the same array in gpu::particle
    velocity.clear();
    particle.get_velocity(velocity);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        velocity.begin()
      , velocity.end()
      , lattice_iterator<velocity_type>(particle.nparticle(), 0)
      , lattice_iterator<velocity_type>(particle.nparticle(), particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of particle forces.
 */
template <typename particle_type>
void particle_force(particle_type& particle)
{
    typedef typename particle_type::force_type force_type;

    // check that species are initialised to zero
    vector<force_type> force;
    particle.get_force(force);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        force.begin()
      , force.end()
      , constant_iterator<force_type>(0, 0)
      , constant_iterator<force_type>(0, particle.nparticle())
    );

    // check that setter method validates input array size
    BOOST_CHECK_THROW(
        particle.set_force(vector<force_type>(force.begin(), force.end() - 1))
      , invalid_argument
    );

    // assign square/cubic lattice vectors
    copy(
         lattice_iterator<force_type>(particle.nparticle(), 0)
       , lattice_iterator<force_type>(particle.nparticle(), particle.nparticle())
       , force.begin()
    );
    particle.set_force(force);

    // zero array, and resize array to mismatching, non-zero size
    // to ensure that getter method properly resizes input array
    fill(force.begin(), force.end(), 0);
    force.resize(particle.nparticle() / 2);

    particle.get_force(force);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        force.begin()
      , force.end()
      , lattice_iterator<force_type>(particle.nparticle(), 0)
      , lattice_iterator<force_type>(particle.nparticle(), particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of potential energy per particle.
 */
template <typename particle_type>
void particle_en_pot(particle_type& particle)
{
    typedef typename particle_type::en_pot_type en_pot_type;

    // check that potential energies are initialised to zero
    vector<en_pot_type> en_pot;
    particle.get_en_pot(en_pot);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        en_pot.begin()
      , en_pot.end()
      , constant_iterator<en_pot_type>(0, 0)
      , constant_iterator<en_pot_type>(0, particle.nparticle())
    );

    // check that setter method validates input array size
    BOOST_CHECK_THROW(
        particle.set_en_pot(vector<en_pot_type>(en_pot.begin(), en_pot.end() - 1))
      , invalid_argument
    );

    // set potential energies to ascending sequence of integers starting at 1 ≠ 0
    copy(
        counting_iterator<en_pot_type>(1)
      , counting_iterator<en_pot_type>(particle.nparticle() + 1)
      , en_pot.begin()
    );
    particle.set_en_pot(en_pot);

    // zero array, and resize array to mismatching, non-zero size
    // to ensure that getter method properly resizes input array
    fill(en_pot.begin(), en_pot.end(), 0);
    en_pot.resize(particle.nparticle() / 2);

    particle.get_en_pot(en_pot);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        en_pot.begin()
      , en_pot.end()
      , counting_iterator<en_pot_type>(1)
      , counting_iterator<en_pot_type>(particle.nparticle() + 1)
    );
}

/**
 * Test initialisation, getter and setter of potential part of stress tensor per particle.
 */
template <typename particle_type>
void particle_stress_pot(particle_type& particle)
{
    typedef typename particle_type::stress_pot_type stress_pot_type;

    // check that stress tensors are initialised to zero
    vector<stress_pot_type> stress_pot;
    particle.get_stress_pot(stress_pot);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        stress_pot.begin()
      , stress_pot.end()
      , constant_iterator<stress_pot_type>(0, 0)
      , constant_iterator<stress_pot_type>(0, particle.nparticle())
    );

    // check that setter method validates input array size
    BOOST_CHECK_THROW(
        particle.set_stress_pot(vector<stress_pot_type>(stress_pot.begin(), stress_pot.end() - 1))
      , invalid_argument
    );

    // assign square/four-dimensional cubic lattice vectors
    copy(
        lattice_iterator<stress_pot_type>(particle.nparticle(), 0)
      , lattice_iterator<stress_pot_type>(particle.nparticle(), particle.nparticle())
      , stress_pot.begin()
    );
    particle.set_stress_pot(stress_pot);

    // zero array, and resize array to mismatching, non-zero size
    // to ensure that getter method properly resizes input array
    fill(stress_pot.begin(), stress_pot.end(), 0);
    stress_pot.resize(particle.nparticle() / 2);

    particle.get_stress_pot(stress_pot);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        stress_pot.begin()
      , stress_pot.end()
      , lattice_iterator<stress_pot_type>(particle.nparticle(), 0)
      , lattice_iterator<stress_pot_type>(particle.nparticle(), particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of hypervirial per particle.
 */
template <typename particle_type>
void particle_hypervirial(particle_type& particle)
{
    typedef typename particle_type::hypervirial_type hypervirial_type;

    // check that hypervirials are initialised to zero
    vector<hypervirial_type> hypervirial;
    particle.get_hypervirial(hypervirial);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        hypervirial.begin()
      , hypervirial.end()
      , constant_iterator<hypervirial_type>(0, 0)
      , constant_iterator<hypervirial_type>(0, particle.nparticle())
    );

    // check that setter method validates input array size
    BOOST_CHECK_THROW(
        particle.set_hypervirial(vector<hypervirial_type>(hypervirial.begin(), hypervirial.end() - 1))
      , invalid_argument
    );

    // set hypervirials to ascending sequence of integers starting at 1 ≠ 0
    copy(
        counting_iterator<hypervirial_type>(1)
      , counting_iterator<hypervirial_type>(particle.nparticle() + 1)
      , hypervirial.begin()
    );
    particle.set_hypervirial(hypervirial);

    // zero array, and resize array to mismatching, non-zero size
    // to ensure that getter method properly resizes input array
    fill(hypervirial.begin(), hypervirial.end(), 0);
    hypervirial.resize(particle.nparticle() / 2);

    particle.get_hypervirial(hypervirial);
    BOOST_CHECK_EQUAL_COLLECTIONS(
        hypervirial.begin()
      , hypervirial.end()
      , counting_iterator<hypervirial_type>(1)
      , counting_iterator<hypervirial_type>(particle.nparticle() + 1)
    );
}

BOOST_AUTO_TEST_SUITE( host )

/**
 * Fixture that constructs instance of host::particle.
 */
template <unsigned int dimension>
struct particle_fixture
{
    particle_fixture() : particle(12345) {}

#ifdef USE_HOST_SINGLE_PRECISION
    halmd::mdsim::host::particle<dimension, float> particle;
#else
    halmd::mdsim::host::particle<dimension, double> particle;
#endif
};

BOOST_FIXTURE_TEST_SUITE( two, particle_fixture<2> )

BOOST_AUTO_TEST_CASE( position )
{
    particle_position(particle);
}

BOOST_AUTO_TEST_CASE( image )
{
    particle_image(particle);
}

BOOST_AUTO_TEST_CASE( velocity )
{
    particle_velocity(particle);
}

BOOST_AUTO_TEST_CASE( tag )
{
    particle_tag(particle);
}

BOOST_AUTO_TEST_CASE( reverse_tag )
{
    particle_reverse_tag(particle);
}

BOOST_AUTO_TEST_CASE( species )
{
    particle_species(particle);
}

BOOST_AUTO_TEST_CASE( mass )
{
    particle_mass(particle);
}

BOOST_AUTO_TEST_CASE( force )
{
    particle_force(particle);
}

BOOST_AUTO_TEST_CASE( en_pot )
{
    particle_en_pot(particle);
}

BOOST_AUTO_TEST_CASE( stress_pot )
{
    particle_stress_pot(particle);
}

BOOST_AUTO_TEST_CASE( hypervirial )
{
    particle_hypervirial(particle);
}

BOOST_AUTO_TEST_SUITE_END() // two

BOOST_FIXTURE_TEST_SUITE( three, particle_fixture<3> )

BOOST_AUTO_TEST_CASE( position )
{
    particle_position(particle);
}

BOOST_AUTO_TEST_CASE( image )
{
    particle_image(particle);
}

BOOST_AUTO_TEST_CASE( velocity )
{
    particle_velocity(particle);
}

BOOST_AUTO_TEST_CASE( tag )
{
    particle_tag(particle);
}

BOOST_AUTO_TEST_CASE( reverse_tag )
{
    particle_reverse_tag(particle);
}

BOOST_AUTO_TEST_CASE( species )
{
    particle_species(particle);
}

BOOST_AUTO_TEST_CASE( mass )
{
    particle_mass(particle);
}

BOOST_AUTO_TEST_CASE( force )
{
    particle_force(particle);
}

BOOST_AUTO_TEST_CASE( en_pot )
{
    particle_en_pot(particle);
}

BOOST_AUTO_TEST_CASE( stress_pot )
{
    particle_stress_pot(particle);
}

BOOST_AUTO_TEST_CASE( hypervirial )
{
    particle_hypervirial(particle);
}

BOOST_AUTO_TEST_SUITE_END() // three

BOOST_AUTO_TEST_SUITE_END() // host

#ifdef HALMD_WITH_GPU
BOOST_AUTO_TEST_SUITE( gpu )

/**
 * Fixture that constructs instance of gpu::particle.
 */
template <unsigned int dimension>
struct particle_fixture : set_cuda_device
{
    particle_fixture() : particle(12345) {}

    halmd::mdsim::gpu::particle<dimension, float> particle;
};

BOOST_FIXTURE_TEST_SUITE( two, particle_fixture<2> )

BOOST_AUTO_TEST_CASE( position )
{
    particle_position(particle);
}

BOOST_AUTO_TEST_CASE( image )
{
    particle_image(particle);
}

BOOST_AUTO_TEST_CASE( velocity )
{
    particle_velocity(particle);
}

BOOST_AUTO_TEST_CASE( tag )
{
    particle_tag(particle);
}

BOOST_AUTO_TEST_CASE( reverse_tag )
{
    particle_reverse_tag(particle);
}

BOOST_AUTO_TEST_CASE( species )
{
    particle_species(particle);
}

BOOST_AUTO_TEST_CASE( mass )
{
    particle_mass(particle);
}

BOOST_AUTO_TEST_CASE( force )
{
    particle_force(particle);
}

BOOST_AUTO_TEST_CASE( en_pot )
{
    particle_en_pot(particle);
}

BOOST_AUTO_TEST_CASE( stress_pot )
{
    particle_stress_pot(particle);
}

BOOST_AUTO_TEST_CASE( hypervirial )
{
    particle_hypervirial(particle);
}

BOOST_AUTO_TEST_SUITE_END() // two

BOOST_FIXTURE_TEST_SUITE( three, particle_fixture<3> )

BOOST_AUTO_TEST_CASE( position )
{
    particle_position(particle);
}

BOOST_AUTO_TEST_CASE( image )
{
    particle_image(particle);
}

BOOST_AUTO_TEST_CASE( velocity )
{
    particle_velocity(particle);
}

BOOST_AUTO_TEST_CASE( tag )
{
    particle_tag(particle);
}

BOOST_AUTO_TEST_CASE( reverse_tag )
{
    particle_reverse_tag(particle);
}

BOOST_AUTO_TEST_CASE( species )
{
    particle_species(particle);
}

BOOST_AUTO_TEST_CASE( mass )
{
    particle_mass(particle);
}

BOOST_AUTO_TEST_CASE( force )
{
    particle_force(particle);
}

BOOST_AUTO_TEST_CASE( en_pot )
{
    particle_en_pot(particle);
}

BOOST_AUTO_TEST_CASE( stress_pot )
{
    particle_stress_pot(particle);
}

BOOST_AUTO_TEST_CASE( hypervirial )
{
    particle_hypervirial(particle);
}

BOOST_AUTO_TEST_SUITE_END() // three

BOOST_AUTO_TEST_SUITE_END() // gpu
#endif /* HALMD_WITH_GPU */
