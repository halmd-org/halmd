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
    particle.set_species(
        counting_iterator<species_type>(1)
      , counting_iterator<species_type>(particle.nparticle() + 1)
    );

    // check that positions are initialised to zero
    vector<position_type> position;
    position.reserve(particle.nparticle());
    particle.get_position(back_inserter(position));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        position.begin()
      , position.end()
      , constant_iterator<position_type>(0, 0)
      , constant_iterator<position_type>(0, particle.nparticle())
    );

    // assign square/cubic lattice vectors
    particle.set_position(
         lattice_iterator<position_type>(particle.nparticle(), 0)
       , lattice_iterator<position_type>(particle.nparticle(), particle.nparticle())
    );

    position.clear();
    particle.get_position(back_inserter(position));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        position.begin()
      , position.end()
      , lattice_iterator<position_type>(particle.nparticle(), 0)
      , lattice_iterator<position_type>(particle.nparticle(), particle.nparticle())
    );

    // check that particle species are preserved, since positions
    // and species are stored in the same array in gpu::particle
    vector<species_type> species;
    species.reserve(particle.nparticle());
    particle.get_species(back_inserter(species));
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
    image.reserve(particle.nparticle());
    particle.get_image(back_inserter(image));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        image.begin()
      , image.end()
      , constant_iterator<image_type>(0, 0)
      , constant_iterator<image_type>(0, particle.nparticle())
    );

    // assign square/cubic lattice vectors
    particle.set_image(
         lattice_iterator<image_type>(particle.nparticle(), 0)
       , lattice_iterator<image_type>(particle.nparticle(), particle.nparticle())
    );

    image.clear();
    particle.get_image(back_inserter(image));
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
    particle.set_mass(
        counting_iterator<mass_type>(2)
      , counting_iterator<mass_type>(particle.nparticle() + 2)
    );

    // check that velocities are initialised to zero
    vector<velocity_type> velocity;
    velocity.reserve(particle.nparticle());
    particle.get_velocity(back_inserter(velocity));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        velocity.begin()
      , velocity.end()
      , constant_iterator<velocity_type>(0, 0)
      , constant_iterator<velocity_type>(0, particle.nparticle())
    );

    // assign square/cubic lattice vectors
    particle.set_velocity(
         lattice_iterator<velocity_type>(particle.nparticle(), 0)
       , lattice_iterator<velocity_type>(particle.nparticle(), particle.nparticle())
    );

    velocity.clear();
    particle.get_velocity(back_inserter(velocity));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        velocity.begin()
      , velocity.end()
      , lattice_iterator<velocity_type>(particle.nparticle(), 0)
      , lattice_iterator<velocity_type>(particle.nparticle(), particle.nparticle())
    );

    // check that particle masses are preserved, since velocities
    // and masses are stored in the same array in gpu::particle
    vector<mass_type> mass;
    mass.reserve(particle.nparticle());
    particle.get_mass(back_inserter(mass));
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
    tag.reserve(particle.nparticle());
    particle.get_tag(back_inserter(tag));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        tag.begin()
      , tag.end()
      , counting_iterator<tag_type>(0)
      , counting_iterator<tag_type>(particle.nparticle())
    );

    // reverse order of particle tags
    particle.set_tag(tag.rbegin(), tag.rend());

    tag.clear();
    particle.get_tag(back_inserter(tag));
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
    reverse_tag.reserve(particle.nparticle());
    particle.get_reverse_tag(back_inserter(reverse_tag));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        reverse_tag.begin()
      , reverse_tag.end()
      , counting_iterator<reverse_tag_type>(0)
      , counting_iterator<reverse_tag_type>(particle.nparticle())
    );

    // reverse order of reverse particle tags
    particle.set_reverse_tag(reverse_tag.rbegin(), reverse_tag.rend());

    // zero memory before dropping all elements
    // zero array, and resize array to mismatching, non-zero size
    // to ensure that getter method properly resizes input array
    fill(reverse_tag.begin(), reverse_tag.end(), 0);
    reverse_tag.resize(particle.nparticle() / 2);

    reverse_tag.clear();
    particle.get_reverse_tag(back_inserter(reverse_tag));
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
    BOOST_CHECK_EQUAL( particle.nspecies(), 1u );

    // assign square/cubic lattice vectors
    particle.set_position(
         lattice_iterator<position_type>(particle.nparticle(), 0)
       , lattice_iterator<position_type>(particle.nparticle(), particle.nparticle())
    );

    // check that species are initialised to zero
    vector<species_type> species;
    species.reserve(particle.nparticle());
    particle.get_species(back_inserter(species));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        species.begin()
      , species.end()
      , constant_iterator<species_type>(0, 0)
      , constant_iterator<species_type>(0, particle.nparticle())
    );

    // set species to ascending sequence of integers starting at 1 ≠ 0
    particle.set_species(
        counting_iterator<species_type>(1)
      , counting_iterator<species_type>(particle.nparticle() + 1)
    );

    species.clear();
    particle.get_species(back_inserter(species));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        species.begin()
      , species.end()
      , counting_iterator<species_type>(1)
      , counting_iterator<species_type>(particle.nparticle() + 1)
    );

    // check that particle positions are preserved, since positions
    // and species are stored in the same array in gpu::particle
    vector<position_type> position;
    position.reserve(particle.nparticle());
    particle.get_position(back_inserter(position));
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
    particle.set_velocity(
         lattice_iterator<velocity_type>(particle.nparticle(), 0)
       , lattice_iterator<velocity_type>(particle.nparticle(), particle.nparticle())
    );

    // check that masses are initialised to unit mass
    vector<mass_type> mass;
    mass.reserve(particle.nparticle());
    particle.get_mass(back_inserter(mass));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        mass.begin()
      , mass.end()
      , constant_iterator<mass_type>(1, 0)
      , constant_iterator<mass_type>(1, particle.nparticle())
    );

    // set masses to ascending sequence of integers starting at 2 ≠ 1
    particle.set_mass(
        counting_iterator<mass_type>(2)
      , counting_iterator<mass_type>(particle.nparticle() + 2)
    );

    mass.clear();
    particle.get_mass(back_inserter(mass));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        mass.begin()
      , mass.end()
      , counting_iterator<mass_type>(2)
      , counting_iterator<mass_type>(particle.nparticle() + 2)
    );

    // check that particle velocities are preserved, since velocities
    // and masses are stored in the same array in gpu::particle
    vector<velocity_type> velocity;
    velocity.reserve(particle.nparticle());
    particle.get_velocity(back_inserter(velocity));
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

    // check that forces are initialised to zero
    vector<force_type> force;
    force.reserve(particle.nparticle());
    particle.get_force(back_inserter(force));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        force.begin()
      , force.end()
      , constant_iterator<force_type>(0, 0)
      , constant_iterator<force_type>(0, particle.nparticle())
    );

    // assign square/cubic lattice vectors
    particle.set_force(
         lattice_iterator<force_type>(particle.nparticle(), 0)
       , lattice_iterator<force_type>(particle.nparticle(), particle.nparticle())
    );

    force.clear();
    particle.get_force(back_inserter(force));
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
    en_pot.reserve(particle.nparticle());
    particle.get_en_pot(back_inserter(en_pot));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        en_pot.begin()
      , en_pot.end()
      , constant_iterator<en_pot_type>(0, 0)
      , constant_iterator<en_pot_type>(0, particle.nparticle())
    );

    // set potential energies to ascending sequence of integers starting at 1 ≠ 0
    particle.set_en_pot(
        counting_iterator<en_pot_type>(1)
      , counting_iterator<en_pot_type>(particle.nparticle() + 1)
    );

    en_pot.clear();
    particle.get_en_pot(back_inserter(en_pot));
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
    stress_pot.reserve(particle.nparticle());
    particle.get_stress_pot(back_inserter(stress_pot));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        stress_pot.begin()
      , stress_pot.end()
      , constant_iterator<stress_pot_type>(0, 0)
      , constant_iterator<stress_pot_type>(0, particle.nparticle())
    );

    // assign square/four-dimensional cubic lattice vectors
    particle.set_stress_pot(
        lattice_iterator<stress_pot_type>(particle.nparticle(), 0)
      , lattice_iterator<stress_pot_type>(particle.nparticle(), particle.nparticle())
    );

    stress_pot.clear();
    particle.get_stress_pot(back_inserter(stress_pot));
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
    hypervirial.reserve(particle.nparticle());
    particle.get_hypervirial(back_inserter(hypervirial));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        hypervirial.begin()
      , hypervirial.end()
      , constant_iterator<hypervirial_type>(0, 0)
      , constant_iterator<hypervirial_type>(0, particle.nparticle())
    );

    // set hypervirials to ascending sequence of integers starting at 1 ≠ 0
    particle.set_hypervirial(
        counting_iterator<hypervirial_type>(1)
      , counting_iterator<hypervirial_type>(particle.nparticle() + 1)
    );

    hypervirial.clear();
    particle.get_hypervirial(back_inserter(hypervirial));
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
