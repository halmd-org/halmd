/*
 * Copyright © 2012 Peter Colberg
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

#include <algorithm>
#include <boost/array.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <cmath>
#include <iterator>

#include <halmd/config.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/positions/lattice_primitive.hpp>
#include <test/tools/constant_iterator.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/init.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/particle.hpp>
# include <test/tools/cuda.hpp>
#endif

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
 * Test initialisation, getter and setter of particle positions.
 */
template <typename particle_type>
static void test_position(particle_type& particle)
{
    typedef typename particle_type::position_type position_type;
    typedef typename particle_type::species_type species_type;
    particle_type const& const_particle = particle;

    // set species to ascending sequence of integers starting at 1 ≠ 0
    particle.set_species(
        boost::counting_iterator<species_type>(1)
      , boost::counting_iterator<species_type>(particle.nparticle() + 1)
    );

    // check that positions are initialised to zero
    std::vector<position_type> position;
    position.reserve(particle.nparticle());
    const_particle.get_position(back_inserter(position));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        position.begin()
      , position.end()
      , constant_iterator<position_type>(0, 0)
      , constant_iterator<position_type>(0, particle.nparticle())
    );

    // assign square/cubic lattice vectors
    equilateral_lattice<position_type> lattice(particle.nparticle());
    particle.set_position(
         make_lattice_iterator(lattice, 0)
       , make_lattice_iterator(lattice, particle.nparticle())
    );

    position.clear();
    const_particle.get_position(back_inserter(position));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        position.begin()
      , position.end()
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, particle.nparticle())
    );

    // check that particle species are preserved, since positions
    // and species are stored in the same array in gpu::particle
    std::vector<species_type> species;
    species.reserve(particle.nparticle());
    const_particle.get_species(back_inserter(species));
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
    std::vector<image_type> image;
    image.reserve(particle.nparticle());
    const_particle.get_image(back_inserter(image));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        image.begin()
      , image.end()
      , constant_iterator<image_type>(0, 0)
      , constant_iterator<image_type>(0, particle.nparticle())
    );

    // assign square/cubic lattice vectors
    equilateral_lattice<image_type> lattice(particle.nparticle());
    particle.set_image(
         make_lattice_iterator(lattice, 0)
       , make_lattice_iterator(lattice, particle.nparticle())
    );

    image.clear();
    const_particle.get_image(back_inserter(image));
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
    particle.set_mass(
        boost::counting_iterator<mass_type>(2)
      , boost::counting_iterator<mass_type>(particle.nparticle() + 2)
    );

    // check that velocities are initialised to zero
    std::vector<velocity_type> velocity;
    velocity.reserve(particle.nparticle());
    const_particle.get_velocity(back_inserter(velocity));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        velocity.begin()
      , velocity.end()
      , constant_iterator<velocity_type>(0, 0)
      , constant_iterator<velocity_type>(0, particle.nparticle())
    );

    // assign square/cubic lattice vectors
    equilateral_lattice<velocity_type> lattice(particle.nparticle());
    particle.set_velocity(
         make_lattice_iterator(lattice, 0)
       , make_lattice_iterator(lattice, particle.nparticle())
    );

    velocity.clear();
    const_particle.get_velocity(back_inserter(velocity));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        velocity.begin()
      , velocity.end()
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, particle.nparticle())
    );

    // check that particle masses are preserved, since velocities
    // and masses are stored in the same array in gpu::particle
    std::vector<mass_type> mass;
    mass.reserve(particle.nparticle());
    const_particle.get_mass(back_inserter(mass));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        mass.begin()
      , mass.end()
      , boost::counting_iterator<mass_type>(2)
      , boost::counting_iterator<mass_type>(particle.nparticle() + 2)
    );
}

/**
 * Test initialisation, getter and setter of particle tags.
 */
template <typename particle_type>
static void test_tag(particle_type& particle)
{
    typedef typename particle_type::tag_type tag_type;
    particle_type const& const_particle = particle;

    // check that tags default to ascending sequence of integers
    std::vector<tag_type> tag;
    tag.reserve(particle.nparticle());
    const_particle.get_tag(back_inserter(tag));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        tag.begin()
      , tag.end()
      , boost::counting_iterator<tag_type>(0)
      , boost::counting_iterator<tag_type>(particle.nparticle())
    );

    // reverse order of particle tags
    particle.set_tag(tag.rbegin(), tag.rend());

    tag.clear();
    const_particle.get_tag(back_inserter(tag));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        tag.rbegin()
      , tag.rend()
      , boost::counting_iterator<tag_type>(0)
      , boost::counting_iterator<tag_type>(particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of particle reverse tags.
 */
template <typename particle_type>
static void test_reverse_tag(particle_type& particle)
{
    typedef typename particle_type::reverse_tag_type reverse_tag_type;
    particle_type const& const_particle = particle;

    // check that reverse tags default to ascending sequence of integers
    std::vector<reverse_tag_type> reverse_tag;
    reverse_tag.reserve(particle.nparticle());
    const_particle.get_reverse_tag(back_inserter(reverse_tag));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        reverse_tag.begin()
      , reverse_tag.end()
      , boost::counting_iterator<reverse_tag_type>(0)
      , boost::counting_iterator<reverse_tag_type>(particle.nparticle())
    );

    // reverse order of reverse particle tags
    particle.set_reverse_tag(reverse_tag.rbegin(), reverse_tag.rend());

    // zero memory before dropping all elements
    // zero array, and resize array to mismatching, non-zero size
    // to ensure that getter method properly resizes input array
    fill(reverse_tag.begin(), reverse_tag.end(), 0);
    reverse_tag.resize(particle.nparticle() / 2);

    reverse_tag.clear();
    const_particle.get_reverse_tag(back_inserter(reverse_tag));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        reverse_tag.rbegin()
      , reverse_tag.rend()
      , boost::counting_iterator<reverse_tag_type>(0)
      , boost::counting_iterator<reverse_tag_type>(particle.nparticle())
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
    particle.set_position(
         make_lattice_iterator(lattice, 0)
       , make_lattice_iterator(lattice, particle.nparticle())
    );

    // check that species are initialised to zero
    std::vector<species_type> species;
    species.reserve(particle.nparticle());
    const_particle.get_species(back_inserter(species));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        species.begin()
      , species.end()
      , constant_iterator<species_type>(0, 0)
      , constant_iterator<species_type>(0, particle.nparticle())
    );

    // set species to ascending sequence of integers starting at 1 ≠ 0
    particle.set_species(
        boost::counting_iterator<species_type>(1)
      , boost::counting_iterator<species_type>(particle.nparticle() + 1)
    );

    species.clear();
    const_particle.get_species(back_inserter(species));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        species.begin()
      , species.end()
      , boost::counting_iterator<species_type>(1)
      , boost::counting_iterator<species_type>(particle.nparticle() + 1)
    );

    // check that particle positions are preserved, since positions
    // and species are stored in the same array in gpu::particle
    std::vector<position_type> position;
    position.reserve(particle.nparticle());
    const_particle.get_position(back_inserter(position));
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
    particle.set_velocity(
         make_lattice_iterator(lattice, 0)
       , make_lattice_iterator(lattice, particle.nparticle())
    );

    // check that masses are initialised to unit mass
    std::vector<mass_type> mass;
    mass.reserve(particle.nparticle());
    const_particle.get_mass(back_inserter(mass));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        mass.begin()
      , mass.end()
      , constant_iterator<mass_type>(1, 0)
      , constant_iterator<mass_type>(1, particle.nparticle())
    );

    // set masses to ascending sequence of integers starting at 2 ≠ 1
    particle.set_mass(
        boost::counting_iterator<mass_type>(2)
      , boost::counting_iterator<mass_type>(particle.nparticle() + 2)
    );

    mass.clear();
    const_particle.get_mass(back_inserter(mass));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        mass.begin()
      , mass.end()
      , boost::counting_iterator<mass_type>(2)
      , boost::counting_iterator<mass_type>(particle.nparticle() + 2)
    );

    // check that particle velocities are preserved, since velocities
    // and masses are stored in the same array in gpu::particle
    std::vector<velocity_type> velocity;
    velocity.reserve(particle.nparticle());
    const_particle.get_velocity(back_inserter(velocity));
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
    particle_type const& const_particle = particle;

    // check that forces are initialised to zero
    std::vector<force_type> force;
    force.reserve(particle.nparticle());
    const_particle.get_force(back_inserter(force));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        force.begin()
      , force.end()
      , constant_iterator<force_type>(0, 0)
      , constant_iterator<force_type>(0, particle.nparticle())
    );

    // assign square/cubic lattice vectors
    equilateral_lattice<force_type> lattice(particle.nparticle());
    particle.set_force(
         make_lattice_iterator(lattice, 0)
       , make_lattice_iterator(lattice, particle.nparticle())
    );

    force.clear();
    const_particle.get_force(back_inserter(force));
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
    particle_type const& const_particle = particle;

    // check that potential energies are initialised to zero
    std::vector<en_pot_type> en_pot;
    en_pot.reserve(particle.nparticle());
    const_particle.get_en_pot(back_inserter(en_pot));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        en_pot.begin()
      , en_pot.end()
      , constant_iterator<en_pot_type>(0, 0)
      , constant_iterator<en_pot_type>(0, particle.nparticle())
    );

    // set potential energies to ascending sequence of integers starting at 1 ≠ 0
    particle.set_en_pot(
        boost::counting_iterator<en_pot_type>(1)
      , boost::counting_iterator<en_pot_type>(particle.nparticle() + 1)
    );

    en_pot.clear();
    const_particle.get_en_pot(back_inserter(en_pot));
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
    particle_type const& const_particle = particle;

    // check that stress tensors are initialised to zero
    std::vector<stress_pot_type> stress_pot;
    stress_pot.reserve(particle.nparticle());
    const_particle.get_stress_pot(back_inserter(stress_pot));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        stress_pot.begin()
      , stress_pot.end()
      , constant_iterator<stress_pot_type>(0, 0)
      , constant_iterator<stress_pot_type>(0, particle.nparticle())
    );

    // assign square/four-dimensional cubic lattice vectors
    equilateral_lattice<stress_pot_type> lattice(particle.nparticle());
    particle.set_stress_pot(
        make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, particle.nparticle())
    );

    stress_pot.clear();
    const_particle.get_stress_pot(back_inserter(stress_pot));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        stress_pot.begin()
      , stress_pot.end()
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, particle.nparticle())
    );
}

/**
 * Test initialisation, getter and setter of hypervirial per particle.
 */
template <typename particle_type>
static void test_hypervirial(particle_type& particle)
{
    typedef typename particle_type::hypervirial_type hypervirial_type;
    particle_type const& const_particle = particle;

    // check that hypervirials are initialised to zero
    std::vector<hypervirial_type> hypervirial;
    hypervirial.reserve(particle.nparticle());
    const_particle.get_hypervirial(back_inserter(hypervirial));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        hypervirial.begin()
      , hypervirial.end()
      , constant_iterator<hypervirial_type>(0, 0)
      , constant_iterator<hypervirial_type>(0, particle.nparticle())
    );

    // set hypervirials to ascending sequence of integers starting at 1 ≠ 0
    particle.set_hypervirial(
        boost::counting_iterator<hypervirial_type>(1)
      , boost::counting_iterator<hypervirial_type>(particle.nparticle() + 1)
    );

    hypervirial.clear();
    const_particle.get_hypervirial(back_inserter(hypervirial));
    BOOST_CHECK_EQUAL_COLLECTIONS(
        hypervirial.begin()
      , hypervirial.end()
      , boost::counting_iterator<hypervirial_type>(1)
      , boost::counting_iterator<hypervirial_type>(particle.nparticle() + 1)
    );
}

template <typename particle_type>
static void
test_suite_host(std::size_t nparticle, boost::unit_test::test_suite* ts)
{
    auto position = [=]() {
        particle_type particle(nparticle);
        test_position(particle);
    };
    ts->add(BOOST_TEST_CASE( position ));

    auto image = [=]() {
        particle_type particle(nparticle);
        test_image(particle);
    };
    ts->add(BOOST_TEST_CASE( image ));

    auto velocity = [=]() {
        particle_type particle(nparticle);
        test_velocity(particle);
    };
    ts->add(BOOST_TEST_CASE( velocity ));

    auto tag = [=]() {
        particle_type particle(nparticle);
        test_tag(particle);
    };
    ts->add(BOOST_TEST_CASE( tag ));

    auto reverse_tag = [=]() {
        particle_type particle(nparticle);
        test_reverse_tag(particle);
    };
    ts->add(BOOST_TEST_CASE( reverse_tag ));

    auto species = [=]() {
        particle_type particle(nparticle);
        test_species(particle);
    };
    ts->add(BOOST_TEST_CASE( species ));

    auto mass = [=]() {
        particle_type particle(nparticle);
        test_mass(particle);
    };
    ts->add(BOOST_TEST_CASE( mass ));

    auto force = [=]() {
        particle_type particle(nparticle);
        test_force(particle);
    };
    ts->add(BOOST_TEST_CASE( force ));

    auto en_pot = [=]() {
        particle_type particle(nparticle);
        test_en_pot(particle);
    };
    ts->add(BOOST_TEST_CASE( en_pot ));

    auto stress_pot = [=]() {
        particle_type particle(nparticle);
        test_stress_pot(particle);
    };
    ts->add(BOOST_TEST_CASE( stress_pot ));

    auto hypervirial = [=]() {
        particle_type particle(nparticle);
        test_hypervirial(particle);
    };
    ts->add(BOOST_TEST_CASE( hypervirial ));
}

template <typename particle_type>
static void
test_suite_gpu(std::size_t nparticle, boost::unit_test::test_suite* ts)
{
    auto position = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle);
        test_position(particle);
    };
    ts->add(BOOST_TEST_CASE( position ));

    auto image = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle);
        test_image(particle);
    };
    ts->add(BOOST_TEST_CASE( image ));

    auto velocity = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle);
        test_velocity(particle);
    };
    ts->add(BOOST_TEST_CASE( velocity ));

    auto tag = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle);
        test_tag(particle);
    };
    ts->add(BOOST_TEST_CASE( tag ));

    auto reverse_tag = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle);
        test_reverse_tag(particle);
    };
    ts->add(BOOST_TEST_CASE( reverse_tag ));

    auto species = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle);
        test_species(particle);
    };
    ts->add(BOOST_TEST_CASE( species ));

    auto mass = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle);
        test_mass(particle);
    };
    ts->add(BOOST_TEST_CASE( mass ));

    auto force = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle);
        test_force(particle);
    };
    ts->add(BOOST_TEST_CASE( force ));

    auto en_pot = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle);
        test_en_pot(particle);
    };
    ts->add(BOOST_TEST_CASE( en_pot ));

    auto stress_pot = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle);
        test_stress_pot(particle);
    };
    ts->add(BOOST_TEST_CASE( stress_pot ));

    auto hypervirial = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle);
        test_hypervirial(particle);
    };
    ts->add(BOOST_TEST_CASE( hypervirial ));
}

HALMD_TEST_INIT( particle )
{
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

    for (unsigned int nparticle : {109, 4789, 42589}) {
#ifdef USE_HOST_SINGLE_PRECISION
        test_suite_host<halmd::mdsim::host::particle<3, float> >(nparticle, ts_host_three);
        test_suite_host<halmd::mdsim::host::particle<2, float> >(nparticle, ts_host_two);
#else
        test_suite_host<halmd::mdsim::host::particle<3, double> >(nparticle, ts_host_three);
        test_suite_host<halmd::mdsim::host::particle<2, double> >(nparticle, ts_host_two);
#endif
#ifdef HALMD_WITH_GPU
        test_suite_gpu<halmd::mdsim::gpu::particle<3, float> >(nparticle, ts_gpu_three);
        test_suite_gpu<halmd::mdsim::gpu::particle<2, float> >(nparticle, ts_gpu_two);
#endif
    }
}
