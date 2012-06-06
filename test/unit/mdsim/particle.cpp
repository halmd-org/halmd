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
      , constant_iterator<position_type>(0, 0)
      , constant_iterator<position_type>(0, particle.nparticle())
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
      , constant_iterator<image_type>(0, 0)
      , constant_iterator<image_type>(0, particle.nparticle())
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
      , constant_iterator<velocity_type>(0, 0)
      , constant_iterator<velocity_type>(0, particle.nparticle())
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
 * Test initialisation, getter and setter of particle tags.
 */
template <typename particle_type>
static void test_tag(particle_type& particle)
{
    typedef typename particle_type::tag_type tag_type;
    particle_type const& const_particle = particle;

    // check that tags default to ascending sequence of integers
    std::vector<tag_type> tag(particle.nparticle());
    BOOST_CHECK(
        get_tag(const_particle, tag.begin()) == tag.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        tag.begin()
      , tag.end()
      , boost::counting_iterator<tag_type>(0)
      , boost::counting_iterator<tag_type>(particle.nparticle())
    );

    // reverse order of particle tags
    BOOST_CHECK(
        set_tag(particle, tag.rbegin()) == tag.rend()
    );
    BOOST_CHECK(
        get_tag(const_particle, tag.begin()) == tag.end()
    );
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
    std::vector<reverse_tag_type> reverse_tag(particle.nparticle());
    BOOST_CHECK(
        get_reverse_tag(const_particle, reverse_tag.begin()) == reverse_tag.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        reverse_tag.begin()
      , reverse_tag.end()
      , boost::counting_iterator<reverse_tag_type>(0)
      , boost::counting_iterator<reverse_tag_type>(particle.nparticle())
    );

    // reverse order of reverse particle tags
    BOOST_CHECK(
        set_reverse_tag(particle, reverse_tag.rbegin()) == reverse_tag.rend()
    );
    BOOST_CHECK(
        get_reverse_tag(const_particle, reverse_tag.begin()) == reverse_tag.end()
    );
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
      , constant_iterator<species_type>(0, 0)
      , constant_iterator<species_type>(0, particle.nparticle())
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
      , constant_iterator<mass_type>(1, 0)
      , constant_iterator<mass_type>(1, particle.nparticle())
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
    particle_type const& const_particle = particle;

    // check that forces are initialised to zero
    std::vector<force_type> force(particle.nparticle());
    BOOST_CHECK(
        get_force(const_particle, force.begin()) == force.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        force.begin()
      , force.end()
      , constant_iterator<force_type>(0, 0)
      , constant_iterator<force_type>(0, particle.nparticle())
    );

    // assign square/cubic lattice vectors
    equilateral_lattice<force_type> lattice(particle.nparticle());
    BOOST_CHECK(
        set_force(particle, make_lattice_iterator(lattice, 0))
            == make_lattice_iterator(lattice, particle.nparticle())
    );
    BOOST_CHECK(
        get_force(const_particle, force.begin()) == force.end()
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
    particle_type const& const_particle = particle;

    // check that potential energies are initialised to zero
    std::vector<en_pot_type> en_pot(particle.nparticle());
    BOOST_CHECK(
        get_en_pot(const_particle, en_pot.begin()) == en_pot.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        en_pot.begin()
      , en_pot.end()
      , constant_iterator<en_pot_type>(0, 0)
      , constant_iterator<en_pot_type>(0, particle.nparticle())
    );

    // set potential energies to ascending sequence of integers starting at 1 ≠ 0
    BOOST_CHECK(
        set_en_pot(particle, boost::counting_iterator<en_pot_type>(1))
            == boost::counting_iterator<en_pot_type>(particle.nparticle() + 1)
    );
    BOOST_CHECK(
        get_en_pot(const_particle, en_pot.begin()) == en_pot.end()
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
    particle_type const& const_particle = particle;

    // check that stress tensors are initialised to zero
    std::vector<stress_pot_type> stress_pot(particle.nparticle());
    BOOST_CHECK(
        get_stress_pot(const_particle, stress_pot.begin()) == stress_pot.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        stress_pot.begin()
      , stress_pot.end()
      , constant_iterator<stress_pot_type>(0, 0)
      , constant_iterator<stress_pot_type>(0, particle.nparticle())
    );

    // assign square/four-dimensional cubic lattice vectors
    equilateral_lattice<stress_pot_type> lattice(particle.nparticle());
    BOOST_CHECK(
        set_stress_pot(particle, make_lattice_iterator(lattice, 0))
            == make_lattice_iterator(lattice, particle.nparticle())
    );
    BOOST_CHECK(
        get_stress_pot(const_particle, stress_pot.begin()) == stress_pot.end()
    );
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
    std::vector<hypervirial_type> hypervirial(particle.nparticle());
    BOOST_CHECK(
        get_hypervirial(const_particle, hypervirial.begin()) == hypervirial.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        hypervirial.begin()
      , hypervirial.end()
      , constant_iterator<hypervirial_type>(0, 0)
      , constant_iterator<hypervirial_type>(0, particle.nparticle())
    );

    // set hypervirials to ascending sequence of integers starting at 1 ≠ 0
    BOOST_CHECK(
        set_hypervirial(particle, boost::counting_iterator<hypervirial_type>(1))
            == boost::counting_iterator<hypervirial_type>(particle.nparticle() + 1)
    );
    BOOST_CHECK(
        get_hypervirial(const_particle, hypervirial.begin()) == hypervirial.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        hypervirial.begin()
      , hypervirial.end()
      , boost::counting_iterator<hypervirial_type>(1)
      , boost::counting_iterator<hypervirial_type>(particle.nparticle() + 1)
    );
}

template <typename particle_type>
static void
test_suite_host(std::size_t nparticle, unsigned int nspecies, boost::unit_test::test_suite* ts)
{
    auto position = [=]() {
        particle_type particle(nparticle, nspecies);
        test_position(particle);
    };
    ts->add(BOOST_TEST_CASE( position ));

    auto image = [=]() {
        particle_type particle(nparticle, nspecies);
        test_image(particle);
    };
    ts->add(BOOST_TEST_CASE( image ));

    auto velocity = [=]() {
        particle_type particle(nparticle, nspecies);
        test_velocity(particle);
    };
    ts->add(BOOST_TEST_CASE( velocity ));

    auto tag = [=]() {
        particle_type particle(nparticle, nspecies);
        test_tag(particle);
    };
    ts->add(BOOST_TEST_CASE( tag ));

    auto reverse_tag = [=]() {
        particle_type particle(nparticle, nspecies);
        test_reverse_tag(particle);
    };
    ts->add(BOOST_TEST_CASE( reverse_tag ));

    auto species = [=]() {
        particle_type particle(nparticle, nspecies);
        test_species(particle);
    };
    ts->add(BOOST_TEST_CASE( species ));

    auto mass = [=]() {
        particle_type particle(nparticle, nspecies);
        test_mass(particle);
    };
    ts->add(BOOST_TEST_CASE( mass ));

    auto force = [=]() {
        particle_type particle(nparticle, nspecies);
        test_force(particle);
    };
    ts->add(BOOST_TEST_CASE( force ));

    auto en_pot = [=]() {
        particle_type particle(nparticle, nspecies);
        test_en_pot(particle);
    };
    ts->add(BOOST_TEST_CASE( en_pot ));

    auto stress_pot = [=]() {
        particle_type particle(nparticle, nspecies);
        test_stress_pot(particle);
    };
    ts->add(BOOST_TEST_CASE( stress_pot ));

    auto hypervirial = [=]() {
        particle_type particle(nparticle, nspecies);
        test_hypervirial(particle);
    };
    ts->add(BOOST_TEST_CASE( hypervirial ));
}

#ifdef HALMD_WITH_GPU
template <typename particle_type>
static void
test_suite_gpu(std::size_t nparticle, unsigned int nspecies, boost::unit_test::test_suite* ts)
{
    auto position = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle, nspecies);
        test_position(particle);
    };
    ts->add(BOOST_TEST_CASE( position ));

    auto image = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle, nspecies);
        test_image(particle);
    };
    ts->add(BOOST_TEST_CASE( image ));

    auto velocity = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle, nspecies);
        test_velocity(particle);
    };
    ts->add(BOOST_TEST_CASE( velocity ));

    auto tag = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle, nspecies);
        test_tag(particle);
    };
    ts->add(BOOST_TEST_CASE( tag ));

    auto reverse_tag = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle, nspecies);
        test_reverse_tag(particle);
    };
    ts->add(BOOST_TEST_CASE( reverse_tag ));

    auto species = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle, nspecies);
        test_species(particle);
    };
    ts->add(BOOST_TEST_CASE( species ));

    auto mass = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle, nspecies);
        test_mass(particle);
    };
    ts->add(BOOST_TEST_CASE( mass ));

    auto force = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle, nspecies);
        test_force(particle);
    };
    ts->add(BOOST_TEST_CASE( force ));

    auto en_pot = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle, nspecies);
        test_en_pot(particle);
    };
    ts->add(BOOST_TEST_CASE( en_pot ));

    auto stress_pot = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle, nspecies);
        test_stress_pot(particle);
    };
    ts->add(BOOST_TEST_CASE( stress_pot ));

    auto hypervirial = [=]() {
        set_cuda_device device;
        particle_type particle(nparticle, nspecies);
        test_hypervirial(particle);
    };
    ts->add(BOOST_TEST_CASE( hypervirial ));
}
#endif

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

    unsigned int const nspecies = 1;

    for (unsigned int nparticle : {109, 4789, 42589}) {
#ifdef USE_HOST_SINGLE_PRECISION
        test_suite_host<halmd::mdsim::host::particle<3, float> >(nparticle, nspecies, ts_host_three);
        test_suite_host<halmd::mdsim::host::particle<2, float> >(nparticle, nspecies, ts_host_two);
#else
        test_suite_host<halmd::mdsim::host::particle<3, double> >(nparticle, nspecies, ts_host_three);
        test_suite_host<halmd::mdsim::host::particle<2, double> >(nparticle, nspecies, ts_host_two);
#endif
#ifdef HALMD_WITH_GPU
        test_suite_gpu<halmd::mdsim::gpu::particle<3, float> >(nparticle, nspecies, ts_gpu_three);
        test_suite_gpu<halmd::mdsim::gpu::particle<2, float> >(nparticle, nspecies, ts_gpu_two);
#endif
    }
}
