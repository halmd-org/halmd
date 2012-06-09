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

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE force
#include <boost/test/unit_test.hpp>

#include <halmd/mdsim/host/force.hpp>
#include <halmd/mdsim/positions/lattice_primitive.hpp>
#include <test/tools/constant_iterator.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/init.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/force.hpp>
# include <test/tools/cuda.hpp>
#endif

#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#ifdef HALMD_WITH_GPU
# include <cuda_wrapper/cuda_wrapper.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <type_traits>

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
    cuda::host::vector<typename array_type::value_type> h_output(g_output.size());
    std::copy(first, last, h_output.begin());
    cuda::copy(h_output.begin(), h_output.end(), g_output.begin());
    array = std::move(g_output);
}
#endif

/**
 * Force with members initialised from iterator range.
 */
template <typename force_type>
class force_from_iterator_range
  : public force_type
{
public:
    typedef typename force_type::net_force_array_type net_force_array_type;
    typedef typename force_type::en_pot_array_type en_pot_array_type;
    typedef typename force_type::stress_pot_array_type stress_pot_array_type;
    typedef typename force_type::hypervirial_array_type hypervirial_array_type;

    virtual halmd::cache<net_force_array_type> const& net_force()
    {
        return net_force_;
    }

    template <typename iterator_type>
    void set_net_force(iterator_type const& first, iterator_type const& last)
    {
        halmd::cache_proxy<net_force_array_type> net_force = net_force_;
        make_array_from_iterator_range(*net_force, first, last);
    }

    virtual halmd::cache<en_pot_array_type> const& en_pot()
    {
        return en_pot_;
    }

    template <typename iterator_type>
    void set_en_pot(iterator_type const& first, iterator_type const& last)
    {
        halmd::cache_proxy<en_pot_array_type> en_pot = en_pot_;
        make_array_from_iterator_range(*en_pot, first, last);
    }

    virtual halmd::cache<stress_pot_array_type> const& stress_pot()
    {
        return stress_pot_;
    }

    template <typename iterator_type>
    void set_stress_pot(iterator_type const& first, iterator_type const& last)
    {
        halmd::cache_proxy<stress_pot_array_type> stress_pot = stress_pot_;
        make_array_from_iterator_range(*stress_pot, first, last);
    }

    virtual halmd::cache<hypervirial_array_type> const& hypervirial()
    {
        return hypervirial_;
    }

    template <typename iterator_type>
    void set_hypervirial(iterator_type const& first, iterator_type const& last)
    {
        halmd::cache_proxy<hypervirial_array_type> hypervirial = hypervirial_;
        make_array_from_iterator_range(*hypervirial, first, last);
    }

private:
    halmd::cache<net_force_array_type> net_force_;
    halmd::cache<en_pot_array_type> en_pot_;
    halmd::cache<stress_pot_array_type> stress_pot_;
    halmd::cache<hypervirial_array_type> hypervirial_;
};

/**
 * Test initialisation, getter and setter of particle forces.
 */
template <typename force_type>
static void test_net_force(unsigned int nparticle)
{
    typedef typename force_type::net_force_type net_force_type;
    force_from_iterator_range<force_type> force;

    // initialise net forces to zero
    force.set_net_force(
        make_constant_iterator(net_force_type(0), 0)
      , make_constant_iterator(net_force_type(0), nparticle)
    );
    std::vector<net_force_type> net_force(nparticle);
    BOOST_CHECK( get_net_force(
        force
      , net_force.begin()) == net_force.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        net_force.begin()
      , net_force.end()
      , make_constant_iterator(net_force_type(0), 0)
      , make_constant_iterator(net_force_type(0), nparticle)
    );

    // assign square/cubic lattice vectors
    equilateral_lattice<net_force_type> lattice(nparticle);
    force.set_net_force(
        make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, nparticle)
    );
    BOOST_CHECK( get_net_force(
        force
      , net_force.begin()) == net_force.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        net_force.begin()
      , net_force.end()
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, nparticle)
    );
}

/**
 * Test initialisation, getter and setter of potential energy per particle.
 */
template <typename force_type>
static void test_en_pot(unsigned int nparticle)
{
    typedef typename force_type::en_pot_type en_pot_type;
    force_from_iterator_range<force_type> force;

    // initialise potential energies to zero
    force.set_en_pot(
        make_constant_iterator(en_pot_type(0), 0)
      , make_constant_iterator(en_pot_type(0), nparticle)
    );
    std::vector<en_pot_type> en_pot(nparticle);
    BOOST_CHECK( get_en_pot(
        force
      , en_pot.begin()) == en_pot.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        en_pot.begin()
      , en_pot.end()
      , make_constant_iterator(en_pot_type(0), 0)
      , make_constant_iterator(en_pot_type(0), nparticle)
    );

    // set potential energies to ascending sequence of integers starting at 1 ≠ 0
    force.set_en_pot(
        boost::counting_iterator<en_pot_type>(1)
      , boost::counting_iterator<en_pot_type>(nparticle + 1)
    );
    BOOST_CHECK( get_en_pot(
        force
      , en_pot.begin()) == en_pot.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        en_pot.begin()
      , en_pot.end()
      , boost::counting_iterator<en_pot_type>(1)
      , boost::counting_iterator<en_pot_type>(nparticle + 1)
    );
}

/**
 * Test initialisation, getter and setter of potential part of stress tensor per particle.
 */
template <typename force_type>
static void test_stress_pot(unsigned int nparticle)
{
    typedef typename force_type::stress_pot_type stress_pot_type;
    force_from_iterator_range<force_type> force;

    // initialise stress tensors to zero
    force.set_stress_pot(
        make_constant_iterator(stress_pot_type(0), 0)
      , make_constant_iterator(stress_pot_type(0), nparticle)
    );
    std::vector<stress_pot_type> stress_pot(nparticle);
    BOOST_CHECK( get_stress_pot(
        force
      , stress_pot.begin()) == stress_pot.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        stress_pot.begin()
      , stress_pot.end()
      , make_constant_iterator(stress_pot_type(0), 0)
      , make_constant_iterator(stress_pot_type(0), nparticle)
    );

    // assign square/four-dimensional cubic lattice vectors
    equilateral_lattice<stress_pot_type> lattice(nparticle);
    force.set_stress_pot(
        make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, nparticle)
    );
    BOOST_CHECK( get_stress_pot(
        force
      , stress_pot.begin()) == stress_pot.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        stress_pot.begin()
      , stress_pot.end()
      , make_lattice_iterator(lattice, 0)
      , make_lattice_iterator(lattice, nparticle)
    );
}

/**
 * Test initialisation, getter and setter of hypervirial per particle.
 */
template <typename force_type>
static void test_hypervirial(unsigned int nparticle)
{
    typedef typename force_type::hypervirial_type hypervirial_type;
    force_from_iterator_range<force_type> force;

    // initialise hypervirials to zero
    force.set_hypervirial(
        make_constant_iterator(hypervirial_type(0), 0)
      , make_constant_iterator(hypervirial_type(0), nparticle)
    );
    std::vector<hypervirial_type> hypervirial(nparticle);
    BOOST_CHECK( get_hypervirial(
        force
      , hypervirial.begin()) == hypervirial.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        hypervirial.begin()
      , hypervirial.end()
      , make_constant_iterator(hypervirial_type(0), 0)
      , make_constant_iterator(hypervirial_type(0), nparticle)
    );

    // set hypervirials to ascending sequence of integers starting at 1 ≠ 0
    force.set_hypervirial(
        boost::counting_iterator<hypervirial_type>(1)
      , boost::counting_iterator<hypervirial_type>(nparticle + 1)
    );
    BOOST_CHECK( get_hypervirial(
        force
      , hypervirial.begin()) == hypervirial.end()
    );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        hypervirial.begin()
      , hypervirial.end()
      , boost::counting_iterator<hypervirial_type>(1)
      , boost::counting_iterator<hypervirial_type>(nparticle + 1)
    );
}

template <typename force_type>
static void
test_suite_host(unsigned int nparticle, boost::unit_test::test_suite* ts)
{
    auto net_force = [=]() {
        test_net_force<force_type>(nparticle);
    };
    ts->add(BOOST_TEST_CASE( net_force ));

    auto en_pot = [=]() {
        test_en_pot<force_type>(nparticle);
    };
    ts->add(BOOST_TEST_CASE( en_pot ));

    auto stress_pot = [=]() {
        test_stress_pot<force_type>(nparticle);
    };
    ts->add(BOOST_TEST_CASE( stress_pot ));

    auto hypervirial = [=]() {
        test_hypervirial<force_type>(nparticle);
    };
    ts->add(BOOST_TEST_CASE( hypervirial ));
}

#ifdef HALMD_WITH_GPU
template <typename force_type>
static void
test_suite_gpu(unsigned int nparticle, boost::unit_test::test_suite* ts)
{
    auto net_force = [=]() {
        set_cuda_device device;
        test_net_force<force_type>(nparticle);
    };
    ts->add(BOOST_TEST_CASE( net_force ));

    auto en_pot = [=]() {
        set_cuda_device device;
        test_en_pot<force_type>(nparticle);
    };
    ts->add(BOOST_TEST_CASE( en_pot ));

    auto stress_pot = [=]() {
        set_cuda_device device;
        test_stress_pot<force_type>(nparticle);
    };
    ts->add(BOOST_TEST_CASE( stress_pot ));

    auto hypervirial = [=]() {
        set_cuda_device device;
        test_hypervirial<force_type>(nparticle);
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

    for (unsigned int nparticle : {109, 4789, 42589}) {
#ifdef USE_HOST_SINGLE_PRECISION
        test_suite_host<halmd::mdsim::host::force<3, float>>(nparticle, ts_host_three);
        test_suite_host<halmd::mdsim::host::force<2, float>>(nparticle, ts_host_two);
#else
        test_suite_host<halmd::mdsim::host::force<3, double>>(nparticle, ts_host_three);
        test_suite_host<halmd::mdsim::host::force<2, double>>(nparticle, ts_host_two);
#endif
#ifdef HALMD_WITH_GPU
        test_suite_gpu<halmd::mdsim::gpu::force<3, float>>(nparticle, ts_gpu_three);
        test_suite_gpu<halmd::mdsim::gpu::force<2, float>>(nparticle, ts_gpu_two);
#endif
    }
}
