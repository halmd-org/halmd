/*
 * Copyright © 2015      Sutapa Roy
 * Copyright © 2011-2015 Felix Höfling
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

#define BOOST_TEST_MODULE planar_wall
#include <boost/test/unit_test.hpp>

#include <boost/numeric/ublas/assignment.hpp> // <<=
#include <boost/numeric/ublas/banded.hpp> // diagonal_matrix
#include <boost/numeric/ublas/io.hpp>
#include <cmath> // std::pow
#include <limits>
#include <numeric> // std::accumulate

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/potentials/external/planar_wall.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/forces/external.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/potentials/external/planar_wall.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif
#include <test/tools/ctest.hpp>
#include <test/tools/dsfloat.hpp>

using namespace halmd;

/** test planar_wall potential
 *
 *  The host module is a conventional functor which can be tested directly. For
 *  the GPU module, we use the external force module to compute some values of
 *  the potential which are compared against the host module. This requires a
 *  special neighbour list module with only one defined neighbour per particle.
 */

template <typename vector_type>
std::tuple<vector_type, double> compute_wall(
    vector_type const& r
  , vector_type const& surface_normal, double offset
  , double epsilon, double sigma, double wetting
  , double cutoff
  , double smoothing
)
{
    // distance to the planar wall
    double d = inner_prod(r, surface_normal) - offset;
    double sign = (d > 0) ? 1 : -1;
    d = std::abs(d);

    if (d > cutoff) {
        return std::make_tuple(vector_type(0), 0);
    }

    // compute magnitude of force and potential energy (shifted at the cutoff)
    double fval = 3 * sign * epsilon / sigma * ((2. / 5) * std::pow(sigma / d, 10) - wetting * std::pow(sigma / d, 4));
    double en_pot = epsilon * ((2. / 15) * std::pow(sigma / d, 9)- wetting * std::pow(sigma / d, 3));
    en_pot -= epsilon * ((2. / 15) * std::pow(sigma / cutoff, 9) - wetting * std::pow(sigma / cutoff, 3));

    // apply smooth truncation
    double xi = (d - cutoff) / smoothing;
    double w = std::pow(xi, 4) / (1 + std::pow(xi, 4));
    double w1 = 4 * std::pow(xi, 3) / std::pow((1 + std::pow(xi, 4)), 2);
    fval = fval * w - w1 * en_pot / smoothing;
    en_pot *= w;

    return std::make_tuple(fval * surface_normal, en_pot);
}

template <typename vector_type>
std::vector<vector_type> make_positions()
{
    // positions of test particles
    std::vector<vector_type> positions;
    positions.push_back(vector_type({ 0.0,  0.0,  0.0}));
    positions.push_back(vector_type({ 5.0,  0.0,  1.0}));
    positions.push_back(vector_type({-4.0,  0.0,  2.0}));
    positions.push_back(vector_type({-5.7,  0.0,  3.0}));
    positions.push_back(vector_type({ 0.0, -3.0,  4.0}));
    positions.push_back(vector_type({-5.2, -8.0,  5.0}));
    positions.push_back(vector_type({ 0.0,  0.0, 20.0}));
    return positions;
}

template <typename float_type>
std::shared_ptr<mdsim::host::potentials::external::planar_wall<3, float_type>>
make_host_potential()
{
    typedef typename mdsim::host::potentials::external::planar_wall<3, float_type> potential_type;
    typedef typename potential_type::vector_type vector_type;

    // define interaction parameters
    unsigned int ntype = 2;  // test a binary mixture
    unsigned int nwall = 3;  // 3 walls

    typename potential_type::matrix_container_type cutoff(nwall, ntype);
    cutoff <<=
        1., 2.
      , 10., 10.
      , 1e2, 1e2;

    typename potential_type::matrix_container_type epsilon(nwall, ntype);
    epsilon <<=
        1., 0.5
      , 2., 0.
      , 0., 1.;

    typename potential_type::matrix_container_type sigma(nwall, ntype);
    sigma <<=
        1., 2.
      , 3., 1.
      , 2., 1.;

    typename potential_type::matrix_container_type wetting(nwall, ntype);
    wetting <<=
        0.4, 0.
      , -1, 1.
      , 0., 2.;

    typename potential_type::vector_container_type surface_normal(nwall);
    surface_normal <<=
        vector_type({1,  0, 0})
      , vector_type({1, -1, 0})
      , vector_type({-1,  0, 2});

    typename potential_type::scalar_container_type offset(nwall);
    offset <<=
        -5, 5, 10;

    float_type smoothing = 0.01;

    // construct module
    return std::make_shared<potential_type>(offset, surface_normal, epsilon, sigma, wetting, cutoff, smoothing);
}

BOOST_AUTO_TEST_CASE( planar_wall_host )
{
#ifndef USE_HOST_SINGLE_PRECISION
    typedef double float_type;
#else
    typedef float float_type;
#endif

    // construct host module with fixed parameters, test in three dimensions only
    auto potential = make_host_potential<float_type>();
    unsigned int nwall = potential->epsilon().size1();
    unsigned int nspecies = potential->epsilon().size2();

    // positions of test particles
    typedef decltype(potential)::element_type::vector_type vector_type;
    auto positions = make_positions<vector_type>();

    for (unsigned int species = 0; species < nspecies; ++species) {
        for (vector_type const& r : positions) {
            BOOST_TEST_MESSAGE("test particle of species " << species << " at " << r);

            // compute force and potential energy using the planar_wall module
            vector_type force;
            float_type en_pot;
            std::tie(force, en_pot) = (*potential)(r, species);  // particle species A

            // compute reference values: sum up contributions from each wall
            vector_type force2 = 0;
            float_type en_pot2 = 0;
            for (unsigned int i = 0; i < nwall; ++i) {
                vector_type f;
                float_type en;
                std::tie(f, en) = compute_wall(
                    r, potential->surface_normal()(i), potential->offset()(i)
                  , potential->epsilon()(i, species), potential->sigma()(i, species), potential->wetting()(i, species)
                  , potential->cutoff()(i, species), potential->smoothing()
                );
                force2 += f;
                en_pot2 += en;
            }

            // compare output from planar_wall module with reference
            const float_type tolerance = 5 * std::numeric_limits<float_type>::epsilon();
            for (unsigned int i = 0; i < force.size(); ++i) {
                BOOST_CHECK_CLOSE_FRACTION(force[i], force2[i], tolerance);
            }
            BOOST_CHECK_CLOSE_FRACTION(en_pot, en_pot2, tolerance);
        }
    }
}

#ifdef HALMD_WITH_GPU

template <typename float_type>
struct planar_wall
{
    enum { dimension = 3 };

    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::potentials::external::planar_wall<dimension, float> potential_type;
#ifndef USE_HOST_SINGLE_PRECISION
    typedef mdsim::host::potentials::external::planar_wall<dimension, double> host_potential_type;
#else
    typedef mdsim::host::potentials::external::planar_wall<dimension, float> host_potential_type;
#endif
    typedef mdsim::gpu::forces::external<dimension, float_type, potential_type> force_type;

    typedef typename particle_type::vector_type vector_type;

    std::shared_ptr<box_type> box;
    std::shared_ptr<potential_type> potential;
    std::shared_ptr<force_type> force;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<host_potential_type> host_potential;

    std::vector<vector_type> positions;

    planar_wall();
    void test();
};

template <typename float_type>
void planar_wall<float_type>::test()
{
    // alternating species of test particles
    std::vector<unsigned int> species(particle->nparticle());
    for (unsigned int i = 0; i < species.size(); ++i) {
        species[i] = i % particle->nspecies();
    }

    // assign positions and species
    BOOST_CHECK( set_position(*particle, positions.begin()) == positions.end() );
    BOOST_CHECK( set_species(*particle, species.begin()) == species.end() );

    // compute and read forces and other stuff from device
    std::vector<float> en_pot(particle->nparticle());
    BOOST_CHECK( get_potential_energy(*particle, en_pot.begin()) == en_pot.end() );

    std::vector<vector_type> f_list(particle->nparticle());
    BOOST_CHECK( get_force(*particle, f_list.begin()) == f_list.end() );

    // on the GPU, the potential is always evaluated in single precision
    const float_type tolerance = 10 * std::numeric_limits<float>::epsilon();

    for (unsigned int i = 0; i < positions.size(); ++i) {
        // reference values from host module
        typedef typename host_potential_type::vector_type host_vector_type;
        host_vector_type r = static_cast<host_vector_type>(positions[i]);
        BOOST_TEST_MESSAGE("test particle #" << i+1 << " of species " << species[i] << " at " << r);

        host_vector_type force2;
        float_type en_pot2;
        std::tie(force2, en_pot2) = (*host_potential)(r, species[i]);

        // compare output from host and GPU implementations
        vector_type force = f_list[i];
        for (unsigned int j = 0; j < dimension; ++j) {
            BOOST_CHECK_CLOSE_FRACTION(force[j], force2[j], tolerance);
        }
    }
}

template <typename float_type>
planar_wall<float_type>::planar_wall()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    positions = make_positions<vector_type>();

    // set module parameters
    float box_length = 100;
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = box_length;
    }

    // create host module for reference
#ifndef USE_HOST_SINGLE_PRECISION
    host_potential = make_host_potential<double>();
#else
    host_potential = make_host_potential<float>();
#endif

    // create modules
    particle = std::make_shared<particle_type>(positions.size(), host_potential->epsilon().size2()); // 2nd argument: #species
    box = std::make_shared<box_type>(edges);

    // fixed_vector<double, N> is not implicitly convertible to fixed_vector<float, N>
    //
    // Therefore, we have to convert the surface normals by hand.
    unsigned int nwall = host_potential->surface_normal().size();
    typename potential_type::vector_container_type surface_normal(nwall);
    for (unsigned int i = 0; i < nwall; ++i) {
        surface_normal[i] = static_cast<typename potential_type::vector_type>(host_potential->surface_normal()[i]);
    }

    potential = std::make_shared<potential_type>(
            host_potential->offset()
          , surface_normal
          , host_potential->epsilon()
          , host_potential->sigma()
          , host_potential->wetting()
          , host_potential->cutoff()
          , host_potential->smoothing()
    );

    force = std::make_shared<force_type>(potential, particle, box);

    particle->on_prepend_force([=](){force->check_cache();});
    particle->on_force([=](){force->apply();});
}

# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE( planar_wall_gpu_dsfloat, device ) {
    planar_wall<dsfloat>().test();
}
#endif
# ifdef USE_GPU_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE( planar_wall_gpu_float, device ) {
    planar_wall<float>().test();
}
#endif

#endif // HALMD_WITH_GPU
