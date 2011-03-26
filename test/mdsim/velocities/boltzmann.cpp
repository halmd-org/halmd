/*
 * Copyright © 2011  Felix Höfling and Peter Colberg
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

#define BOOST_TEST_MODULE velocity_boltzmann
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <boost/assign.hpp>
#include <limits>
#include <string>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/utility/init.hpp>
#include <test/modules.hpp>

using namespace boost;
using namespace halmd;
using namespace halmd::test;
using namespace std;

const double eps = numeric_limits<double>::epsilon();
const float eps_float = numeric_limits<float>::epsilon();

/**
 * test initialisation of particle velocities: boltzmann module
 */

template <int dimension>
void boltzmann(string const& backend)
{
    typedef typename mdsim::type_traits<dimension, double>::vector_type vector_type;
    typedef typename mdsim::type_traits<dimension, float>::vector_type gpu_vector_type;

    unsigned npart = (backend == "gpu") ? 10000 : 3000;
    double temp = 2.0;
    double density = 0.3;
    char const* random_file = "/dev/urandom";

    // enable logging to console
    shared_ptr<logger> log(new logger);
    log->log_to_console(
#ifdef NDEBUG
        logger::warning
#else
        logger::debug
#endif
    );

    // init modules
    BOOST_TEST_MESSAGE("initialise modules");
#ifdef WITH_CUDA
    shared_ptr<utility::gpu::device> device = make_device(backend);
#endif /* WITH_CUDA */

    shared_ptr<halmd::random::random> random = make_random(backend, random_file);

    shared_ptr<mdsim::particle<dimension> > particle =
        make_particle<dimension>(backend, npart);

    shared_ptr<mdsim::velocity<dimension> > velocity =
        make_boltzmann(backend, particle, random, 0);        // actual temperature is set below

    shared_ptr<mdsim::box<dimension> > box = make_box(particle, density);

    // measure velocity distribution via thermodynamics module
    shared_ptr<observables::thermodynamics<dimension> > thermodynamics =
        make_thermodynamics(
            backend, particle, box
          , make_zero_force<dimension>(backend, particle)
        );

    // destroy and reconstruct module
    velocity = make_boltzmann(backend, particle, random, temp);

    // generate velocity distribution
    LOG_DEBUG("set particle tags");
    particle->set();
    BOOST_TEST_MESSAGE("generate Maxwell-Boltzmann distribution");
    velocity->set();

    //
    // test velocity distribution of final state
    //
    // centre-of-mass velocity ⇒ mean of velocity distribution
    // each particle is an independent "measurement"
    // limit is 3σ, σ = √(<v_x²> / (N - 1)) where <v_x²> = k T
    double vcm_limit = 3 * sqrt(temp / (npart - 1));
    BOOST_TEST_MESSAGE("Absolute limit on instantaneous centre-of-mass velocity: " << vcm_limit);
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);  //< norm_inf tests the max. value

    // temperature ⇒ variance of velocity distribution
    // we have only one measurement of the variance
    // limit is 3σ, σ = √<ΔT²> where <ΔT²> / T² = 2 / (dimension × N)
    double rel_temp_limit = 3 * sqrt(2. / (dimension * npart)) / temp;
    BOOST_TEST_MESSAGE("Relative limit on instantaneous temperature: " << rel_temp_limit);
    BOOST_CHECK_CLOSE_FRACTION(thermodynamics->temp(), temp, rel_temp_limit);

    //
    // test shifting and rescaling
    //
    // multiplication of the velocities by a constant factor
    double scale = 1.5;
    velocity->rescale(scale);
    BOOST_CHECK_CLOSE_FRACTION(thermodynamics->temp(), scale * scale * temp, rel_temp_limit);

    // shift mean velocity to zero
    vector_type v_cm = thermodynamics->v_cm();
    velocity->shift(-v_cm);
    vcm_limit = (backend == "gpu") ? 0.1 * eps_float : 2 * eps;
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    // first shift, then rescale in one step
    velocity->shift_rescale(v_cm, 1 / scale);
    BOOST_CHECK_CLOSE_FRACTION(thermodynamics->temp(), temp, rel_temp_limit);
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm() - v_cm), vcm_limit);
}

HALMD_INIT( init_unit_test_suite )
{
    using namespace boost::assign;
    using namespace boost::unit_test;
    using namespace boost::unit_test::framework;

    // parametrize specific program options
    vector<string> backend = list_of
        ("host")
#ifdef WITH_CUDA
        ("gpu")
#endif /* WITH_CUDA */
        ;

    test_suite* ts1 = BOOST_TEST_SUITE( "boltzmann" );

    test_suite* ts11 = BOOST_TEST_SUITE( "host" );

    test_suite* ts111 = BOOST_TEST_SUITE( "2d" );
    ts111->add( BOOST_PARAM_TEST_CASE( &boltzmann<2>, backend.begin(), backend.begin() + 1 ) );

    test_suite* ts112 = BOOST_TEST_SUITE( "3d" );
    ts112->add( BOOST_PARAM_TEST_CASE( &boltzmann<3>, backend.begin(), backend.begin() + 1 ) );

    ts11->add( ts111 );
    ts11->add( ts112 );
    ts1->add( ts11 );

#ifdef WITH_CUDA
    test_suite* ts12 = BOOST_TEST_SUITE( "gpu" );

    test_suite* ts121 = BOOST_TEST_SUITE( "2d" );
    ts121->add( BOOST_PARAM_TEST_CASE( &boltzmann<2>, backend.begin() + 1, backend.end() ) );

    test_suite* ts122 = BOOST_TEST_SUITE( "3d" );
    ts122->add( BOOST_PARAM_TEST_CASE( &boltzmann<3>, backend.begin() + 1, backend.end() ) );

    ts12->add( ts121 );
    ts12->add( ts122 );
    ts1->add( ts12 );
#endif

    master_test_suite().add( ts1 );
}
