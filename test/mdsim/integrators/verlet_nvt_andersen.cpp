/*
 * Copyright © 2010  Felix Höfling and Peter Colberg
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

#define BOOST_TEST_MODULE integrator_verlet_nvt_andersen
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <utility>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <test/modules.hpp>

using namespace boost;
using namespace boost::assign;
using namespace halmd;
using namespace halmd::test;
using namespace std;

/**
 * test NVT verlet integrator with stochastic Andersen thermostat
 */

const double eps = numeric_limits<double>::epsilon();
const float eps_float = numeric_limits<float>::epsilon();

template <int dimension>
void verlet_nvt_andersen(string const& backend)
{
    typedef typename mdsim::type_traits<dimension, double>::vector_type vector_type;
    typedef typename mdsim::type_traits<dimension, float>::vector_type gpu_vector_type;

    float temp = 1.;
    float density = 0.3;
    unsigned npart = (backend == "gpu") ? 5000 : 1500;
    double timestep = 0.01;
    double coll_rate = 10;
    char const* random_file = "/dev/urandom";
    fixed_vector<double, dimension> box_ratios =
        (dimension == 3) ? list_of(1.)(2.)(1.01) : list_of(1.)(2.);

    // enable logging to console
    shared_ptr<logger> log(new logger);
    log->log_to_console(
#ifdef NDEBUG
        logger::warning
#else
        logger::debug
#endif
    );

    BOOST_TEST_MESSAGE("using backend '" << backend << "' in " <<
                       dimension << " dimensions");

#ifdef WITH_CUDA
    shared_ptr<utility::gpu::device> device = make_device(backend);
#endif /* WITH_CUDA */
    shared_ptr<halmd::random::random> random = make_random(backend, random_file);

    BOOST_TEST_MESSAGE("initialise simulation modules");
    // init core module and all dependencies
    shared_ptr<mdsim::core<dimension> > core(new mdsim::core<dimension>);

    core->particle = make_particle<dimension>(backend, npart);

    core->box = make_box<dimension>(core->particle, density, box_ratios);

    core->integrator = make_verlet_nvt_andersen_integrator<dimension>(
        backend, core->particle, core->box, random, timestep, temp, coll_rate
    );

    core->force = make_zero_force<dimension>(backend, core->particle);

    core->position = make_lattice(backend, core->particle, core->box, random);

    core->velocity = make_boltzmann(backend, core->particle, random, temp);

    // use thermodynamics module to measure temperature (velocity distribution)
    shared_ptr<observables::thermodynamics<dimension> > thermodynamics =
        make_thermodynamics(backend, core->particle, core->box, core->force);

    // run for Δt*=500
    uint64_t steps = static_cast<uint64_t>(ceil(500 / timestep));
    // ensure that sampling period is sufficiently large such that
    // the samples can be considered independent
    uint64_t period = static_cast<uint64_t>(round(3. / (coll_rate * timestep)));
    accumulator<double> temp_;
    array<accumulator<double>, dimension> v_cm;   //< accumulate velocity component-wise

    core->force->aux_disable();                   //< we don't need potential energy, pressure, etc.
    core->prepare();
    BOOST_TEST_MESSAGE("run NVT integrator over " << steps << " steps");
    for (uint64_t i = 0; i < steps; ++i) {
        core->mdstep();
        if(i % period == 0) {
            temp_(thermodynamics->temp());
            fixed_vector<double, dimension> v(thermodynamics->v_cm());
            for (unsigned int i = 0; i < dimension; ++i) {
                v_cm[i](v[i]);
            }
        }
    }

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
    // test velocity distribution averaged over the whole simulation run
    //
    // centre-of-mass velocity ⇒ mean of velocity distribution
    // #measurements = #particles × #samples
    // limit is 3σ, σ = √(<v_x²> / (N × C - 1)) where <v_x²> = k T
    vcm_limit = 3 * sqrt(temp / (npart * count(v_cm[0]) - 1));
    BOOST_TEST_MESSAGE("Absolute limit on centre-of-mass velocity: " << vcm_limit);
    for (unsigned int i = 0; i < dimension; ++i) {
        BOOST_CHECK_SMALL(mean(v_cm[i]), 3 * vcm_limit);
        BOOST_CHECK_SMALL(error_of_mean(v_cm[i]), vcm_limit);
    }

    // mean temperature ⇒ variance of velocity distribution
    // each sample should constitute an independent measurement
    // limit is 3σ, σ = √(<ΔT²> / (C - 1)) where <ΔT²> / T² = 2 / (dimension × N)
    rel_temp_limit = 3 * sqrt(2. / (dimension * npart * (count(temp_) - 1))) / temp;
    BOOST_TEST_MESSAGE("Relative limit on temperature: " << rel_temp_limit);
    BOOST_CHECK_CLOSE_FRACTION(mean(temp_), temp, rel_temp_limit);

    // specific heat per particle ⇒ temperature fluctuations
    // c_V = k × (dimension × N / 2)² <ΔT²> / T² / N = k × dimension / 2
    // where we have used <ΔT²> / T² = 2 / (dimension × N)
    // limit is 3σ, with the approximation
    // σ² = Var[ΔE² / (k T²)] / C → (dimension / 2) × (dimension + 6 / N) / C
    // (one measurement only from the average over C samples)
    double cv = pow(.5 * dimension, 2.) * npart * variance(temp_);
    double cv_variance=  (.5 * dimension) * (dimension + 6. / npart) / count(temp_);
    double rel_cv_limit = 3 * sqrt(cv_variance) / (.5 * dimension);
    BOOST_TEST_MESSAGE("Relative limit on specific heat: " << rel_cv_limit);
    BOOST_CHECK_CLOSE_FRACTION(cv, .5 * dimension, rel_cv_limit);
}

static void __attribute__((constructor)) init_unit_test_suite()
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

    test_suite* ts1 = BOOST_TEST_SUITE( "verlet_nvt_andersen" );

    test_suite* ts11 = BOOST_TEST_SUITE( "host" );

    test_suite* ts111 = BOOST_TEST_SUITE( "2d" );
    ts111->add( BOOST_PARAM_TEST_CASE( &verlet_nvt_andersen<2>, backend.begin(), backend.begin() + 1 ) );

    test_suite* ts112 = BOOST_TEST_SUITE( "3d" );
    ts112->add( BOOST_PARAM_TEST_CASE( &verlet_nvt_andersen<3>, backend.begin(), backend.begin() + 1 ) );

    ts11->add( ts111 );
    ts11->add( ts112 );
    ts1->add( ts11 );

#ifdef WITH_CUDA
    test_suite* ts12 = BOOST_TEST_SUITE( "gpu" );

    test_suite* ts121 = BOOST_TEST_SUITE( "2d" );
    ts121->add( BOOST_PARAM_TEST_CASE( &verlet_nvt_andersen<2>, backend.begin() + 1, backend.end() ) );

    test_suite* ts122 = BOOST_TEST_SUITE( "3d" );
    ts122->add( BOOST_PARAM_TEST_CASE( &verlet_nvt_andersen<3>, backend.begin() + 1, backend.end() ) );

    ts12->add( ts121 );
    ts12->add( ts122 );
    ts1->add( ts12 );
#endif

    master_test_suite().add( ts1 );
}
