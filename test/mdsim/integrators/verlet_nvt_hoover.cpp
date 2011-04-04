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

#define BOOST_TEST_MODULE integrator_verlet_nvt_hoover
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

#include <iomanip>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <test/modules.hpp>
#include <test/tools/init.hpp>

using namespace boost;
using namespace boost::assign;
using namespace halmd;
using namespace halmd::test;
using namespace std;

/**
 * test NVT Verlet integrator with Nosé-Hoover chain thermostat
 */

const double eps = numeric_limits<double>::epsilon();
const float eps_float = numeric_limits<float>::epsilon();

/**
 * heat capacity from canonical fluctuations (variance) of potential and kinetic energy
 */
inline double heat_capacity_nvt(double var_en_pot, double var_en_kin, double temperature, unsigned npart)
{
    return npart * (var_en_pot + var_en_kin) / (temperature * temperature);
}

template <int dimension>
void verlet_nvt_hoover(string const& backend)
{
    typedef typename mdsim::type_traits<dimension, double>::vector_type vector_type;
    typedef typename mdsim::type_traits<dimension, float>::vector_type gpu_vector_type;

    float temp = 1.;
    float start_temp = 3.;
    float density = 0.1;
    unsigned npart = (backend == "gpu") ? 1500 : 1500;
    double timestep = 0.002;
    double resonance_frequency = 5.;
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

    core->integrator = make_verlet_nvt_hoover_integrator<dimension, double>(
        backend, core->particle, core->box, timestep, temp, resonance_frequency
    );

    core->force = make_lennard_jones_force<dimension>(
        backend, core->particle, core->box
      , list_of(pow(2.f, 1.f/6))(0.f)(0.f)     /* cutoff, WCA potential */
      , list_of(1.f)(0.f)(0.f)                 /* epsilon */
      , list_of(1.f)(0.f)(0.f)                 /* sigma */
    );

    core->neighbour = make_neighbour(backend, core->particle, core->box, core->force);

    core->position = make_lattice(backend, core->particle, core->box, random);

    core->velocity = make_boltzmann(backend, core->particle, random, start_temp);

    // use thermodynamics module to measure temperature (velocity distribution)
    shared_ptr<observables::thermodynamics<dimension> > thermodynamics =
        make_thermodynamics(backend, core->particle, core->box, core->force);

    // run for Δt*=500
    uint64_t steps = static_cast<uint64_t>(ceil(500 / timestep));
    // skip Δt*=50 for equilibration
    uint64_t skip = static_cast<uint64_t>(ceil(50 / timestep));
    // ensure that sampling period is sufficiently large such that
    // the samples can be considered independent
    uint64_t period = static_cast<uint64_t>(round(3 / (resonance_frequency * timestep)));
    accumulator<double> temp_;
    array<accumulator<double>, dimension> v_cm;   //< accumulate velocity component-wise
    double max_en_diff = 0;                       // integral of motion: Hamiltonian extended by NHC terms

    BOOST_TEST_MESSAGE("prepare system");
    core->force->aux_enable();                    //< enable computation of potential energy
    core->prepare();

    // equilibrate the system,
    // this avoids a jump in the conserved energy at the very beginning
    BOOST_TEST_MESSAGE("equilibrate over " << steps / 20 << " steps");
    for (uint64_t i = 0; i < steps / 20; ++i) {
        core->mdstep();
    }

    // compute modified Hamiltonian
    double en_nhc0 = thermodynamics->en_tot();
    if (backend == "gpu") {
#ifdef WITH_CUDA
        typedef mdsim::gpu::integrators::verlet_nvt_hoover<dimension, double> integrator_type;
        shared_ptr<integrator_type> integrator = dynamic_pointer_cast<integrator_type>(core->integrator);
        en_nhc0 += integrator->en_nhc();
#endif
    }
    else if (backend == "host") {
        typedef mdsim::host::integrators::verlet_nvt_hoover<dimension, double> integrator_type;
        shared_ptr<integrator_type> integrator = dynamic_pointer_cast<integrator_type>(core->integrator);
        en_nhc0 += integrator->en_nhc();
    }

    BOOST_TEST_MESSAGE("run NVT integrator over " << steps << " steps");
    core->force->aux_disable();
    for (uint64_t i = 0; i < steps; ++i) {
        // enable auxiliary variables in force module
        if(i % period == 0) {
            core->force->aux_enable();
        }

        // perform MD step
        core->mdstep();

        // measurement
        if(i % period == 0) {
            // measure temperature after thermalisation
            if (i >= skip) {
                temp_(thermodynamics->temp());
            }
            // track centre-of-mass velocity over the whole run
            fixed_vector<double, dimension> v(thermodynamics->v_cm());
            for (unsigned int k = 0; k < dimension; ++k) {
                v_cm[k](v[k]);
            }

            // compute modified Hamiltonian
            double en_nhc_;
            fixed_vector<double, 2> xi(0), v_xi(0);
            en_nhc_ = thermodynamics->en_tot();
            if (backend == "gpu") {
#ifdef WITH_CUDA
                typedef mdsim::gpu::integrators::verlet_nvt_hoover<dimension, double> integrator_type;
                shared_ptr<integrator_type> integrator = dynamic_pointer_cast<integrator_type>(core->integrator);
                xi = integrator->xi;
                v_xi = integrator->v_xi;
                en_nhc_ += integrator->en_nhc();
#endif
            }
            else if (backend == "host") {
                typedef mdsim::host::integrators::verlet_nvt_hoover<dimension, double> integrator_type;
                shared_ptr<integrator_type> integrator = dynamic_pointer_cast<integrator_type>(core->integrator);
                xi = integrator->xi;
                v_xi = integrator->v_xi;
                en_nhc_ += integrator->en_nhc();
            }
            LOG_TRACE(setprecision(12)
                << "en_nhc: " << i * timestep
                << " " << en_nhc_ << " " << thermodynamics->temp()
                << " " << xi[0] << " " << xi[1] << " " << v_xi[0] << " " << v_xi[1]
                << setprecision(6)
            );
            max_en_diff = max(abs(en_nhc_ - en_nhc0), max_en_diff);
            core->force->aux_disable();
        }
    }

    //
    // test conservation of pseudo-Hamiltonian
    //
    const double en_limit = max(2e-5, steps * 1e-12);
    BOOST_CHECK_SMALL(max_en_diff / fabs(en_nhc0), en_limit);

    //
    // test conservation of total momentum
    //
    double vcm_limit = (backend == "gpu") ? 0.1 * eps_float : 20 * eps;
    BOOST_TEST_MESSAGE("Absolute limit on centre-of-mass velocity: " << vcm_limit);
    for (unsigned int i = 0; i < dimension; ++i) {
        BOOST_CHECK_SMALL(mean(v_cm[i]), vcm_limit);
        BOOST_CHECK_SMALL(error_of_mean(v_cm[i]), vcm_limit);
    }

    //
    // test final distribution of velocities
    //
    // mean (total momentum) should be zero
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

    // kinetic part of specific heat per particle ⇒ temperature fluctuations
    // c_V = k × (dimension × N / 2)² <ΔT²> / T² / N = k × dimension / 2
    // where we have used <ΔT²> / T² = 2 / (dimension × N)
    // limit is 3σ, with the approximation
    // σ² = Var[ΔE² / (k T²)] / C → (dimension / 2) × (dimension + 6 / N) / C
    // (one measurement only from the average over C samples)
    double cv = pow(.5 * dimension, 2.) * npart * variance(temp_);
    double cv_variance=  (.5 * dimension) * (dimension + 6. / npart) / count(temp_);
    double rel_cv_limit = 3 * sqrt(cv_variance) / (.5 * dimension);
    BOOST_TEST_MESSAGE("Relative limit on kinetic part of specific heat: " << rel_cv_limit);
    BOOST_CHECK_CLOSE_FRACTION(cv, .5 * dimension, rel_cv_limit);
}

HALMD_TEST_INIT( init_unit_test_suite )
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

    test_suite* ts1 = BOOST_TEST_SUITE( "verlet_nvt_hoover" );

    test_suite* ts11 = BOOST_TEST_SUITE( "host" );

    test_suite* ts111 = BOOST_TEST_SUITE( "2d" );
    ts111->add( BOOST_PARAM_TEST_CASE( &verlet_nvt_hoover<2>, backend.begin(), backend.begin() + 1 ) );

    test_suite* ts112 = BOOST_TEST_SUITE( "3d" );
    ts112->add( BOOST_PARAM_TEST_CASE( &verlet_nvt_hoover<3>, backend.begin(), backend.begin() + 1 ) );

    ts11->add( ts111 );
    ts11->add( ts112 );
    ts1->add( ts11 );

#ifdef WITH_CUDA
    test_suite* ts12 = BOOST_TEST_SUITE( "gpu" );

    test_suite* ts121 = BOOST_TEST_SUITE( "2d" );
    ts121->add( BOOST_PARAM_TEST_CASE( &verlet_nvt_hoover<2>, backend.begin() + 1, backend.end() ) );

    test_suite* ts122 = BOOST_TEST_SUITE( "3d" );
    ts122->add( BOOST_PARAM_TEST_CASE( &verlet_nvt_hoover<3>, backend.begin() + 1, backend.end() ) );

    ts12->add( ts121 );
    ts12->add( ts122 );
    ts1->add( ts12 );
#endif

    master_test_suite().add( ts1 );
}
