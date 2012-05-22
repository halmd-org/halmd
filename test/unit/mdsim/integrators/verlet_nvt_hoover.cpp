/*
 * Copyright © 2010-2012  Felix Höfling
 * Copyright © 2010-2011  Peter Colberg
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
#include <boost/make_shared.hpp>
#include <limits>
#include <iomanip>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/core.hpp>
#include <halmd/mdsim/host/forces/pair_trunc.hpp>
#include <halmd/mdsim/host/integrators/verlet_nvt_hoover.hpp>
#include <halmd/mdsim/host/maximum_squared_displacement.hpp>
#include <halmd/mdsim/host/neighbour.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/mdsim/host/potentials/lennard_jones.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/predicates/greater.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
# include <halmd/mdsim/gpu/integrators/verlet_nvt_hoover.hpp>
# include <halmd/mdsim/gpu/maximum_squared_displacement.hpp>
# include <halmd/mdsim/gpu/neighbours/from_binning.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/positions/lattice.hpp>
# include <halmd/mdsim/gpu/potentials/lennard_jones.hpp>
# include <halmd/mdsim/gpu/velocities/boltzmann.hpp>
# include <halmd/observables/gpu/thermodynamics.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif
#include <test/tools/ctest.hpp>

using namespace boost;
using namespace boost::assign; // list_of
using namespace halmd;
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

/** test Verlet integrator: 'ideal' gas without interactions (setting ε=0) */

template <typename modules_type>
struct verlet_nvt_hoover
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::potential_type potential_type;
    typedef typename modules_type::force_type force_type;
    typedef typename modules_type::binning_type binning_type;
    typedef typename modules_type::neighbour_type neighbour_type;
    typedef typename modules_type::max_displacement_type max_displacement_type;
    typedef typename modules_type::integrator_type integrator_type;
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::position_type position_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::thermodynamics_type thermodynamics_type;
    typedef typename modules_type::velocity_type velocity_type;
    static bool const gpu = modules_type::gpu;

    typedef mdsim::clock clock_type;
    typedef typename clock_type::time_type time_type;
    typedef typename clock_type::step_type step_type;
    typedef mdsim::core core_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    typedef predicates::greater<float_type> greater_type;
    static unsigned int const dimension = vector_type::static_size;

    float temp;
    float start_temp;
    float density;
    unsigned int npart;
    double timestep;
    double resonance_frequency;
    fixed_vector<double, dimension> box_ratios;
    float skin;

    shared_ptr<box_type> box;
    shared_ptr<clock_type> clock;
    shared_ptr<core_type> core;
    shared_ptr<potential_type> potential;
    shared_ptr<force_type> force;
    shared_ptr<binning_type> binning;
    shared_ptr<neighbour_type> neighbour;
    shared_ptr<max_displacement_type> max_displacement;
    shared_ptr<integrator_type> integrator;
    shared_ptr<particle_type> particle;
    shared_ptr<position_type> position;
    shared_ptr<random_type> random;
    shared_ptr<thermodynamics_type> thermodynamics;
    shared_ptr<velocity_type> velocity;

    void test();
    verlet_nvt_hoover();
    void connect();
};

template <typename modules_type>
void verlet_nvt_hoover<modules_type>::test()
{
    typedef typename mdsim::type_traits<dimension, double>::vector_type vector_type;
    typedef typename mdsim::type_traits<dimension, float>::vector_type gpu_vector_type;

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
    core->setup();

    // equilibrate the system,
    // this avoids a jump in the conserved energy at the very beginning
    BOOST_TEST_MESSAGE("equilibrate over " << steps / 20 << " steps");
    for (uint64_t i = 0; i < steps / 20; ++i) {
        if (i + 1 == steps / 20) {
            force->aux_enable();                    //< enable computation of potential energy
        }
        core->mdstep();
    }

    // compute modified Hamiltonian
    double en_nhc0 = thermodynamics->en_tot() + integrator->en_nhc();

    BOOST_TEST_MESSAGE("run NVT integrator over " << steps << " steps");
    for (uint64_t i = 0; i < steps; ++i) {
        // enable auxiliary variables in force module
        if(i % period == 0) {
            force->aux_enable();
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
            double en_nhc_ = thermodynamics->en_tot() + integrator->en_nhc();
            LOG_TRACE(setprecision(12)
                << "en_nhc: " << i * timestep
                << " " << en_nhc_ << " " << thermodynamics->temp()
                << " " << integrator->xi << " " << integrator->v_xi
                << setprecision(6)
            );
            max_en_diff = max(abs(en_nhc_ - en_nhc0), max_en_diff);
        }
    }

    //
    // test conservation of pseudo-Hamiltonian
    //
    const double en_tolerance = max(3e-5, steps * 1e-12);
    BOOST_CHECK_SMALL(max_en_diff / fabs(en_nhc0), en_tolerance);

    //
    // test conservation of total momentum
    //
    double vcm_tolerance = gpu ? 0.1 * eps_float : 20 * eps;
    BOOST_TEST_MESSAGE("Absolute tolerance on centre-of-mass velocity: " << vcm_tolerance);
    for (unsigned int i = 0; i < dimension; ++i) {
        BOOST_CHECK_SMALL(mean(v_cm[i]), vcm_tolerance);
        BOOST_CHECK_SMALL(error_of_mean(v_cm[i]), vcm_tolerance);
    }

    //
    // test final distribution of velocities
    //
    // mean (total momentum) should be zero
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_tolerance);  //< norm_inf tests the max. value

    // temperature ⇒ variance of velocity distribution
    // we have only one measurement of the variance
    // tolerance is 4.5σ, σ = √<ΔT²> where <ΔT²> / T² = 2 / (dimension × N),
    // with this choice, a single test passes with 99.999% probability
    double rel_temp_tolerance = 4.5 * sqrt(2. / (dimension * npart)) / temp;
    BOOST_TEST_MESSAGE("Relative tolerance on instantaneous temperature: " << rel_temp_tolerance);
    BOOST_CHECK_CLOSE_FRACTION(thermodynamics->temp(), temp, rel_temp_tolerance);

    //
    // test velocity distribution averaged over the whole simulation run
    //
    // centre-of-mass velocity ⇒ mean of velocity distribution
    // #measurements = #particles × #samples,
    // tolerance is 4.5σ, σ = √(<v_x²> / (N × C - 1)) where <v_x²> = k T
    vcm_tolerance = 4.5 * sqrt(temp / (npart * count(v_cm[0]) - 1));
    BOOST_TEST_MESSAGE("Absolute tolerance on centre-of-mass velocity: " << vcm_tolerance);
    for (unsigned int i = 0; i < dimension; ++i) {
        BOOST_CHECK_SMALL(mean(v_cm[i]), vcm_tolerance);
        BOOST_CHECK_SMALL(error_of_mean(v_cm[i]), vcm_tolerance);
    }

    // mean temperature ⇒ variance of velocity distribution
    // each sample should constitute an independent measurement
    // tolerance is 4.5σ, σ = √(<ΔT²> / (C - 1)) where <ΔT²> / T² = 2 / (dimension × N)
    rel_temp_tolerance = 4.5 * sqrt(2. / (dimension * npart * (count(temp_) - 1))) / temp;
    BOOST_TEST_MESSAGE("Relative tolerance on temperature: " << rel_temp_tolerance);
    BOOST_CHECK_CLOSE_FRACTION(mean(temp_), temp, rel_temp_tolerance);

    // kinetic part of specific heat per particle ⇒ temperature fluctuations
    // c_V = k × (dimension × N / 2)² <ΔT²> / T² / N = k × dimension / 2
    // where we have used <ΔT²> / T² = 2 / (dimension × N)
    // tolerance is 4.5σ, with the approximation
    // σ² = Var[ΔE² / (k T²)] / C → (dimension / 2) × (dimension + 6 / N) / C
    // (one measurement only from the average over C samples)
    double cv = pow(.5 * dimension, 2.) * npart * variance(temp_);
    double cv_variance=  (.5 * dimension) * (dimension + 6. / npart) / count(temp_);
    double rel_cv_tolerance = 4.5 * sqrt(cv_variance) / (.5 * dimension);
    BOOST_TEST_MESSAGE("Relative tolerance on kinetic part of specific heat: " << rel_cv_tolerance);
    BOOST_CHECK_CLOSE_FRACTION(cv, .5 * dimension, rel_cv_tolerance);
}

template <typename modules_type>
verlet_nvt_hoover<modules_type>::verlet_nvt_hoover()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    temp = 1.;
    start_temp = 3.;
    density = 0.1;
    npart = 1500;
    timestep = 0.002;
    resonance_frequency = 5.;
    box_ratios = (dimension == 3) ? list_of(1.)(2.)(1.01) : list_of(1.)(2.);
    skin = 0.5;

    vector<unsigned int> npart_vector = list_of(npart);

    // create modules
    particle = make_shared<particle_type>(npart_vector);
    box = make_shared<box_type>(npart, density, box_ratios);
    random = make_shared<random_type>();
    integrator = make_shared<integrator_type>(particle, box, timestep, temp, resonance_frequency);
    potential = make_shared<potential_type>(
        particle->ntype
      , list_of(pow(2.f, 1.f/6))(0.f)(0.f)     /* cutoff, WCA potential */
      , list_of(1.f)(0.f)(0.f)                 /* epsilon */
      , list_of(1.f)(0.f)(0.f)                 /* sigma */
    );
    binning = make_shared<binning_type>(particle, box, potential->r_cut(), skin);
    neighbour = make_shared<neighbour_type>(particle, box, binning, potential->r_cut(), skin);
    force = make_shared<force_type>(potential, particle, box, neighbour);
    position = make_shared<position_type>(particle, box, random, 1);
    velocity = make_shared<velocity_type>(particle, random, start_temp);
    clock = make_shared<clock_type>(timestep);
    thermodynamics = make_shared<thermodynamics_type>(particle, box, clock, force);
    max_displacement = make_shared<max_displacement_type>(particle, box);

    // create core and connect module slots to core signals
    this->connect();
}

template <typename modules_type>
void verlet_nvt_hoover<modules_type>::connect()
{
    core = make_shared<core_type>(clock);
    // system preparation
    core->on_prepend_setup( bind(&particle_type::set, particle) );
    core->on_setup( bind(&position_type::set, position) );
    core->on_setup( bind(&velocity_type::set, velocity) );
    core->on_append_setup( bind(&max_displacement_type::zero, max_displacement) );
    core->on_append_setup( bind(&binning_type::update, binning) );
    core->on_append_setup( bind(&neighbour_type::update, neighbour) );
    core->on_append_setup( bind(&force_type::compute, force) );

    // integration step
    core->on_integrate( bind(&integrator_type::integrate, integrator) );
    core->on_force( bind(&force_type::compute, force) );
    core->on_finalize( bind(&integrator_type::finalize, integrator) );

    // update neighbour lists if maximum squared displacement is greater than (skin/2)²
    float_type limit = pow(neighbour->r_skin() / 2, 2);
    shared_ptr<greater_type> greater =
        make_shared<greater_type>(bind(&max_displacement_type::compute, max_displacement), limit);
    greater->on_greater( bind(&max_displacement_type::zero, max_displacement) );
    greater->on_greater( bind(&binning_type::update, binning) );
    greater->on_greater( bind(&neighbour_type::update, neighbour) );
    core->on_prepend_force( bind(&greater_type::evaluate, greater) );
}

template <int dimension, typename float_type>
struct host_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::potentials::lennard_jones<float_type> potential_type;
    typedef mdsim::host::forces::pair_trunc<dimension, float_type, potential_type> force_type;
    typedef mdsim::host::binning<dimension, float_type> binning_type;
    typedef mdsim::host::neighbour<dimension, float_type> neighbour_type;
    typedef mdsim::host::maximum_squared_displacement<dimension, float_type> max_displacement_type;
    typedef mdsim::host::integrators::verlet_nvt_hoover<dimension, float_type> integrator_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::positions::lattice<dimension, float_type> position_type;
    typedef halmd::random::host::random random_type;
    typedef mdsim::host::velocities::boltzmann<dimension, float_type> velocity_type;
    typedef observables::host::thermodynamics<dimension, float_type> thermodynamics_type;
    static bool const gpu = false;
};

BOOST_AUTO_TEST_CASE( verlet_nvt_hoover_host_2d ) {
    verlet_nvt_hoover<host_modules<2, double> >().test();
}
BOOST_AUTO_TEST_CASE( verlet_nvt_hoover_host_3d ) {
    verlet_nvt_hoover<host_modules<3, double> >().test();
}

#ifdef WITH_CUDA
template <int dimension, typename float_type>
struct gpu_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::potentials::lennard_jones<float_type> potential_type;
    typedef mdsim::gpu::forces::pair_trunc<dimension, float_type, potential_type> force_type;
    typedef mdsim::gpu::binning<dimension, float_type> binning_type;
    typedef mdsim::gpu::neighbours::from_binning<dimension, float_type> neighbour_type;
    typedef mdsim::gpu::maximum_squared_displacement<dimension, float_type> max_displacement_type;
    typedef mdsim::gpu::integrators::verlet_nvt_hoover<dimension, double> integrator_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::positions::lattice<dimension, float_type, halmd::random::gpu::rand48> position_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    typedef observables::gpu::thermodynamics<dimension, float_type> thermodynamics_type;
    typedef mdsim::gpu::velocities::boltzmann<dimension, float_type, halmd::random::gpu::rand48> velocity_type;
    static bool const gpu = true;
};

BOOST_FIXTURE_TEST_CASE( verlet_nvt_hoover_gpu_2d, device ) {
    verlet_nvt_hoover<gpu_modules<2, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( verlet_nvt_hoover_gpu_3d, device ) {
    verlet_nvt_hoover<gpu_modules<3, float> >().test();
}
#endif // WITH_CUDA
