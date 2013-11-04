/*
 * Copyright © 2010-2012 Felix Höfling
 * Copyright © 2010-2012 Peter Colberg
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

#define BOOST_TEST_MODULE thermodynamics
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <boost/numeric/ublas/banded.hpp>
#include <limits>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/forces/pair_trunc.hpp>
#include <halmd/mdsim/host/integrators/verlet.hpp>
#include <halmd/mdsim/host/integrators/verlet_nvt_andersen.hpp>
#include <halmd/mdsim/host/max_displacement.hpp>
#include <halmd/mdsim/host/neighbours/from_binning.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/all.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/mdsim/host/potentials/lennard_jones.hpp>
#include <halmd/mdsim/host/velocities/boltzmann.hpp>
#include <halmd/mdsim/host/velocity.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/thermodynamics.hpp>
#include <halmd/random/host/random.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/forces/pair_trunc.hpp>
# include <halmd/mdsim/gpu/integrators/verlet.hpp>
# include <halmd/mdsim/gpu/integrators/verlet_nvt_andersen.hpp>
# include <halmd/mdsim/gpu/max_displacement.hpp>
# include <halmd/mdsim/gpu/neighbours/from_binning.hpp>
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_groups/all.hpp>
# include <halmd/mdsim/gpu/positions/lattice.hpp>
# include <halmd/mdsim/gpu/potentials/lennard_jones.hpp>
# include <halmd/mdsim/gpu/velocities/boltzmann.hpp>
# include <halmd/mdsim/gpu/velocity.hpp>
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
 * test simulation results for thermodynamic variables in the dilute limit
 *
 * reference values for the state variables of a 3-dim. LJ fluid can be found in
 * L. Verlet, Phys. Rev. 159, 98 (1967)
 * and in Hansen & McDonald, Table 4.2
 *
 * a more detailed analysis and more accurate values are found in
 * Johnson, Zollweg, and Gubbins, Mol. Phys. 78, 591 (1993).
 *
 * and we refer to results from integral equations theory:
 * Ayadim, Oettel, Amokrane, J. Phys.: Condens. Matter 21, 115103 (2009).
 */

const double eps = numeric_limits<double>::epsilon();
const float eps_float = numeric_limits<float>::epsilon();

/**
 * heat capacity from microcanonical fluctuations of kinetic energy
 * see Lebowitz, Percus, and Verlet, Phys. Rev. 153, 250 (1967) for details
 */
inline double heat_capacity_nve(double en_kin, double variance, unsigned npart)
{
    return 1 / (2./3 - npart * variance / (en_kin * en_kin));
}

template <typename modules_type>
struct lennard_jones_fluid
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::potential_type potential_type;
    typedef typename modules_type::force_type force_type;
    typedef typename modules_type::binning_type binning_type;
    typedef typename modules_type::neighbour_type neighbour_type;
    typedef typename modules_type::msd_type msd_type;
    typedef typename modules_type::nve_integrator_type nve_integrator_type;
    typedef typename modules_type::nvt_integrator_type nvt_integrator_type;
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::particle_group_type particle_group_type;
    typedef typename modules_type::position_type position_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::thermodynamics_type thermodynamics_type;
    typedef typename modules_type::velocity_type velocity_type;
    static bool const gpu = modules_type::gpu;

    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    static unsigned int const dimension = vector_type::static_size;

    float density;
    float temp;
    float rc;
    float epsilon;
    float sigma;
    float skin;
    double timestep;
    unsigned npart;
    fixed_vector<double, dimension> box_ratios;
    fixed_vector<double, dimension> slab;

    std::shared_ptr<box_type> box;
    std::shared_ptr<potential_type> potential;
    std::shared_ptr<force_type> force;
    std::shared_ptr<binning_type> binning;
    std::shared_ptr<neighbour_type> neighbour;
    std::shared_ptr<msd_type> msd;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<position_type> position;
    std::shared_ptr<random_type> random;
    std::shared_ptr<thermodynamics_type> thermodynamics;
    std::shared_ptr<velocity_type> velocity;

    void test();
    lennard_jones_fluid();
};

template <typename modules_type>
void lennard_jones_fluid<modules_type>::test()
{
    // create NVT integrator
    std::shared_ptr<nvt_integrator_type> nvt_integrator = std::make_shared<nvt_integrator_type>(
        particle
      , force
      , box
      , random
      , 0.005 /* time step */
      , temp
      , 1 /* collision rate */
    );

    // relax configuration and thermalise at given temperature, run for t*=30
    BOOST_TEST_MESSAGE("thermalise initial state at T=" << temp);
    position->set();
    velocity->set();
    unsigned int steps = static_cast<unsigned int>(ceil(30 / nvt_integrator->timestep()));
    for (unsigned int i = 0; i < steps; ++i) {
        nvt_integrator->integrate();
        nvt_integrator->finalize();
    }

    // set different timestep and choose NVE integrator
    std::shared_ptr<nve_integrator_type> nve_integrator = std::make_shared<nve_integrator_type>(particle, force, box, timestep);

    // stochastic thermostat => centre particle velocities around zero
    fixed_vector<double, dimension> v_cm = thermodynamics->v_cm();
    BOOST_TEST_MESSAGE("shift particle velocities by" << v_cm);
    shift_velocity(*particle, -v_cm);

    const double vcm_limit = gpu ? 0.5 * eps_float : 30 * eps;
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    // take averages of fluctuating quantities,
    accumulator<double> temp_, press, en_pot;

    // equilibration run, measure temperature in second half
    BOOST_TEST_MESSAGE("equilibrate initial state");
    steps = static_cast<unsigned int>(ceil(30 / timestep));
    unsigned int period = static_cast<unsigned int>(round(0.01 / timestep));
    for (unsigned int i = 0; i < steps; ++i) {
        nve_integrator->integrate();
        nve_integrator->finalize();
        if(i > steps / 2 && i % period == 0) {
            temp_(thermodynamics->temp());
        }
    }
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    // rescale velocities to match the exact temperature
    double v_scale = sqrt(temp / mean(temp_));
    BOOST_TEST_MESSAGE("rescale velocities by factor " << v_scale);
    rescale_velocity(*particle, v_scale);

    double en_tot = thermodynamics->en_tot();
    double max_en_diff = 0; // maximal absolut deviation from initial total energy
    temp_.reset();

    // microcanonical simulation run
    BOOST_TEST_MESSAGE("run NVE simulation");
    steps = static_cast<unsigned int>(ceil(60 / nve_integrator->timestep()));
    period = static_cast<unsigned int>(round(
        ((dimension == 3) ? 0.32 : 0.45) / nve_integrator->timestep()
    )); // relaxation time (from VACF)
    for (unsigned int i = 0; i < steps; ++i) {
        // perform MD step
        nve_integrator->integrate();
        nve_integrator->finalize();

        // measurement
        if(i % period == 0) {
            temp_(thermodynamics->temp());
            press(thermodynamics->pressure());
            en_pot(thermodynamics->en_pot());
            max_en_diff = max(abs(thermodynamics->en_tot() - en_tot), max_en_diff);
        }
    }
    BOOST_CHECK_SMALL(norm_inf(thermodynamics->v_cm()), vcm_limit);

    // with the first released version of HAL's MD package (commit f5283a2),
    // an energy drift of less than 5e-6 ε was obtained over 2e8 MD steps
    // using a potential with smooth cutoff (dt*=0.001, h=0.005)
    const double en_limit = max(3e-5, steps * 1e-12);
    BOOST_CHECK_SMALL(max_en_diff / fabs(en_tot), en_limit);

    // use tolerance of 4.5σ, see below;
    // empirically determined standard deviation from many test runs:
    // σ(T) = 0.004 for N=4000 and σ(T) = 0.007 for N=1500
    BOOST_CHECK_CLOSE_FRACTION(temp, mean(temp_), 4.5 * (gpu ? 0.004 : 0.007) / temp);
    BOOST_CHECK_CLOSE_FRACTION(density, (float)thermodynamics->density(), eps_float);

    // compute response coefficients from fluctuations
    double cV = heat_capacity_nve(mean(temp_), variance(temp_), npart);

    // long-tail corrections for trunctated LJ potential,
    // see e.g. book by Allen & Tildesley or Johnsen et al. (1993)
    //
    // \f$ U_\mathrm{corr} = 0.5 \rho S_d \int_{r_c}^\infty r^{d-1} U(r) dr \f$
    // where \f$ S_d \f$ denotes the surface of the unit sphere
    double en_corr = (dimension == 3) ?
        8./9 * M_PI * density * (pow(rc, -6) - 3) * pow(rc, -3)
      : 2./5 * M_PI * density * (pow(rc, -6) - 2.5) * pow(rc, -4);

    // \f$ P_\mathrm{corr} = - \rho^2 S_d / (2 * d) \int_{r_c}^\infty r^{d-1} r U'(r) dr \f$
    // where \f$ S_d \f$ denotes the surface of the unit sphere
    double press_corr = (dimension == 3) ?
        32./9 * M_PI * pow(density, 2) * (pow(rc, -6) - 1.5) * pow(rc, -3)
      : 12./5 * M_PI * pow(density, 2) * (pow(rc, -6) - 1.25) * pow(rc, -4);

    BOOST_TEST_MESSAGE("Density: " << density);
    BOOST_TEST_MESSAGE("Temperature: " << mean(temp_) << " ± " << error_of_mean(temp_));
    BOOST_TEST_MESSAGE("Pressure: " << mean(press) << " ± " << error_of_mean(press));
    BOOST_TEST_MESSAGE("Pressure (corrected): " << mean(press) + press_corr);
    BOOST_TEST_MESSAGE("Potential energy: " << mean(en_pot) << " ± " << error_of_mean(en_pot));
    BOOST_TEST_MESSAGE("Potential energy (corrected): " << mean(en_pot) + en_corr);
    BOOST_TEST_MESSAGE("β P / ρ = " << mean(press) / mean(temp_) / density);
    BOOST_TEST_MESSAGE("β U / N = " << mean(en_pot) / mean(temp_));
    BOOST_TEST_MESSAGE("Heat capacity c_V = " << cV);

    // tolerances are 4.5σ, where σ is the standard deviation of the test
    // results, with this choice, the test should pass with 99.999% probability
    //
    // quoted values for σ were obtained empirically from at least 100 test runs
    if (dimension == 3) {
        // values from Johnson et al.: P = 1.023, Epot = -1.673  (Npart = 864)
        // values from RFA theory (Ayadim et al.): P = 1.0245, Epot = -1.6717,
        // standard deviations: σ(P) = 0.003, σ(Epot) = 0.002 for N = 4000,
        // σ(P) = 0.005, σ(Epot) = 0.003 for N = 1500
        BOOST_CHECK_CLOSE_FRACTION(mean(press), 1.023, 4.5 * (gpu ? 0.003 : 0.005) / 1.023);
        BOOST_CHECK_CLOSE_FRACTION(mean(en_pot), -1.673, 4.5 * (gpu ? 0.002 : 0.003) / 1.673);
        // our own measurements using HAL's MD package FIXME find reference values:
        // c_V = 1.648, σ(c_V) = 0.02 for GPU and host test cases
        BOOST_CHECK_CLOSE_FRACTION(cV, 1.648, 4.5 * 0.02 / 1.648);
    }
    else if (dimension == 2) {
        // our own measurements using HAL's MD package FIXME find reference values
        // P = 1.235, σ(P) = 0.004 for N = 4000, σ(P) = 0.006 for N = 1500
        // Epot = -0.589, σ(Epot) = 0.002 for N = 4000, σ(Epot) = 0.003 for N = 1500
        BOOST_CHECK_CLOSE_FRACTION(mean(press), 1.235, 4.5 * (gpu ? 0.004 : 0.006) / 1.235);
        BOOST_CHECK_CLOSE_FRACTION(mean(en_pot), -0.589, 4.5 * (gpu ? 0.002 : 0.003) / 0.589);
        // c_V = 1.68, σ(c_V) = 0.03 for GPU and host test cases
        BOOST_CHECK_CLOSE_FRACTION(cV, 1.68, 4.5 * 0.03 / 1.68);
    }
}

template <typename modules_type>
lennard_jones_fluid<modules_type>::lennard_jones_fluid()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    // set module parameters
    density = 0.3;
    temp = 3;
    rc = 4;
    epsilon = 1;
    sigma = 1;
    skin = 0.5;
    timestep = 0.001;
    npart = gpu ? 4000 : 1500;
    box_ratios = (dimension == 3) ? list_of(1.)(2.)(1.01) : list_of(1.)(2.);
    double det = accumulate(box_ratios.begin(), box_ratios.end(), 1., multiplies<double>());
    double volume = npart / density;
    double edge_length = pow(volume / det, 1. / dimension);
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = edge_length * box_ratios[i];
    }
    slab = 1;

    typename potential_type::matrix_type rc_mat(1, 1);
    typename potential_type::matrix_type epsilon_mat(1, 1);
    typename potential_type::matrix_type sigma_mat(1, 1);
    rc_mat(0, 0) = rc;
    epsilon_mat(0, 0) = epsilon;
    sigma_mat(0, 0) = sigma;

    // create modules
    random = std::make_shared<random_type>();
    particle = std::make_shared<particle_type>(npart, 1);
    box = std::make_shared<box_type>(edges);
    potential = std::make_shared<potential_type>(particle->nspecies(), particle->nspecies(), rc_mat, epsilon_mat, sigma_mat);
    binning = std::make_shared<binning_type>(particle, box, potential->r_cut(), skin);
    msd = std::make_shared<msd_type>(particle, box);
    neighbour = std::make_shared<neighbour_type>(std::make_pair(particle, particle), std::make_pair(binning, binning), msd, box, potential->r_cut(), skin);
    position = std::make_shared<position_type>(particle, box, slab);
    velocity = std::make_shared<velocity_type>(particle, random, temp);
    force = std::make_shared<force_type>(potential, particle, particle, box, neighbour);
    std::shared_ptr<particle_group_type> group = std::make_shared<particle_group_type>(particle);
    thermodynamics = std::make_shared<thermodynamics_type>(particle, force, group, box);
}

template <int dimension, typename float_type>
struct host_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::potentials::lennard_jones<float_type> potential_type;
    typedef mdsim::host::forces::pair_trunc<dimension, float_type, potential_type> force_type;
    typedef mdsim::host::binning<dimension, float_type> binning_type;
    typedef mdsim::host::neighbours::from_binning<dimension, float_type> neighbour_type;
    typedef mdsim::host::max_displacement<dimension, float_type> msd_type;
    typedef mdsim::host::integrators::verlet<dimension, float_type> nve_integrator_type;
    typedef mdsim::host::integrators::verlet_nvt_andersen<dimension, float_type> nvt_integrator_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::particle_groups::all<particle_type> particle_group_type;
    typedef mdsim::host::positions::lattice<dimension, float_type> position_type;
    typedef halmd::random::host::random random_type;
    typedef mdsim::host::velocities::boltzmann<dimension, float_type> velocity_type;
    typedef observables::host::thermodynamics<dimension, float_type> thermodynamics_type;
    static bool const gpu = false;
};

BOOST_AUTO_TEST_CASE( lennard_jones_fluid_host_2d ) {
    lennard_jones_fluid<host_modules<2, double> >().test();
}
BOOST_AUTO_TEST_CASE( lennard_jones_fluid_host_3d ) {
    lennard_jones_fluid<host_modules<3, double> >().test();
}

#ifdef HALMD_WITH_GPU
template <int dimension, typename float_type>
struct gpu_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::potentials::lennard_jones<float_type> potential_type;
    typedef mdsim::gpu::forces::pair_trunc<dimension, float_type, potential_type> force_type;
    typedef mdsim::gpu::binning<dimension, float_type> binning_type;
    typedef mdsim::gpu::neighbours::from_binning<dimension, float_type> neighbour_type;
    typedef mdsim::gpu::max_displacement<dimension, float_type> msd_type;
    typedef mdsim::gpu::integrators::verlet<dimension, float_type> nve_integrator_type;
    typedef mdsim::gpu::integrators::verlet_nvt_andersen<dimension, float_type, halmd::random::gpu::rand48> nvt_integrator_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::particle_groups::all<particle_type> particle_group_type;
    typedef mdsim::gpu::positions::lattice<dimension, float_type> position_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    typedef observables::gpu::thermodynamics<dimension, float_type> thermodynamics_type;
    typedef mdsim::gpu::velocities::boltzmann<dimension, float_type, halmd::random::gpu::rand48> velocity_type;
    static bool const gpu = true;
};

BOOST_FIXTURE_TEST_CASE( lennard_jones_fluid_gpu_2d, device ) {
    lennard_jones_fluid<gpu_modules<2, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( lennard_jones_fluid_gpu_3d, device ) {
    lennard_jones_fluid<gpu_modules<3, float> >().test();
}
#endif // HALMD_WITH_GPU
