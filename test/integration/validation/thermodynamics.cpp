/*
 * Copyright © 2010-2014 Felix Höfling
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

#include <boost/lexical_cast.hpp>
#include <h5xx/h5xx.hpp>
#include <limits>

#include <halmd/numeric/accumulator.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <test/tools/ctest.hpp>

// using namespace boost;
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

BOOST_AUTO_TEST_CASE( validation )
{
    // parse command-line arguments:
    // program_name filename USE_HOST
    using namespace boost::unit_test::framework;

    BOOST_REQUIRE_EQUAL( master_test_suite().argc, 3 );
    char** argv = master_test_suite().argv;

    bool gpu = boost::lexical_cast<bool>(argv[2]);

    // open file from production run
    H5::H5File file(argv[1], H5F_ACC_RDONLY);
    H5::Group observables = file.openGroup("observables");

    // read parameters
    unsigned int dimension = h5xx::read_attribute<unsigned int>(observables, "dimension");
    unsigned int npart;
    h5xx::read_dataset(observables.openDataSet("particle_number"), npart);
    double density;
    h5xx::read_dataset(observables.openDataSet("density"), density);

    BOOST_TEST_MESSAGE("space dimension: " << dimension);
    BOOST_TEST_MESSAGE("number of particles: " << npart);
    BOOST_TEST_MESSAGE("particle density: " << density);

    double temperature = 3;
    double rc = 4;
    BOOST_CHECK_CLOSE_FRACTION(density, 0.3, eps_float);

    // CM velocity at first step
    const double vcm_limit = gpu ? 0.5 * eps_float : 30 * eps;
    fixed_vector<double, 3> v_cm;
    h5xx::detail::read_chunked_dataset<double, 1>(
        observables.openDataSet("center_of_mass_velocity/value")
      , reinterpret_cast<double*>(&v_cm)
      , 0
    );
//    h5xx::read_chunked_dataset(
//        observables.openDataSet("center_of_mass_velocity/value")
//      , v_cm
//      , 0
//    );
// FIXME we use the generic method on raw pointers since halmd::fixed_vector
// inherits from std::array, which is not yet supported by h5xx. Second,
// "reinterpret_cast<boost::array<double, 3>*>(&v_cm)" would violate
// strict-aliasing rules

    BOOST_CHECK_SMALL(norm_inf(v_cm), vcm_limit);

    // CM velocity at last step
    h5xx::detail::read_chunked_dataset<double, 1>(  // FIXME see above
        observables.openDataSet("center_of_mass_velocity/value")
      , reinterpret_cast<double*>(&v_cm)
      , -1
    );
    BOOST_CHECK_SMALL(norm_inf(v_cm), vcm_limit);

    // number of integration steps
    size_t steps;
    h5xx::read_chunked_dataset(observables.openDataSet("temperature/step"), steps, ssize_t(-1));

    // read time series data for thermodynamic observables
    std::vector<double> en_tot_data, en_pot_data, temp_data, press_data;
    h5xx::read_dataset(observables.openDataSet("internal_energy/value"), en_tot_data);
    h5xx::read_dataset(observables.openDataSet("potential_energy/value"), en_pot_data);
    h5xx::read_dataset(observables.openDataSet("pressure/value"), press_data);
    h5xx::read_dataset(observables.openDataSet("temperature/value"), temp_data);

    // total internal energy at first step
    double en_tot = en_tot_data[0];
    double max_en_diff = 0; // maximal absolut deviation from initial total energy

    // take averages of fluctuating quantities,
    accumulator<double> temp, press, en_pot;
    for (unsigned int i = 0; i < temp_data.size(); ++i) {
        temp(temp_data[i]);
        press(press_data[i]);
        en_pot(en_pot_data[i]);
        max_en_diff = max(abs(en_tot_data[i] - en_tot), max_en_diff);
    }

    // with the first released version of HAL's MD package (commit f5283a2),
    // an energy drift of less than 5e-6 ε was obtained over 2e8 MD steps
    // using a potential with smooth cutoff (dt*=0.001, h=0.005)
    const double en_limit = max(3e-5, steps * 1e-12);
    BOOST_CHECK_SMALL(max_en_diff / fabs(en_tot), en_limit);

    // use tolerance of 4.5σ, see below;
    // empirically determined standard deviation from many test runs:
    // σ(T) = 0.004 for N=4000 and σ(T) = 0.007 for N=1500
    BOOST_CHECK_CLOSE_FRACTION(temperature, mean(temp), 4.5 * (gpu ? 0.004 : 0.007) / temperature);

    // compute response coefficients from fluctuations
    double cV = heat_capacity_nve(mean(temp), variance(temp), npart);

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
    BOOST_TEST_MESSAGE("Temperature: " << mean(temp) << " ± " << error_of_mean(temp));
    BOOST_TEST_MESSAGE("Pressure: " << mean(press) << " ± " << error_of_mean(press));
    BOOST_TEST_MESSAGE("Pressure (corrected): " << mean(press) + press_corr);
    BOOST_TEST_MESSAGE("Potential energy: " << mean(en_pot) << " ± " << error_of_mean(en_pot));
    BOOST_TEST_MESSAGE("Potential energy (corrected): " << mean(en_pot) + en_corr);
    BOOST_TEST_MESSAGE("β P / ρ = " << mean(press) / mean(temp) / density);
    BOOST_TEST_MESSAGE("β U / N = " << mean(en_pot) / mean(temp));
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
        // our own measurements using HAL's MD package:
        // c_V = 1.648, σ(c_V) = 0.02 for GPU and host test cases
        // FIXME find reference values from independent sources
        BOOST_CHECK_CLOSE_FRACTION(cV, 1.648, 4.5 * 0.02 / 1.648);
    }
    else if (dimension == 2) {
        // our own measurements using HAL's MD package:
        // P = 1.235, σ(P) = 0.004 for N = 4000, σ(P) = 0.006 for N = 1500
        // Epot = -0.589, σ(Epot) = 0.002 for N = 4000, σ(Epot) = 0.003 for N = 1500
        // FIXME find reference values from independent sources
        BOOST_CHECK_CLOSE_FRACTION(mean(press), 1.235, 4.5 * (gpu ? 0.004 : 0.006) / 1.235);
        BOOST_CHECK_CLOSE_FRACTION(mean(en_pot), -0.589, 4.5 * (gpu ? 0.002 : 0.003) / 0.589);
        // c_V = 1.68, σ(c_V) = 0.03 for GPU and host test cases
        BOOST_CHECK_CLOSE_FRACTION(cV, 1.68, 4.5 * 0.03 / 1.68);
    }
}

