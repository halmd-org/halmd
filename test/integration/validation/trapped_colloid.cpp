/*
 * Copyright © 2024 Max Orteu
 * Copyright © 2024 Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE trapped_colloid
#include <boost/test/unit_test.hpp>

#include <halmd/io/utility/hdf5.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <test/tools/ctest.hpp>

#include <h5xx/h5xx.hpp>
#include <limits>
#include <cmath>

using namespace halmd;

/**
 * test simulation results for overdamped Brownian integrator with an external
 * harmonic potential:
 * @f[ \dot{\vec r}(t) = - K \vec r(t) + \sqrt{2D} \vec\xi(t), @f]
 * where @f$ \vec\xi(t) @f$ is unit white noise.
 *
 * reference values come from the analytical solution of the Ornstein-Uhlenbeck
 * process
 */

#ifndef USE_HOST_SINGLE_PRECISION
const double eps = std::numeric_limits<double>::epsilon();
#else
const double eps = std::numeric_limits<float>::epsilon();
#endif
const float eps_float = std::numeric_limits<float>::epsilon();

// simulation parameters
struct
{
    double T = 3;       // temperature
    double K = 2;       // potential stiffness
    double D = 0.3;     // diffusion constant
    double gamma = D * K / T;
    double sigma2 = 2 * D;
} const params;

template <int dimension>
void check_cm(
    std::vector<double> const& en_pot_time
  , std::vector<fixed_vector<double, dimension>> const& cm_data
)
{
    double x0 = cm_data[0][0];      // initial position

    for (unsigned int i = 0; i < en_pot_time.size(); ++i) {
        double t = en_pot_time[i];
        auto const& cm = cm_data[i];

        // reference value from analytical solution
        double ref_cm = std::exp(-t * params.gamma) * x0;
        double var_cm = (1 - std::exp(-2 * t * params.gamma)) * params.sigma2 / (2 * params.gamma);

        // absolute tolerance, using 4.5-sigma rule
        double tol = 4.5 * std::sqrt(var_cm);
        for (unsigned int j = 0; j < cm.size(); ++j) {
            BOOST_CHECK_SMALL(cm[j] - ref_cm, tol);
        }
    }
}

void check_en_pot(
    std::vector<double> const& en_pot_time
  , std::vector<double> const& en_pot_data
  , double x0   // initial position
  , unsigned int dimension
)
{
    for (unsigned int i = 0; i < en_pot_data.size(); ++i) {
        double t = en_pot_time[i];
        auto const& en_pot = en_pot_data[i];

        // reference value from analytical solution,
        // use mean and variance of centre-of-mass position
        double cm = std::exp(-t * params.gamma) * x0;
        double var_cm = (1 - std::exp(-2 * t * params.gamma)) * params.sigma2 / (2 * params.gamma);
        double ref_en_pot = (params.K / 2) * dimension * (var_cm +  cm * cm);

        // analytical solutions for the variance
        double var_en_pot = std::pow(params.sigma2, 2) / (2 * std::pow(params.gamma, 2));

        // relative tolerance, using 4.5-sigma rule
        double tol = 4.5 * std::sqrt(var_en_pot) / ref_en_pot;
        BOOST_CHECK_CLOSE_FRACTION(en_pot, ref_en_pot, tol);
    }
}

void check_msd(
    boost::multi_array<double, 2> const& msd_time
  , boost::multi_array<double, 2> const& msd_value
  , boost::multi_array<double, 2> const& msd_error
  , unsigned int dimension
)
{
    // iterate over the three arrays simultaneously,
    // the outer loop is over all blocks
    for (unsigned int i = 0; i < msd_time.size(); ++i) {            // number of blocks
        auto block_time = msd_time[i];
        auto block_value = msd_value[i];
        auto block_error = msd_error[i];
        for (unsigned int j = 0; j < block_time.size(); ++j) {      // block size
            double time = block_time[j];
            double value = block_value[j];
            double error = block_error[j];

            // reference: analytical solution for one Cartesian component,
            // MSD(t) = (1 - exp(- γ t)) * (σ² / γ)
            double ref_msd = dimension * (1 - std::exp(-time * params.gamma)) * params.sigma2 / params.gamma;

            // relative tolerance, using 4.5-sigma rule
            double tol = 4.5 * std::max(error, 1e-3) / ref_msd;               // enforce minimum on standard error
            BOOST_CHECK_CLOSE_FRACTION(value, ref_msd, tol);
        }
    }
}

BOOST_AUTO_TEST_CASE( validation )
{
    // parse command-line arguments:
    // program_name filename USE_HOST
    using namespace boost::unit_test::framework;

    BOOST_REQUIRE_EQUAL( master_test_suite().argc, 2 );
    char** argv = master_test_suite().argv;

    // open file from production run
    H5::H5File file(argv[1], H5F_ACC_RDONLY);
    H5::Group observables = file.openGroup("observables");
    H5::Group dynamics = file.openGroup("dynamics/all");

    // read parameters
    unsigned int npart;
    h5xx::read_dataset(observables.openDataSet("particle_number"), npart);
    unsigned int dimension = h5xx::read_attribute<unsigned int>(observables, "dimension");

    BOOST_TEST_MESSAGE("space dimension: " << dimension);
    BOOST_TEST_MESSAGE("number of particles: " << npart);

    // read and check time series data for thermodynamic observables:
    // centre-of-mass position
    std::vector<double> cm_time;
    h5xx::read_dataset(observables.openDataSet("center_of_mass/time"), cm_time);

    double x0;      // initial position
    if (dimension == 2) {
        std::vector<fixed_vector<double, 2>> cm_data;
        h5xx::read_dataset(observables.openDataSet("center_of_mass/value"), cm_data);
        x0 = cm_data[0][0];

        check_cm<2>(cm_time, cm_data);
    }
    else {
        std::vector<fixed_vector<double, 3>> cm_data;
        h5xx::read_dataset(observables.openDataSet("center_of_mass/value"), cm_data);
        x0 = cm_data[0][0];

        check_cm<3>(cm_time, cm_data);
    }

    // read and check potential energy
    std::vector<double> en_pot_data, en_pot_time;
    h5xx::read_dataset(observables.openDataSet("potential_energy/value"), en_pot_data);
    h5xx::read_dataset(observables.openDataSet("potential_energy/time"), en_pot_time);

    check_en_pot(en_pot_time, en_pot_data, x0, dimension);

    // read mean-square displacement
    boost::multi_array<double, 2> msd_value, msd_time, msd_error;
    h5xx::read_dataset(dynamics.openDataSet("mean_square_displacement/value"), msd_value);
    h5xx::read_dataset(dynamics.openDataSet("mean_square_displacement/time"), msd_time);
    h5xx::read_dataset(dynamics.openDataSet("mean_square_displacement/error"), msd_error);

    check_msd(msd_time, msd_value, msd_error, dimension);
}

