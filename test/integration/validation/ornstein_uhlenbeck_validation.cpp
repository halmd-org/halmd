/*
 * Copyright © 2010-2014 Felix Höfling
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

#define BOOST_TEST_MODULE ornstein_uhlenbeck_validation // need to write appropriate name
#include <boost/test/unit_test.hpp>

#include <halmd/io/utility/hdf5.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <test/tools/ctest.hpp>

#include <boost/lexical_cast.hpp>
#include <h5xx/h5xx.hpp>
#include <limits>
#include <cmath>

// using namespace boost;
using namespace halmd;

/**
 * test simulation results for brownian integrator with harmonic potential (colloid?)
 *
 * reference values come from analytical solution for Ornstein-Uhlenbeck process
 */

#ifndef USE_HOST_SINGLE_PRECISION
const double eps = std::numeric_limits<double>::epsilon();
#else
const double eps = std::numeric_limits<float>::epsilon();
#endif
const float eps_float = std::numeric_limits<float>::epsilon();

struct simulation_parameters
{
    double T;
    double K;
    double D;
    double gamma;
    double sigma_sq;

    simulation_parameters() 
        : T(3)
        , K(2)
        , D(0.3)
        , gamma(D * K / T)
        , sigma_sq(2 * D)
    {
    }
};

simulation_parameters params;

template <int dimension> 
void parse_cm_data(
    H5::Group& observables
  , std::vector<halmd::fixed_vector<double, dimension>>& cm_data
)
{
    static_assert(dimension == 2 || dimension == 3, "Dimension must be either 2 or 3");

    std::vector<halmd::fixed_vector<double, dimension>> cm_data_raw;
    h5xx::read_dataset(observables.openDataSet("center_of_mass/value"), cm_data_raw);
    cm_data.resize(cm_data_raw.size());

    for (size_t i = 0; i < cm_data_raw.size(); ++i) {
        for (int j = 0; j < dimension; ++j) {
            cm_data[i][j] = cm_data_raw[i][j];
        }
    }
}

template <int dimension> 
void check_cm(
    std::vector<double>& en_pot_time
  , std::vector<halmd::fixed_vector<double, dimension>> cm_data
) 
{
    double ti, exp_cm, var_cm;
    double x0 = cm_data[0][0];

    for (unsigned int i = 0; i < en_pot_time.size(); ++i) {
        ti = en_pot_time[i];

        exp_cm = std::exp(-ti*params.gamma)*x0;
        var_cm = (1 - std::exp(-2*ti*params.gamma))*params.sigma_sq/(2*params.gamma);

        for (int j = 0; j < dimension; ++j) {
            BOOST_CHECK_SMALL(cm_data[i][j] - exp_cm, 4.5*std::sqrt(var_cm));
        }
    }
}

template <int dimension>
void check_en_pot(
    int npart
  , std::vector<double>& en_pot_time
  , std::vector<double>& en_pot_data
  , std::vector<halmd::fixed_vector<double, dimension>> cm_data
) 
{
    double ti, var_cm, var_en_pot, r_cm;

    for (unsigned int i = 0; i < en_pot_data.size(); ++i) {
        ti = en_pot_time[i];

        r_cm = 0;
        for (int j = 0; j < dimension; ++j) {
            r_cm += std::pow(cm_data[i][j], 2);
        }
        r_cm = std::sqrt(r_cm);

        var_cm = (1 - std::exp(-2*ti*params.gamma))*params.sigma_sq/(2*params.gamma);
        var_en_pot = std::pow(params.sigma_sq,2) / (2*std::pow(params.gamma,2));
        BOOST_CHECK_CLOSE(en_pot_data[i], (params.K/2)*(dimension*var_cm + std::pow(r_cm,2)), 4.5*std::sqrt(var_en_pot));
    }
}

template <int dimension>
void check_msd(
    std::vector<double>& msd_time_flat
  , std::vector<double>& msd_data_flat
  , std::vector<double>& msd_error_flat
)
{
    double ti, exp_msd;

    // MSD values are taken only from equilibration onwards
    for (unsigned int i = 0; i < msd_data_flat.size(); ++i) {
        ti = msd_time_flat[i];
        exp_msd = (1 - std::exp(-ti*params.gamma))*params.sigma_sq/params.gamma;

        BOOST_CHECK_SMALL(std::abs(msd_data_flat[i]/dimension - exp_msd), 4.5*std::max(msd_error_flat[i], 1e-3)); // std error is corrected
    }
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
    H5::Group dynamics = file.openGroup("dynamics/all");

    // read parameters
    unsigned int npart;
    h5xx::read_dataset(observables.openDataSet("particle_number"), npart);
    unsigned int dimension = h5xx::read_attribute<unsigned int>(observables, "dimension");

    BOOST_TEST_MESSAGE("space dimension: " << dimension);
    BOOST_TEST_MESSAGE("number of particles: " << npart);

    size_t steps;
    h5xx::read_chunked_dataset(observables.openDataSet("potential_energy/step"), steps, ssize_t(-1));

    // read time series data for thermodynamic observables
    std::vector<double> en_pot_data, en_pot_time;

    h5xx::read_dataset(observables.openDataSet("potential_energy/value"), en_pot_data);
    h5xx::read_dataset(observables.openDataSet("potential_energy/time"), en_pot_time);

    std::vector<halmd::fixed_vector<double, 10>> msd_data;
    std::vector<halmd::fixed_vector<double, 10>> msd_time;
    std::vector<halmd::fixed_vector<double, 10>> msd_error;
    h5xx::read_dataset(dynamics.openDataSet("mean_square_displacement/value"), msd_data); 
    h5xx::read_dataset(dynamics.openDataSet("mean_square_displacement/time"), msd_time);
    h5xx::read_dataset(dynamics.openDataSet("mean_square_displacement/error"), msd_error);

    std::vector<double> msd_data_flat(40), msd_time_flat(40), msd_error_flat(40);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 10; j++) {
            msd_data_flat[10*i+j] = msd_data[i][j];
            msd_time_flat[10*i+j] = msd_time[i][j];
            msd_error_flat[10*i+j] = msd_error[i][j];    
        }
    }

    if (dimension == 2) {
        std::vector<halmd::fixed_vector<double, 2>> cm_data;
        parse_cm_data<2>(observables, cm_data);

        check_cm<2>(en_pot_time, cm_data);
        check_en_pot<2>(npart, en_pot_time, en_pot_data, cm_data);
        check_msd<2>(msd_time_flat, msd_data_flat, msd_error_flat);
    }
    else {
        std::vector<halmd::fixed_vector<double, 3>> cm_data;
        parse_cm_data<3>(observables, cm_data);

        check_cm<3>(en_pot_time, cm_data);
        check_en_pot<3>(npart, en_pot_time, en_pot_data, cm_data);
        check_msd<3>(msd_time_flat, msd_data_flat, msd_error_flat);
    } 

    
}

