/*
 * Copyright © 2011  Felix Höfling
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

#define BOOST_TEST_MODULE ssf
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

#include <halmd/io/logger.hpp>
#include <halmd/observables/utility/wavevector.hpp>
#include <test/unit/modules.hpp>
#include <test/tools/init.hpp>

using namespace boost;
using namespace boost::assign;
using namespace halmd;
using namespace halmd::test;
using namespace std;

/**
 * test computation of static structure factor
 *
 * The test is analogous to the one for mdsim/positions/lattice, where the structure
 * factor was computed manually to check the generation of an fcc lattice.
 */
template <int dimension>
void ssf(string const& backend)
{
    typedef fixed_vector<double, dimension> vector_type;

    fixed_vector<unsigned, dimension> ncell =
        (dimension == 3) ? list_of(6)(12)(12) : list_of(4)(1024);
    if (dimension == 3 && backend == "gpu") {
        ncell[0] *= 19; // prime
    }
    unsigned nunit_cell = (dimension == 3) ? 4 : 2;  //< number of particles per unit cell
    unsigned npart = nunit_cell * accumulate(ncell.begin(), ncell.end(), 1, multiplies<unsigned>());
    float density = 0.3;
    float lattice_constant = pow(nunit_cell / density, 1.f / dimension);
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

    BOOST_TEST_MESSAGE("#particles: " << npart << ", #unit cells: " << ncell <<
                       ", lattice constant: " << lattice_constant);

    // init modules
    BOOST_TEST_MESSAGE("initialise modules");
#ifdef WITH_CUDA
    shared_ptr<utility::gpu::device> device = make_device(backend);
#endif /* WITH_CUDA */

    shared_ptr<halmd::random::random> random = make_random(backend, random_file);

    shared_ptr<mdsim::particle<dimension> > particle =
        make_particle<dimension>(backend, npart);

    shared_ptr<mdsim::box<dimension> > box =
        make_box<dimension>(particle, density, static_cast<vector_type>(ncell));

    shared_ptr<mdsim::position<dimension> > position =
        make_lattice(backend, particle, box, random);

    // construct phase space sampler
    shared_ptr<observables::host::samples::phase_space<dimension, double> > sample_host;
#ifdef WITH_CUDA
    shared_ptr<observables::gpu::samples::phase_space<dimension, float> > sample_gpu;
#endif
    shared_ptr<observables::phase_space<dimension> > phase_space;
    if (backend == "host") {
        sample_host = make_shared<observables::host::samples::phase_space<dimension, double> >(
            particle->ntypes
        );
        phase_space = make_phase_space_host(sample_host, particle, box);
    }
#ifdef WITH_CUDA
    else if (backend == "gpu") {
        sample_gpu = make_shared<observables::gpu::samples::phase_space<dimension, float> >(
            particle->ntypes
        );
        phase_space = make_phase_space_gpu(sample_gpu, particle, box);
    }
#endif

    // list of reference results: (wavenumber, ssf, count)
    // S_q = {N if h,k,l all even or odd; 0 if h,k,l mixed parity}
    // see e.g., http://en.wikipedia.org/wiki/Structure_factor
    //
    // entries must be ordered with ascending wavenumber
    //
    // FIXME for boxes with a high aspect ratio. some of the wavevectors only fit approximately,
    // reducing the tolerance on the wavevector magnitude, however, would discard many valid ones as well
    vector<tuple<double, double, unsigned> > ssf_ref;
    double q_lat = 2 * M_PI / lattice_constant;
    if (dimension == 2) {
        ssf_ref.push_back(make_tuple(q_lat / ncell[1], 0, 1));          // hkl = (0, 1/1024)
        ssf_ref.push_back(make_tuple(q_lat / ncell[0], 0, 2));          // hkl = (1/4, 0), (0, 256/1024)
        ssf_ref.push_back(make_tuple(sqrt(  1.) * q_lat, 0, 3));         // hkl = (0,1), (1,0), approx.: (4, 1) × (1/4, 1/1024)
        ssf_ref.push_back(make_tuple(sqrt(  4.) * q_lat, npart / 2., 4));// hkl = (0,2), (2,0), approx.: (4, 1), (8, 1)
        ssf_ref.push_back(make_tuple(sqrt(  9.) * q_lat, 0, 4));         // hkl = (0,3), (3,0), approx.: (3, 1), (4, 1)
        ssf_ref.push_back(make_tuple(sqrt( 16.) * q_lat, npart / 2., 4));// hkl = (0,4), (4,0), approx.: (4, 1), (8, 1)
        ssf_ref.push_back(make_tuple(sqrt(256.) * q_lat, npart / 2., 4));// hkl = (0,16), (16,0), approx.: (4, 1), (8, 1)
    }
    else if (dimension == 3) {
        if (backend == "host") {
            ssf_ref.push_back(make_tuple(q_lat / ncell[1], 0, 2));         // hkl = (0, 0, 1/12), (0, 1/12, 0)
            ssf_ref.push_back(make_tuple(q_lat / ncell[0], 0, 3));         // hkl = (1/6, 0, 0), (0, 2/12, 0), (0, 0, 2/12)
        }
        else if (backend == "gpu") {
            ssf_ref.push_back(make_tuple(q_lat / ncell[0], 0, 1));         // hkl = (1/(6*19), 0, 0)
            ssf_ref.push_back(make_tuple(q_lat / ncell[1], 0, 2));         // hkl = (0, 1/12, 0), (0, 0, 1/12)
        }
        ssf_ref.push_back(make_tuple(sqrt( 1.) * q_lat, 0, 6));         // hkl = (0,0,1), (1/3, 2/3, 2/3) and permutations
        ssf_ref.push_back(make_tuple(sqrt( 2.) * q_lat, 0, 6));         // hkl = (0,1,1), (1/3, 1/3, 4/3), ...
        ssf_ref.push_back(make_tuple(sqrt( 3.) * q_lat, npart / 4., 4));// hkl = (1,1,1), (1/3, 1/3, 5/3), ..., only the 1st one contributes
        ssf_ref.push_back(make_tuple(sqrt( 4.) * q_lat, npart / 2., 6));// hkl = (0,0,2), (2/3, 4/3, 4/3), ...
        ssf_ref.push_back(make_tuple(sqrt( 5.) * q_lat, 0, 6));         // hkl = (0,1,2), ...
        ssf_ref.push_back(make_tuple(sqrt( 6.) * q_lat, 0, 6));         // hkl = (1,1,2), (1/3, 2/3, 7/3), ...
        ssf_ref.push_back(make_tuple(sqrt( 8.) * q_lat, npart / 2., 6));// hkl = (0,2,2), (2/3, 2/3, 8/3), ...
        ssf_ref.push_back(make_tuple(sqrt( 9.) * q_lat, 0, 6));         // hkl = (1,2,2), (0,0,3), ...
        ssf_ref.push_back(make_tuple(sqrt(10.) * q_lat, 0, 6));         // hkl = (0,1,3), ...
        ssf_ref.push_back(make_tuple(sqrt(12.) * q_lat, npart / 4., 4));// hkl = (2,2,2), (2, 2, 10), ...
    }

    // setup wavevectors
    vector<double> wavenumber(ssf_ref.size());
    double const& (*get0)(tuple<double, double, unsigned>::inherited const&) = &boost::get<0>;
    transform(ssf_ref.begin(), ssf_ref.end(), wavenumber.begin(), bind(get0, _1));

    typedef observables::utility::wavevector<dimension> wavevector_type;
    shared_ptr<wavevector_type> wavevector =
        make_shared<wavevector_type>(wavenumber, box->length(), 1e-6, 2 * dimension); // FIXME tolerance, see above

    // construct modules for density modes and static structure factor
    shared_ptr<observables::density_mode<dimension> > density_mode =
        make_density_mode(backend, phase_space, wavevector);

    typedef observables::ssf<dimension> ssf_type;
    shared_ptr<ssf_type> ssf = make_shared<ssf_type>(density_mode, particle->nbox);

    // generate lattices
    LOG_DEBUG("set particle tags");
    particle->set();
    BOOST_TEST_MESSAGE("generate fcc lattice");
    position->set();

    // acquire phase space sample
    BOOST_TEST_MESSAGE("acquire phase space sample");
    phase_space->acquire(0);

    // compute density modes
    BOOST_TEST_MESSAGE("compute density modes");
    density_mode->acquire(0);

    // compute static structure factor
    BOOST_TEST_MESSAGE("compute static structure factor");
    ssf->sample(0);
    vector<typename ssf_type::result_type> const& result = ssf->value()[0]; // particle type 0
    BOOST_CHECK(result.size() == ssf_ref.size());

    // compare with reference values
    double eps = (backend == "gpu") ? 100 * numeric_limits<float>::epsilon() : numeric_limits<double>::epsilon();
    for (unsigned i = 0; i < result.size(); ++i) {
        // check wavenumber
        double q = ssf->wavevector().wavenumber()[i];
        BOOST_CHECK_CLOSE_FRACTION(q, get<0>(ssf_ref[i]), 1e-6);

#ifndef NDEBUG
        // convert boost::array to fixed_vector for output
        fixed_vector<double, 3> S_q;
        copy(result[i].begin(), result[i].end(), S_q.begin());
        LOG_TRACE(
            "S(q = " <<  q << "): " << S_q
         << "; expected: S_q = " << get<1>(ssf_ref[i]) << " (count: " << get<2>(ssf_ref[i]) << ")"
        );
#endif

        // check structure factor
        BOOST_CHECK_SMALL(fabs(result[i][0] - get<1>(ssf_ref[i])) / npart, eps);
        // check error estimate on structure factor
        // from N contributions assume n times the value npart, zero otherwise
        // then variance(S_q) = npart² × (n/N) × (1-n/N) and
        // S_q_err = sqrt(variance(S_q) / (N-1))
        unsigned count = get<2>(ssf_ref[i]);
        if (count > 1) {
            double S_q_err = sqrt(get<1>(ssf_ref[i]) * (npart - get<1>(ssf_ref[i])) / (count - 1));
            BOOST_CHECK_SMALL(fabs(result[i][1] - S_q_err) / max(S_q_err, 1.), eps);
        }
        else {
            BOOST_CHECK(result[i][1] == 0);
        }
        // check count, i.e., number of wavevectors
        BOOST_CHECK(result[i][2] == count);
    }
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

    test_suite* ts1 = BOOST_TEST_SUITE( "ssf" );

    test_suite* ts11 = BOOST_TEST_SUITE( "host" );

    test_suite* ts111 = BOOST_TEST_SUITE( "2d" );
    ts111->add( BOOST_PARAM_TEST_CASE( &ssf<2>, backend.begin(), backend.begin() + 1 ) );

    test_suite* ts112 = BOOST_TEST_SUITE( "3d" );
    ts112->add( BOOST_PARAM_TEST_CASE( &ssf<3>, backend.begin(), backend.begin() + 1 ) );

    ts11->add( ts111 );
    ts11->add( ts112 );
    ts1->add( ts11 );

#ifdef WITH_CUDA
    test_suite* ts12 = BOOST_TEST_SUITE( "gpu" );

    test_suite* ts121 = BOOST_TEST_SUITE( "2d" );
    ts121->add( BOOST_PARAM_TEST_CASE( &ssf<2>, backend.begin() + 1, backend.end() ) );

    test_suite* ts122 = BOOST_TEST_SUITE( "3d" );
    ts122->add( BOOST_PARAM_TEST_CASE( &ssf<3>, backend.begin() + 1, backend.end() ) );

    ts12->add( ts121 );
    ts12->add( ts122 );
    ts1->add( ts12 );
#endif

    master_test_suite().add( ts1 );
}
