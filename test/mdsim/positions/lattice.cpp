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

#define BOOST_TEST_MODULE position_lattice
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
#include <halmd/numeric/accumulator.hpp>
#include <test/modules.hpp>
#include <test/tools/init.hpp>

using namespace boost;
using namespace boost::assign;
using namespace halmd;
using namespace halmd::test;
using namespace std;

/**
 * test initialisation of particle positions: lattice, ...
 */

/** compute static structure factor of phase space sample for some wavevectors */
template <typename sample_type, typename vector_type>
vector<double> compute_ssf(
    shared_ptr<sample_type> sample
  , vector<vector_type> const& wavevector
)
{
    unsigned nq = wavevector.size();
    vector<double> ssf(nq);
    vector<double> cos_(nq, 0);
    vector<double> sin_(nq, 0);
    size_t npart = 0;
    for (size_t i = 0; i < sample->r.size(); ++i) {
        BOOST_FOREACH(typename sample_type::vector_type const& r, *sample->r[i]) {
            for (unsigned j = 0; j < nq; ++j) {
                double qr = inner_prod(wavevector[j], static_cast<vector_type>(r));
                cos_[j] += cos(qr);
                sin_[j] += sin(qr);
            }
        }
        npart += sample->r[i]->size();
    }
    for (size_t j = 0; j < nq; ++j) {
        ssf[j] = (cos_[j] * cos_[j] + sin_[j] * sin_[j]) / npart;
    }
    return ssf;
}

/** similar as std::accumulate, but iterate over all particles of the sample */
template <typename sample_type, typename vector_type, typename BinaryFunction>
vector_type accumulate(
    shared_ptr<sample_type> sample
  , vector_type initial_value
  , BinaryFunction fct
)
{
    typedef typename sample_type::sample_vector_ptr sample_vector_ptr;

    vector_type result(initial_value);
    BOOST_FOREACH(sample_vector_ptr const r, sample->r) {
        result = std::accumulate(r->begin(), r->end(), result, fct);
    }
    return result;
}

template <int dimension>
void lattice(string const& backend)
{
    typedef typename mdsim::type_traits<dimension, double>::vector_type vector_type;
    typedef typename mdsim::type_traits<dimension, float>::vector_type gpu_vector_type;

    fixed_vector<unsigned, dimension> ncell =
        (dimension == 3) ? list_of(3)(6)(6) : list_of(4)(1024);
    unsigned nunit_cell = (dimension == 3) ? 4 : 2;  //< number of particles per unit cell
    unsigned npart = nunit_cell * accumulate(ncell.begin(), ncell.end(), 1, multiplies<unsigned>());
    float density = 0.3;
    float lattice_constant = pow(nunit_cell / density, 1.f / dimension);
    char const* random_file = "/dev/urandom";

    vector_type slab = (dimension == 3) ? list_of(1.)(.5)(1.) : list_of(1.)(1.);
    double slab_vol_frac = accumulate(slab.begin(), slab.end(), 1., multiplies<double>());
    npart *= slab_vol_frac;
    // adjust density to make sure that the slab can accomodate an fcc lattice with the
    // same lattice spacing (a mismatch is a likely reason for failure of the test)
    density *= slab_vol_frac;

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
                       ", lattice constant: " << lattice_constant << ", slab extents: " << slab);

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
        make_lattice(backend, particle, box, random, slab);

    // generate lattices
    LOG_DEBUG("set particle tags");
    particle->set();
    BOOST_TEST_MESSAGE("generate fcc lattice");
    position->set();

    // acquire phase space samples
    LOG_DEBUG("acquire phase space sample");
    shared_ptr<observables::host::samples::phase_space<dimension, double> > sample_host;
#ifdef WITH_CUDA
    shared_ptr<observables::host::samples::phase_space<dimension, float> > sample_gpu;
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
        sample_gpu = make_shared<observables::host::samples::phase_space<dimension, float> >(
            particle->ntypes
        );
        phase_space = make_phase_space_gpu(sample_gpu, particle, box);
    }
#endif
    phase_space->acquire(0);

    // compute static structure factors for a set of wavenumbers
    // which are points of the reciprocal lattice
    BOOST_TEST_MESSAGE("compute static structure factors from particle positions");
    vector<vector_type> q;
    double qlat = 2 * M_PI / lattice_constant;
    if (dimension == 3) {
        for (unsigned i = 7; i > 0; --i) {
            vector_type q_ = list_of((i >> 2) & 1)((i >> 1) & 1)(i & 1);
            q.push_back(qlat * q_);
        }
    }
    else if (dimension == 2) {
        for (unsigned i = 3; i > 0; --i) {
            vector_type q_ = list_of((i >> 1) & 1)(i & 1);
            q.push_back(qlat * q_);
        }
    }
    q.push_back(.5 * q[0]);
    q.push_back(2 * q[0]);

    vector<double> ssf;
    vector_type r_cm, r_min, r_max;
    if (backend == "host") {
        // compute structure factor
        ssf = compute_ssf(sample_host, q);
        // centre of mass
        r_cm  = accumulate(sample_host, vector_type(0), plus<vector_type>()) / npart;
        // minimal and maximal coordinates
        using namespace halmd::detail::numeric::blas;
        r_min = accumulate(sample_host, vector_type(0), bind(element_min<double, dimension>, _1, _2));
        r_max = accumulate(sample_host, vector_type(0), bind(element_max<double, dimension>, _1, _2));
    }
#ifdef WITH_CUDA
    else if (backend == "gpu") {
        // compute structure factor
        ssf = compute_ssf(sample_gpu, q);
        // centre of mass
        r_cm  = static_cast<vector_type>(  //< from fixed_vector<float, N> to fixed_vector<double, N>
                    accumulate(sample_gpu, gpu_vector_type(0), plus<gpu_vector_type>())
                ) / npart;
        // minimal and maximal coordinates
        using namespace halmd::detail::numeric::blas;
        r_min = static_cast<vector_type>(
                    accumulate(sample_gpu, gpu_vector_type(0), bind(element_min<float, dimension>, _1, _2))
                );
        r_max = static_cast<vector_type>(
                    accumulate(sample_gpu, gpu_vector_type(0), bind(element_max<float, dimension>, _1, _2))
                );
    }
#endif
    else {
        return;
    }

    double eps = (double)numeric_limits<float>::epsilon();
    BOOST_CHECK_CLOSE_FRACTION(ssf.front(), npart, eps);
    BOOST_CHECK_CLOSE_FRACTION(ssf.back(), npart, eps);
    for (unsigned i = 1; i < ssf.size() - 1; ++i) {
        BOOST_CHECK_SMALL(ssf[i] / npart, eps);
    }

    // check centre and corners
    vector_type corner = .5 * element_prod(box->length(), slab);  //< upper right corner
    vector_type offset(lattice_constant);                               //< diagonal of the unit cell
    BOOST_CHECK_SMALL(norm_1(r_cm + offset / 4) / norm_1(corner), 2 * eps);
    BOOST_CHECK_SMALL(norm_1(r_min + corner) / norm_1(corner), eps);
    BOOST_CHECK_SMALL(norm_1(r_max - corner + offset / 2) / norm_1(corner), eps);
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

    test_suite* ts1 = BOOST_TEST_SUITE( "lattice" );

    test_suite* ts11 = BOOST_TEST_SUITE( "host" );

    test_suite* ts111 = BOOST_TEST_SUITE( "2d" );
    ts111->add( BOOST_PARAM_TEST_CASE( &lattice<2>, backend.begin(), backend.begin() + 1 ) );

    test_suite* ts112 = BOOST_TEST_SUITE( "3d" );
    ts112->add( BOOST_PARAM_TEST_CASE( &lattice<3>, backend.begin(), backend.begin() + 1 ) );

    ts11->add( ts111 );
    ts11->add( ts112 );
    ts1->add( ts11 );

#ifdef WITH_CUDA
    test_suite* ts12 = BOOST_TEST_SUITE( "gpu" );

    test_suite* ts121 = BOOST_TEST_SUITE( "2d" );
    ts121->add( BOOST_PARAM_TEST_CASE( &lattice<2>, backend.begin() + 1, backend.end() ) );

    test_suite* ts122 = BOOST_TEST_SUITE( "3d" );
    ts122->add( BOOST_PARAM_TEST_CASE( &lattice<3>, backend.begin() + 1, backend.end() ) );

    ts12->add( ts121 );
    ts12->add( ts122 );
    ts1->add( ts12 );
#endif

    master_test_suite().add( ts1 );
}
