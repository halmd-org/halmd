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
#include <boost/foreach.hpp>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <utility>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/trajectory.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/read_integer.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/positions/lattice.hpp>
# include <halmd/observables/gpu/trajectory.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif

using namespace boost;
using namespace boost::assign;
using namespace halmd;
using namespace std;

/**
 * test initialisation of particle positions: lattice, ...
 */

#ifdef WITH_CUDA
shared_ptr<utility::gpu::device> make_device(string const& backend)
{
    if (backend == "gpu") {
        static weak_ptr<utility::gpu::device> device;
        shared_ptr<utility::gpu::device> device_(device.lock());
        if (!device_) {
            device_ = make_shared<utility::gpu::device>(
                vector<int>()   // devices
              , 128             // threads
            );
            device = device_;
        }
        return device_;
    }
    if (backend == "host") {
        return shared_ptr<utility::gpu::device>(); // null pointer
    }
    throw runtime_error("unknown backend: " + backend);
}
#endif /* WITH_CUDA */

template <int dimension>
shared_ptr<mdsim::particle<dimension> > make_particle(
    string const& backend
  , unsigned int npart
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return make_shared<mdsim::gpu::particle<dimension, float> >(
            make_device(backend)
          , vector<unsigned int>(1, npart)
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return make_shared<mdsim::host::particle<dimension, double> >(
            vector<unsigned int>(1, npart)
        );
    }
    throw runtime_error("unknown backend: " + backend);
}

shared_ptr<halmd::random::random> make_random(
    string const& backend
  , unsigned int seed
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
        return make_shared<random_type>(
            make_device(backend)
          , seed
          , 32                  // blocks
          , 32 << DEVICE_SCALE  // threads
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return make_shared<halmd::random::host::random>(
            seed
        );
    }
    throw runtime_error("unknown backend: " + backend);
}

shared_ptr<halmd::random::random> make_random(
    string const& backend
  , string const& filename
)
{
    return make_random(backend, read_integer<unsigned int>(filename));
}

template <int dimension>
shared_ptr<mdsim::position<dimension> > make_lattice(
    string const& backend
  , shared_ptr<mdsim::particle<dimension> > particle
  , shared_ptr<mdsim::box<dimension> > box
  , shared_ptr<halmd::random::random> random
)
{
#ifdef WITH_CUDA
    if (backend == "gpu") {
        return make_shared<mdsim::gpu::positions::lattice<dimension, float, halmd::random::gpu::rand48> >(
            dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
          , dynamic_pointer_cast<halmd::random::gpu::random<halmd::random::gpu::rand48> >(random)
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return make_shared<mdsim::host::positions::lattice<dimension, double> >(
            dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
          , dynamic_pointer_cast<halmd::random::host::random>(random)
        );
    }
    throw runtime_error("unknown backend: " + backend);
}

template <int dimension, typename float_type>
shared_ptr<observables::trajectory<dimension> > make_trajectory_host(
    shared_ptr<observables::host::samples::trajectory<dimension, float_type> > sample
  , shared_ptr<mdsim::particle<dimension> > particle
  , shared_ptr<mdsim::box<dimension> > box
)
{
    return make_shared<observables::host::trajectory<dimension, float_type> >(
        sample
      , dynamic_pointer_cast<mdsim::host::particle<dimension, float_type> >(particle)
      , box
    );
}

template <int dimension, typename float_type>
shared_ptr<observables::trajectory<dimension> > make_trajectory_gpu(
    shared_ptr<observables::host::samples::trajectory<dimension, float_type> > sample
  , shared_ptr<mdsim::particle<dimension> > particle
  , shared_ptr<mdsim::box<dimension> > box
)
{
    return make_shared<observables::gpu::trajectory<observables::host::samples::trajectory<dimension, float_type> > >(
        sample
      , dynamic_pointer_cast<mdsim::gpu::particle<dimension, float_type> >(particle)
      , box
    );
}

/** compute static structure factor of trajectory sample for some wavevectors */
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

template <int dimension>
void lattice(string const& backend)
{
    typedef typename mdsim::type_traits<dimension, double>::vector_type vector_type;
    typedef typename mdsim::type_traits<dimension, float>::vector_type gpu_vector_type;

    fixed_vector<unsigned, dimension> ncell =
        (dimension == 3) ? list_of(3)(3)(6) : list_of(4)(1024);
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
        make_shared<mdsim::box<dimension> >(particle, density, static_cast<vector_type>(ncell));

    shared_ptr<mdsim::position<dimension> > position =
        make_lattice(backend, particle, box, random);

    // generate lattices
    LOG_DEBUG("set particle tags");
    particle->set();
    BOOST_TEST_MESSAGE("generate fcc lattice");
    position->set();

    // acquire trajectory samples
    LOG_DEBUG("acquire trajectory sample");
    shared_ptr<observables::host::samples::trajectory<dimension, double> > sample_host;
    shared_ptr<observables::host::samples::trajectory<dimension, float> > sample_gpu;
    shared_ptr<observables::trajectory<dimension> > trajectory;
    if (backend == "host") {
        sample_host = make_shared<observables::host::samples::trajectory<dimension, double> >(
            particle->ntypes
        );
        trajectory = make_trajectory_host<dimension, double>(sample_host, particle, box);
    }
#ifdef WITH_CUDA
    else if (backend == "gpu") {
        sample_gpu = make_shared<observables::host::samples::trajectory<dimension, float> >(
            particle->ntypes
        );
        trajectory = make_trajectory_gpu<dimension, float>(sample_gpu, particle, box);
    }
#endif
    trajectory->acquire(0.);

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
    if (backend == "host") {
        ssf = compute_ssf(sample_host, q);
    }
    else if (backend == "gpu") {
        ssf = compute_ssf(sample_gpu, q);
    }

    double eps = (double)numeric_limits<float>::epsilon();
    BOOST_CHECK_CLOSE_FRACTION(ssf.front(), npart, eps);
    BOOST_CHECK_CLOSE_FRACTION(ssf.back(), npart, eps);
    for (unsigned i = 1; i < ssf.size() - 1; ++i) {
        BOOST_CHECK_SMALL(ssf[i] / npart, eps);
    }
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
