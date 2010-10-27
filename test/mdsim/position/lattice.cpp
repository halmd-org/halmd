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

#define BOOST_TEST_MODULE position
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <boost/program_options.hpp>
#include <limits>
#include <map>
#include <string>
#include <utility>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/position/lattice.hpp>
#include <halmd/mdsim/host/sampler/trajectory.hpp>
#include <halmd/mdsim/particle.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/options.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/read_integer.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/position/lattice.hpp>
# include <halmd/mdsim/gpu/sampler/trajectory.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif

using namespace boost;
using namespace halmd;
using namespace std;

/**
 * test initialisation of particle positions: lattice, ...
 */

#ifdef WITH_CUDA
shared_ptr<utility::gpu::device> make_device(
    string const& backend
  , vector<int> devices
  , unsigned int threads
)
{
    if (backend == "gpu") {
        static weak_ptr<utility::gpu::device> device;
        shared_ptr<utility::gpu::device> device_(device.lock());
        if (!device_) {
            device_ = make_shared<utility::gpu::device>(
                devices
              , threads
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
            make_device(
                backend
              , utility::gpu::device::default_devices()
              , utility::gpu::device::default_threads()
            )
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
            make_device(
                backend
              , utility::gpu::device::default_devices()
              , utility::gpu::device::default_threads()
            )
          , seed
          , random_type::default_blocks()
          , random_type::default_threads()
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
        return make_shared<mdsim::gpu::position::lattice<dimension, float, halmd::random::gpu::rand48> >(
            dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
          , box
          , dynamic_pointer_cast<halmd::random::gpu::random<halmd::random::gpu::rand48> >(random)
        );
    }
#endif /* WITH_CUDA */
    if (backend == "host") {
        return make_shared<mdsim::host::position::lattice<dimension, double> >(
            dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
          , box
          , dynamic_pointer_cast<halmd::random::host::random>(random)
        );
    }
    throw runtime_error("unknown backend: " + backend);
}

template <int dimension>
shared_ptr<mdsim::samples::host::trajectory<dimension, double> > make_sample_host(
    shared_ptr<mdsim::particle<dimension> > particle
  , shared_ptr<mdsim::box<dimension> > box
)
{
    typedef mdsim::host::sampler::trajectory<dimension, double> sampler_type;
    return make_shared<sampler_type>(
        dynamic_pointer_cast<mdsim::host::particle<dimension, double> >(particle)
      , box
    );
}

template <int dimension>
shared_ptr<mdsim::samples::host::trajectory<dimension, float> > make_sample_gpu(
    shared_ptr<mdsim::particle<dimension> > particle
  , shared_ptr<mdsim::box<dimension> > box
)
{
#ifdef WITH_CUDA
    typedef mdsim::samples::host::trajectory<dimension, float> sample_type;
    typedef mdsim::gpu::sampler::trajectory<sample_type> sampler_type;
    return make_shared<sampler_type>(
        dynamic_pointer_cast<mdsim::gpu::particle<dimension, float> >(particle)
      , box
    );
#endif /* WITH_CUDA */
    throw runtime_error("unknown backend: gpu");
}

/** test GPU and host implementation separately */
template <int dimension>
void lattice(string const& backend)
{
    // FIXME
    BOOST_CHECK( true );
}

/** test whether GPU and host implementation yield the same results */
template <int dimension>
void lattice_both()
{
#ifdef WITH_CUDA
    float density = 0.3;
    unsigned npart = 997;       //< the prime nearest to 1000
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

    BOOST_TEST_MESSAGE("compare host and GPU backend in " <<
                       dimension << " dimensions");

    // init modules
    BOOST_TEST_MESSAGE("initialise modules");
    shared_ptr<utility::gpu::device> device = make_device(
        "gpu"
      , utility::gpu::device::default_devices()
      , utility::gpu::device::default_threads()
    );

    shared_ptr<halmd::random::random> random_host = make_random("host", random_file);
    shared_ptr<halmd::random::random> random_gpu = make_random("gpu", random_file);

    shared_ptr<mdsim::particle<dimension> > particle_host =
        make_particle<dimension>("host", npart);
    shared_ptr<mdsim::particle<dimension> > particle_gpu =
        make_particle<dimension>("gpu", npart);

    shared_ptr<mdsim::box<dimension> > box =
        make_shared<mdsim::box<dimension> >(particle_host, density, 1 /*< cube aspect ratios */);

    shared_ptr<mdsim::position<dimension> > position_host =
        make_lattice("host", particle_host, box, random_host);
    shared_ptr<mdsim::position<dimension> > position_gpu =
        make_lattice("gpu", particle_gpu, box, random_gpu);

    shared_ptr<mdsim::samples::host::trajectory<dimension, double> > sample_host =
        make_sample_host(particle_host, box);
    shared_ptr<mdsim::samples::host::trajectory<dimension, float> > sample_gpu =
        make_sample_gpu(particle_gpu, box);

    // generate lattices
    BOOST_TEST_MESSAGE("set particle tags at host");
    particle_host->set();
    BOOST_TEST_MESSAGE("set particle tags at GPU");
    particle_gpu->set();
    BOOST_TEST_MESSAGE("generate fcc lattice at host");
    position_host->set();
    BOOST_TEST_MESSAGE("generate fcc lattice at GPU");
    position_gpu->set();

    // acquire trajectory samples
    BOOST_TEST_MESSAGE("acquire sample from host");
    sample_host->acquire(0);
    BOOST_TEST_MESSAGE("acquire sample from GPU");
    sample_gpu->acquire(0);

    // compute deviations between particle positions from both implementations
    BOOST_TEST_MESSAGE("compare particle positions");
    double diff = 0;
    for (size_t i = 0; i < particle_host->ntype; i++) {
        typedef fixed_vector<double, dimension> vector_type;
        for (size_t j = 0; j < particle_host->ntypes[j]; j++) {
            vector_type dr((*sample_host->r[i])[j] - static_cast<vector_type>((*sample_gpu->r[i])[j]));
            diff += inner_prod(dr, dr);
        }
    }
    diff /= particle_host->nbox;

    BOOST_CHECK_SMALL(diff, (double)numeric_limits<float>::epsilon());
#endif /* WITH_CUDA */
}

static void __attribute__((constructor)) init_unit_test_suite()
{
    typedef boost::program_options::variable_value variable_value;
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

    test_suite* ts11 = BOOST_TEST_SUITE( "2d" );
    ts11->add( BOOST_PARAM_TEST_CASE( &lattice<2>, backend.begin(), backend.begin() + 1 ) );
#ifdef WITH_CUDA
    ts11->add( BOOST_PARAM_TEST_CASE( &lattice<2>, backend.begin() + 1, backend.end() ) );
    ts11->add( BOOST_TEST_CASE( &lattice_both<2> ) );
#endif

    test_suite* ts12 = BOOST_TEST_SUITE( "3d" );
    ts12->add( BOOST_PARAM_TEST_CASE( &lattice<3>, backend.begin(), backend.begin() + 1 ) );
#ifdef WITH_CUDA
    ts12->add( BOOST_PARAM_TEST_CASE( &lattice<3>, backend.begin() + 1, backend.end() ) );
    ts12->add( BOOST_TEST_CASE( &lattice_both<3> ) );
#endif

    ts1->add( ts11 );
    ts1->add( ts12 );

    master_test_suite().add( ts1 );
}
