/*
 * Copyright © 2011  Felix Höfling and Peter Colberg
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

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <cmath>
#include <functional> // std::multiplies
#include <limits>
#include <numeric> // std::accumulate

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/density_mode.hpp>
#include <halmd/observables/host/phase_space.hpp>
#include <halmd/observables/ssf.hpp>
#include <halmd/observables/utility/wavevector.hpp>
#include <halmd/random/host/random.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/positions/lattice.hpp>
# include <halmd/observables/gpu/density_mode.hpp>
# include <halmd/observables/gpu/phase_space.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif

using namespace boost;
using namespace boost::assign; // list_of
using namespace halmd;
using namespace std;

/**
 * test computation of static structure factor
 *
 * The test is analogous to the one for mdsim/positions/lattice, where the structure
 * factor was computed manually to check the generation of an fcc lattice.
 */

template <typename modules_type>
struct lattice
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::position_type position_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::sample_type sample_type;
    typedef typename modules_type::phase_space_type phase_space_type;
    typedef typename modules_type::density_mode_type density_mode_type;
    static bool const gpu = modules_type::gpu;
    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    static unsigned int const dimension = vector_type::static_size;
    typedef observables::utility::wavevector<dimension> wavevector_type;
    typedef observables::ssf<dimension> ssf_type;
    typedef mdsim::clock clock_type;

    fixed_vector<unsigned, dimension> ncell;
    unsigned nunit_cell;
    unsigned npart;
    float density;
    float lattice_constant;
    fixed_vector<double, dimension> slab;

    shared_ptr<box_type> box;
    shared_ptr<particle_type> particle;
    shared_ptr<position_type> position;
    shared_ptr<random_type> random;
    shared_ptr<sample_type> sample;
    shared_ptr<phase_space_type> phase_space;
    shared_ptr<wavevector_type> wavevector;
    shared_ptr<density_mode_type> density_mode;
    shared_ptr<ssf_type> ssf;
    shared_ptr<clock_type> clock;

    void test();
    lattice();
};

template <typename modules_type>
void lattice<modules_type>::test()
{
    BOOST_TEST_MESSAGE("#particles: " << npart << ", #unit cells: " << ncell <<
                       ", lattice constant: " << lattice_constant);

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
        if (!gpu) {
            ssf_ref.push_back(make_tuple(q_lat / ncell[1], 0, 2));         // hkl = (0, 0, 1/12), (0, 1/12, 0)
            ssf_ref.push_back(make_tuple(q_lat / ncell[0], 0, 3));         // hkl = (1/6, 0, 0), (0, 2/12, 0), (0, 0, 2/12)
        }
        else {
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

    wavevector = make_shared<wavevector_type>(wavenumber, box->length(), 1e-6, 2 * dimension); // FIXME tolerance, see above

    // construct modules for density modes and static structure factor
    density_mode = make_shared<density_mode_type>(phase_space, wavevector, clock);
    ssf = make_shared<ssf_type>(density_mode, clock, particle->nbox);

    // generate lattices
    BOOST_TEST_MESSAGE("set particle tags");
    particle->set();
    BOOST_TEST_MESSAGE("generate fcc lattice");
    position->set();

    // acquire phase space sample
    BOOST_TEST_MESSAGE("acquire phase space sample");
    phase_space->acquire();

    // compute density modes
    BOOST_TEST_MESSAGE("compute density modes");
    density_mode->acquire();

    // compute static structure factor
    BOOST_TEST_MESSAGE("compute static structure factor");
    ssf->sample();
    vector<typename ssf_type::result_type> const& result = ssf->value()[0]; // particle type 0
    BOOST_CHECK(result.size() == ssf_ref.size());

    // compare with reference values
    double eps = gpu ? 100 * numeric_limits<float>::epsilon() : numeric_limits<double>::epsilon();
    for (unsigned i = 0; i < result.size(); ++i) {
        // check wavenumber
        double q = ssf->wavevector().wavenumber()[i];
        BOOST_CHECK_CLOSE_FRACTION(q, get<0>(ssf_ref[i]), 1e-6);

#ifndef NDEBUG
        // convert boost::array to fixed_vector for output
        fixed_vector<double, 3> S_q;
        copy(result[i].begin(), result[i].end(), S_q.begin());
        BOOST_TEST_MESSAGE(
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

template <typename modules_type>
lattice<modules_type>::lattice()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    ncell = (dimension == 3) ? list_of(6)(12)(12) : list_of(4)(1024);
    if (dimension == 3 && gpu) {
        ncell[0] *= 19; // prime
    }
    nunit_cell = (dimension == 3) ? 4 : 2;  //< number of particles per unit cell
    npart = nunit_cell * accumulate(ncell.begin(), ncell.end(), 1, multiplies<unsigned>());
    density = 0.3;
    lattice_constant = pow(nunit_cell / density, 1.f / dimension);
    slab = 1;

    vector<unsigned int> npart_vector = list_of(npart);

    particle = make_shared<particle_type>(npart_vector);
    box = make_shared<box_type>(npart, density, fixed_vector<double, dimension>(ncell));
    random = make_shared<random_type>();
    position = make_shared<position_type>(particle, box, random, slab);
    sample = make_shared<sample_type>(particle->ntypes);
    clock = make_shared<clock_type>(0); // bogus time-step
    phase_space = make_shared<phase_space_type>(sample, particle, box, clock);
}

template <int dimension, typename float_type>
struct host_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::positions::lattice<dimension, float_type> position_type;
    typedef halmd::random::host::random random_type;
    typedef observables::host::samples::phase_space<dimension, float_type> sample_type;
    typedef observables::host::phase_space<dimension, float_type> phase_space_type;
    typedef observables::host::density_mode<dimension, float_type> density_mode_type;
    static bool const gpu = false;
};

BOOST_AUTO_TEST_CASE( ssf_host_2d ) {
    lattice<host_modules<2, double> >().test();
}
BOOST_AUTO_TEST_CASE( ssf_host_3d ) {
    lattice<host_modules<3, double> >().test();
}

#ifdef WITH_CUDA
template <int dimension, typename float_type>
struct gpu_modules
{
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::positions::lattice<dimension, float_type, halmd::random::gpu::rand48> position_type;
    typedef halmd::random::gpu::random<halmd::random::gpu::rand48> random_type;
    typedef observables::gpu::samples::phase_space<dimension, float_type> sample_type;
    typedef observables::gpu::phase_space<sample_type> phase_space_type;
    typedef observables::gpu::density_mode<dimension, float_type> density_mode_type;
    static bool const gpu = true;
};

BOOST_FIXTURE_TEST_CASE( ssf_gpu_2d, device ) {
    lattice<gpu_modules<2, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( ssf_gpu_3d, device ) {
    lattice<gpu_modules<3, float> >().test();
}
#endif // WITH_CUDA
