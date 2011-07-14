/*
 * Copyright © 2010-2011  Felix Höfling and Peter Colberg
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

#define BOOST_TEST_MODULE lattice
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <cmath>
#include <limits>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/observables/host/phase_space.hpp>
#include <halmd/random/host/random.hpp>
#ifdef WITH_CUDA
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/positions/lattice.hpp>
# include <halmd/observables/gpu/phase_space.hpp>
# include <halmd/random/gpu/random.hpp>
# include <halmd/utility/gpu/device.hpp>
#endif

using namespace boost;
using namespace boost::assign; // list_of
using namespace halmd;
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

template <typename modules_type>
struct lattice
{
    typedef typename modules_type::box_type box_type;
    typedef typename modules_type::particle_type particle_type;
    typedef typename modules_type::position_type position_type;
    typedef typename modules_type::random_type random_type;
    typedef typename modules_type::sample_type sample_type;
    typedef typename modules_type::phase_space_type phase_space_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename vector_type::value_type float_type;
    static unsigned int const dimension = vector_type::static_size;
    static bool const gpu = modules_type::gpu;

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

    void test();
    lattice();
};

template <typename modules_type>
void lattice<modules_type>::test()
{
    BOOST_TEST_MESSAGE("#particles: " << npart << ", #unit cells: " << ncell <<
                       ", lattice constant: " << lattice_constant << ", slab extents: " << slab);

    // generate lattices
    BOOST_TEST_MESSAGE("set particle tags");
    particle->set();
    BOOST_TEST_MESSAGE("generate fcc lattice");
    position->set();

    // acquire phase space samples
    BOOST_TEST_MESSAGE("acquire phase space sample");
    phase_space->acquire(0);

    // compute static structure factors for a set of wavenumbers
    // which are points of the reciprocal lattice
    BOOST_TEST_MESSAGE("compute static structure factors from particle positions");
    vector<fixed_vector<double, dimension> > q;
    double qlat = 2 * M_PI / lattice_constant;
    if (dimension == 3) {
        for (unsigned i = 7; i > 0; --i) {
            fixed_vector<double, dimension> q_ = list_of((i >> 2) & 1)((i >> 1) & 1)(i & 1);
            q.push_back(qlat * q_);
        }
    }
    else if (dimension == 2) {
        for (unsigned i = 3; i > 0; --i) {
            fixed_vector<double, dimension> q_ = list_of((i >> 1) & 1)(i & 1);
            q.push_back(qlat * q_);
        }
    }
    q.push_back(.5 * q[0]);
    q.push_back(2 * q[0]);

    // compute structure factor
    vector<double> ssf = compute_ssf(sample, q);
    // centre of mass
    fixed_vector<double, dimension> r_cm(
        accumulate(
            sample->r[0]->begin(), sample->r[0]->end(), vector_type(0)
          , plus<vector_type>()
        ) / npart
    );
    // minimal and maximal coordinates
    using namespace halmd::detail::numeric::blas;
    fixed_vector<double, dimension> r_min(
        accumulate(
            sample->r[0]->begin(), sample->r[0]->end(), vector_type(0)
          , bind(element_min<float_type, dimension>, _1, _2)
        )
    );
    fixed_vector<double, dimension> r_max(
        accumulate(
            sample->r[0]->begin(), sample->r[0]->end(), vector_type(0)
          , bind(element_max<float_type, dimension>, _1, _2)
        )
    );

    double eps = numeric_limits<float>::epsilon();
    BOOST_CHECK_CLOSE_FRACTION(ssf.front(), npart, eps);
    BOOST_CHECK_CLOSE_FRACTION(ssf.back(), npart, eps);
    for (unsigned i = 1; i < ssf.size() - 1; ++i) {
        BOOST_CHECK_SMALL(ssf[i] / npart, eps);
    }

    // check centre and corners
    fixed_vector<double, dimension> corner = .5 * element_prod(box->length(), slab);  //< upper right corner
    fixed_vector<double, dimension> offset = lattice_constant;                               //< diagonal of the unit cell
    BOOST_CHECK_SMALL(norm_1(r_cm + offset / 4) / norm_1(corner), 2 * eps);
    BOOST_CHECK_SMALL(norm_1(r_min + corner) / norm_1(corner), eps);
    BOOST_CHECK_SMALL(norm_1(r_max - corner + offset / 2) / norm_1(corner), eps);
}

template <typename modules_type>
lattice<modules_type>::lattice()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");

    ncell = (dimension == 3) ? list_of(3)(6)(6) : list_of(4)(1024);
    nunit_cell = (dimension == 3) ? 4 : 2;  //< number of particles per unit cell
    npart = nunit_cell * accumulate(ncell.begin(), ncell.end(), 1, multiplies<unsigned int>());
    density = 0.3;
    lattice_constant = pow(nunit_cell / density, 1.f / dimension);

    slab = (dimension == 3) ? list_of(1.)(.5)(1.) : list_of(1.)(1.);
    double slab_vol_frac = accumulate(slab.begin(), slab.end(), 1., multiplies<double>());
    npart *= slab_vol_frac;
    // adjust density to make sure that the slab can accomodate an fcc lattice with the
    // same lattice spacing (a mismatch is a likely reason for failure of the test)
    density *= slab_vol_frac;

    vector<unsigned int> npart_vector = list_of(npart);

    particle = make_shared<particle_type>(npart_vector);
    box = make_shared<box_type>(npart, density, fixed_vector<double, dimension>(ncell));
    random = make_shared<random_type>();
    position = make_shared<position_type>(particle, box, random, slab);
    sample = make_shared<sample_type>(particle->ntypes);
    phase_space = make_shared<phase_space_type>(sample, particle, box);
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
    static bool const gpu = false;
};

BOOST_AUTO_TEST_CASE( lattice_host_2d ) {
    lattice<host_modules<2, double> >().test();
}
BOOST_AUTO_TEST_CASE( lattice_host_3d ) {
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
    typedef observables::host::samples::phase_space<dimension, float_type> sample_type;
    typedef observables::gpu::phase_space<sample_type> phase_space_type;
    static bool const gpu = true;
};

BOOST_FIXTURE_TEST_CASE( lattice_gpu_2d, device ) {
    lattice<gpu_modules<2, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( lattice_gpu_3d, device ) {
    lattice<gpu_modules<3, float> >().test();
}
#endif // WITH_CUDA
