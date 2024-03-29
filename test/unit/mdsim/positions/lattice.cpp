/*
 * Copyright © 2010-2011  Felix Höfling and Peter Colberg
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

#define BOOST_TEST_MODULE lattice
#include <boost/test/unit_test.hpp>

#include <boost/bind/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/particle_groups/all.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/observables/host/phase_space.hpp>
#ifdef HALMD_WITH_GPU
# include <halmd/mdsim/gpu/particle.hpp>
# include <halmd/mdsim/gpu/particle_groups/all.hpp>
# include <halmd/mdsim/gpu/positions/lattice.hpp>
# include <halmd/observables/gpu/phase_space.hpp>
# include <halmd/utility/gpu/device.hpp>
# include <test/tools/cuda.hpp>
#endif
#include <test/tools/ctest.hpp>

using namespace boost;
using namespace halmd;
using namespace std;

/**
 * test initialisation of particle positions: lattice, ...
 */

/** compute static structure factor of phase space sample for some wavevectors */
template <typename position_sample_type, typename vector_type>
vector<double> compute_ssf(
    std::shared_ptr<position_sample_type> position_sample
  , vector<vector_type> const& wavevector
)
{
    unsigned nq = wavevector.size();
    vector<double> ssf(nq);
    vector<double> cos_(nq, 0);
    vector<double> sin_(nq, 0);
    size_t npart = position_sample->data().size();
    BOOST_FOREACH(typename position_sample_type::data_type const& r, position_sample->data()) {
        for (unsigned j = 0; j < nq; ++j) {
            double qr = inner_prod(wavevector[j], static_cast<vector_type>(r));
            cos_[j] += cos(qr);
            sin_[j] += sin(qr);
        }
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
    typedef typename modules_type::particle_group_type particle_group_type;
    typedef typename modules_type::position_type position_type;
    typedef typename modules_type::position_sample_type position_sample_type;
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
    typename modules_type::slab_type slab;

    std::shared_ptr<box_type> box;
    std::shared_ptr<particle_type> particle;
    std::shared_ptr<position_type> position;
    std::shared_ptr<phase_space_type> phase_space;

    void test();
    lattice();
};

template <typename vector_type>
static vector_type wrap_element_max(vector_type const& v, vector_type const& w)
{
    return element_max(v, w);
}

template <typename vector_type>
static vector_type wrap_element_min(vector_type const& v, vector_type const& w)
{
    return element_min(v, w);
}

template <typename modules_type>
void lattice<modules_type>::test()
{
    using namespace boost::placeholders;

    typedef typename modules_type::vector_type vector_type;
    typedef typename modules_type::slab_type slab_type;
    BOOST_TEST_MESSAGE("#particles: " << npart << ", #unit cells: " << ncell <<
                       ", lattice constant: " << lattice_constant << ", slab extents: " << slab);

    // generate lattices
    BOOST_TEST_MESSAGE("generate fcc lattice");
    position->set();

    // acquire phase space samples
    BOOST_TEST_MESSAGE("acquire phase space sample");
    auto position_sample = phase_space->template acquire<position_sample_type>("position");

    // compute static structure factors for a set of wavenumbers
    // which are points of the reciprocal lattice
    BOOST_TEST_MESSAGE("compute static structure factors from particle positions");
    vector<vector_type> q;
    double qlat = 2 * M_PI / lattice_constant;
    if (dimension == 3) {
        for (unsigned i = 7; i > 0; --i) {
            vector_type q_{(i >> 2) & 1, (i >> 1) & 1, i & 1};
            q.push_back(qlat * q_);
        }
    }
    else if (dimension == 2) {
        for (unsigned i = 3; i > 0; --i) {
            vector_type q_{(i >> 1) & 1, i & 1};
            q.push_back(qlat * q_);
        }
    }
    q.push_back(.5 * q[0]);
    q.push_back(2 * q[0]);

    // compute structure factor
    vector<double> ssf = compute_ssf(position_sample, q);
    // centre of mass
    slab_type r_cm(
        accumulate(
            position_sample->data().begin(), position_sample->data().end(), vector_type(0)
          , plus<vector_type>()
        ) / npart
    );
    // minimal and maximal coordinates
    slab_type r_min(
        accumulate(
            position_sample->data().begin(), position_sample->data().end(), vector_type(0)
          , bind(wrap_element_min<vector_type>, _1, _2)
        )
    );
    slab_type r_max(
        accumulate(
            position_sample->data().begin(), position_sample->data().end(), vector_type(0)
          , bind(wrap_element_max<vector_type>, _1, _2)
        )
    );

    double eps = numeric_limits<float>::epsilon();
    BOOST_CHECK_CLOSE_FRACTION(ssf.front(), npart, eps);
    BOOST_CHECK_CLOSE_FRACTION(ssf.back(), npart, eps);
    for (unsigned i = 1; i < ssf.size() - 1; ++i) {
        BOOST_CHECK_SMALL(ssf[i] / npart, eps);
    }

    // check centre and corners
    slab_type corner = .5 * element_prod(static_cast<slab_type>(box->length()), slab);  //< upper right corner
    slab_type offset = lattice_constant;                               //< diagonal of the unit cell
    BOOST_CHECK_SMALL(norm_1(r_cm) / norm_1(corner), 2 * static_cast<typename slab_type::value_type>(eps));
    BOOST_CHECK_SMALL(norm_1(r_min + corner - offset / 4) / norm_1(corner), static_cast<typename slab_type::value_type>(eps));
    BOOST_CHECK_SMALL(norm_1(r_max - corner + offset / 4) / norm_1(corner), static_cast<typename slab_type::value_type>(eps));
}

template <typename modules_type>
lattice<modules_type>::lattice()
{
    BOOST_TEST_MESSAGE("initialise simulation modules");
    typedef fixed_vector<unsigned, dimension> cell_vector;
    typedef typename modules_type::slab_type slab_type;

    ncell = (dimension == 3) ? cell_vector{3, 6, 6} : cell_vector{4, 1024};
    nunit_cell = (dimension == 3) ? 4 : 2;  //< number of particles per unit cell
    npart = nunit_cell * accumulate(ncell.begin(), ncell.end(), 1, multiplies<unsigned int>());
    density = 0.3;
    lattice_constant = pow(nunit_cell / density, 1.f / dimension);
    typename box_type::vector_type box_ratios(ncell);
    boost::numeric::ublas::diagonal_matrix<typename box_type::matrix_type::value_type> edges(dimension);
    for (unsigned int i = 0; i < dimension; ++i) {
        edges(i, i) = lattice_constant * box_ratios[i];
    }

    slab = (dimension == 3) ? slab_type{1., .5, 1.} : slab_type{1., 1.};
    double slab_vol_frac = accumulate(slab.begin(), slab.end(), 1., multiplies<double>());
    // adjust density to make sure that the slab can accomodate an fcc lattice with the
    // same lattice spacing (a mismatch is a likely reason for failure of the test)
    npart *= slab_vol_frac;

    particle = std::make_shared<particle_type>(npart, 1);
    box = std::make_shared<box_type>(edges);
    position = std::make_shared<position_type>(particle, box, slab);
    std::shared_ptr<particle_group_type> particle_group = std::make_shared<particle_group_type>(particle);
    phase_space = std::make_shared<phase_space_type>(particle, particle_group, box);
}

template <int dimension, typename float_type>
struct host_modules
{
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef fixed_vector<float_type, dimension> slab_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::particle_groups::all<particle_type> particle_group_type;
    typedef mdsim::host::positions::lattice<dimension, float_type> position_type;
    typedef observables::host::samples::sample<dimension, float_type> position_sample_type;
    typedef observables::host::phase_space<dimension, float_type> phase_space_type;
    static bool const gpu = false;
};


#ifndef USE_HOST_SINGLE_PRECISION
BOOST_AUTO_TEST_CASE( lattice_host_2d ) {
    lattice<host_modules<2, double> >().test();
}
BOOST_AUTO_TEST_CASE( lattice_host_3d ) {
    lattice<host_modules<3, double> >().test();
}
#else
BOOST_AUTO_TEST_CASE( lattice_host_2d ) {
    lattice<host_modules<2, float> >().test();
}
BOOST_AUTO_TEST_CASE( lattice_host_3d ) {
    lattice<host_modules<3, float> >().test();
}
#endif

#ifdef HALMD_WITH_GPU
template <int dimension, typename float_type>
struct gpu_modules
{
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef fixed_vector<double, dimension> slab_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::particle_groups::all<particle_type> particle_group_type;
    typedef mdsim::gpu::positions::lattice<dimension, float_type> position_type;
    typedef observables::host::samples::sample<dimension, float> position_sample_type;
    typedef observables::gpu::phase_space<dimension, float_type> phase_space_type;
    static bool const gpu = true;
};

# ifdef USE_GPU_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE( lattice_gpu_float_2d, set_cuda_device ) {
    lattice<gpu_modules<2, float> >().test();
}
BOOST_FIXTURE_TEST_CASE( lattice_gpu_float_3d, set_cuda_device ) {
    lattice<gpu_modules<3, float> >().test();
}
# endif
# ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
BOOST_FIXTURE_TEST_CASE( lattice_gpu_dsfloat_2d, set_cuda_device ) {
    lattice<gpu_modules<2, dsfloat> >().test();
}
BOOST_FIXTURE_TEST_CASE( lattice_gpu_dsfloat_3d, set_cuda_device ) {
    lattice<gpu_modules<3, dsfloat> >().test();
}
# endif
#endif // HALMD_WITH_GPU
