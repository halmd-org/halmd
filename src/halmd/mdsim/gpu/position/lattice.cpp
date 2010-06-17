/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <algorithm>
#include <boost/array.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <cmath>
#include <limits>
#include <numeric>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/position/lattice_kernel.cuh>
#include <halmd/mdsim/gpu/position/lattice.hpp>
#include <halmd/util/exception.hpp>
#include <halmd/utility/module.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace position
{

// using namespace boost;
// using namespace numeric::gpu::blas;
using namespace halmd;
using namespace std;

/**
 * Resolve module dependencies
 */
template <int dimension, typename float_type>
void lattice<dimension, float_type>::depends()
{
    modules::required<_Self, particle_type>();
    modules::required<_Self, box_type>();
    modules::required<_Self, random_type>();
}

template <int dimension, typename float_type>
void lattice<dimension, float_type>::select(po::options const& vm)
{
    if (vm["position"].as<string>() != "lattice") {
        throw unsuitable_module("mismatching option position");
    }
}

template <int dimension, typename float_type>
lattice<dimension, float_type>::lattice(po::options const& vm)
  : _Base(vm)
  // dependency injection
  , particle(modules::fetch<particle_type>(vm))
  , box(modules::fetch<box_type>(vm))
  , random(modules::fetch<random_type>(vm))
{}

/**
 * Place particles on a face-centered cubic (fcc) lattice
 *
 * The task is to determine the minimum lattice distance of an fcc lattice
 * that fits into a rectangular parallelepiped (the simulation box) and has
 * equally many or more lattice sites than the number of particles.
 *
 * The number of lattice unit cells is defined as
 *
 *    s = floor(L_x / a) floor(L_y / a) floor(L_z / a) ≥ ceil(N / u)
 *
 *    N       number of particles
 *    L_i     box edge lengths for i ∈ {x, y, z}
 *    a       lattice distance
 *    u       number of particle per unit cell (4 in 3D, 2 in 2D)
 *
 * The upper bound for the lattice distance is given by
 *
 *    a ≤ (L_x L_y L_z / ceil(N / u))^(1/3)
 *
 * This yields lower bounds for the number of unit cells per dimension
 *
 *    n_i ≥ floor(L_i / a)
 *
 * The minimum a is then determined by iteratively increasing any of the
 * n_i to yield the nearest smaller value for the lattice distance until
 *
 *    s * u ≥ N
 *
 * is satisfied.
 */
template <int dimension, typename float_type>
void lattice<dimension, float_type>::set()
{
    // TODO: handle non-cubic boxes
    for (unsigned i=1; i < dimension; i++) {
        assert(box->length()[0] == box->length()[i]);
    }

    LOG("placing particles on face-centered cubic (fcc) lattice");

    // particles per 2- or 3-dimensional unit cell
    unsigned int const m = 2 * (dimension - 1);
    // lower boundary for number of particles per lattice dimension
    unsigned int n = static_cast<unsigned int>(pow(particle->nbox / m, 1. / dimension));
    // lower boundary for total number of lattice sites
    unsigned int N = m * static_cast<unsigned int>(pow(n, dimension));

    if (N < particle->nbox) {
        n += 1;
        N = m * static_cast<unsigned int>(pow(n, dimension));
    }
    if (N > particle->nbox) {
        LOG_WARNING("lattice not fully occupied (" << N << " sites)");
    }

    // minimum distance in 2- or 3-dimensional fcc lattice
    LOG("minimum lattice distance: " << (box->length()[0] / n) / sqrt(2));

#ifdef USE_VERLET_DSFUN
    cuda::memset(particle->g_r, 0, particle->g_r.capacity());
#endif

//     boost::array<high_resolution_timer, 2> timer;
    cuda::thread::synchronize();
    try {
//         timer[0].record();
//         cuda::configure(particle->dim_.grid, particle->dim_.block);
        lattice_wrapper<dimension>::fcc(particle->g_r, n, box->length()[0]);
        cuda::thread::synchronize();
//         timer[1].record();
    }
    catch (cuda::error const& e) {
        LOG_ERROR("CUDA: " << e.what());
        throw exception("failed to compute particle lattice positions on GPU");
    }
//     m_times["lattice"] += timer[1] - timer[0];

    // randomly permute particle coordinates for binary mixture
//     algorith::gpu::permute(particle->g_r);
//     particle->rearrange(particle->g_r);
//     LOG("randomly permuting particle coordinates");
//
//     cuda::vector<unsigned int> g_sort_index(npart);
//     g_sort_index.reserve(dim_.threads());
//
//     try {
//         rng_.get(g_sort_index);
//         radix_sort_(g_sort_index, g_r);
//         cuda::thread::synchronize();
//     }
//     catch (cuda::error const& e) {
//         LOG_ERROR("CUDA: " << e.what());
//         throw exception("failed to randomly permute particle coordinates on GPU");
//     }

    // ??? shift particle positions to range (-L/2, L/2)
//    box->reduce_periodic(r);


//    ??? assign_positions();

/*

    // assign particle types in random order
    for (size_t i = 0, j = 0, k = 0; j < particle->ntype; ++j) {
        k += particle->ntypes[j];
        for (; i < k; ++i) {
            particle->type[i] = j;
        }
    }
    random->shuffle(particle->type.begin(), particle->type.end());

    // assign particle tags
    copy(counting_iterator<size_t>(0), counting_iterator<size_t>(particle->nbox), particle->tag.begin());

    // assign lattice coordinates
    vector_type L = box->length();
    double u = (dimension == 3) ? 4 : 2;
    double V = accumulate(L.begin(), L.end(), 1. / ceil(particle->nbox / u), multiplies<double>());
    double a = pow(V, 1. / dimension);
    vector<unsigned int, dimension> n(L / a);
    while (particle->nbox > u * accumulate(n.begin(), n.end(), 1, multiplies<unsigned int>())) {
        vector_type t;
        for (size_t i = 0; i < dimension; ++i) {
            t[i] = L[i] / (n[i] + 1);
        }
        a = *max_element(t.begin(), t.end());
        for (size_t i = 0; i < dimension; ++i) {
            if (t[i] == a) {
                n[i]++;
            }
        }
    }
    for (size_t i = 0; i < particle->nbox; ++i) {
        vector_type& r = particle->r[i] = a;
        if (dimension == 3) {
            r[0] *= ((i >> 2) % n[0]) + ((i ^ (i >> 1)) & 1) / 2.;
            r[1] *= ((i >> 2) / n[0] % n[1]) + (i & 1) / 2.;
            r[2] *= ((i >> 2) / n[0] / n[1]) + (i & 2) / 4.;
        }
        else {
            r[0] *= ((i >> 1) % n[0]) + (i & 1) / 2.;
            r[1] *= ((i >> 1) / n[0]) + (i & 1) / 2.;
        }
        // shift particle positions to range (-L/2, L/2)
        box->reduce_periodic(r);
    }

    // assign particle image vectors
    fill(particle->image.begin(), particle->image.end(), 0);*/

//     LOG("particles placed on fcc lattice: a = " << a);
}

// explicit instantiation
template class lattice<3, float>;
template class lattice<2, float>;

}}} // namespace mdsim::gpu::position

template class module<mdsim::gpu::position::lattice<3, float> >;
template class module<mdsim::gpu::position::lattice<2, float> >;

} // namespace halmd
