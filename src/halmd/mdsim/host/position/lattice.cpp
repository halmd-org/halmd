/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <halmd/mdsim/host/position/lattice.hpp>
#include <halmd/util/logger.hpp>

namespace halmd
{
namespace mdsim { namespace host { namespace position
{

using namespace boost;
using namespace numeric::host::blas;
using namespace std;

template <int dimension, typename float_type>
lattice<dimension, float_type>::lattice(options const& vm)
  : _Base(vm)
  // dependency injection
  , particle(module<particle_type>::fetch(vm))
  , box(module<box_type>::fetch(vm))
  , random(module<random_type>::fetch(vm))
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
    // assign particle types in random order
    for (size_t i = 0, j = 0, k = particle->ntypes[j]; j < particle->ntype; ++j, k += particle->ntypes[j]) {
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
    fill(particle->image.begin(), particle->image.end(), 0);

    LOG("placed particles on fcc lattice: a = " << a);
}

template <int dimension, typename float_type>
typename lattice<dimension, float_type>::pointer
lattice<dimension, float_type>::create(options const& vm)
{
    if (module<particle_type>::fetch(vm)) {
        return pointer(new lattice<dimension, float_type>(vm));
    }
    return pointer();
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class lattice<3, double>;
template class lattice<2, double>;
#else
template class lattice<3, float>;
template class lattice<2, float>;
#endif

}}} // namespace mdsim::host::position

#ifndef USE_HOST_SINGLE_PRECISION
template class module<mdsim::host::position::lattice<3, double> >;
template class module<mdsim::host::position::lattice<2, double> >;
#else
template class module<mdsim::host::position::lattice<3, float> >;
template class module<mdsim::host::position::lattice<2, float> >;
#endif

} // namespace halmd
