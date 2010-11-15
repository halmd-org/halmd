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
#include <cmath>
#include <limits>
#include <numeric>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/position/lattice.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace position
{

template <int dimension, typename float_type>
lattice<dimension, float_type>::lattice(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<random_type> random
)
  // dependency injection
  : particle(particle)
  , box(box)
  , random(random)
{
}

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
    // randomise particle types if there are more than 1
    if (particle->ntypes.size() > 1) {
        LOG("randomly permuting particle types");
        random->shuffle(particle->type.begin(), particle->type.end());
    }

    // assign lattice coordinates
    vector_type L = box->length();
    double u = (dimension == 3) ? 4 : 2;
    double V = box->volume() / ceil(particle->nbox / u);
    double a = pow(V, 1. / dimension);
    fixed_vector<unsigned int, dimension> n(L / a);
    while (particle->nbox > u * accumulate(n.begin(), n.end(), 1, multiplies<unsigned int>())) {
        vector_type t;
        for (size_t i = 0; i < dimension; ++i) {
            t[i] = L[i] / (n[i] + 1);
        }
        typename vector_type::iterator it = max_element(t.begin(), t.end());
        a = *it;
        ++n[it - t.begin()];
    }
    LOG("placing particles on fcc lattice: a = " << a);

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
}

template <int dimension, typename float_type>
void lattice<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("lattice_" + lexical_cast<string>(dimension) + "_");
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("host")
                [
                    namespace_("position")
                    [
                        class_<lattice, shared_ptr<_Base>, _Base>(class_name.c_str())
                            .def(constructor<
                                 shared_ptr<particle_type>
                               , shared_ptr<box_type>
                               , shared_ptr<random_type>
                             >())
                    ]
                ]
            ]
        ]
    ];
}

static __attribute__((constructor)) void register_lua()
{
    lua_wrapper::register_(1) //< distance of derived to base class
#ifndef USE_HOST_SINGLE_PRECISION
    [
        &lattice<3, double>::luaopen
    ]
    [
        &lattice<2, double>::luaopen
    ];
#else
    [
        &lattice<3, float>::luaopen
    ]
    [
        &lattice<2, float>::luaopen
    ];
#endif
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

} // namespace halmd
