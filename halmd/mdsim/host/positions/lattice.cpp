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
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace host { namespace positions
{

template <int dimension, typename float_type>
lattice<dimension, float_type>::lattice(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<random_type> random
  , vector_type const& slab
)
  // dependency injection
  : particle(particle)
  , box(box)
  , random(random)
  , slab_(slab)
{
    if (*min_element(slab_.begin(), slab_.end()) <= 0 ||
        *max_element(slab_.begin(), slab_.end()) > 1
       ) {
        throw std::logic_error("slab extents must be a fraction between 0 and 1");
    }

    if (*min_element(slab_.begin(), slab_.end()) < 1) {
        LOG("restrict initial particle positions to slab: " << slab_);
    }
}

template <int dimension, typename float_type>
void lattice<dimension, float_type>::set()
{
    assert(particle->r.size() == particle->nbox);

    // assign fcc lattice points to a fraction of the particles in a slab at the centre
    vector_type length = element_prod(box->length(), slab_);
    vector_type offset = -length / 2;
    fcc(particle->r.begin(), particle->r.end(), length, offset);

    // randomise particle positions if there is more than 1 particle type
    // FIXME this requires a subsequent sort
    // FIXME this will fail greatly once we support polymers
    if (particle->ntypes.size() > 1) {
        LOG("randomly permuting particle positions");
        random->shuffle(particle->r.begin(), particle->r.end());
    }

    // assign particle image vectors
    fill(particle->image.begin(), particle->image.end(), 0);
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
template <int dimension, typename float_type> template <typename position_iterator>
void lattice<dimension, float_type>::fcc(
    position_iterator first, position_iterator last
  , vector_type const& length, vector_type const& offset
)
{
    typedef fixed_vector<unsigned int, dimension> index_type;

    LOG_TRACE("generating fcc lattice for " << last - first << " particles, box: " << length << ", offset: " << offset);
    size_t npart = last - first;
    double u = (dimension == 3) ? 4 : 2;
    double V = accumulate(length.begin(), length.end(), 1., multiplies<double>()) / ceil(npart / u);
    double a = pow(V, 1. / dimension);
    index_type n(length / a);
    while (npart > u * accumulate(n.begin(), n.end(), 1, multiplies<unsigned int>())) {
        vector_type t;
        for (size_t i = 0; i < dimension; ++i) {
            t[i] = length[i] / (n[i] + 1);
        }
        typename vector_type::iterator it = max_element(t.begin(), t.end());
        a = *it;
        unsigned int m = n[it - t.begin()];
        n = static_cast<index_type>(length / a); //< recompute to preserve aspect ratios of box
        n[it - t.begin()] = m + 1;               //< ensure increment of at least one component
    }
    LOG("placing particles on fcc lattice: a = " << a);
    LOG_DEBUG("number of fcc unit cells: " << n);

    unsigned int N = static_cast<unsigned int>(
        u * accumulate(n.begin(), n.end(), 1, multiplies<unsigned int>())
    );
    if (N > npart) {
        LOG_WARNING("lattice not fully occupied (" << N << " sites)");
    }

    // insert a vacancy every 'skip' sites
    unsigned int skip = (N - npart) ? static_cast<unsigned int>(ceil(static_cast<double>(N) / (N - npart))) : 0;
    if (skip) {
        LOG_TRACE("insert a vacancy at every " << skip << "th site");
    }

    size_t i = 0;
    for (position_iterator r_it = first; r_it != last; ++r_it, ++i) {
        // skip vacant lattice points
        if (skip && i % skip == skip - 1) {
            ++i;
        }

        vector_type& r = *r_it = a;
        if (dimension == 3) {
            r[0] *= ((i >> 2) % n[0]) + ((i ^ (i >> 1)) & 1) / 2.;
            r[1] *= ((i >> 2) / n[0] % n[1]) + (i & 1) / 2.;
            r[2] *= ((i >> 2) / n[0] / n[1]) + (i & 2) / 4.;
        }
        else {
            r[0] *= ((i >> 1) % n[0]) + (i & 1) / 2.;
            r[1] *= ((i >> 1) / n[0]) + (i & 1) / 2.;
        }
        // shift origin of lattice to offset
        r += offset;
    }
    assert(i <= N);
    LOG_DEBUG("number of particles inserted: " << last - first);
}

template <int dimension, typename float_type>
static char const* module_name_wrapper(lattice<dimension, float_type> const&)
{
    return lattice<dimension, float_type>::module_name();
}

template <int dimension, typename float_type>
void lattice<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name(module_name() + ("_" + lexical_cast<string>(dimension) + "_"));
    module(L)
    [
        namespace_("libhalmd")
        [
            namespace_("mdsim")
            [
                namespace_("host")
                [
                    namespace_("positions")
                    [
                        class_<lattice, shared_ptr<_Base>, _Base>(class_name.c_str())
                            .def(constructor<
                                 shared_ptr<particle_type>
                               , shared_ptr<box_type>
                               , shared_ptr<random_type>
                               , vector_type const&
                            >())
                            .property("slab", &lattice::slab)
                            .property("module_name", &module_name_wrapper<dimension, float_type>)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_host_positions_lattice(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    lattice<3, double>::luaopen(L);
    lattice<2, double>::luaopen(L);
#else
    lattice<3, float>::luaopen(L);
    lattice<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class lattice<3, double>;
template class lattice<2, double>;
#else
template class lattice<3, float>;
template class lattice<2, float>;
#endif

}}} // namespace mdsim::host::positions

} // namespace halmd
