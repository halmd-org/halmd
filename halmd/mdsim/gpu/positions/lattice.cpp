/*
 * Copyright © 2008-2011  Felix Höfling
 * Copyright © 2013       Nicolas Höft
 * Copyright © 2008-2011  Peter Colberg
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

#include <algorithm>
#include <boost/array.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <cmath>
#include <limits>
#include <numeric>

#include <halmd/mdsim/gpu/positions/lattice_kernel.hpp>
#include <halmd/mdsim/gpu/positions/lattice.hpp>
#include <halmd/mdsim/positions/lattice_primitive.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/gpu/configure_kernel.hpp>

using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace positions {

template <int dimension, typename float_type>
lattice<dimension, float_type>::lattice(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<box_type const> box
  , typename box_type::vector_type const& slab
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , logger_(logger)
  , slab_(slab)
{
    if (*min_element(slab_.begin(), slab_.end()) <= 0 ||
        *max_element(slab_.begin(), slab_.end()) > 1
       ) {
        throw std::logic_error("slab extents must be a fraction between 0 and 1");
    }

    if (*min_element(slab_.begin(), slab_.end()) < 1) {
        LOG("restrict lattice positions to centred slab: " << slab_);
    }
}

template <int dimension, typename float_type>
void lattice<dimension, float_type>::set()
{
    auto position = make_cache_mutable(particle_->position());
    auto image = make_cache_mutable(particle_->image());

    // assign fcc lattice points to a fraction of the particles in a slab at the centre
    gpu_vector_type length = static_cast<gpu_vector_type>(element_prod(box_->length(), slab_));
    gpu_vector_type offset = -length / 2;

    fcc(position->data(), particle_->nparticle(), length, offset);

    // reset particle image vectors
    cuda::memset(image->begin(), image->begin() + image->capacity(), 0);
}

/**
 * Place particles on a face-centreed cubic (fcc) lattice
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
template <typename position_iterator>
void lattice<dimension, float_type>::fcc(
    position_iterator first, size_t npart
  , gpu_vector_type const& length, gpu_vector_type const& offset
)
{
    typedef close_packed_lattice<vector_type, index_type> lattice_type;

    scoped_timer_type timer(runtime_.set);
    // determine maximal lattice constant
    // use the same floating point precision as the CUDA device,
    // assign lattice coordinates to (sub-)volume of the box
    LOG_DEBUG("generating fcc lattice for " << npart << " particles, box: " << length << ", offset: " << offset);
    float u = lattice_type(1).size();
    float V = accumulate(
        length.begin(), length.end()
      , 1.f / ceil(npart / u)
      , multiplies<float>()
    );
    float a = pow(V, 1.f / dimension);
    index_type n(length / a);
    while (npart > u * accumulate(n.begin(), n.end(), 1, multiplies<unsigned int>())) {
        gpu_vector_type t;
        for (size_t i = 0; i < dimension; ++i) {
            t[i] = length[i] / (n[i] + 1);
        }
        typename gpu_vector_type::iterator it = max_element(t.begin(), t.end());
        a = *it;
        // recompute n to preserve aspect ratios of box, ensure that
        // no compoment decreases and that at least one component
        // is incremented
        index_type m = n;
        n = element_max(m, static_cast<index_type>(length / a));
        if (m == n) {
            n += index_type(1);
        }
    }
    LOG("placing particles on fcc lattice: a = " << a);
    LOG("number of fcc unit cells: " << n);

    unsigned int N = static_cast<unsigned int>(
        u * accumulate(n.begin(), n.end(), 1, multiplies<unsigned int>())
    );
    if (N > npart) {
        LOG_WARNING("lattice not fully occupied (" << N << " sites)");
    }

    // insert a vacancy every 'skip' sites
    unsigned int skip = (N - npart) ? static_cast<unsigned int>(ceil(static_cast<double>(N) / (N - npart))) : 0;
    if (skip) {
        LOG_INFO("insert a vacancy after every " << skip << " sites");
    }

    try {
        auto& lattice_kernel = get_lattice_kernel<float_type, lattice_type>().lattice;
        configure_kernel(lattice_kernel, particle_->dim(), false);
        lattice_kernel(first, npart, a, skip, offset, n);
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to generate particle lattice on GPU");
        throw;
    }
    LOG_DEBUG("number of particles inserted: " << npart);
}

template <int dimension, typename float_type>
void lattice<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("positions")
            [
                class_<lattice>()
                    .property("slab", &lattice::slab)
                    .def("set", &lattice::set)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("set", &runtime::set)
                    ]
                    .def_readonly("runtime", &lattice::runtime_)
              , def("lattice", &std::make_shared<lattice
                  , std::shared_ptr<particle_type>
                  , std::shared_ptr<box_type const>
                  , typename box_type::vector_type const&
                  , std::shared_ptr<logger>
                  >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_positions_lattice(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    lattice<3, float>::luaopen(L);
    lattice<2, float>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    lattice<3, dsfloat>::luaopen(L);
    lattice<2, dsfloat>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class lattice<3, float>;
template class lattice<2, float>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class lattice<3, dsfloat>;
template class lattice<2, dsfloat>;
#endif

} // namespace mdsim
} // namespace gpu
} // namespace positions
} // namespace halmd
