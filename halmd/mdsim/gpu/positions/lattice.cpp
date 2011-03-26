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
#include <halmd/mdsim/gpu/positions/lattice_kernel.hpp>
#include <halmd/mdsim/gpu/positions/lattice.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu { namespace positions
{

template <int dimension, typename float_type, typename RandomNumberGenerator>
lattice<dimension, float_type, RandomNumberGenerator>::lattice(
    shared_ptr<particle_type> particle
  , shared_ptr<box_type> box
  , shared_ptr<random_type> random
  , typename box_type::vector_type const& slab
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

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type, typename RandomNumberGenerator>
void lattice<dimension, float_type, RandomNumberGenerator>::register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.set, "set", "setting particle positions on lattice");
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void lattice<dimension, float_type, RandomNumberGenerator>::set()
{
#ifdef USE_VERLET_DSFUN
    // set hi parts of dsfloat values to zero
    cuda::memset(particle->g_r, 0, particle->g_r.capacity());
#endif

    assert(particle->g_r.size() == particle->nbox);

    // assign fcc lattice points to a fraction of the particles in a slab at the centre
    gpu_vector_type length = static_cast<gpu_vector_type>(element_prod(box->length(), slab_));
    gpu_vector_type offset = -length / 2;

    float4* r_it = particle->g_r.data(); // use pointer as substitute for missing iterator
    fcc(r_it, r_it + particle->nbox, length, offset);

    // randomise particle positions if there is more than 1 particle type
    // FIXME this requires a subsequent sort
    // FIXME this will fail greatly once we support polymers
    if (particle->ntypes.size() > 1) {
        LOG("randomly permuting particle positions");
        random->shuffle(particle->g_r);
    }

    // reset particle image vectors
    cuda::memset(particle->g_image, 0, particle->g_image.capacity());
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
template <int dimension, typename float_type, typename RandomNumberGenerator>
template <typename position_iterator>
void lattice<dimension, float_type, RandomNumberGenerator>::fcc(
    position_iterator first, position_iterator last
  , gpu_vector_type const& length, gpu_vector_type const& offset
)
{
    // determine maximal lattice constant
    // use the same floating point precision as the CUDA device,
    // assign lattice coordinates to (sub-)volume of the box
    LOG_TRACE("generating fcc lattice for " << last - first << " particles, box: " << length << ", offset: " << offset);
    size_t npart = last - first;
    float_type u = (dimension == 3) ? 4 : 2;
    float_type V = accumulate(
        length.begin(), length.end()
      , float_type(1) / ceil(npart / u)
      , multiplies<float_type>()
    );
    float_type a = pow(V, float_type(1) / dimension);
    index_type n(length / a);
    while (npart > u * accumulate(n.begin(), n.end(), 1, multiplies<unsigned int>())) {
        gpu_vector_type t;
        for (size_t i = 0; i < dimension; ++i) {
            t[i] = length[i] / (n[i] + 1);
        }
        typename gpu_vector_type::iterator it = max_element(t.begin(), t.end());
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
        LOG_TRACE("insert a vacancy after every " << skip << " sites");
    }

    // set kernel globals in constant memory
    lattice_wrapper<dimension> const& kernel = get_lattice_kernel<dimension>();
    cuda::copy(offset, kernel.offset);
    cuda::copy(n, kernel.ncell);

    cuda::thread::synchronize();
    try {
        scoped_timer<timer> timer_(runtime_.set);
        cuda::configure(particle->dim.grid, particle->dim.block);
        kernel.fcc(first, npart, a, skip);
        cuda::thread::synchronize();
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to generate particle lattice on GPU");
        throw;
    }
    LOG_DEBUG("number of particles inserted: " << npart);
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
static char const* module_name_wrapper(lattice<dimension, float_type, RandomNumberGenerator> const&)
{
    return lattice<dimension, float_type, RandomNumberGenerator>::module_name();
}

template <int dimension, typename float_type, typename RandomNumberGenerator>
void lattice<dimension, float_type, RandomNumberGenerator>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name(module_name() + ("_" + lexical_cast<string>(dimension) + "_"));
    module(L)
    [
        namespace_("halmd_wrapper")
        [
            namespace_("mdsim")
            [
                namespace_("gpu")
                [
                    namespace_("positions")
                    [
                        class_<lattice, shared_ptr<_Base>, bases<_Base> >(class_name.c_str())
                            .def(constructor<
                                 shared_ptr<particle_type>
                               , shared_ptr<box_type>
                               , shared_ptr<random_type>
                               , typename box_type::vector_type const&
                             >())
                            .property("slab", &lattice::slab)
                            .def("register_runtimes", &lattice::register_runtimes)
                            .property("module_name", &module_name_wrapper<dimension, float_type, RandomNumberGenerator>)
                    ]
                ]
            ]
        ]
    ];
}

HALMD_INIT( register_luaopen )
{
    lua_wrapper::register_(1) //< distance of derived to base class
    [
        &lattice<3, float, random::gpu::rand48>::luaopen
    ]
    [
        &lattice<2, float, random::gpu::rand48>::luaopen
    ];
}

// explicit instantiation
template class lattice<3, float, random::gpu::rand48>;
template class lattice<2, float, random::gpu::rand48>;

}}} // namespace mdsim::gpu::positions

} // namespace halmd
